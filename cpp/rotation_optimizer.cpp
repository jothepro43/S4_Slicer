#include "rotation_optimizer.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <queue>
#include <thread>
#include <utility>

namespace s4 {
namespace {

constexpr double kPi = 3.14159265358979323846;

template <typename Func>
void parallel_for(std::size_t count, Func&& func) {
    if (count == 0) {
        return;
    }
    const std::size_t hardware_threads = std::max<std::size_t>(1, std::thread::hardware_concurrency());
    const std::size_t chunk = (count + hardware_threads - 1) / hardware_threads;
    std::vector<std::thread> workers;
    workers.reserve(hardware_threads);
    for (std::size_t index = 0; index < hardware_threads; ++index) {
        const std::size_t start = index * chunk;
        if (start >= count) {
            break;
        }
        const std::size_t end = std::min(count, start + chunk);
        workers.emplace_back([start, end, &func]() {
            for (std::size_t i = start; i < end; ++i) {
                func(i);
            }
        });
    }
    for (auto& worker : workers) {
        worker.join();
    }
}

struct QueueEntry {
    double distance;
    int node;
    int bottom;
};

inline double clamp(double value, double min_value, double max_value) {
    if (value < min_value) {
        return min_value;
    }
    if (value > max_value) {
        return max_value;
    }
    return value;
}

std::array<double, 2> xy_components(const std::array<double, 3>& value) {
    return {value[0], value[1]};
}

double euclidean_distance(const std::array<double, 3>& a, const std::array<double, 3>& b) {
    const double dx = a[0] - b[0];
    const double dy = a[1] - b[1];
    const double dz = a[2] - b[2];
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

double norm2d(const std::array<double, 2>& v) {
    return std::hypot(v[0], v[1]);
}

std::array<double, 2> normalize2d(const std::array<double, 2>& v) {
    const double n = norm2d(v);
    if (n <= 1e-12) {
        return {0.0, 0.0};
    }
    return {v[0] / n, v[1] / n};
}

std::array<double, 3> normalize3d(const std::array<double, 3>& v) {
    const double n = std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
    if (n <= 1e-12) {
        return {0.0, 0.0, 0.0};
    }
    return {v[0] / n, v[1] / n, v[2] / n};
}

double dot2d(const std::array<double, 2>& a, const std::array<double, 2>& b) {
    return a[0] * b[0] + a[1] * b[1];
}

struct PlaneParameters {
    double a;
    double b;
    double c;
    bool valid;
};

PlaneParameters fit_plane(const std::vector<int>& indices,
                          const std::vector<std::array<double, 3>>& cell_centers,
                          const std::vector<double>& distances) {
    const std::size_t count = indices.size();
    if (count < 3) {
        return {0.0, 0.0, 0.0, false};
    }

    double sum_x = 0.0;
    double sum_y = 0.0;
    double sum_z = 0.0;
    double sum_xx = 0.0;
    double sum_xy = 0.0;
    double sum_yy = 0.0;
    double sum_xz = 0.0;
    double sum_yz = 0.0;

    for (const int index : indices) {
        const double x = cell_centers[index][0];
        const double y = cell_centers[index][1];
        const double z = distances[index];

        sum_x += x;
        sum_y += y;
        sum_z += z;
        sum_xx += x * x;
        sum_xy += x * y;
        sum_yy += y * y;
        sum_xz += x * z;
        sum_yz += y * z;
    }

    double matrix[3][4] = {
        {sum_xx, sum_xy, sum_x, sum_xz},
        {sum_xy, sum_yy, sum_y, sum_yz},
        {sum_x, sum_y, static_cast<double>(count), sum_z},
    };

    for (int pivot = 0; pivot < 3; ++pivot) {
        int best_row = pivot;
        for (int row = pivot + 1; row < 3; ++row) {
            if (std::fabs(matrix[row][pivot]) > std::fabs(matrix[best_row][pivot])) {
                best_row = row;
            }
        }
        if (std::fabs(matrix[best_row][pivot]) <= 1e-12) {
            return {0.0, 0.0, 0.0, false};
        }
        if (best_row != pivot) {
            for (int col = pivot; col < 4; ++col) {
                std::swap(matrix[pivot][col], matrix[best_row][col]);
            }
        }

        const double pivot_value = matrix[pivot][pivot];
        for (int col = pivot; col < 4; ++col) {
            matrix[pivot][col] /= pivot_value;
        }
        for (int row = 0; row < 3; ++row) {
            if (row == pivot) {
                continue;
            }
            const double factor = matrix[row][pivot];
            for (int col = pivot; col < 4; ++col) {
                matrix[row][col] -= factor * matrix[pivot][col];
            }
        }
    }

    const double a = matrix[0][3];
    const double b = matrix[1][3];
    const double c = matrix[2][3];
    if (!std::isfinite(a) || !std::isfinite(b) || !std::isfinite(c)) {
        return {0.0, 0.0, 0.0, false};
    }
    return {a, b, c, true};
}

struct DijkstraResult {
    std::vector<double> distance;
    std::vector<int> closest_bottom;
};

DijkstraResult multi_source_dijkstra(const PathGradientInput& input) {
    const std::size_t count = input.cell_centers.size();
    DijkstraResult result;
    result.distance.assign(count, std::numeric_limits<double>::infinity());
    result.closest_bottom.assign(count, -1);

    auto cmp = [](const QueueEntry& lhs, const QueueEntry& rhs) {
        return lhs.distance > rhs.distance;
    };
    std::priority_queue<QueueEntry, std::vector<QueueEntry>, decltype(cmp)> queue(cmp);

    for (int bottom : input.bottom_cells) {
        if (bottom < 0 || static_cast<std::size_t>(bottom) >= count) {
            continue;
        }
        result.distance[bottom] = 0.0;
        result.closest_bottom[bottom] = bottom;
        queue.push({0.0, bottom, bottom});
    }

    while (!queue.empty()) {
        const QueueEntry entry = queue.top();
        queue.pop();
        if (entry.distance > result.distance[entry.node] + 1e-12) {
            continue;
        }
        const auto& neighbours = input.point_neighbors[entry.node];
        for (int neighbour : neighbours) {
            if (neighbour < 0 || static_cast<std::size_t>(neighbour) >= count) {
                continue;
            }
            const double weight = euclidean_distance(input.cell_centers[entry.node],
                                                     input.cell_centers[neighbour]);
            const double candidate = entry.distance + weight;
            if (candidate + 1e-12 < result.distance[neighbour]) {
                result.distance[neighbour] = candidate;
                result.closest_bottom[neighbour] = entry.bottom;
                queue.push({candidate, neighbour, entry.bottom});
            }
        }
    }

    return result;
}

std::vector<int> build_bottom_mask(const std::vector<int>& bottom_cells, std::size_t count) {
    std::vector<int> mask(count, 0);
    for (int bottom : bottom_cells) {
        if (bottom >= 0 && static_cast<std::size_t>(bottom) < count) {
            mask[bottom] = 1;
        }
    }
    return mask;
}

std::vector<int> unique_indices(std::vector<int> values) {
    std::sort(values.begin(), values.end());
    values.erase(std::unique(values.begin(), values.end()), values.end());
    return values;
}

}  // namespace

PathGradientResult compute_path_length_to_base_gradient(const PathGradientInput& input) {
    PathGradientResult output;
    const std::size_t count = input.cell_centers.size();
    output.gradient.assign(count, 0.0);
    output.distance_to_bottom.assign(count, std::numeric_limits<double>::quiet_NaN());
    output.closest_bottom_index.assign(count, -1);

    if (count == 0) {
        return output;
    }

    const auto bottom_mask = build_bottom_mask(input.bottom_cells, count);
    const DijkstraResult dijkstra_result = multi_source_dijkstra(input);

    const double max_overhang_radians = (input.max_overhang_degrees + 90.0) * kPi / 180.0;

    for (std::size_t index = 0; index < count; ++index) {
        const double dot_with_up = input.face_normals[index][2];
        const double angle_cosine = clamp(dot_with_up, -1.0, 1.0);
        const double angle = std::acos(angle_cosine);
        const bool is_overhang = angle > max_overhang_radians;
        if (is_overhang && !bottom_mask[index] && dijkstra_result.closest_bottom[index] != -1) {
            output.distance_to_bottom[index] = dijkstra_result.distance[index];
            output.closest_bottom_index[index] = dijkstra_result.closest_bottom[index];
        }
    }

    for (std::size_t index = 0; index < count; ++index) {
        if (!std::isfinite(output.distance_to_bottom[index])) {
            continue;
        }

        std::vector<int> local_cells = input.edge_neighbors[index];
        local_cells.push_back(static_cast<int>(index));
        local_cells = unique_indices(std::move(local_cells));

        std::vector<int> valid_cells;
        valid_cells.reserve(local_cells.size());
        for (int cell : local_cells) {
            if (cell >= 0 && static_cast<std::size_t>(cell) < count &&
                std::isfinite(output.distance_to_bottom[cell])) {
                valid_cells.push_back(cell);
            }
        }

        if (valid_cells.size() < 3) {
            const int closest_bottom = output.closest_bottom_index[index];
            if (closest_bottom < 0 || static_cast<std::size_t>(closest_bottom) >= count) {
                continue;
            }
            const auto location_to_roll_to = xy_components(input.cell_centers[closest_bottom]);
            auto direction_to_bottom = normalize2d({location_to_roll_to[0] - input.cell_centers[index][0],
                                                    location_to_roll_to[1] - input.cell_centers[index][1]});
            auto cell_center = normalize2d(xy_components(input.cell_centers[index]));
            const double dot_value = dot2d(cell_center, direction_to_bottom);
            if (std::fabs(dot_value) <= 1e-12 || std::isnan(dot_value)) {
                output.gradient[index] = 0.0;
            } else {
                output.gradient[index] = dot_value / std::fabs(dot_value);
            }
            continue;
        }

        const PlaneParameters plane = fit_plane(valid_cells, input.cell_centers, output.distance_to_bottom);
        if (!plane.valid) {
            double sum = 0.0;
            int count_valid = 0;
            for (int cell : valid_cells) {
                const double value = output.gradient[cell];
                if (!std::isnan(value)) {
                    sum += value;
                    ++count_valid;
                }
            }
            output.gradient[index] = (count_valid > 0) ? sum / count_valid : 0.0;
            continue;
        }

        const auto plane_normal = normalize3d({plane.a, plane.b, -1.0});
        const auto cell_center_direction = normalize2d(xy_components(input.cell_centers[index]));
        const double gradient_value = dot2d(cell_center_direction, {plane_normal[0], plane_normal[1]});
        if (std::isnan(gradient_value)) {
            double sum = 0.0;
            int count_valid = 0;
            for (int cell : valid_cells) {
                const double value = output.gradient[cell];
                if (!std::isnan(value)) {
                    sum += value;
                    ++count_valid;
                }
            }
            output.gradient[index] = (count_valid > 0) ? sum / count_valid : 0.0;
        } else {
            output.gradient[index] = gradient_value;
        }
    }

    if (input.initial_rotation_field_smoothing > 0) {
        std::vector<double> smoothed(count, 0.0);
        for (int iteration = 0; iteration < input.initial_rotation_field_smoothing; ++iteration) {
            std::fill(smoothed.begin(), smoothed.end(), 0.0);
            parallel_for(count, [&](std::size_t index) {
                const double value = output.gradient[index];
                if (value == 0.0 || std::isnan(value)) {
                    smoothed[index] = value;
                    return;
                }
                std::vector<int> local_cells = input.point_neighbors[index];
                for (int neighbour : input.point_neighbors[index]) {
                    const auto& neighbour_neighbours = input.point_neighbors[neighbour];
                    local_cells.insert(local_cells.end(), neighbour_neighbours.begin(), neighbour_neighbours.end());
                }
                local_cells.push_back(static_cast<int>(index));
                local_cells = unique_indices(std::move(local_cells));

                double sum = 0.0;
                int count_valid = 0;
                for (int cell : local_cells) {
                    if (cell < 0 || static_cast<std::size_t>(cell) >= count) {
                        continue;
                    }
                    const double neighbour_value = output.gradient[cell];
                    if (neighbour_value != 0.0 && !std::isnan(neighbour_value)) {
                        sum += neighbour_value;
                        ++count_valid;
                    }
                }
                smoothed[index] = (count_valid > 0) ? sum / count_valid : value;
            });
            output.gradient = smoothed;
        }
    }

    if (!input.set_initial_rotation_to_zero) {
        for (double& value : output.gradient) {
            if (std::fabs(value) <= 1e-12) {
                value = std::numeric_limits<double>::quiet_NaN();
            }
        }
    }

    return output;
}

}  // namespace s4

