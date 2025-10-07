#include "stl_mesh.hpp"

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iterator>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>

namespace s4 {
namespace {

struct VertexKey {
    long long x;
    long long y;
    long long z;

    bool operator==(const VertexKey& other) const {
        return x == other.x && y == other.y && z == other.z;
    }
};

struct VertexKeyHash {
    std::size_t operator()(const VertexKey& key) const noexcept {
        const std::size_t hx = std::hash<long long>{}(key.x);
        const std::size_t hy = std::hash<long long>{}(key.y);
        const std::size_t hz = std::hash<long long>{}(key.z);
        return hx ^ (hy << 1) ^ (hz << 2);
    }
};

VertexKey make_key(const std::array<double, 3>& vertex) {
    constexpr double kScale = 1'000'000.0;  // micrometre precision
    return {static_cast<long long>(std::llround(vertex[0] * kScale)),
            static_cast<long long>(std::llround(vertex[1] * kScale)),
            static_cast<long long>(std::llround(vertex[2] * kScale))};
}

std::array<double, 3> compute_face_normal(const std::array<double, 3>& a,
                                          const std::array<double, 3>& b,
                                          const std::array<double, 3>& c) {
    const double ux = b[0] - a[0];
    const double uy = b[1] - a[1];
    const double uz = b[2] - a[2];
    const double vx = c[0] - a[0];
    const double vy = c[1] - a[1];
    const double vz = c[2] - a[2];

    std::array<double, 3> normal = {
        uy * vz - uz * vy,
        uz * vx - ux * vz,
        ux * vy - uy * vx,
    };

    const double length = std::sqrt(normal[0] * normal[0] + normal[1] * normal[1] +
                                    normal[2] * normal[2]);
    if (length <= 1e-12) {
        return {0.0, 0.0, 0.0};
    }
    normal[0] /= length;
    normal[1] /= length;
    normal[2] /= length;
    return normal;
}

std::array<double, 3> compute_centroid(const std::array<double, 3>& a,
                                       const std::array<double, 3>& b,
                                       const std::array<double, 3>& c) {
    return {(a[0] + b[0] + c[0]) / 3.0,
            (a[1] + b[1] + c[1]) / 3.0,
            (a[2] + b[2] + c[2]) / 3.0};
}

TriangleMesh mesh_from_raw_faces(
    const std::vector<std::array<double, 3>>& raw_vertices,
    const std::vector<std::array<double, 3>>& raw_normals) {
    TriangleMesh mesh;
    std::unordered_map<VertexKey, int, VertexKeyHash> vertex_map;
    mesh.vertices.reserve(raw_vertices.size());
    mesh.faces.reserve(raw_vertices.size() / 3);
    mesh.face_normals.reserve(raw_vertices.size() / 3);
    mesh.face_centroids.reserve(raw_vertices.size() / 3);

    for (std::size_t i = 0; i + 2 < raw_vertices.size(); i += 3) {
        const auto& a = raw_vertices[i];
        const auto& b = raw_vertices[i + 1];
        const auto& c = raw_vertices[i + 2];
        std::array<int, 3> face{};
        const std::array<std::array<double, 3>, 3> verts = {a, b, c};
        for (int v = 0; v < 3; ++v) {
            const VertexKey key = make_key(verts[v]);
            auto it = vertex_map.find(key);
            if (it == vertex_map.end()) {
                const int index = static_cast<int>(mesh.vertices.size());
                vertex_map.emplace(key, index);
                mesh.vertices.push_back(verts[v]);
                face[v] = index;
            } else {
                face[v] = it->second;
            }
        }
        mesh.faces.push_back(face);

        std::array<double, 3> normal = raw_normals[i / 3];
        const double length = std::sqrt(normal[0] * normal[0] + normal[1] * normal[1] +
                                        normal[2] * normal[2]);
        if (length <= 1e-6) {
            normal = compute_face_normal(a, b, c);
        } else {
            normal[0] /= length;
            normal[1] /= length;
            normal[2] /= length;
        }
        mesh.face_normals.push_back(normal);
        mesh.face_centroids.push_back(compute_centroid(a, b, c));
    }

    return mesh;
}

bool looks_like_ascii(const std::vector<char>& buffer) {
    const std::string header(buffer.begin(), buffer.begin() + std::min<std::size_t>(buffer.size(), 80));
    if (header.rfind("solid", 0) != 0) {
        return false;
    }
    for (char c : buffer) {
        if (c == '\0') {
            return false;
        }
    }
    return true;
}

TriangleMesh load_binary_stl(const std::vector<char>& buffer) {
    if (buffer.size() < 84) {
        throw std::runtime_error("STL file too small to contain header");
    }
    std::uint32_t triangle_count = 0;
    std::memcpy(&triangle_count, buffer.data() + 80, sizeof(std::uint32_t));
    const std::size_t expected_size = 84 + static_cast<std::size_t>(triangle_count) * 50;
    if (expected_size > buffer.size()) {
        throw std::runtime_error("Binary STL reports more triangles than bytes in file");
    }

    std::vector<std::array<double, 3>> vertices;
    std::vector<std::array<double, 3>> normals;
    vertices.reserve(triangle_count * 3);
    normals.reserve(triangle_count);

    std::size_t offset = 84;
    for (std::uint32_t i = 0; i < triangle_count; ++i) {
        std::array<float, 12> record{};
        std::memcpy(record.data(), buffer.data() + offset, sizeof(float) * 12);
        offset += 48;
        offset += 2;  // attribute byte count

        std::array<double, 3> normal = {record[0], record[1], record[2]};
        normals.push_back({normal[0], normal[1], normal[2]});
        for (int v = 0; v < 3; ++v) {
            vertices.push_back({record[3 + v * 3], record[4 + v * 3], record[5 + v * 3]});
        }
    }
    return mesh_from_raw_faces(vertices, normals);
}

TriangleMesh load_ascii_stl(const std::vector<char>& buffer) {
    std::string text(buffer.begin(), buffer.end());
    std::istringstream stream(text);
    std::string token;
    std::vector<std::array<double, 3>> vertices;
    std::vector<std::array<double, 3>> normals;

    std::array<double, 3> current_normal = {0.0, 0.0, 0.0};
    while (stream >> token) {
        for (char& c : token) {
            c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
        }
        if (token == "facet") {
            std::string normal_token;
            stream >> normal_token;
            if (normal_token != "normal") {
                throw std::runtime_error("Malformed ASCII STL: expected normal after facet");
            }
            stream >> current_normal[0] >> current_normal[1] >> current_normal[2];
            normals.push_back(current_normal);
        } else if (token == "vertex") {
            std::array<double, 3> vertex{};
            stream >> vertex[0] >> vertex[1] >> vertex[2];
            vertices.push_back(vertex);
        }
    }

    if (vertices.empty() || vertices.size() % 3 != 0 || normals.size() * 3 != vertices.size()) {
        throw std::runtime_error("Malformed ASCII STL: inconsistent triangle data");
    }

    return mesh_from_raw_faces(vertices, normals);
}

}  // namespace

TriangleMesh load_stl(const std::string& path) {
    std::ifstream file(path, std::ios::binary);
    if (!file) {
        throw std::runtime_error("Failed to open STL file: " + path);
    }
    const std::vector<char> buffer((std::istreambuf_iterator<char>(file)), std::istreambuf_iterator<char>());
    if (buffer.size() < 84) {
        throw std::runtime_error("STL file is too small: " + path);
    }

    if (!looks_like_ascii(buffer)) {
        return load_binary_stl(buffer);
    }

    try {
        return load_ascii_stl(buffer);
    } catch (const std::exception&) {
        return load_binary_stl(buffer);
    }
}

void write_binary_stl(const TriangleMesh& mesh, const std::string& path) {
    std::ofstream file(path, std::ios::binary);
    if (!file) {
        throw std::runtime_error("Failed to open output STL file: " + path);
    }

    char header[80] = {};
    std::snprintf(header, sizeof(header), "Generated by S4 mesh processor");
    file.write(header, sizeof(header));

    const std::uint32_t triangle_count = static_cast<std::uint32_t>(mesh.faces.size());
    file.write(reinterpret_cast<const char*>(&triangle_count), sizeof(triangle_count));

    for (std::size_t i = 0; i < mesh.faces.size(); ++i) {
        const auto& face = mesh.faces[i];
        const auto& normal = (i < mesh.face_normals.size()) ? mesh.face_normals[i] : std::array<double, 3>{0.0, 0.0, 1.0};
        const auto write_float = [&file](double value) {
            const float f = static_cast<float>(value);
            file.write(reinterpret_cast<const char*>(&f), sizeof(float));
        };

        write_float(normal[0]);
        write_float(normal[1]);
        write_float(normal[2]);

        for (int v = 0; v < 3; ++v) {
            const auto& vertex = mesh.vertices[face[v]];
            write_float(vertex[0]);
            write_float(vertex[1]);
            write_float(vertex[2]);
        }

        const std::uint16_t attribute_byte_count = 0;
        file.write(reinterpret_cast<const char*>(&attribute_byte_count), sizeof(attribute_byte_count));
    }
}

void write_colored_ply(const TriangleMesh& mesh,
                       const std::vector<double>& face_values,
                       const std::string& path) {
    if (mesh.faces.size() != face_values.size()) {
        throw std::runtime_error("Face values do not match mesh triangle count");
    }

    std::vector<std::array<double, 3>> vertex_colors(mesh.vertices.size(), {0.6, 0.6, 0.6});
    std::vector<int> contributions(mesh.vertices.size(), 0);

    double min_value = std::numeric_limits<double>::infinity();
    double max_value = -std::numeric_limits<double>::infinity();
    for (double value : face_values) {
        if (std::isfinite(value)) {
            min_value = std::min(min_value, value);
            max_value = std::max(max_value, value);
        }
    }
    if (!std::isfinite(min_value) || !std::isfinite(max_value) || min_value == max_value) {
        min_value = -1.0;
        max_value = 1.0;
    }

    auto value_to_color = [&](double value) {
        if (!std::isfinite(value)) {
            return std::array<double, 3>{0.5, 0.5, 0.5};
        }
        const double t = std::clamp((value - min_value) / (max_value - min_value), 0.0, 1.0);
        return std::array<double, 3>{t, 0.2, 1.0 - t};
    };

    for (std::size_t face_index = 0; face_index < mesh.faces.size(); ++face_index) {
        const auto color = value_to_color(face_values[face_index]);
        for (int vertex : mesh.faces[face_index]) {
            auto& target = vertex_colors[vertex];
            target[0] += color[0];
            target[1] += color[1];
            target[2] += color[2];
            contributions[vertex] += 1;
        }
    }

    for (std::size_t i = 0; i < mesh.vertices.size(); ++i) {
        if (contributions[i] > 0) {
            vertex_colors[i][0] /= static_cast<double>(contributions[i] + 1);
            vertex_colors[i][1] /= static_cast<double>(contributions[i] + 1);
            vertex_colors[i][2] /= static_cast<double>(contributions[i] + 1);
        }
    }

    std::ofstream file(path, std::ios::binary);
    if (!file) {
        throw std::runtime_error("Failed to open output PLY file: " + path);
    }

    file << "ply\n";
    file << "format ascii 1.0\n";
    file << "comment Generated by S4 mesh processor\n";
    file << "element vertex " << mesh.vertices.size() << "\n";
    file << "property float x\n";
    file << "property float y\n";
    file << "property float z\n";
    file << "property uchar red\n";
    file << "property uchar green\n";
    file << "property uchar blue\n";
    file << "element face " << mesh.faces.size() << "\n";
    file << "property list uchar int vertex_indices\n";
    file << "end_header\n";

    for (std::size_t i = 0; i < mesh.vertices.size(); ++i) {
        const auto& vertex = mesh.vertices[i];
        const auto& color = vertex_colors[i];
        const int r = static_cast<int>(std::clamp(color[0], 0.0, 1.0) * 255.0);
        const int g = static_cast<int>(std::clamp(color[1], 0.0, 1.0) * 255.0);
        const int b = static_cast<int>(std::clamp(color[2], 0.0, 1.0) * 255.0);
        file << static_cast<float>(vertex[0]) << ' ' << static_cast<float>(vertex[1]) << ' '
             << static_cast<float>(vertex[2]) << ' ' << r << ' ' << g << ' ' << b << "\n";
    }

    for (const auto& face : mesh.faces) {
        file << "3 " << face[0] << ' ' << face[1] << ' ' << face[2] << "\n";
    }
}

PathGradientInput build_path_gradient_input(const TriangleMesh& mesh,
                                            double bottom_slice_fraction,
                                            double overhang_angle_degrees,
                                            int smoothing_iterations) {
    PathGradientInput input;
    const std::size_t face_count = mesh.faces.size();
    input.cell_centers = mesh.face_centroids;
    input.face_normals = mesh.face_normals;
    input.point_neighbors.assign(face_count, {});
    input.edge_neighbors.assign(face_count, {});
    input.max_overhang_degrees = overhang_angle_degrees;
    input.initial_rotation_field_smoothing = smoothing_iterations;
    input.set_initial_rotation_to_zero = false;

    if (face_count == 0) {
        return input;
    }

    std::vector<std::vector<int>> vertex_to_faces(mesh.vertices.size());
    for (std::size_t face_index = 0; face_index < face_count; ++face_index) {
        for (int vertex : mesh.faces[face_index]) {
            vertex_to_faces[vertex].push_back(static_cast<int>(face_index));
        }
    }

    for (std::size_t face_index = 0; face_index < face_count; ++face_index) {
        std::vector<int> neighbours;
        const auto& face = mesh.faces[face_index];
        for (int vertex : face) {
            const auto& adjacent = vertex_to_faces[vertex];
            neighbours.insert(neighbours.end(), adjacent.begin(), adjacent.end());
        }
        std::sort(neighbours.begin(), neighbours.end());
        neighbours.erase(std::unique(neighbours.begin(), neighbours.end()), neighbours.end());
        neighbours.erase(std::remove(neighbours.begin(), neighbours.end(), static_cast<int>(face_index)),
                         neighbours.end());
        input.point_neighbors[face_index] = std::move(neighbours);
    }

    std::unordered_map<std::uint64_t, std::vector<int>> edge_map;
    edge_map.reserve(face_count * 3);
    auto make_edge_key = [](int a, int b) {
        const std::uint64_t min_v = static_cast<std::uint64_t>(std::min(a, b));
        const std::uint64_t max_v = static_cast<std::uint64_t>(std::max(a, b));
        return (min_v << 32) | max_v;
    };

    for (std::size_t face_index = 0; face_index < face_count; ++face_index) {
        const auto& face = mesh.faces[face_index];
        for (int e = 0; e < 3; ++e) {
            const int a = face[e];
            const int b = face[(e + 1) % 3];
            edge_map[make_edge_key(a, b)].push_back(static_cast<int>(face_index));
        }
    }

    for (const auto& [_, faces] : edge_map) {
        if (faces.size() < 2) {
            continue;
        }
        for (int f : faces) {
            auto& list = input.edge_neighbors[f];
            list.insert(list.end(), faces.begin(), faces.end());
        }
    }

    for (auto& neighbours : input.edge_neighbors) {
        std::sort(neighbours.begin(), neighbours.end());
        neighbours.erase(std::unique(neighbours.begin(), neighbours.end()), neighbours.end());
    }

    double min_z = std::numeric_limits<double>::infinity();
    double max_z = -std::numeric_limits<double>::infinity();
    for (const auto& centroid : mesh.face_centroids) {
        min_z = std::min(min_z, centroid[2]);
        max_z = std::max(max_z, centroid[2]);
    }
    const double range_z = std::max(max_z - min_z, 1e-6);
    const double threshold = min_z + range_z * bottom_slice_fraction;

    for (std::size_t face_index = 0; face_index < face_count; ++face_index) {
        if (mesh.face_centroids[face_index][2] <= threshold) {
            input.bottom_cells.push_back(static_cast<int>(face_index));
        }
    }

    if (input.bottom_cells.empty()) {
        const auto it = std::min_element(mesh.face_centroids.begin(), mesh.face_centroids.end(),
                                         [](const auto& a, const auto& b) { return a[2] < b[2]; });
        if (it != mesh.face_centroids.end()) {
            const int index = static_cast<int>(std::distance(mesh.face_centroids.begin(), it));
            input.bottom_cells.push_back(index);
        }
    }

    return input;
}

}  // namespace s4

