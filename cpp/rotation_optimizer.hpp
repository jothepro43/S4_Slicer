#pragma once

#include <array>
#include <cstddef>
#include <vector>

namespace s4 {

struct PathGradientInput {
    std::vector<std::array<double, 3>> cell_centers;
    std::vector<std::array<double, 3>> face_normals;
    std::vector<std::vector<int>> point_neighbors;
    std::vector<std::vector<int>> edge_neighbors;
    std::vector<int> bottom_cells;
    double max_overhang_degrees = 0.0;
    int initial_rotation_field_smoothing = 0;
    bool set_initial_rotation_to_zero = true;
};

struct PathGradientResult {
    std::vector<double> gradient;
    std::vector<double> distance_to_bottom;
    std::vector<int> closest_bottom_index;
};

PathGradientResult compute_path_length_to_base_gradient(const PathGradientInput& input);

}  // namespace s4

