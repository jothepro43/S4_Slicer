#include "rotation_optimizer.hpp"

#include <iomanip>
#include <iostream>

int main() {
    s4::PathGradientInput input;
    input.cell_centers = {
        {0.0, 0.0, 0.0},
        {1.0, 0.0, 1.0},
        {0.0, 1.0, 1.2},
        {-1.0, 0.0, 1.1},
        {0.0, -1.0, 0.9},
    };

    input.face_normals = {
        {0.0, 0.0, 1.0},
        {0.0, 0.0, -1.0},
        {0.1, 0.0, -0.99},
        {-0.1, 0.0, -0.99},
        {0.0, 0.2, -0.98},
    };

    input.point_neighbors = {
        {1, 2, 3, 4},
        {0, 2, 4},
        {0, 1, 3},
        {0, 2, 4},
        {0, 1, 3},
    };

    input.edge_neighbors = {
        {1, 2, 3, 4},
        {0, 2, 3, 4},
        {0, 1, 3},
        {0, 1, 2, 4},
        {0, 1, 3},
    };

    input.bottom_cells = {0};
    input.max_overhang_degrees = 30.0;
    input.initial_rotation_field_smoothing = 2;
    input.set_initial_rotation_to_zero = false;

    const s4::PathGradientResult result = s4::compute_path_length_to_base_gradient(input);

    std::cout << std::fixed << std::setprecision(6);
    std::cout << "Cell gradients:\n";
    for (std::size_t i = 0; i < result.gradient.size(); ++i) {
        std::cout << "  " << i << ": " << result.gradient[i] << '\n';
    }

    std::cout << "\nDistances to bottom:\n";
    for (std::size_t i = 0; i < result.distance_to_bottom.size(); ++i) {
        std::cout << "  " << i << ": " << result.distance_to_bottom[i] << '\n';
    }

    std::cout << "\nClosest bottom indices:\n";
    for (std::size_t i = 0; i < result.closest_bottom_index.size(); ++i) {
        std::cout << "  " << i << ": " << result.closest_bottom_index[i] << '\n';
    }

    return 0;
}

