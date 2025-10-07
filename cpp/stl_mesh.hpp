#pragma once

#include <array>
#include <string>
#include <vector>

#include "rotation_optimizer.hpp"

namespace s4 {

struct TriangleMesh {
    std::vector<std::array<double, 3>> vertices;
    std::vector<std::array<int, 3>> faces;
    std::vector<std::array<double, 3>> face_normals;
    std::vector<std::array<double, 3>> face_centroids;
};

TriangleMesh load_stl(const std::string& path);
void write_binary_stl(const TriangleMesh& mesh, const std::string& path);
void write_colored_ply(const TriangleMesh& mesh,
                       const std::vector<double>& face_values,
                       const std::string& path);

PathGradientInput build_path_gradient_input(const TriangleMesh& mesh,
                                            double bottom_slice_fraction = 0.02,
                                            double overhang_angle_degrees = 30.0,
                                            int smoothing_iterations = 2);

}  // namespace s4

