#include "rotation_optimizer.hpp"
#include "stl_mesh.hpp"

#include <algorithm>
#include <cmath>
#include <exception>
#include <iomanip>
#include <iostream>
#include <string>

namespace {

void print_usage(const char* executable) {
    std::cerr << "Usage: " << executable
              << " <input.stl> <output.stl> <visualisation.ply> [overhang_degrees]" << std::endl;
}

}  // namespace

int main(int argc, char** argv) {
    if (argc < 4) {
        print_usage(argv[0]);
        return 1;
    }

    const std::string input_path = argv[1];
    const std::string output_stl_path = argv[2];
    const std::string output_ply_path = argv[3];
    double overhang_degrees = 30.0;
    if (argc >= 5) {
        overhang_degrees = std::stod(argv[4]);
    }

    try {
        const s4::TriangleMesh mesh = s4::load_stl(input_path);
        const s4::PathGradientInput gradient_input =
            s4::build_path_gradient_input(mesh, 0.02, overhang_degrees, 3);
        const s4::PathGradientResult result = s4::compute_path_length_to_base_gradient(gradient_input);

        s4::write_binary_stl(mesh, output_stl_path);
        s4::write_colored_ply(mesh, result.gradient, output_ply_path);

        std::cout << "Processed mesh with " << mesh.faces.size() << " triangles and "
                  << mesh.vertices.size() << " unique vertices." << std::endl;
        std::size_t reachable = 0;
        double max_distance = 0.0;
        for (double value : result.distance_to_bottom) {
            if (std::isfinite(value)) {
                ++reachable;
                max_distance = std::max(max_distance, value);
            }
        }
        std::cout << "Reachable faces: " << reachable << "/" << mesh.faces.size()
                  << ", max distance: " << std::fixed << std::setprecision(3) << max_distance << " units"
                  << std::endl;
        std::cout << "Visualisation written to " << output_ply_path << std::endl;
        std::cout << "Binary STL written to " << output_stl_path << std::endl;
    } catch (const std::exception& error) {
        std::cerr << "Error: " << error.what() << std::endl;
        return 1;
    }

    return 0;
}

