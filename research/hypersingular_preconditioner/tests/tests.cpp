//
// Created by evgen on 08.08.2025.
//

#include "../DiscreteLaplacian.hpp"
#include "../OperatorT_simple.hpp"

#include "mesh/Parser.hpp"
#include "mesh/SurfaceMesh.hpp"

#include "math/integration/gauss_quadrature/GaussLegenderPoints.hpp"

#include "../matrix/MatrixPreconditioned.hpp"
#include "../matrix/MatrixTraits.hpp"

#include "visualisation/VTKFunctions.hpp"

#include "../Solve.hpp"
#include "../matrix/SpecialPreconditioner.hpp"

#include "unsupported/Eigen/IterativeSolvers"

#include <random>

#include <string>

using namespace EMW;

int main() {
    long unsigned int N_POINTS = 50;
    long unsigned int N_INNER_CELLS = (N_POINTS - 3) * (N_POINTS - 3);
    const std::string file_nodes =
        "/home/evgen/Education/Schools/Sirius2025/meshes/" + std::to_string(N_POINTS) + "/nodes.csv";
    const std::string file_cells =
        "/home/evgen/Education/Schools/Sirius2025/meshes/" + std::to_string(N_POINTS) + "/cells.csv";
    const std::string path_to_res = "/home/evgen/Education/Schools/Sirius2025/results/";

    // делаем сетку
    const auto [nodes, cells] = EMW::Parser::parseMesh(file_nodes, file_cells);
    auto mesh = EMW::Mesh::SurfaceMesh(nodes, cells);
    mesh.setName("Test mesh");
    const Types::index n_cells = mesh.getCells().size();
    // ищем граничные узлы
    auto boundary_cells_indexes = DiscreteLaplacian::detail::get_boundary_cells_indexes(N_POINTS);
    for (int i = 0; i < mesh.getCells().size(); i++) {
        if (boundary_cells_indexes.contains(i))
            mesh.getCells()[i].tag = Mesh::IndexedCell::Tag::BOUNDARY_CELL;
        else
            mesh.getCells()[i].tag = Mesh::IndexedCell::Tag::INTERNAL_CELL;
    }

    // TEST 1: PROJECTORS
    const auto P = DiscreteLaplacian::getProjector(N_POINTS, N_POINTS);
    const auto P_i = DiscreteLaplacian::getInverseProjector(N_POINTS, N_POINTS);

    std::cout << "|| P * P_i - I ||^2 = " << (P * P_i - Types::MatrixXd::Identity(N_INNER_CELLS, N_INNER_CELLS)).squaredNorm() << std::endl;
    std::cout << "|| P_i * P - I ||^2 = " << (P_i * P - Types::MatrixXd::Identity(n_cells, n_cells)).squaredNorm() << std::endl;
    std::cout << "|| P_i * P - I ||^2 must be " << boundary_cells_indexes.size() << std::endl;

    // генирируем рандомное поле и смотрим на корректность проекции на него
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dis(0.5, 1.0);
    Math::SurfaceScalarField<Types::scalar> random_field{
        mesh, [&](const Mesh::IndexedCell& c) -> Types::scalar {return c.tag == Mesh::IndexedCell::Tag::BOUNDARY_CELL ? 0 : dis(gen);}};
    random_field.setName("Random field");
    const auto full_vector = random_field.formVector();
    const auto inner_cells_vector = random_field.formVector(Mesh::IndexedCell::Tag::INTERNAL_CELL);
    const Types::VectorXd P_full_vector = P * full_vector;
    std::cout << "|| P * f - real_values || = " << (inner_cells_vector - P_full_vector).norm() << std::endl;
    std::cout << "|| P_i * P * f - full_vector || = " << (full_vector - P_i * P_full_vector).norm() << std::endl;

    // TEST 2: правильное нахождение граничных узлов
    auto tag_field = Math::SurfaceScalarField<Types::scalar>::tagField(mesh);
    tag_field.setName("Tag field");
    auto reconstructed_field = Math::SurfaceScalarField<Types::scalar>::fromSLAESolution(mesh, P_i * P_full_vector);
    reconstructed_field.setName("Reconstructed field");
    VTK::united_snapshot(Containers::vector{tag_field, random_field, reconstructed_field}, {}, mesh, path_to_res + "paraview_snapshots/");
}
