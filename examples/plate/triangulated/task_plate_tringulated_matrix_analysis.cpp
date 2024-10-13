//
// Created by evgen on 05.10.2024.
//

#include "experiment/PhysicalCondition.hpp"
#include "math/MathConstants.hpp"
#include "mesh/Parser.hpp"
#include "mesh/VolumeMesh.hpp"
#include "meshes/plate/PlateGrid.hpp"
#include "slae_generation/MatrixGeneration.hpp"
#include "visualisation/VTKFunctions.hpp"

#include <Eigen/Core>
#include <Eigen/IterativeLinearSolvers>
#include <unsupported/Eigen/IterativeSolvers>
#include <Eigen/Eigenvalues>

#include <iostream>

using namespace EMW;
using namespace EMW::Types;

std::array<Types::Vector3d, 3> equalLocalBasis_1(const Mesh::IndexedCell &cell) {
    return {Vector3d{1, 0, 0}, Vector3d{0, 1, 0}, cell.normal};
}

std::array<Types::Vector3d, 3> equalLocalBasis_2(const Mesh::IndexedCell &cell) {
    return {Vector3d{1, 1, 0}.normalized(), Vector3d{-1, 1, 0}.normalized(), cell.normal};
}

int main() {
    // сетка на пластинке
    const std::string nodesFile = "/home/evgen/Education/MasterDegree/thesis/Electromagnetic-Waves-Scattering/meshes/"
                                  "plate/triangulated/1191_nodes.csv";
    const std::string cellsFile = "/home/evgen/Education/MasterDegree/thesis/Electromagnetic-Waves-Scattering/meshes/"
                                  "plate/triangulated/2256_cells.csv";
    const EMW::Types::index nNodes = 1191;
    const EMW::Types::index nCells = 2256;
    auto surfaceMesh = Mesh::SurfaceMesh{EMW::Parser::parseMesh(nodesFile, cellsFile, nNodes, nCells)};

    surfaceMesh.setName("surface_mesh_triangular_basis_1" + std::to_string(nNodes));
    surfaceMesh.customLocalBasis(equalLocalBasis_1);

    // физика
    EMW::Physics::planeWaveCase physics{Vector3d{0, 1, 0}.normalized(), 4 * Math::Constants::PI<scalar>(),
                                        Vector3d{0, 0, 1}.normalized()};
    const auto initial_field_function = [physics](const Mesh::point_t &point) -> Vector3c {
        return physics.value(point);
    };

    // след падающего поля на расчетной поверхности
    const Math::SurfaceVectorField incidentField(surfaceMesh, initial_field_function);

    // генерируем СЛАУ
    const VectorXc b3 = incidentField.asSLAERHS();
    const MatrixXc A3 = Matrix::getMatrix(physics.k, surfaceMesh);

    std::cout << "Matrix has been obtained" << std::endl;

    // анализ собственных чисел матрицы
    Eigen::ComplexEigenSolver<MatrixXc> eigenSolver;
    eigenSolver.compute(A3);
    std::cout << eigenSolver.info() << std::endl;
    std::cout << eigenSolver.eigenvalues() << std::endl;
}