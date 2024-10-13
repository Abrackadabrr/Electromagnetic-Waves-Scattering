//
// Created by evgen on 27.06.24.
//

#include "examples/pathes.hpp"
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

#include <iostream>

using namespace EMW;
using namespace EMW::Types;

std::array<Types::Vector3d, 3> equalLocalBasis_1(const Mesh::IndexedCell &cell) {
    return {Vector3d{1, 0, 0}, Vector3d{0, 1, 0}, cell.normal};
}

std::array<Types::Vector3d, 3> equalLocalBasis_2(const Mesh::IndexedCell &cell) {
    return {Vector3d{1, 1, 0}.normalized(), Vector3d{-1, 1, 0}.normalized(), cell.normal};
}

Types::Vector3d collocationPoint(const Mesh::IndexedCell &cell) {
    const auto vertex = cell.getVertexAsArray();
    return (1./3.) * (vertex[0] + vertex[1] + vertex[2]);
}

Types::Vector3d collocationPointBad(const Mesh::IndexedCell &cell) {
    const auto vertex = cell.getVertexAsArray();
    return (1./ 45) * vertex[0] + (2./ 45) * vertex[1] + (43./45) * vertex[2];
}


Math::SurfaceVectorField solve(const Math::SurfaceVectorField &incidentField, const Physics::planeWaveCase &physics, const scalar tolerance) {
    const Mesh::SurfaceMesh &surfaceMesh = incidentField.getManifold();
    const VectorXc b3 = incidentField.asSLAERHS();
    const MatrixXc A3 = Matrix::getMatrix(physics.k, surfaceMesh);

    // умное предобуславливание
    const VectorXc JacobiAux = A3.diagonal().cwiseInverse();
    const MatrixXc newMatrix = A3 * JacobiAux.asDiagonal();

    // полулчили уравнение APy = b, дальше решаем с матрицей AP
    // и потом y умножаем на P, чтобы получить x = Py (A(x) = b)

    auto method = Eigen::GMRES<MatrixXc>{};
    method.setMaxIterations(30000);
    std::cout << method.maxIterations() << std::endl;
    method.setTolerance(tolerance);
    method.set_restart(1000);
    method.compute(newMatrix);
    const VectorXc y = VectorXc{method.solve(b3)};
    const VectorXc j_vec = JacobiAux.cwiseProduct(y);

    auto j = Math::SurfaceVectorField::TangentField(surfaceMesh, j_vec);
    j.setName("j");
    std::cout << "total iterations: " << method.iterations() << std::endl;
    std::cout << "total error: " << method.error() << std::endl;
    std::cout << "Info: " << static_cast<int>(method.info()) << std::endl;
    return j;
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

    surfaceMesh.setName("triangular_basis_1_collPoints_bad");
    surfaceMesh.customLocalBasis(equalLocalBasis_1);
    surfaceMesh.customCollocationPpoints(collocationPointBad);

    // физика
    EMW::Physics::planeWaveCase physics{Vector3d{0, 1, 0}.normalized(), 4 * Math::Constants::PI<scalar>(),
                                        Vector3d{0, 0, 1}.normalized()};
    const auto initial_field_function = [physics](const Mesh::point_t &point) -> Vector3c {
        return physics.value(point);
    };

    // след падающего поля на расчетной поверхности
    const Math::SurfaceVectorField incidentField(surfaceMesh, initial_field_function);
    Math::SurfaceVectorField j = solve(incidentField, physics, 1e-3);
    j.setName("j_1");
    Math::SurfaceVectorField j_2 = solve(incidentField, physics, 1e-5);
    j_2.setName("j_2");
    auto diff = j - j_2;
    diff.setName("diff_j");
    VTK::united_snapshot({j, j_2, diff}, {}, surfaceMesh, Pathes::examples + "plane/new_discretization/");

    /*// изменение параметров сетки и повторное решение методом
    surfaceMesh.customLocalBasis(equalLocalBasis_2);
    surfaceMesh.setName("surface_mesh_triangular_basis_2" + std::to_string(nNodes));
    const auto j2 = solve(incidentField, physics);

    auto diff = j - j2;
    auto diff_norm = diff.fieldNorm();
    diff.setName("diff");
    diff_norm.setName("diff_norm");

    VTK::united_snapshot({j2, diff}, {diff_norm}, surfaceMesh, Pathes::examples + "plane/new_discretization/");*/
}
