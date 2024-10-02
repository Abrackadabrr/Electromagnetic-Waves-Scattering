//
// Created by evgen on 05.02.24.
//

#include <vector>
#include <utility>
#include <tuple>
#include <ranges>

#include "gtest/gtest.h"
#include "types/Types.hpp"
#include "mesh/SurfaceMesh.hpp"
#include "mesh/MeshTypes.hpp"
#include "slae_generation/MatrixGeneration.hpp"
#include "experiment/PhysicalCondition.hpp"
#include "math/MathConstants.hpp"

#include <Eigen/Core>
#include <Eigen/IterativeLinearSolvers>
#include <unsupported/Eigen/IterativeSolvers>

using namespace EMW;

const Types::scalar h = 0.1;

template<typename Range1, typename Range2, typename OutputIterator>
void cartesian_product(Range1 const &r1, Range2 const &r2, OutputIterator out) {
    using std::begin;
    using std::end;

    for (auto i = begin(r1); i != end(r1); ++i) {
        for (auto j = begin(r2); j != end(r2); ++j) {
            *out++ = Types::Vector3d{static_cast<Types::scalar>(*j), static_cast<Types::scalar>(*i), 0} * h;
        }
    }
}

TEST(MARTIX, MATRIX_SMALL) {
    int N = 3;
    std::vector<Mesh::point_t> meshgrid;
    meshgrid.reserve(N * N);
    cartesian_product(std::ranges::views::iota(0, N), std::ranges::views::iota(0, N), std::back_inserter(meshgrid));

    const auto cellsView = std::views::iota(0, (N - 1) * (N - 1)) | std::views::transform(
            [N](int index) {
                Types::index i = index + index / (N - 1);
                const auto point = Mesh::IndexedCell::nodes_t{i, i + 1, i + 1 + N, i + N};
                return point;
            }
    );

    const auto cells = Containers::vector<Mesh::IndexedCell::nodes_t>{std::ranges::begin(cellsView),
                                                                      std::ranges::end(cellsView)};

    // физика
    EMW::Physics::physicalConditionsCase physics{
            .E0 = Types::Vector3d{0, 1, 0}.normalized(),
            .k = Types::complex_d{4 * Math::Constants::PI<Types::scalar>(), 0},
            .k_vec = Types::Vector3d{0, 0, 1}.normalized()
    };

    Mesh::SurfaceMesh mesh1{meshgrid, cells};
    Mesh::SurfaceMesh mesh2{meshgrid, cells};
    const auto matrix1 = Matrix::getMatrix(1, mesh1.getCells());
    const Types::VectorXc rhs1 = Types::VectorXc{
            Matrix::getRHS(physics.E0, physics.k.real() * physics.k_vec, mesh1.getCells())};

    mesh2.basisHack();
    const auto matrix2 = Matrix::getMatrix(1, mesh2.getCells());
    const Types::VectorXc rhs2 = Types::VectorXc{
            Matrix::getRHS(physics.E0, physics.k.real() * physics.k_vec, mesh2.getCells())};

    // решаем обе системы и смотрим на получившиеся токи
    auto method = Eigen::GMRES<Types::MatrixXc>{};
    method.setTolerance(1e-5);

    method.compute(matrix1);
    const auto j1 = Types::VectorXc{method.solve(rhs1)};

    method.compute(matrix2);
    const auto j2 = Types::VectorXc{method.solve(rhs2)};

    mesh1.fillJ(j1);
    mesh2.fillJ(j2);
    int a = 0;
}
