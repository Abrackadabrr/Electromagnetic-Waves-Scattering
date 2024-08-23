//
// Created by evgen on 04.02.24.
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

using namespace EMW;
using namespace EMW::Types;

template<typename Range1, typename Range2, typename OutputIterator>
void cartesian_product(Range1 const &r1, Range2 const &r2, OutputIterator out, scalar h) {
    using std::begin;
    using std::end;

    for (auto i = begin(r1); i != end(r1); ++i) {
        for (auto j = begin(r2); j != end(r2); ++j) {
            *out++ = Types::Vector3d{static_cast<Types::scalar>(*j) - 1. / 2, static_cast<Types::scalar>(*i) - 1. / 2,
                                     0} * h;
        }
    }
}

TEST(MARTIX, MATRIX_COEFFICIENTS) {
    int N = 2;
    std::vector<Mesh::point_t> meshgrid;
    meshgrid.reserve(N * N);
    cartesian_product(std::ranges::views::iota(0, N), std::ranges::views::iota(0, N), std::back_inserter(meshgrid), 2);

    const auto cellsView = std::views::iota(0, (N - 1) * (N - 1)) | std::views::transform(
            [N](int index) {
                Types::index i = index + index / (N - 1);
                const auto point = Mesh::IndexedCell::nodes_t{i, i + 1, i + 1 + N, i + N};
                return point;
            }
    );

    const auto cells = Containers::vector<Mesh::IndexedCell::nodes_t>{std::ranges::begin(cellsView),
                                                                      std::ranges::end(cellsView)};

    const Mesh::SurfaceMesh mesh{meshgrid, cells};

    // Тест на расчет K0
    const auto k0_value = Matrix::getZeroPartIntegral(0, 0, 1, mesh.getCells());
    const scalar inner_error = norm(k0_value(1, 1) - k0_value(0, 0));
    ASSERT_NEAR(inner_error, 0, 1e-15);
    std::cout << "K0 value: " << k0_value(1, 1) << std::endl;
    const scalar analytical_error = norm(
            k0_value(1, 1) - complex_d{-2 * 0.1615145224880715150268185498864132, -2 * 0.04632231448567410781500901437237300});
    ASSERT_NEAR(analytical_error, 0, 1e-6);

    // Teст на расчет К1
    const auto k1_value = Matrix::getFirstPartIntegral(0, 0, 1, mesh.getCells());
    std::cout << "K1 value: " << k1_value << std::endl;
}
# if 0
TEST(MARTIX, MATRIX) {
    int N = 100;
    std::vector<Mesh::point_t> meshgrid;
    meshgrid.reserve(N * N);
    cartesian_product(std::ranges::views::iota(0, N), std::ranges::views::iota(0, N), std::back_inserter(meshgrid), 0.05);

    const auto cellsView = std::views::iota(0, (N - 1) * (N - 1)) | std::views::transform(
            [N](int index) {
                Types::index i = index + index / (N - 1);
                const auto point = Mesh::IndexedCell::nodes_t{i, i + 1, i + 1 + N, i + N};
                return point;
            }
    );

    const auto cells = Containers::vector<Mesh::IndexedCell::nodes_t>{std::ranges::begin(cellsView),
                                                                      std::ranges::end(cellsView)};

    const Mesh::SurfaceMesh mesh{meshgrid, cells};

    const auto k1_value = Matrix::getMatrix(Types::complex_d{1, 0}, mesh.getCells());

//    const scalar error = norm(k1_value - complex_d{1, 1});
//    ASSERT_NEAR(error, 0, 1e-12);
//    std::cout << k1_value.norm() * k1_value.inverse().norm() << std::endl;
}
#endif