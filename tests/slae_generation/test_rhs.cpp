//
// Created by evgen on 06.02.24.
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

TEST(SLAE, RHS) {
    int N = 3;
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

   const Vector3d E0{0, 1, 0};

   const auto b = EMW::Matrix::getRHS(E0, Vector3d{1, 0, 0}, mesh.getCells());

   std::cout << b << std::endl;
}
