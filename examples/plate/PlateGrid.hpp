//
// Created by evgen on 30.01.24.
//

#ifndef ELECTROMAGNETIC_WAVES_SCATTERING_PLATEGRID_HPP
#define ELECTROMAGNETIC_WAVES_SCATTERING_PLATEGRID_HPP

#include <utility>
#include <tuple>
#include <ranges>
#include "types/Types.hpp"
#include "mesh/Mesh.hpp"

namespace EMW::Examples {
    template<typename Range1, typename Range2, typename OutputIterator>
    void cartesian_product(Range1 const &r1, Range2 const &r2, OutputIterator out, Types::index N, Types::scalar h) {
        using std::begin;
        using std::end;

        for (auto i = begin(r1); i != end(r1); ++i) {
            for (auto j = begin(r2); j != end(r2); ++j) {
                *out++ = Types::Vector3d{0, static_cast<Types::scalar>(*j) - static_cast<Types::scalar>(N / 2),
                                         static_cast<Types::scalar>(*i) - static_cast<Types::scalar>(N / 2)} * h;
            }
        }
    }

    template<typename Range1, typename Range2, typename OutputIterator>
    void
    cartesian_product_uneven(Range1 const &r1, Range2 const &r2, OutputIterator out, Types::index N1, Types::index N2,
                             Types::scalar h1, Types::scalar h2) {
        using std::begin;
        using std::end;

        for (auto i = begin(r2); i != end(r2); ++i) {
            for (auto j = begin(r1); j != end(r1); ++j) {
                *out++ = Types::Vector3d{0, (static_cast<Types::scalar>(*j) - static_cast<Types::scalar>(N1 / 2)) * h1,
                                         (static_cast<Types::scalar>(*i) - static_cast<Types::scalar>(N2 / 2)) * h2};
            }
        }
    }


    namespace Plate {
        EMW::Mesh::SurfaceMesh generatePlatePrimaryMesh(int N, Types::scalar h) {
            std::vector<Mesh::Point> meshgrid;
            meshgrid.reserve(N * N);
            cartesian_product(std::ranges::views::iota(0, N), std::ranges::views::iota(0, N),
                              std::back_inserter(meshgrid), N, h);

            const auto cellsView = std::views::iota(0, (N - 1) * (N - 1)) | std::views::transform(
                    [N](int index) {
                        Types::index i = index + index / (N - 1);
                        const auto point = Mesh::IndexedCell::nodes_t{i, i + 1, i + 1 + N, i + N};
                        return point;
                    }
            );

            const auto cells = Containers::vector<Mesh::IndexedCell::nodes_t>{std::ranges::begin(cellsView),
                                                                              std::ranges::end(cellsView)};

            return Mesh::SurfaceMesh{meshgrid, cells};
        }

        EMW::Mesh::SurfaceMesh
        generateRectangularMesh(int N1, int N2, Types::scalar h1, Types::scalar h2) {
            std::vector<Mesh::Point> meshgrid;
            meshgrid.reserve(N1 * N2);
            cartesian_product_uneven(std::ranges::views::iota(0, N1), std::ranges::views::iota(0, N2),
                                     std::back_inserter(meshgrid), N1, N2, h1, h2);

            const auto cellsView = std::views::iota(0, (N1 - 1) * (N2 - 1)) | std::views::transform(
                    [N1](int index) {
                        Types::index i = index + index / (N1 - 1);  // левый нижний индекс зависит от номера СТРОКИ
                                                                    // в которой располагается ячейка, поэтому
                                                                    // здесь стоит N1
                        const auto point = Mesh::IndexedCell::nodes_t{i, i + 1, i + 1 + N1, i + N1};
                        return point;
                    }
            );

            const auto cells = Containers::vector<Mesh::IndexedCell::nodes_t>{std::ranges::begin(cellsView),
                                                                              std::ranges::end(cellsView)};

            return Mesh::SurfaceMesh{meshgrid, cells};
        }
    }
}

#endif //ELECTROMAGNETIC_WAVES_SCATTERING_PLATEGRID_HPP
