//
// Created by evgen on 06.07.24.
//

#ifndef ELECTROMAGNETIC_WAVES_SCATTERING_UTILS_HPP
#define ELECTROMAGNETIC_WAVES_SCATTERING_UTILS_HPP


#include <utility>
#include <tuple>
#include <ranges>
#include "types/Types.hpp"
#include "mesh/SurfaceMesh.hpp"

namespace EMW::Utils {
    template<typename Range1, typename Range2, typename Range3, typename OutputIterator>
    void cartesian_product(Range1 const &r1, Range2 const &r2, Range3 const &r3,
                           Types::index N1, Types::index N2, Types::index N3,
                           Types::scalar h1, Types::scalar h2, Types::scalar h3,
                           OutputIterator out) {
        using std::begin;
        using std::end;

        for (auto i = begin(r3); i != end(r3); ++i) {
            for (auto j = begin(r2); j != end(r2); ++j) {
                for (auto k = begin(r1); k != end(r1); ++k) {
                    *out++ = Types::Vector3d{h1 * (static_cast<Types::scalar>(*k) - static_cast<Types::scalar>(N1 / 2)),
                                             h2 * (static_cast<Types::scalar>(*j) - static_cast<Types::scalar>(N2 / 2)),
                                             h1 * (static_cast<Types::scalar>(*i)
                                                   - static_cast<Types::scalar>(N3 / 2))};
                }
            }
        }
    }

    template<typename Range1, typename Range2, typename OutputIterator>
    void cartesian_productYZ(Range1 const &r1, Range2 const &r2, OutputIterator out, Types::index N, Types::scalar h) {
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
    cartesian_product_unevenXY(Range1 const &r1, Range2 const &r2, OutputIterator out, Types::index N1, Types::index N2,
                               Types::scalar h1, Types::scalar h2) {
        using std::begin;
        using std::end;

        for (auto i = begin(r2); i != end(r2); ++i) {
            for (auto j = begin(r1); j != end(r1); ++j) {
                *out++ = Types::Vector3d{ (static_cast<Types::scalar>(*j) - static_cast<Types::scalar>(N1 / 2) +
                                             static_cast<Types::scalar>((N1 - 1) % 2) / 2.) * h1,
                                         (static_cast<Types::scalar>(*i) - static_cast<Types::scalar>(N2 / 2) +
                                          static_cast<Types::scalar>((N2 - 1) % 2) / 2.) * h2, 0};
            }
        }
    }
}

#endif //ELECTROMAGNETIC_WAVES_SCATTERING_UTILS_HPP
