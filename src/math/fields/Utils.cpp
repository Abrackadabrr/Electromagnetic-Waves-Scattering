//
// Created by evgen on 12.10.2024.
//

#include "Utils.hpp"
#include <ranges>

namespace EMW::Math::FieldUtils {

Math::SurfaceScalarField<EMW::Types::complex_d> relativeError(const SurfaceVectorField &v1, const SurfaceVectorField &v2) {
    using std::begin;
    using std::end;
    const auto data_view = std::views::zip_transform(
        [](const auto &f1, const auto &f2) -> Types::complex_d {
            if (f1.norm() < 1e-14)
                return {0, 0};
            return {std::abs(f1.real().norm() - f2.real().norm()) / f1.real().norm(),
                    std::abs(f1.imag().norm() - f2.imag().norm()) / f1.imag().norm()};
        },
        v1.getField(), v2.getField());
    return {v1.getManifold(), {begin(data_view), end(data_view)}};
}

} // namespace EMW::Math::FieldUtils