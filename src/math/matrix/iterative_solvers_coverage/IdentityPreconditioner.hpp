//
// Created by evgen on 04.03.2025.
//

#ifndef IDENTITYPRECONDITIONER_HPP
#define IDENTITYPRECONDITIONER_HPP

#include "types/Types.hpp"

namespace EMW::Math::LinAgl::Matrix::Preconditioning {
template <typename scalar_t, typename matrix_t> struct IdentityPreconditioner {
    using vector_t = Types::VectorX<scalar_t>;

    IdentityPreconditioner(const matrix_t &matrix) {}
    // Решается система уравнений Dx = rsh
    // D -- это наша диагонлаьная матрица
    inline vector_t solve(const vector_t &rhs) const noexcept { return rhs;}
};
}

#endif //IDENTITYPRECONDITIONER_HPP
