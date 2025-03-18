//
// Created by evgen on 04.03.2025.
//

#ifndef DIAGONALPRECONDITIONER_HPP
#define DIAGONALPRECONDITIONER_HPP

#include "types/Types.hpp"

namespace EMW::Math::LinAgl::Matrix::Preconditioning {
template <typename scalar_t, typename matrix_t> class DiagonalPreconditioner {
    using vector_t = Types::VectorX<scalar_t>;
    vector_t invdiag;
    Types::scalar tolerance = 1e-10;

  public:
    DiagonalPreconditioner(const matrix_t &matrix) {
        invdiag.resize(matrix.cols());
        for (Types::index j = 0; j < matrix.cols(); ++j)
            invdiag(j) = std::abs(matrix(j, j)) > tolerance ? static_cast<Types::scalar>(1) / matrix(j, j) : scalar_t(1);
    }
    // Решается система уравнений Dx = rsh
    // D -- это наша диагонлаьная матрица
    vector_t solve(const vector_t &rhs) const {
        return rhs.cwiseProduct(invdiag);
    }
};

} // namespace EMW::EMW::Math::LinAgl::Matrix::Preconditioning

#endif // DIAGONALPRECONDITIONER_HPP
