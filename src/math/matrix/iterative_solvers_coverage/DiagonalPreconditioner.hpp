//
// Created by evgen on 04.03.2025.
//

#ifndef DIAGONALPRECONDITIONER_HPP
#define DIAGONALPRECONDITIONER_HPP

#include "math/matrix/MatrixTraits.hpp"

#include "types/Types.hpp"

#include <iostream>

namespace EMW::Math::LinAgl::Matrix::Preconditioning {
template <typename scalar_t, typename matrix_t> class DiagonalPreconditioner {
    using matrix_traits = MatrixTraits<matrix_t>;
    using vector_t = Types::VectorX<scalar_t>;

    // Храним только обратные элементы к диагональным
    vector_t invdiag;
    Types::scalar tolerance = 1e-10;

  public:
    /** ну тут проблема с делением на ноль но я забил */
    DiagonalPreconditioner(const matrix_t &matrix) : invdiag(matrix.diagonal().cwiseInverse()) {}

    // Решается система уравнений Dx = rsh
    // D -- это наша диагонлаьная матрица
    vector_t solve(const vector_t &rhs) const { return rhs.cwiseProduct(invdiag); }
};

template <typename scalar_t, typename matrix_t> class BlockDiagonalPreconditioner {
    using matrix_traits = MatrixTraits<matrix_t>;
    using vector_t = Types::VectorX<scalar_t>;

    // Храним только обратные элементы к диагональным
    Types::MatrixX<scalar_t> invdiag_block;
    Types::index block_rows;
    Types::index n_blocks;
    Types::scalar tolerance = 1e-10;

  public:
    /** ну тут проблема с делением на ноль но я забил */
    BlockDiagonalPreconditioner(const matrix_t &matrix)
        : invdiag_block(matrix.get_block(0, 0).get_block(0, 0).template get<0>()
          .inverse()), block_rows(invdiag_block.rows()),
          n_blocks(matrix.rows() / block_rows) {
        std::cout << "BlockDiagonalPreconditioner is constructed" << std::endl;
        std::cout << "with " << n_blocks << "equal diagonal blocks" << std::endl;
    }

    // Решается система уравнений Dx = rsh
    // D -- это наша диагонлаьная матрица
    vector_t solve(const vector_t &rhs) const {
        vector_t result = vector_t::Zero(block_rows * n_blocks);
        for (Types::index i = 0; i < n_blocks; i++) {
          result.block(i * block_rows, 0, block_rows, 1) = invdiag_block * rhs.block(i * block_rows, 0, block_rows, 1);
      }
      return result;
  }
};


} // namespace EMW::EMW::Math::LinAgl::Matrix::Preconditioning

#endif // DIAGONALPRECONDITIONER_HPP
