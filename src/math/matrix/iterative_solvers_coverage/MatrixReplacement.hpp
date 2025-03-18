//
// Created by evgen on 27.02.2025.
//

#ifndef MATRIXREPLACEMENT_HPP
#define MATRIXREPLACEMENT_HPP

#include "IdentityPreconditioner.hpp"
#include "types/Types.hpp"
#include <Eigen/Core>

namespace EMW::Math::LinAgl::Matrix::Wrappers {
/**
 * Класс-обёртка комплекной матрицы для работы с солверами внутри Eigen
 * @tparam MatrixType -- собственно та самая матрица
 */
template <typename MatrixType,
          typename PreconditionerType = Preconditioning::IdentityPreconditioner<Types::complex_d, MatrixType>>
class MatrixReplacement : public Eigen::EigenBase<MatrixReplacement<MatrixType, PreconditionerType>> {
  public:
    // Required typedefs, constants, and method:
    typedef Types::complex_d Scalar;
    typedef Eigen::NumTraits<Scalar>::Real RealScalar;
    typedef Types::integer StorageIndex;
    constexpr static Types::integer Options_ = 0;
    enum { ColsAtCompileTime = Eigen::Dynamic, MaxColsAtCompileTime = Eigen::Dynamic, IsRowMajor = false };

    decltype(auto) rows() const { return mp_mat->rows(); }
    decltype(auto) cols() const { return mp_mat->cols(); }

    template <typename Rhs>
    Eigen::Product<MatrixReplacement, Rhs, Eigen::AliasFreeProduct> operator*(const Eigen::MatrixBase<Rhs> &x) const {
        return Eigen::Product<MatrixReplacement, Rhs, Eigen::AliasFreeProduct>(*this, x.derived());
    }

    // Custom API:
    MatrixReplacement() = default;
    MatrixReplacement(const MatrixType &matrix) : mp_mat(&matrix), precond(matrix) {}
    void attachMyMatrix(const MatrixType &mat) {
        mp_mat = &mat;
        precond = PreconditionerType(mat);
    }
    const MatrixType &my_matrix() const { return *mp_mat; }
    const PreconditionerType &get_preconditioner() const { return precond; }

  private:
    const MatrixType *mp_mat;
    const PreconditionerType precond;
};
} // namespace EMW::EMW::Math::LinAgl::Matrix::Wrappers

// Написание обёртки для умножения матрицы на вектор через специальную структуру внутри Eigen:
// MatrixReplacement * Eigen::DenseVector though a specialization of internal::generic_product_impl:
namespace Eigen::internal {
template <typename MatrixType, typename PreconditionerType>
using matrix_replacement = EMW::Math::LinAgl::Matrix::Wrappers::MatrixReplacement<MatrixType, PreconditionerType>;

template <typename MatrixType, typename PreconditionerType, typename Rhs>
struct generic_product_impl<matrix_replacement<MatrixType, PreconditionerType>, Rhs, SparseShape, DenseShape,
                            GemvProduct> // GEMV stands for matrix-vector
    : generic_product_impl_base<matrix_replacement<MatrixType, PreconditionerType>, Rhs,
                                generic_product_impl<matrix_replacement<MatrixType, PreconditionerType>, Rhs>> {
    typedef typename Product<matrix_replacement<MatrixType, PreconditionerType>, Rhs>::Scalar Scalar;

    template <typename Dest>
    static void scaleAndAddTo(Dest &dst, const matrix_replacement<MatrixType, PreconditionerType> &lhs, const Rhs &rhs,
                              const Scalar &alpha) {
        // This method should implement "dst += alpha * lhs * rhs" inplace,
        // however, for iterative solvers, alpha is always equal to 1, so let's not bother about it.
#if DEBUG
        assert(alpha == Scalar(1) && "scaling is not implemented");
        EIGEN_ONLY_USED_FOR_DEBUG(alpha);
#endif
        dst.noalias() += alpha * (lhs.my_matrix() * (lhs.get_preconditioner().solve(rhs)));
    }
};
}

#endif //MATRIXREPLACEMENT_HPP
