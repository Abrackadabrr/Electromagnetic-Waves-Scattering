//
// Created by evgen on 27.02.2025.
//

#ifndef MATRIXREPLACEMENT_HPP
#define MATRIXREPLACEMENT_HPP

#include "types/Types.hpp"
#include <Eigen/Core>

namespace EMW::Wrappers {
/**
 * Класс-обёртка комплекной матрицы для работы с солверами внутри Eigen
 * @tparam MatrixType -- собственно та самая матрица
 */
template <typename MatrixType> class MatrixReplacement : public Eigen::EigenBase<MatrixReplacement<MatrixType>> {
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
    explicit MatrixReplacement(const MatrixType & matrix) : mp_mat(&matrix) {}
    void attachMyMatrix(const MatrixType &mat) { mp_mat = &mat; }
    const MatrixType &my_matrix() const { return *mp_mat; }

  private:
    const MatrixType *mp_mat;
};
} // namespace EMW::Wrappers

// Написание обёртки для умножения матрицы на вектор через специальную структуру внутри Eigen:
// MatrixReplacement * Eigen::DenseVector though a specialization of internal::generic_product_impl:
namespace Eigen::internal {
template <typename MatrixType> using matrix_replacement = EMW::Wrappers::MatrixReplacement<MatrixType>;

template <typename MatrixType, typename Rhs>
struct generic_product_impl<matrix_replacement<MatrixType>, Rhs, SparseShape, DenseShape,
                            GemvProduct> // GEMV stands for matrix-vector
    : generic_product_impl_base<matrix_replacement<MatrixType>, Rhs,
                                generic_product_impl<matrix_replacement<MatrixType>, Rhs>> {
    typedef typename Product<matrix_replacement<MatrixType>, Rhs>::Scalar Scalar;

    template <typename Dest>
    static void scaleAndAddTo(Dest &dst, const matrix_replacement<MatrixType> &lhs, const Rhs &rhs,
                              const Scalar &alpha) {
        // This method should implement "dst += alpha * lhs * rhs" inplace,
        // however, for iterative solvers, alpha is always equal to 1, so let's not bother about it.
#if DEBUG
        assert(alpha == Scalar(1) && "scaling is not implemented");
        EIGEN_ONLY_USED_FOR_DEBUG(alpha);
#endif
        // Here we could simply call dst.noalias() += lhs.my_matrix() * rhs,
        // but let's do something fancier (and less efficient):
        dst.noalias() += alpha * (lhs.my_matrix() * rhs);
    }
};
}

#endif //MATRIXREPLACEMENT_HPP
