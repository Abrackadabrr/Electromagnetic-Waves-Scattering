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
    MatrixReplacement(const MatrixType &matrix)
    : mp_mat(&matrix), precond(matrix), blocks_to_ignore({}), equal_blocks_n(1) {}
    MatrixReplacement(const MatrixType &matrix, const Containers::set<Types::index> &to_ignore, Types::index equal_bl_n)
        : mp_mat(&matrix), precond(matrix), blocks_to_ignore(to_ignore), equal_blocks_n(equal_bl_n) {}
    void attachMyMatrix(const MatrixType &mat) {
        mp_mat = &mat;
        precond = PreconditionerType(mat);
    }
    const MatrixType &my_matrix() const { return *mp_mat; }
    const PreconditionerType &get_preconditioner() const { return precond; }
    const Containers::set<Types::index> &get_blocks_to_ignore() const { return blocks_to_ignore; }
    const Types::index &get_equal_blocks_n() const { return equal_blocks_n; }

  private:
    const MatrixType *mp_mat;
    const PreconditionerType precond;
    const Containers::set<Types::index> blocks_to_ignore;
    const Types::index equal_blocks_n;
};
} // namespace EMW::Math::LinAgl::Matrix::Wrappers

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
#define MEGA_HACK 1
#if MEGA_HACK  // это сейчас основной сценарий использования
        // хак с занулением после умонжения
        // иммитация расчета по хакнутой матрице
        EMW::Types::VectorXc result = lhs.my_matrix() * (lhs.get_preconditioner().solve(rhs));
        const EMW::Types::index size_of_one_j = rhs.size() / lhs.get_equal_blocks_n();
        for (auto &&i : lhs.get_blocks_to_ignore())
            result.block(i * size_of_one_j, 0, size_of_one_j, 1) = EMW::Types::VectorXc::Zero(size_of_one_j);
        dst.noalias() += alpha * result;
#else
        dst.noalias() += alpha * (lhs.my_matrix() * lhs.get_preconditioner().solve(rhs));
#endif
    }
};
}

#endif //MATRIXREPLACEMENT_HPP
