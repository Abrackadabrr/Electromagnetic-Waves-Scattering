//
// Created by evgen on 25.03.2026.
//

#ifndef RESEARCH_MATRUXREPLACEMENT_HPP
#define RESEARCH_MATRUXREPLACEMENT_HPP


#include "math/matrix/iterative_solvers_coverage/IdentityPreconditioner.hpp"
#include "types/Types.hpp"

#include <Eigen/Core>

namespace EMW::Math::LinAgl::Matrix::Wrappers
{
    /**
     * Класс-обёртка для факторизованный комплексной матрицы из двух матриц
     * @tparam MatrixType -- собственно та самая матрица
     */
    template <typename MatrixType1, typename MatrixType2,
              typename PreconditionerType = Preconditioning::IdentityPreconditioner<Types::complex_d, MatrixType1>>
    class VolumeOperatorMatrixReplacement : public Eigen::EigenBase<VolumeOperatorMatrixReplacement<
            MatrixType1, MatrixType2, PreconditionerType>>
    {
    public:
        // Required typedefs, constants, and method:
        typedef Types::complex_d Scalar;
        typedef Eigen::NumTraits<Scalar>::Real RealScalar;
        typedef Types::integer StorageIndex;
        constexpr static Types::integer Options_ = 0;

        enum { ColsAtCompileTime = Eigen::Dynamic, MaxColsAtCompileTime = Eigen::Dynamic, IsRowMajor = false };

        decltype(auto) rows() const { return mat1.rows(); }
        decltype(auto) cols() const { return mat2.cols(); }

        template <typename Rhs>
        Eigen::Product<VolumeOperatorMatrixReplacement, Rhs, Eigen::AliasFreeProduct> operator*(
            const Eigen::MatrixBase<Rhs>& x) const
        {
            return Eigen::Product<VolumeOperatorMatrixReplacement, Rhs, Eigen::AliasFreeProduct>(*this, x.derived());
        }

        // Custom API:
        VolumeOperatorMatrixReplacement() = default;

        VolumeOperatorMatrixReplacement(const MatrixType1& mat1_, const MatrixType2& mat2_): mat1(mat1_),
                                                                           mat2(mat2_),
                                                                           precond(mat1)
        {
        }

        const PreconditionerType& get_preconditioner() const { return precond; }
        const MatrixType1& get_mat1() const { return mat1; }
        const MatrixType2& get_mat2() const { return mat2; }

    private:
        const MatrixType1& mat1;
        const MatrixType2& mat2;
        PreconditionerType precond;
    };
} // namespace EMW::Math::LinAgl::Matrix::Wrappers

// Написание обёртки для умножения матрицы на вектор через специальную структуру внутри Eigen:
// MatrixReplacement * Eigen::DenseVector though a specialization of internal::generic_product_impl:
namespace Eigen::internal
{
    template <typename MatrixType1, typename MatrixType2, typename PreconditionerType>
    using factored_matrix_replacement = EMW::Math::LinAgl::Matrix::Wrappers::VolumeOperatorMatrixReplacement<
        MatrixType1, MatrixType2, PreconditionerType>;

    template <typename MatrixType1, typename MatrixType2, typename PreconditionerType, typename Rhs>
    struct generic_product_impl<factored_matrix_replacement<MatrixType1, MatrixType2, PreconditionerType>, Rhs,
                                SparseShape, DenseShape,
                                GemvProduct> // GEMV stands for matrix-vector
        : generic_product_impl_base<factored_matrix_replacement<MatrixType1, MatrixType2, PreconditionerType>, Rhs,
                                    generic_product_impl<
                                        factored_matrix_replacement<MatrixType1, MatrixType2, PreconditionerType>, Rhs>>
    {
        typedef typename Product<factored_matrix_replacement<MatrixType1, MatrixType2, PreconditionerType>, Rhs>::Scalar
        Scalar;

        template <typename Dest>
        static void scaleAndAddTo(
            Dest& dst, const factored_matrix_replacement<MatrixType1, MatrixType2, PreconditionerType>& lhs,
            const Rhs& rhs,
            const Scalar& alpha)
        {
            // This method should implement "dst += alpha * lhs * rhs" inplace,
            // however, for iterative solvers, alpha is always equal to 1, so let's not bother about it.
            // Для оператора объемного уравнения пишем dst += (I - K * eps) P^{-1} * rhs
#if DEBUG
        assert(alpha == Scalar(1) && "scaling is not implemented");
        EIGEN_ONLY_USED_FOR_DEBUG(alpha);
#endif
            Rhs preconditioned = lhs.get_mat2() * lhs.get_preconditioner().solve(rhs);
            dst.noalias() += lhs.get_mat2() * (preconditioned - (lhs.get_mat1() * preconditioned));
        }
    };
}

#endif //RESEARCH_MATRUXREPLACEMENT_HPP
