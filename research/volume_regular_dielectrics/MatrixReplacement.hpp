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
    template <typename MatrixType1,
              typename PreconditionerType = Preconditioning::IdentityPreconditioner<Types::complex_d, MatrixType1>>
    class VolumeOperatorMatrixReplacement : public Eigen::EigenBase<VolumeOperatorMatrixReplacement<
            MatrixType1, PreconditionerType>>
    {
    public:
        // Required typedefs, constants, and method:
        typedef Types::complex_d Scalar;
        typedef Eigen::NumTraits<Scalar>::Real RealScalar;
        typedef Types::integer StorageIndex;
        constexpr static Types::integer Options_ = 0;
        constexpr static auto zeroLike = std::numeric_limits<RealScalar>::epsilon();

        enum { ColsAtCompileTime = Eigen::Dynamic, MaxColsAtCompileTime = Eigen::Dynamic, IsRowMajor = false };

        decltype(auto) rows() const { return mat_.rows(); }
        decltype(auto) cols() const { return mat_.cols(); }

        template <typename Rhs>
        Eigen::Product<VolumeOperatorMatrixReplacement, Rhs, Eigen::AliasFreeProduct> operator*(
            const Eigen::MatrixBase<Rhs>& x) const
        {
            return Eigen::Product<VolumeOperatorMatrixReplacement, Rhs, Eigen::AliasFreeProduct>(*this, x.derived());
        }

        // Custom API:
        VolumeOperatorMatrixReplacement() = default;

        VolumeOperatorMatrixReplacement(const MatrixType1& mat, const Types::VectorX<Scalar>& eps_vec): mat_(mat),
            epsilon_vec_(eps_vec),
            precond_(mat_), mask_(Types::VectorX<RealScalar>::Zero(mat_.rows())),
            inverse_epsilon_vec_(Types::VectorX<RealScalar>::Zero(mat_.rows()))
        {
            if (epsilon_vec_.size() != mat.cols())
            {
                throw std::invalid_argument("Epsilon is not set up correctly");
            }
            for (size_t i = 0; i < mat.rows(); ++i)
            {
                mask_[i] = (std::abs(eps_vec[i]) > zeroLike);
                inverse_epsilon_vec_[i] = std::abs(eps_vec[i]) < zeroLike ? 1 : 1. / epsilon_vec_[i];
            }
        }

        const PreconditionerType& get_preconditioner() const { return precond_; }
        const MatrixType1& get_mat() const { return mat_; }
        [[nodiscard]] const Types::VectorX<Scalar>& get_epsilon_vec() const { return epsilon_vec_; }
        [[nodiscard]] const Types::VectorX<Scalar>& get_inverse_epsilon_vec() const { return inverse_epsilon_vec_; }
        [[nodiscard]] const Types::VectorX<RealScalar>& get_mask() const { return mask_; }

        // Изменение правой части согласно насчитанной маске
        //
        // Это нужно для того, чтобы занулить фиктивные компоненты в невязке, которые необходимы
        // только для устраивания дважды тёплицевой структуры.
        //
        template <typename Rhs>
        [[nodiscard]] Types::VectorX<Scalar> modify_rhs_according_to_mask(Rhs&& rhs) const
        {
            return std::forward<Rhs>(rhs).cwiseProduct(mask_);
        }

        template <typename Solution>
        [[nodiscard]] Types::VectorXc modify_vec_according_to_inverse_epsilon_sqrt(Solution&& sol) const
        {
            return std::forward<Solution>(sol).cwiseProduct(inverse_epsilon_vec_.cwiseSqrt());
        }

        template <typename Solution>
        [[nodiscard]] Types::VectorXc modify_vec_according_to_inverse_epsilon(Solution&& sol) const
        {
            return std::forward<Solution>(sol).cwiseProduct(inverse_epsilon_vec_);
        }

    private:
        const MatrixType1& mat_;
        // Вектора, которые отвечают за диагональную матрицу epsilon - 1
        const Types::VectorX<Scalar> epsilon_vec_;
        Types::VectorX<Scalar> inverse_epsilon_vec_;
        // Нужна для отсечения лишних элементов в векторе неизвестных (так как кое-где не нужно решать систему)
        Types::VectorX<RealScalar> mask_;
        PreconditionerType precond_;
    };
} // namespace EMW::Math::LinAgl::Matrix::Wrappers

// Написание обёртки для умножения матрицы на вектор через специальную структуру внутри Eigen:
// MatrixReplacement * Eigen::DenseVector though a specialization of internal::generic_product_impl:
namespace Eigen::internal
{
    template <typename MatrixType1, typename PreconditionerType>
    using factored_matrix_replacement = EMW::Math::LinAgl::Matrix::Wrappers::VolumeOperatorMatrixReplacement<
        MatrixType1, PreconditionerType>;

    template <typename MatrixType1, typename PreconditionerType, typename Rhs>
    struct generic_product_impl<factored_matrix_replacement<MatrixType1, PreconditionerType>, Rhs,
                                SparseShape, DenseShape,
                                GemvProduct> // GEMV stands for matrix-vector
        : generic_product_impl_base<factored_matrix_replacement<MatrixType1, PreconditionerType>, Rhs,
                                    generic_product_impl<
                                        factored_matrix_replacement<MatrixType1, PreconditionerType>, Rhs>>
    {
        typedef typename Product<factored_matrix_replacement<MatrixType1, PreconditionerType>, Rhs>::Scalar
        Scalar;

        template <typename Dest>
        static void scaleAndAddTo(
            Dest& dst, const factored_matrix_replacement<MatrixType1, PreconditionerType>& lhs,
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
            Rhs preconditioned = lhs.get_preconditioner().solve(rhs);
            const auto sqrt_eps_vec = lhs.get_epsilon_vec().cwiseSqrt();
            const auto inverse_sqrt_eps_vec = lhs.get_inverse_epsilon_vec().cwiseSqrt();
            dst.noalias() += (rhs.cwiseProduct(inverse_sqrt_eps_vec) - (lhs.get_mat() *
                Eigen::VectorXcd{rhs.cwiseProduct(sqrt_eps_vec)})).cwiseProduct(lhs.get_mask());
        }
    };
}

#endif //RESEARCH_MATRUXREPLACEMENT_HPP
