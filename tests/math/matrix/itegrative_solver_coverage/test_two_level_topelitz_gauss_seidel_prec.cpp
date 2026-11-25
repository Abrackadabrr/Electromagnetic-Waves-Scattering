#include "EMW/math/matrix/iterative_solvers_coverage/ToeplitzTwoLevelGaussSeidel.hpp"

#include <gtest/gtest.h>

#include <complex>

using namespace EMW::Math::LinAgl::Matrix;
using namespace EMW;

TEST(TwoLevelBlockLowerPreconditioner, ForwardSubstitutionAndInPlace) {
    using scalar_t = std::complex<double>;
    using dense_t = Types::MatrixX<scalar_t>;
    using vector_t = Types::VectorX<scalar_t>;
    using inner_t = ToeplitzStructure<scalar_t, dense_t>;
    using outer_t = ToeplitzStructure<scalar_t, inner_t>;
    using preconditioner_t = TwoLevelBlockLowerPreconditioner<scalar_t, inner_t>;

    static_assert(outer_t::toeplitz_depth == 2);

    constexpr Types::index n_outer = 10;
    constexpr Types::index n_inner = 10;
    constexpr Eigen::Index b = 32;

    Containers::vector<inner_t> outer_blocks;
    outer_blocks.reserve(outer_t::get_size_of_container(n_outer, n_outer));

    // Unique block order: first column bottom-to-top, then first row left-to-right.
    for (Types::index oi = 0; oi < outer_t::get_size_of_container(n_outer, n_outer); ++oi) {
        Containers::vector<dense_t> inner_blocks;
        inner_blocks.reserve(inner_t::get_size_of_container(n_inner, n_inner));
        for (Types::index ii = 0; ii < inner_t::get_size_of_container(n_inner, n_inner); ++ii) {
            dense_t block = dense_t::Random(b, b);
            if (oi == n_outer - 1 && ii == n_inner - 1) {
                block += scalar_t(b) * dense_t::Identity(b, b);
                Eigen::JacobiSVD<Eigen::MatrixXcd> svd(block);

                const auto &singular_values = svd.singularValues();

                double cond = singular_values(0) / singular_values(singular_values.size() - 1);
                std::cout << cond << std::endl;
                ASSERT_LT(cond, 2);
            }
            inner_blocks.push_back(std::move(block));
        }
        outer_blocks.emplace_back(n_inner, n_inner, std::move(inner_blocks));
    }

    const outer_t A(n_outer, n_outer, std::move(outer_blocks));
    auto L = A.to_dense();

    // Drop the upper triangle in innermost block units.
    const Eigen::Index n_blocks = static_cast<Eigen::Index>(n_outer * n_inner);
    for (Eigen::Index i = 0; i < n_blocks; ++i) {
        for (Eigen::Index j = i + 1; j < n_blocks; ++j) {
            L.block(i * b, j * b, b, b).setZero();
        }
    }

    const vector_t x_exact = vector_t::Random(n_blocks * b);
    const vector_t rhs = L * x_exact;
    const preconditioner_t prec(A);
    auto workspace = prec.make_workspace();

    Types::scalar tolerance = std::numeric_limits<Types::scalar>::epsilon() * 1e1;

    vector_t x = vector_t::Zero(n_blocks * b);
    prec.solve(rhs, x, workspace);
    EXPECT_LT((x - x_exact).norm() / x_exact.norm(), tolerance);
    EXPECT_LT((L * x - rhs).norm() / rhs.norm(), tolerance);

    vector_t inplace = rhs;
    prec.solve(inplace, inplace, workspace);
    EXPECT_LT((inplace - x_exact).norm() / x_exact.norm(), tolerance);

    // A subsequent call reuses the same precomputed inverse and workspace.
    prec.solve(rhs, x, workspace);
    EXPECT_LT((x - x_exact).norm() / x_exact.norm(), 1e-12);
};

// Простейший пользовательский класс матрицы, подходящий в теплицев гасс-зейдель прек
class DenseMatrixWrapper {
  public:
    DenseMatrixWrapper() = default;

    explicit DenseMatrixWrapper(Eigen::MatrixXd matrix) : matrix_(std::move(matrix)) {}

    [[nodiscard]] Eigen::Index rows() const noexcept { return matrix_.rows(); }

    [[nodiscard]] Eigen::Index cols() const noexcept { return matrix_.cols(); }

    [[nodiscard]] Eigen::MatrixXd to_dense() const { return matrix_; }

    void matvec(const Eigen::VectorXd &rhs, Eigen::VectorXd &dest) const { dest.noalias() += matrix_ * rhs; }

  private:
    Eigen::MatrixXd matrix_;
};

TEST(TwoLevelBlockLowerPreconditioner, ForwardSubstitutionAndInPlaceForASpecialClass) {
    using scalar_t = double;
    using dense_t = Eigen::MatrixXd;
    using vector_t = Eigen::VectorXd;

    using inner_t = ToeplitzStructure<scalar_t, DenseMatrixWrapper>;
    using outer_t = ToeplitzStructure<scalar_t, inner_t>;
    using preconditioner_t = TwoLevelBlockLowerPreconditioner<scalar_t, inner_t>;

    static_assert(outer_t::toeplitz_depth == 2);
    static_assert(InnermostBlock<DenseMatrixWrapper, scalar_t>);

    constexpr Types::index n_outer = 10;
    constexpr Types::index n_inner = 10;
    constexpr Eigen::Index b = 128;

    Containers::vector<inner_t> outer_blocks;
    outer_blocks.reserve(outer_t::get_size_of_container(n_outer, n_outer));

    // Unique block order: first column bottom-to-top,
    // then first row left-to-right.
    for (Types::index oi = 0; oi < outer_t::get_size_of_container(n_outer, n_outer); ++oi) {

        Containers::vector<DenseMatrixWrapper> inner_blocks;
        inner_blocks.reserve(inner_t::get_size_of_container(n_inner, n_inner));

        for (Types::index ii = 0; ii < inner_t::get_size_of_container(n_inner, n_inner); ++ii) {

            dense_t block = dense_t::Random(b, b);

            if (oi == n_outer - 1 && ii == n_inner - 1) {
                block += scalar_t(b) * dense_t::Identity(b, b);

                Eigen::JacobiSVD<dense_t> svd(block);

                const auto &singular_values = svd.singularValues();

                const double cond = singular_values(0) / singular_values(singular_values.size() - 1);

                ASSERT_LT(cond, 2.0);
            }

            inner_blocks.emplace_back(std::move(block));
        }

        outer_blocks.emplace_back(n_inner, n_inner, std::move(inner_blocks));
    }

    const outer_t A(n_outer, n_outer, std::move(outer_blocks));

    // Convert the entire matrix to dense format.
    dense_t L = A.to_dense();

    // Drop the upper triangle in innermost block units.
    const Eigen::Index n_blocks = static_cast<Eigen::Index>(n_outer * n_inner);

    for (Eigen::Index i = 0; i < n_blocks; ++i) {
        for (Eigen::Index j = i + 1; j < n_blocks; ++j) {
            L.block(i * b, j * b, b, b).setZero();
        }
    }

    const vector_t x_exact = vector_t::Random(n_blocks * b);
    const vector_t rhs = L * x_exact;

    const preconditioner_t prec(A);
    auto workspace = prec.make_workspace();

    constexpr double tolerance = 1e-12;

    // Out-of-place solve.
    vector_t x = vector_t::Zero(n_blocks * b);

    prec.solve(rhs, x, workspace);

    EXPECT_LT((x - x_exact).norm() / x_exact.norm(), tolerance);

    EXPECT_LT((L * x - rhs).norm() / rhs.norm(), tolerance);

    // In-place solve.
    vector_t inplace = rhs;

    prec.solve(inplace, inplace, workspace);

    EXPECT_LT((inplace - x_exact).norm() / x_exact.norm(), tolerance);

    // Reuse the precomputed LU factorization and workspace.
    prec.solve(rhs, x, workspace);

    EXPECT_LT((x - x_exact).norm() / x_exact.norm(), tolerance);
}
