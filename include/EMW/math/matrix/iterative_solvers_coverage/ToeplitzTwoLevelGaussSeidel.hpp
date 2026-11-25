
#ifndef EMW_TOEPLITZ_GAUSS_SEIDEL
#define EMW_TOEPLITZ_GAUSS_SEIDEL

#include "EMW/math/matrix/ToeplitzFullyTemplated.hpp"
#include "third_party/eigen/Eigen/LU"

#include <concepts>
#include <stdexcept>
#include <type_traits>
#include <utility>

namespace EMW::Math::LinAgl::Matrix {

// Either an ordinary Eigen-compatible product, or the library's accumulating
// matvec API (dest += block * rhs).
template <typename Block, typename Scalar>
concept BlockVectorProduct =
    requires(const Block &block, const Types::VectorX<Scalar> &rhs) { Types::VectorX<Scalar>(block * rhs); };

template <typename Block, typename Scalar>
concept BlockMatvec = requires(const Block &block, const Types::VectorX<Scalar> &rhs, Types::VectorX<Scalar> &dest) {
    { block.matvec(rhs, dest) } -> std::same_as<void>;
};

// The conversion must produce a dense matrix of the same scalar type.
// In particular, to_dense() is required for all *custom* block types.
template <typename Block, typename Scalar>
concept BlockToDense = requires(const Block &block) {
    { block.to_dense() } -> std::convertible_to<Types::MatrixX<Scalar>>;
};

// Compatibility with the previous dense-block API: Eigen::MatrixX<Scalar>
// is already dense and has no to_dense() member.
template <typename Block, typename Scalar>
concept NativeDenseBlock = std::same_as<std::remove_cvref_t<Block>, Types::MatrixX<Scalar>>;

template <typename Block, typename Scalar>
concept BlockDimensionsAndProduct = requires(const Block &block) {
    { block.rows() } -> std::convertible_to<Eigen::Index>;
    { block.cols() } -> std::convertible_to<Eigen::Index>;
} && (BlockVectorProduct<Block, Scalar> || BlockMatvec<Block, Scalar>);

// Strict contract for a *custom* innermost block.
template <typename Block, typename Scalar>
concept InnermostBlock = BlockDimensionsAndProduct<Block, Scalar> && BlockToDense<Block, Scalar>;

// Eigen::MatrixX<Scalar> is a backwards-compatible exception to the contract:
// it is dense already and therefore does not require a to_dense() member.
template <typename Block, typename Scalar>
concept SupportedInnermostBlock =
    InnermostBlock<Block, Scalar> || (BlockDimensionsAndProduct<Block, Scalar> && NativeDenseBlock<Block, Scalar>);

// For two levels: ToeplitzStructure<Scalar, ToeplitzStructure<Scalar, Block>>.
// The lower triangular part is taken with respect to the innermost Blocks.
template <typename Scalar, typename InnerToeplitz> class TwoLevelBlockLowerPreconditioner {
  public:
    using matrix_t = ToeplitzStructure<Scalar, InnerToeplitz>;
    using block_t = typename InnerToeplitz::block_type;
    using vector_t = Types::VectorX<Scalar>;
    using dense_matrix_t = Types::MatrixX<Scalar>;
    using lu_t = Eigen::FullPivLU<dense_matrix_t>;

    static_assert(matrix_t::toeplitz_depth == 2,
                  "TwoLevelBlockLowerPreconditioner requires exactly two Toeplitz levels");
    static_assert(std::same_as<InnerToeplitz, ToeplitzStructure<Scalar, block_t>>,
                  "Inner level must be ToeplitzStructure<Scalar, block_t>");
    static_assert(SupportedInnermostBlock<block_t, Scalar>,
                  "Innermost block must expose rows()/cols(), support either operator*(vector) "
                  "or void matvec(const vector_t&, vector_t&), and have "
                  "to_dense() returning a dense matrix (unless it is already Eigen::MatrixX)");

    struct Workspace {
        explicit Workspace(Eigen::Index b) : diagonal_rhs(b), diagonal_result(b), block_rhs(b), block_result(b) {}

        vector_t diagonal_rhs;
        vector_t diagonal_result;
        vector_t block_rhs;
        vector_t block_result;
    };

    explicit TwoLevelBlockLowerPreconditioner(const matrix_t &A) : A_(A) {
        outer_size_ = A_.rows_in_toeplitrz();
        if (outer_size_ == 0 || outer_size_ != A_.cols_in_toeplitrz()) {
            throw std::invalid_argument("Outer Toeplitz level must be nonempty and square");
        }

        const auto &inner = A_.get_block(0, 0);
        inner_size_ = inner.rows_in_toeplitrz();
        if (inner_size_ == 0 || inner_size_ != inner.cols_in_toeplitrz()) {
            throw std::invalid_argument("Inner Toeplitz level must be nonempty and square");
        }

        const block_t &diagonal_block = inner.get_block(0, 0);
        block_size_ = static_cast<Eigen::Index>(diagonal_block.rows());
        if (block_size_ <= 0 || static_cast<Eigen::Index>(diagonal_block.cols()) != block_size_) {
            throw std::invalid_argument("Innermost diagonal block must be nonempty and square");
        }

        const Eigen::Index full_size = static_cast<Eigen::Index>(outer_size_ * inner_size_) * block_size_;
        if (A_.rows() != full_size || A_.cols() != full_size) {
            throw std::invalid_argument("Inconsistent Toeplitz matrix dimensions");
        }

        // Only the common diagonal block is densified; exactly once in setup.
        // A native Eigen::MatrixX block needs no conversion.
        if constexpr (BlockToDense<block_t, Scalar>) {
            dense_matrix_t D = diagonal_block.to_dense();
            if (D.rows() != block_size_ || D.cols() != block_size_) {
                throw std::invalid_argument("to_dense() returned a matrix with incorrect dimensions");
            }
            diagonal_lu_.compute(D);
        } else {
            diagonal_lu_.compute(diagonal_block);
        }

        if (!diagonal_lu_.isInvertible()) {
            throw std::runtime_error("Innermost diagonal block is singular");
        }
    }

    [[nodiscard]] Workspace make_workspace() const { return Workspace(block_size_); }

    // x := block_tril(A)^{-1} rhs. Supports x and rhs referring to the same vector.
    // Work memory is entirely owned by Workspace and can be reused across solves.
    void solve(const vector_t &rhs, vector_t &x, Workspace &workspace) const {
        if (rhs.size() != A_.rows()) {
            throw std::invalid_argument("RHS size does not match matrix size");
        }
        if (workspace.diagonal_rhs.size() != block_size_ || workspace.diagonal_result.size() != block_size_ ||
            workspace.block_rhs.size() != block_size_ || workspace.block_result.size() != block_size_) {
            throw std::invalid_argument("Workspace size does not match innermost block size");
        }

        x = rhs;
        for (Types::index outer_row = 0; outer_row < outer_size_; ++outer_row) {
            for (Types::index inner_row = 0; inner_row < inner_size_; ++inner_row) {
                const Eigen::Index row = static_cast<Eigen::Index>(outer_row * inner_size_ + inner_row);
                auto x_row = x.segment(row * block_size_, block_size_);

                for (Types::index outer_col = 0; outer_col <= outer_row; ++outer_col) {
                    const auto &inner_A = A_.get_block(outer_row, outer_col);
                    const Types::index inner_col_stop = (outer_col == outer_row) ? inner_row : inner_size_;

                    for (Types::index inner_col = 0; inner_col < inner_col_stop; ++inner_col) {
                        const Eigen::Index col = static_cast<Eigen::Index>(outer_col * inner_size_ + inner_col);
                        const auto &Aij = inner_A.get_block(inner_row, inner_col);

                        if constexpr (NativeDenseBlock<block_t, Scalar>) {
                            // Eigen supports multiplying directly by a vector segment.
                            x_row.noalias() -= Aij * x.segment(col * block_size_, block_size_);
                        } else if constexpr (BlockVectorProduct<block_t, Scalar>) {
                            // For arbitrary user-defined operator*, its argument may
                            // need to be an actual vector_t, not an Eigen segment.
                            workspace.block_rhs = x.segment(col * block_size_, block_size_);
                            workspace.block_result = Aij * workspace.block_rhs;
                            x_row -= workspace.block_result;
                        } else {
                            // Supports the accumulating matvec API: dest += Aij * rhs.
                            workspace.block_rhs = x.segment(col * block_size_, block_size_);
                            workspace.block_result.setZero();
                            Aij.matvec(workspace.block_rhs, workspace.block_result);
                            x_row -= workspace.block_result;
                        }
                    }
                }

                workspace.diagonal_rhs = x_row;
                workspace.diagonal_result = diagonal_lu_.solve(workspace.diagonal_rhs);
                x_row = workspace.diagonal_result;
            }
        }
    }

    void solve(const vector_t &rhs, vector_t &x) const {
        auto workspace = make_workspace();
        solve(rhs, x, workspace);
    }

    [[nodiscard]] Eigen::Index rows() const noexcept { return A_.rows(); }
    [[nodiscard]] Eigen::Index cols() const noexcept { return A_.cols(); }
    [[nodiscard]] const lu_t &diagonal_lu() const noexcept { return diagonal_lu_; }

  private:
    const matrix_t &A_; // Must outlive the preconditioner and not be modified.
    Types::index outer_size_ = 0;
    Types::index inner_size_ = 0;
    Eigen::Index block_size_ = 0;
    lu_t diagonal_lu_; // Factored only in the constructor.
};
} // namespace EMW::Math::LinAgl::Matrix

namespace EMW::Math::LinAgl::Matrix::Preconditioning {
template <typename MatrixType> class ToeplitzGaussSeidelEigen {
  public:
    using scalar_t = typename MatrixType::scalar_type;
    using vector_t = Types::VectorX<scalar_t>;

    using preconditioner_t = TwoLevelBlockLowerPreconditioner<scalar_t, typename MatrixType::block_type>;

    explicit ToeplitzGaussSeidelEigen(const MatrixType &matrix)
        : preconditioner_(matrix), workspace_(preconditioner_.make_workspace()) {}

    template <typename Derived> [[nodiscard]] vector_t solve(const Eigen::MatrixBase<Derived> &rhs) const {
        vector_t rhs_vector = rhs;
        vector_t result(rhs.size());

        preconditioner_.solve(rhs_vector, result, workspace_);

        return result;
    }

  private:
    preconditioner_t preconditioner_;
    mutable typename preconditioner_t::Workspace workspace_;
};

} // namespace EMW::Math::LinAgl::Matrix::Preconditioning

#endif
