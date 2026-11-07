#ifndef EMW_MATH_MATRIX_TOEPLITZ_MULTILEVEL_NEW_HPP
#define EMW_MATH_MATRIX_TOEPLITZ_MULTILEVEL_NEW_HPP

#include "TeoplitzTemplateNew.hpp"

#include <functional>
#include <span>
#include <utility>
#include <vector>

namespace EMW::Math::LinAgl::Matrix {

template <ToeplitzMatrixBlock Block> class MultilevelToeplitz;

/** RAII storage for a runtime-depth multilevel Toeplitz matrix.
 *
 * Every level is described by its block-row and block-column count. Dense leaf
 * blocks are stored in mixed-radix diagonal order, with the innermost level
 * varying fastest. For Eigen leaves, toeplitz_new_detail::Storage consolidates
 * all coefficients into one aligned arena.
 */
template <ToeplitzMatrixBlock Block> class MultilevelToeplitzStorage {
  public:
    using block_type = Block;
    using scalar_type = typename ToeplitzBlockTraits<Block>::scalar_type;
    using index_type = Eigen::Index;
    using extent_span = std::span<const index_type>;
    static constexpr bool has_packed_coefficients = toeplitz_new_detail::PackedTree<Block>::value;

  private:
    template <ToeplitzMatrixBlock> friend class MultilevelToeplitz;
    int levels_ = 0;
    std::vector<index_type> block_rows_;
    std::vector<index_type> block_cols_;
    index_type rows_ = 0, cols_ = 0;
    toeplitz_new_detail::Storage<Block> blocks_;
    toeplitz_new_detail::FlatToeplitzIndex plan_;

#ifndef NDEBUG
    static void check_description(int levels, extent_span block_rows, extent_span block_cols) {
        if (levels <= 0)
            throw std::invalid_argument("Multilevel Toeplitz needs at least one level");
        if (std::cmp_not_equal(block_rows.size(), levels) || std::cmp_not_equal(block_cols.size(), levels))
            throw std::invalid_argument("Multilevel Toeplitz extent count must equal levels");
        for (int level = 0; level < levels; ++level)
            if (block_rows[static_cast<std::size_t>(level)] <= 0 ||
                block_cols[static_cast<std::size_t>(level)] <= 0)
                throw std::invalid_argument("Multilevel Toeplitz extents must be positive");
    }
#endif

    static std::vector<index_type> repeated_extents(int levels, index_type value) {
#ifndef NDEBUG
        if (levels <= 0 || value <= 0)
            throw std::invalid_argument("Multilevel Toeplitz levels and extents must be positive");
#endif
        return std::vector<index_type>(static_cast<std::size_t>(levels), value);
    }

    static index_type leaf_count(extent_span block_rows, extent_span block_cols) {
        index_type result = 1;
        for (std::size_t level = 0; level < block_rows.size(); ++level)
            result = toeplitz_new_detail::checked_product(
                result, toeplitz_new_detail::unique_count(block_rows[level], block_cols[level]));
        return result;
    }

    template <typename Factory>
    static std::vector<Block> generate(int levels, extent_span block_rows, extent_span block_cols,
                                       Factory &&factory) {
#ifndef NDEBUG
        check_description(levels, block_rows, block_cols);
#endif
        const index_type count = leaf_count(block_rows, block_cols);
        std::vector<Block> blocks;
        blocks.reserve(static_cast<std::size_t>(count));
        std::vector<index_type> row(static_cast<std::size_t>(levels));
        std::vector<index_type> col(static_cast<std::size_t>(levels));
        for (index_type leaf = 0; leaf < count; ++leaf) {
            index_type encoded = leaf;
            for (int reverse = levels; reverse-- != 0;) {
                const auto position = static_cast<std::size_t>(reverse);
                const index_type unique =
                    toeplitz_new_detail::unique_count(block_rows[position], block_cols[position]);
                const index_type diagonal = encoded % unique;
                encoded /= unique;
                if (diagonal < block_cols[position]) {
                    row[position] = 0;
                    col[position] = diagonal;
                } else {
                    row[position] = diagonal - block_cols[position] + 1;
                    col[position] = 0;
                }
            }
            blocks.emplace_back(std::invoke(factory, extent_span(row), extent_span(col)));
        }
        return blocks;
    }

    void initialize_plan() {
        if (levels_ == 0 || plan_.leaf_blocks == 0)
            return;
        plan_.levels.resize(static_cast<std::size_t>(levels_));
        index_type rows_in_block = plan_.leaf_rows;
        index_type cols_in_block = plan_.leaf_cols;
        for (int reverse = levels_; reverse-- != 0;) {
            const auto position = static_cast<std::size_t>(reverse);
            plan_.levels[position] = toeplitz_new_detail::FlatToeplitzLevel{
                block_rows_[position], block_cols_[position], rows_in_block, cols_in_block,
                toeplitz_new_detail::unique_count(block_rows_[position], block_cols_[position])};
            rows_in_block = toeplitz_new_detail::checked_product(block_rows_[position], rows_in_block);
            cols_in_block = toeplitz_new_detail::checked_product(block_cols_[position], cols_in_block);
        }
        rows_ = rows_in_block;
        cols_ = cols_in_block;
        plan_.enabled = true;
        toeplitz_new_detail::build_flat_runs(plan_);
    }

  public:
    MultilevelToeplitzStorage() = default;

    MultilevelToeplitzStorage(int levels, std::vector<index_type> block_rows,
                             std::vector<index_type> block_cols, std::vector<Block> blocks)
        : levels_(levels), block_rows_(std::move(block_rows)), block_cols_(std::move(block_cols)) {
#ifndef NDEBUG
        check_description(levels_, block_rows_, block_cols_);
        if (std::cmp_not_equal(blocks.size(), leaf_count(block_rows_, block_cols_)))
            throw std::invalid_argument("Multilevel Toeplitz received the wrong number of leaf blocks");
#endif
        if (blocks.empty())
            return;
        plan_.leaf_rows = toeplitz_new_detail::positive_extent(blocks.front().rows());
        plan_.leaf_cols = toeplitz_new_detail::positive_extent(blocks.front().cols());
        plan_.leaf_blocks = leaf_count(block_rows_, block_cols_);
#ifndef NDEBUG
        for (const auto &block : blocks)
            if (std::cmp_not_equal(block.rows(), plan_.leaf_rows) ||
                std::cmp_not_equal(block.cols(), plan_.leaf_cols))
                throw std::invalid_argument("Multilevel Toeplitz leaf blocks must have identical dimensions");
#endif
        blocks_ = toeplitz_new_detail::Storage<Block>(std::move(blocks));
        initialize_plan();
    }

    template <typename Factory>
        requires std::invocable<Factory &, extent_span, extent_span> &&
                 std::convertible_to<std::invoke_result_t<Factory &, extent_span, extent_span>, Block>
    MultilevelToeplitzStorage(int levels, std::vector<index_type> block_rows,
                             std::vector<index_type> block_cols, Factory &&factory)
        : MultilevelToeplitzStorage(
              levels, block_rows, block_cols,
              generate(levels, extent_span(block_rows), extent_span(block_cols), std::forward<Factory>(factory))) {}

    template <typename Factory>
        requires std::invocable<Factory &, extent_span, extent_span> &&
                 std::convertible_to<std::invoke_result_t<Factory &, extent_span, extent_span>, Block>
    MultilevelToeplitzStorage(int levels, index_type block_rows, index_type block_cols, Factory &&factory)
        : MultilevelToeplitzStorage(levels, repeated_extents(levels, block_rows),
                                    repeated_extents(levels, block_cols), std::forward<Factory>(factory)) {}

    MultilevelToeplitzStorage(const MultilevelToeplitzStorage &) = default;
    MultilevelToeplitzStorage(MultilevelToeplitzStorage &&other) noexcept { swap(other); }

    MultilevelToeplitzStorage &operator=(MultilevelToeplitzStorage other) noexcept {
        swap(other);
        return *this;
    }

    void swap(MultilevelToeplitzStorage &other) noexcept {
        using std::swap;
        swap(levels_, other.levels_);
        block_rows_.swap(other.block_rows_);
        block_cols_.swap(other.block_cols_);
        swap(rows_, other.rows_);
        swap(cols_, other.cols_);
        swap(blocks_, other.blocks_);
        plan_.swap(other.plan_);
    }

    [[nodiscard]] int levels() const noexcept { return levels_; }
    [[nodiscard]] index_type rows() const noexcept { return rows_; }
    [[nodiscard]] index_type cols() const noexcept { return cols_; }
    [[nodiscard]] index_type leaf_rows() const noexcept { return plan_.leaf_rows; }
    [[nodiscard]] index_type leaf_cols() const noexcept { return plan_.leaf_cols; }
    [[nodiscard]] index_type unique_leaf_blocks() const noexcept { return plan_.leaf_blocks; }
    [[nodiscard]] extent_span block_rows() const noexcept { return block_rows_; }
    [[nodiscard]] extent_span block_cols() const noexcept { return block_cols_; }

    [[nodiscard]] decltype(auto) unique_leaf_block(index_type leaf) const {
#ifndef NDEBUG
        if (leaf < 0 || leaf >= unique_leaf_blocks())
            throw std::out_of_range("Multilevel Toeplitz leaf index out of range");
#endif
        return blocks_.get(leaf);
    }

    [[nodiscard]] std::span<const scalar_type> packed_coefficients() const requires(has_packed_coefficients) {
        return blocks_.coefficients();
    }
};

/** Runtime-depth multilevel Toeplitz arithmetic facade.
 *
 * The constructor builds the complete compressed gather/scatter plan. matvec
 * subsequently performs a linear traversal of leaf blocks and plan arrays; it
 * neither descends through nested Toeplitz objects nor converts mixed-radix
 * indices. Two- and three-level plan builders are compile-time specializations
 * in TeoplitzTemplateNew.hpp; arbitrary positive runtime depth uses the generic
 * construction fallback and the same branch-free multiplication loop.
 */
template <ToeplitzMatrixBlock Block> class MultilevelToeplitz {
  public:
    using block_type = Block;
    using traits_type = ToeplitzBlockTraits<Block>;
    using scalar_type = typename traits_type::scalar_type;
    using index_type = Eigen::Index;
    using dense_matrix_t = ToeplitzDenseMatrix<scalar_type>;
    using vector_t = Eigen::Matrix<scalar_type, Eigen::Dynamic, 1>;
    using const_matrix_ref = ToeplitzConstMatrixRef<scalar_type>;
    using matrix_ref = ToeplitzMatrixRef<scalar_type>;
    using const_vector_ref = Eigen::Ref<const vector_t>;
    using vector_ref = Eigen::Ref<vector_t>;
    using workspace_type = ToeplitzWorkspace<scalar_type, typename traits_type::workspace_type>;
    using storage_type = MultilevelToeplitzStorage<Block>;
    using extent_span = typename storage_type::extent_span;

  private:
    storage_type storage_;

#ifndef NDEBUG
    void check_product_shapes(const_matrix_ref rhs, matrix_ref dest) const {
        if (rhs.rows() != cols() || dest.rows() != rows() || rhs.cols() != dest.cols())
            throw std::invalid_argument("Multilevel Toeplitz multiplication dimension mismatch");
        if (rhs.outerStride() < 0 || dest.outerStride() < 0)
            throw std::invalid_argument("Multilevel Toeplitz multiplication requires nonnegative column strides");
    }
#endif

    void apply(const_matrix_ref rhs, matrix_ref dest, workspace_type &workspace, bool add) const {
#ifndef NDEBUG
        check_product_shapes(rhs, dest);
#endif
        if (toeplitz_new_detail::overlaps<scalar_type>(rhs, dest)) {
            workspace.alias_copy = rhs;
            if (!add)
                dest.setZero();
            multiply_add_noalias(workspace.alias_copy, dest, workspace);
        } else {
            if (!add)
                dest.setZero();
            multiply_add_noalias(rhs, dest, workspace);
        }
    }

  public:
    MultilevelToeplitz() = default;
    explicit MultilevelToeplitz(storage_type storage) : storage_(std::move(storage)) {}

    MultilevelToeplitz(int levels, std::vector<index_type> block_rows,
                       std::vector<index_type> block_cols, std::vector<Block> blocks)
        : storage_(levels, std::move(block_rows), std::move(block_cols), std::move(blocks)) {}

    template <typename Factory>
        requires std::invocable<Factory &, extent_span, extent_span> &&
                 std::convertible_to<std::invoke_result_t<Factory &, extent_span, extent_span>, Block>
    MultilevelToeplitz(int levels, std::vector<index_type> block_rows,
                       std::vector<index_type> block_cols, Factory &&factory)
        : storage_(levels, std::move(block_rows), std::move(block_cols), std::forward<Factory>(factory)) {}

    template <typename Factory>
        requires std::invocable<Factory &, extent_span, extent_span> &&
                 std::convertible_to<std::invoke_result_t<Factory &, extent_span, extent_span>, Block>
    MultilevelToeplitz(int levels, index_type block_rows, index_type block_cols, Factory &&factory)
        : storage_(levels, block_rows, block_cols, std::forward<Factory>(factory)) {}

    [[nodiscard]] const storage_type &get_storage() const noexcept { return storage_; }
    [[nodiscard]] int levels() const noexcept { return storage_.levels(); }
    [[nodiscard]] index_type rows() const noexcept { return storage_.rows(); }
    [[nodiscard]] index_type cols() const noexcept { return storage_.cols(); }
    [[nodiscard]] index_type leaf_rows() const noexcept { return storage_.leaf_rows(); }
    [[nodiscard]] index_type leaf_cols() const noexcept { return storage_.leaf_cols(); }
    [[nodiscard]] index_type unique_leaf_blocks() const noexcept { return storage_.unique_leaf_blocks(); }
    [[nodiscard]] extent_span block_rows() const noexcept { return storage_.block_rows(); }
    [[nodiscard]] extent_span block_cols() const noexcept { return storage_.block_cols(); }
    [[nodiscard]] decltype(auto) unique_leaf_block(index_type leaf) const {
        return storage_.unique_leaf_block(leaf);
    }

    void multiply_add_noalias(const_matrix_ref rhs, matrix_ref dest, workspace_type &workspace) const {
#ifndef NDEBUG
        check_product_shapes(rhs, dest);
#endif
        if (storage_.plan_.leaf_blocks == 0 || rhs.cols() == 0)
            return;
        const auto &plan = storage_.plan_;
        const index_type q = rhs.cols();
        for (index_type leaf = 0; leaf < plan.leaf_blocks; ++leaf) {
            const auto leaf_position = static_cast<std::size_t>(leaf);
            const index_type occurrences = plan.leaf_occurrences[leaf_position];
            const index_type run_length = plan.leaf_run_length[leaf_position];
            const index_type run_begin = plan.leaf_run_begin[leaf_position];
            const index_type run_end = plan.leaf_run_begin[leaf_position + 1];
            decltype(auto) block = storage_.blocks_.get(leaf);

            if (occurrences == 1 && q == 1) {
                Eigen::Map<const dense_matrix_t> x(
                    rhs.data() + plan.run_source_offset[static_cast<std::size_t>(run_begin)],
                    plan.leaf_cols, 1);
                Eigen::Map<dense_matrix_t> y(
                    dest.data() + plan.run_destination_offset[static_cast<std::size_t>(run_begin)],
                    plan.leaf_rows, 1);
                traits_type::multiply_add(block, x, y, workspace.child);
                continue;
            }

            const index_type width = toeplitz_new_detail::checked_product(occurrences, q);
            toeplitz_new_detail::grow(
                workspace.rhs_panel, toeplitz_new_detail::checked_product(plan.leaf_cols, width));
            toeplitz_new_detail::grow(
                workspace.product_panel, toeplitz_new_detail::checked_product(plan.leaf_rows, width));
            Eigen::Map<dense_matrix_t> x(workspace.rhs_panel.data(), plan.leaf_cols, width);
            Eigen::Map<dense_matrix_t> y(workspace.product_panel.data(), plan.leaf_rows, width);
            for (index_type column = 0; column < q; ++column)
                for (index_type run = run_begin; run < run_end; ++run) {
                    const index_type panel_column = plan.run_panel_offset[static_cast<std::size_t>(run)];
                    Eigen::Map<const dense_matrix_t> source(
                        rhs.data() + column * rhs.outerStride() +
                            plan.run_source_offset[static_cast<std::size_t>(run)],
                        plan.leaf_cols, run_length);
                    x.middleCols(column * occurrences + panel_column, run_length) = source;
                }
            y.setZero();
            traits_type::multiply_add(block, x, y, workspace.child);
            for (index_type column = 0; column < q; ++column)
                for (index_type run = run_begin; run < run_end; ++run) {
                    const index_type panel_column = plan.run_panel_offset[static_cast<std::size_t>(run)];
                    Eigen::Map<dense_matrix_t> target(
                        dest.data() + column * dest.outerStride() +
                            plan.run_destination_offset[static_cast<std::size_t>(run)],
                        plan.leaf_rows, run_length);
                    target += y.middleCols(column * occurrences + panel_column, run_length);
                }
        }
    }

    void matmul(const_matrix_ref rhs, matrix_ref dest, workspace_type &workspace) const {
        apply(rhs, dest, workspace, false);
    }
    void matmul_add(const_matrix_ref rhs, matrix_ref dest, workspace_type &workspace) const {
        apply(rhs, dest, workspace, true);
    }
    void matmul(const_matrix_ref rhs, matrix_ref dest) const {
        workspace_type workspace;
        matmul(rhs, dest, workspace);
    }
    void matmul_add(const_matrix_ref rhs, matrix_ref dest) const {
        workspace_type workspace;
        matmul_add(rhs, dest, workspace);
    }
    [[nodiscard]] dense_matrix_t matmul(const_matrix_ref rhs) const {
#ifndef NDEBUG
        (void)toeplitz_new_detail::checked_product(rows(), rhs.cols());
#endif
        dense_matrix_t result(rows(), rhs.cols());
        workspace_type workspace;
        matmul(rhs, result, workspace);
        return result;
    }

    void matvec(const_vector_ref rhs, vector_ref dest, workspace_type &workspace) const {
        matmul(rhs, dest, workspace);
    }
    void matvec_add(const_vector_ref rhs, vector_ref dest, workspace_type &workspace) const {
        matmul_add(rhs, dest, workspace);
    }
    void matvec(const_vector_ref rhs, vector_ref dest) const {
        workspace_type workspace;
        matvec(rhs, dest, workspace);
    }
    void matvec_add(const_vector_ref rhs, vector_ref dest) const {
        workspace_type workspace;
        matvec_add(rhs, dest, workspace);
    }
    [[nodiscard]] vector_t matvec(const_vector_ref rhs) const {
        vector_t result(rows());
        workspace_type workspace;
        matvec(rhs, result, workspace);
        return result;
    }
    [[nodiscard]] vector_t operator*(const_vector_ref rhs) const { return matvec(rhs); }
};

} // namespace EMW::Math::LinAgl::Matrix

#endif // EMW_MATH_MATRIX_TOEPLITZ_MULTILEVEL_NEW_HPP
