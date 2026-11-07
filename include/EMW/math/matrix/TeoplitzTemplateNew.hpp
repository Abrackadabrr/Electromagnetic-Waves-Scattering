#ifndef EMW_MATH_MATRIX_TEOPLITZ_TEMPLATE_NEW_HPP
#define EMW_MATH_MATRIX_TEOPLITZ_TEMPLATE_NEW_HPP

#include <Eigen/Core>

#include <algorithm>
#include <array>
#include <concepts>
#include <functional>
#include <limits>
#include <memory>
#include <span>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

namespace EMW::Math::LinAgl::Matrix {

template <typename Scalar>
using ToeplitzDenseMatrix = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;
template <typename Scalar> using ToeplitzConstMatrixRef = Eigen::Ref<const ToeplitzDenseMatrix<Scalar>>;
template <typename Scalar> using ToeplitzMatrixRef = Eigen::Ref<ToeplitzDenseMatrix<Scalar>>;

struct ToeplitzEmptyWorkspace {};

/** Static customization point for a block type.
 *
 * A native block supplies scalar_type, workspace_type, integral rows()/cols(),
 * coefficient access operator()(i,j), and multiply_add_noalias(X,Y,workspace).
 * The latter MUST implement Y += B*X for ANY number of columns, without resizing
 * Y or retaining references to X/Y/workspace. X and Y do not overlap. All scalar
 * products are ordinary products, not adjoints. A default-constructed workspace
 * must work. Specialize this trait to adapt an existing block with another API.
 *
 * Eigen owning matrices are supported below; Eigen expressions/maps should be
 * materialized into owning matrices when constructing the Toeplitz object.
 */
template <typename Block, typename = void> struct ToeplitzBlockTraits {};

template <typename Block>
struct ToeplitzBlockTraits<Block, std::void_t<typename Block::scalar_type, typename Block::workspace_type>> {
    using scalar_type = typename Block::scalar_type;
    using workspace_type = typename Block::workspace_type;

    static void multiply_add(const Block &block, ToeplitzConstMatrixRef<scalar_type> rhs,
                             ToeplitzMatrixRef<scalar_type> dest, workspace_type &workspace)
        requires requires { block.multiply_add_noalias(rhs, dest, workspace); }
    {
        block.multiply_add_noalias(rhs, dest, workspace);
    }
};

template <typename Scalar, int Rows, int Cols, int Options, int MaxRows, int MaxCols>
struct ToeplitzBlockTraits<Eigen::Matrix<Scalar, Rows, Cols, Options, MaxRows, MaxCols>> {
    using scalar_type = Scalar;
    using workspace_type = ToeplitzEmptyWorkspace;

    // Also accepts the read-only Map returned by the packed storage below.
    template <typename Derived>
    static void multiply_add(const Eigen::MatrixBase<Derived> &block, ToeplitzConstMatrixRef<Scalar> rhs,
                             ToeplitzMatrixRef<Scalar> dest, workspace_type &) {
        dest.noalias() += block.derived() * rhs;
    }
};

template <typename Block>
concept ToeplitzMatrixBlock = std::copy_constructible<Block> && requires(
    const Block &block, Eigen::Index i, ToeplitzConstMatrixRef<typename ToeplitzBlockTraits<Block>::scalar_type> rhs,
    ToeplitzMatrixRef<typename ToeplitzBlockTraits<Block>::scalar_type> dest,
    typename ToeplitzBlockTraits<Block>::workspace_type &workspace) {
    requires std::default_initializable<typename ToeplitzBlockTraits<Block>::workspace_type>;
    requires std::integral<std::remove_cvref_t<decltype(block.rows())>>;
    requires std::integral<std::remove_cvref_t<decltype(block.cols())>>;
    { block(i, i) } -> std::convertible_to<typename ToeplitzBlockTraits<Block>::scalar_type>;
    { ToeplitzBlockTraits<Block>::multiply_add(block, rhs, dest, workspace) } -> std::same_as<void>;
};

template <ToeplitzMatrixBlock Block> class ToeplitzTemplateNew;

namespace toeplitz_new_detail {

using index = Eigen::Index;
template <typename Scalar> using Buffer = std::vector<Scalar, Eigen::aligned_allocator<Scalar>>;

struct FlatToeplitzLevel {
    index block_rows = 0;
    index block_cols = 0;
    index rows_in_block = 0;
    index cols_in_block = 0;
    index unique_blocks = 0;
};

struct FlatToeplitzDiagonal {
    index row_start = 0;
    index col_start = 0;
    index repeats = 0;
};

struct FlatToeplitzIndex {
    std::vector<FlatToeplitzLevel> levels;
    // Fully precomputed compressed gather/scatter plan. Runs contain adjacent
    // occurrences along the innermost Toeplitz diagonal.
    std::vector<index> leaf_run_begin;
    std::vector<index> leaf_occurrences;
    std::vector<index> leaf_run_length;
    std::vector<index> run_source_offset;
    std::vector<index> run_destination_offset;
    std::vector<index> run_panel_offset;
    // Scalar offset of the persistent gathered panel for every leaf. The last
    // entry is the total prepared-vector storage size.
    std::vector<index> leaf_panel_offset;
    index leaf_rows = 0;
    index leaf_cols = 0;
    index leaf_blocks = 0;
    bool enabled = false;
    bool recurrent_layout = false;

    void reset() noexcept {
        levels.clear();
        leaf_run_begin.clear();
        leaf_occurrences.clear();
        leaf_run_length.clear();
        run_source_offset.clear();
        run_destination_offset.clear();
        run_panel_offset.clear();
        leaf_panel_offset.clear();
        leaf_rows = leaf_cols = leaf_blocks = 0;
        enabled = false;
        recurrent_layout = false;
    }

    void swap(FlatToeplitzIndex &other) noexcept {
        levels.swap(other.levels);
        leaf_run_begin.swap(other.leaf_run_begin);
        leaf_occurrences.swap(other.leaf_occurrences);
        leaf_run_length.swap(other.leaf_run_length);
        run_source_offset.swap(other.run_source_offset);
        run_destination_offset.swap(other.run_destination_offset);
        run_panel_offset.swap(other.run_panel_offset);
        leaf_panel_offset.swap(other.leaf_panel_offset);
        std::swap(leaf_rows, other.leaf_rows);
        std::swap(leaf_cols, other.leaf_cols);
        std::swap(leaf_blocks, other.leaf_blocks);
        std::swap(enabled, other.enabled);
        std::swap(recurrent_layout, other.recurrent_layout);
    }
};

inline index checked_add(index a, index b) {
#ifndef NDEBUG
    if (a < 0 || b < 0 || a > std::numeric_limits<index>::max() - b)
        throw std::length_error("Toeplitz size overflow");
#endif
    return a + b;
}

inline index checked_product(index a, index b) {
#ifndef NDEBUG
    if (a < 0 || b < 0 || (b != 0 && a > std::numeric_limits<index>::max() / b))
        throw std::length_error("Toeplitz size overflow");
#endif
    return a * b;
}

template <std::integral Integer> index positive_extent(Integer value) {
#ifndef NDEBUG
    if (!std::in_range<index>(value) || value <= 0)
        throw std::invalid_argument("Toeplitz blocks must have positive representable dimensions");
#endif
    return static_cast<index>(value);
}

inline index unique_count(index rows, index cols) {
    if (rows == 0 && cols == 0)
        return 0;
#ifndef NDEBUG
    if (rows <= 0 || cols <= 0)
        throw std::invalid_argument("Toeplitz grid must be positive, or 0 by 0");
#endif
    return checked_add(rows, cols - 1);
}

inline FlatToeplitzDiagonal flat_diagonal(const FlatToeplitzLevel &level, index diagonal) {
    FlatToeplitzDiagonal result;
    if (diagonal < level.block_cols) {
        result.col_start = diagonal;
    } else {
        result.row_start = diagonal - level.block_cols + 1;
    }
    result.repeats = std::min(level.block_rows - result.row_start,
                              level.block_cols - result.col_start);
    return result;
}

inline void append_flat_run(FlatToeplitzIndex &plan, index source, index destination) {
    plan.run_source_offset.push_back(source);
    plan.run_destination_offset.push_back(destination);
}

template <std::size_t Levels> struct FlatRunBuilder;

template <> struct FlatRunBuilder<1> {
    static void append(FlatToeplitzIndex &plan, const std::array<FlatToeplitzDiagonal, 1> &diagonals,
                       index source, index destination) {
        (void)diagonals;
        append_flat_run(plan, source, destination);
    }
};

template <> struct FlatRunBuilder<2> {
    static void append(FlatToeplitzIndex &plan, const std::array<FlatToeplitzDiagonal, 2> &diagonals,
                       index source, index destination) {
        const auto &level = plan.levels[0];
        for (index first = 0; first < diagonals[0].repeats; ++first)
            append_flat_run(plan, source + first * level.cols_in_block,
                            destination + first * level.rows_in_block);
    }
};

template <> struct FlatRunBuilder<3> {
    static void append(FlatToeplitzIndex &plan, const std::array<FlatToeplitzDiagonal, 3> &diagonals,
                       index source, index destination) {
        const auto &outer = plan.levels[0];
        const auto &middle = plan.levels[1];
        for (index first = 0; first < diagonals[0].repeats; ++first)
            for (index second = 0; second < diagonals[1].repeats; ++second)
                append_flat_run(plan,
                                source + first * outer.cols_in_block + second * middle.cols_in_block,
                                destination + first * outer.rows_in_block + second * middle.rows_in_block);
    }
};

template <std::size_t Levels> void append_flat_leaf(FlatToeplitzIndex &plan, index leaf) {
    static_assert(Levels > 0);
    std::array<FlatToeplitzDiagonal, Levels> diagonals;
    index encoded = leaf;
    index source = 0, destination = 0, occurrences = 1;
    for (std::size_t reverse = Levels; reverse-- != 0;) {
        const auto &level = plan.levels[reverse];
        diagonals[reverse] = flat_diagonal(level, encoded % level.unique_blocks);
        encoded /= level.unique_blocks;
        source += diagonals[reverse].col_start * level.cols_in_block;
        destination += diagonals[reverse].row_start * level.rows_in_block;
        occurrences = checked_product(occurrences, diagonals[reverse].repeats);
    }
    plan.leaf_occurrences[static_cast<std::size_t>(leaf)] = occurrences;
    plan.leaf_run_length[static_cast<std::size_t>(leaf)] = diagonals.back().repeats;
    FlatRunBuilder<Levels>::append(plan, diagonals, source, destination);
}

inline void append_flat_leaf_dynamic(FlatToeplitzIndex &plan, index leaf,
                                     std::vector<FlatToeplitzDiagonal> &diagonals,
                                     std::vector<index> &position) {
    const std::size_t levels = plan.levels.size();
    index encoded = leaf;
    index source_base = 0, destination_base = 0, occurrences = 1;
    for (std::size_t reverse = levels; reverse-- != 0;) {
        const auto &level = plan.levels[reverse];
        diagonals[reverse] = flat_diagonal(level, encoded % level.unique_blocks);
        encoded /= level.unique_blocks;
        source_base += diagonals[reverse].col_start * level.cols_in_block;
        destination_base += diagonals[reverse].row_start * level.rows_in_block;
        occurrences = checked_product(occurrences, diagonals[reverse].repeats);
    }
    plan.leaf_occurrences[static_cast<std::size_t>(leaf)] = occurrences;
    plan.leaf_run_length[static_cast<std::size_t>(leaf)] = diagonals.back().repeats;

    std::fill(position.begin(), position.end(), index{0});
    const index groups = occurrences / diagonals.back().repeats;
    for (index group = 0; group < groups; ++group) {
        index source = source_base, destination = destination_base;
        for (std::size_t level_index = 0; level_index + 1 < levels; ++level_index) {
            source += position[level_index] * plan.levels[level_index].cols_in_block;
            destination += position[level_index] * plan.levels[level_index].rows_in_block;
        }
        append_flat_run(plan, source, destination);
        for (std::size_t reverse = levels - 1; reverse-- != 0;) {
            if (++position[reverse] < diagonals[reverse].repeats)
                break;
            position[reverse] = 0;
        }
    }
}

inline void build_flat_runs(FlatToeplitzIndex &plan) {
    const std::size_t levels = plan.levels.size();
    plan.leaf_run_begin.assign(static_cast<std::size_t>(plan.leaf_blocks + 1), index{0});
    plan.leaf_occurrences.resize(static_cast<std::size_t>(plan.leaf_blocks));
    plan.leaf_run_length.resize(static_cast<std::size_t>(plan.leaf_blocks));
    plan.run_source_offset.clear();
    plan.run_destination_offset.clear();
    plan.run_panel_offset.clear();

    std::vector<FlatToeplitzDiagonal> diagonals(levels);
    std::vector<index> position(levels > 0 ? levels - 1 : 0);
    for (index leaf = 0; leaf < plan.leaf_blocks; ++leaf) {
        plan.leaf_run_begin[static_cast<std::size_t>(leaf)] =
            static_cast<index>(plan.run_source_offset.size());
        switch (levels) {
        case 1:
            append_flat_leaf<1>(plan, leaf);
            break;
        case 2:
            append_flat_leaf<2>(plan, leaf);
            break;
        case 3:
            append_flat_leaf<3>(plan, leaf);
            break;
        default:
            append_flat_leaf_dynamic(plan, leaf, diagonals, position);
            break;
        }
    }
    plan.leaf_run_begin[static_cast<std::size_t>(plan.leaf_blocks)] =
        static_cast<index>(plan.run_source_offset.size());
    plan.run_panel_offset.resize(plan.run_source_offset.size());
    plan.leaf_panel_offset.resize(static_cast<std::size_t>(plan.leaf_blocks + 1));
    index panel_offset = 0;
    for (index leaf = 0; leaf < plan.leaf_blocks; ++leaf) {
        const auto position = static_cast<std::size_t>(leaf);
        plan.leaf_panel_offset[position] = panel_offset;
        const index begin = plan.leaf_run_begin[position];
        const index end = plan.leaf_run_begin[position + 1];
        const index length = plan.leaf_run_length[position];
        for (index run = begin; run < end; ++run)
            plan.run_panel_offset[static_cast<std::size_t>(run)] = (run - begin) * length;
        panel_offset = checked_add(
            panel_offset,
            checked_product(plan.leaf_cols, plan.leaf_occurrences[position]));
    }
    plan.leaf_panel_offset[static_cast<std::size_t>(plan.leaf_blocks)] = panel_offset;
}

template <typename Scalar> void grow(Buffer<Scalar> &buffer, index size) {
    if (std::cmp_less(buffer.size(), size))
        buffer.resize(static_cast<std::size_t>(size));
}

template <typename Block> struct PackedTree : std::false_type {};
template <typename S, int R, int C, int O, int MR, int MC>
struct PackedTree<Eigen::Matrix<S, R, C, O, MR, MC>> : std::true_type {};
template <ToeplitzMatrixBlock Block>
struct PackedTree<ToeplitzTemplateNew<Block>> : PackedTree<Block> {};

template <typename Block> struct PackedLeafMatrix;
template <typename S, int R, int C, int O, int MR, int MC>
struct PackedLeafMatrix<Eigen::Matrix<S, R, C, O, MR, MC>> {
    // The arena is always column-major. A fixed row vector uses RowMajor
    // because Eigen forbids a column-major Matrix<Scalar,1,N>; both layouts
    // are byte-identical for one row.
    static constexpr int options = R == 1 && C != 1 ? Eigen::RowMajor : Eigen::ColMajor;
    using type = Eigen::Matrix<S, R, C, options, MR, MC>;
};
template <ToeplitzMatrixBlock Block>
struct PackedLeafMatrix<ToeplitzTemplateNew<Block>> : PackedLeafMatrix<Block> {};

template <typename Block, bool = PackedTree<Block>::value> struct PackedArithmeticMatrices {
    using scalar_type = typename ToeplitzBlockTraits<Block>::scalar_type;
    using leaf_type = ToeplitzDenseMatrix<scalar_type>;
    using rhs_panel_type = ToeplitzDenseMatrix<scalar_type>;
    using product_panel_type = ToeplitzDenseMatrix<scalar_type>;
};

template <typename Block> struct PackedArithmeticMatrices<Block, true> {
    using scalar_type = typename ToeplitzBlockTraits<Block>::scalar_type;
    using leaf_type = typename PackedLeafMatrix<Block>::type;
    using rhs_panel_type = Eigen::Matrix<scalar_type, leaf_type::ColsAtCompileTime,
                                         Eigen::Dynamic, Eigen::ColMajor>;
    using product_panel_type = Eigen::Matrix<scalar_type, leaf_type::RowsAtCompileTime,
                                             Eigen::Dynamic, Eigen::ColMajor>;
};

/** Generic blocks: contiguous objects. Nested Toeplitz blocks additionally
 * consolidate all dense leaf coefficients into ONE aligned allocation.
 * Opaque custom blocks retain their own internal allocation policy.
 */
template <typename Block> class Storage {
    using scalar_type = typename ToeplitzBlockTraits<Block>::scalar_type;
    using buffer_type = Buffer<scalar_type>;
    std::vector<Block> blocks_;
    index coefficient_count_ = 0;

  public:
    static constexpr bool packed = PackedTree<Block>::value;
    Storage() = default;

    explicit Storage(std::vector<Block> blocks) : blocks_(std::move(blocks)) {
        if constexpr (packed) {
            for (const auto &block : blocks_)
                coefficient_count_ = checked_add(coefficient_count_, block.storage_.coefficient_count());
            if (coefficient_count_ != 0) {
                auto arena = std::make_shared<buffer_type>(static_cast<std::size_t>(coefficient_count_));
                index offset = 0;
                repack(arena, offset);
            }
        }
    }

    [[nodiscard]] const Block &get(index i) const { return blocks_[static_cast<std::size_t>(i)]; }
    [[nodiscard]] index coefficient_count() const noexcept { return coefficient_count_; }

    [[nodiscard]] std::span<const scalar_type> coefficients() const requires(packed) {
        if (blocks_.empty())
            return {};
        return {blocks_.front().storage_.coefficients().data(), static_cast<std::size_t>(coefficient_count_)};
    }

    void repack(const std::shared_ptr<buffer_type> &arena, index &offset) requires(packed) {
        for (auto &block : blocks_)
            block.storage_.repack(arena, offset);
    }
};

/** Each block is column-major, immediately followed by the next block.
 * No per-block padding and no per-block heap allocation survive construction.
 * The base is Eigen-aligned; individual block starts need not be aligned, so
 * all Maps deliberately use Eigen::Unaligned. Never claim stronger alignment.
 */
template <typename S, int R, int C, int O, int MR, int MC>
class Storage<Eigen::Matrix<S, R, C, O, MR, MC>> {
    using block_type = Eigen::Matrix<S, R, C, O, MR, MC>;
    using matrix_type = typename PackedLeafMatrix<block_type>::type;
    using buffer_type = Buffer<S>;
    std::shared_ptr<const buffer_type> arena_;
    index rows_ = 0, cols_ = 0, offset_ = 0, size_ = 0;

  public:
    static constexpr bool packed = true;
    Storage() = default;

    explicit Storage(std::vector<block_type> blocks) {
        if (blocks.empty())
            return;
        rows_ = blocks.front().rows();
        cols_ = blocks.front().cols();
        const index stride = checked_product(rows_, cols_);
        size_ = checked_product(stride, positive_extent(blocks.size()));
        auto arena = std::make_shared<buffer_type>(static_cast<std::size_t>(size_));
        index offset = 0;
        for (const auto &block : blocks) {
            Eigen::Map<matrix_type>(arena->data() + offset, rows_, cols_) = block;
            offset += stride;
        }
        arena_ = std::move(arena);
    }

    [[nodiscard]] Eigen::Map<const matrix_type> get(index i) const {
        return {arena_->data() + offset_ + i * rows_ * cols_, rows_, cols_};
    }
    [[nodiscard]] index coefficient_count() const noexcept { return arena_ ? size_ : 0; }
    [[nodiscard]] std::span<const S> coefficients() const {
        if (!arena_)
            return {};
        return {arena_->data() + offset_, static_cast<std::size_t>(size_)};
    }

    void repack(const std::shared_ptr<buffer_type> &arena, index &offset) {
        if (size_ == 0)
            return;
        std::copy_n(arena_->data() + offset_, size_, arena->data() + offset);
        arena_ = arena;
        offset_ = offset;
        offset += size_;
    }
};

template <typename Scalar>
bool overlaps(ToeplitzConstMatrixRef<Scalar> rhs, ToeplitzMatrixRef<Scalar> dest) {
    if (rhs.size() == 0 || dest.size() == 0)
        return false;
    // The public API accepts positive column strides (including submatrices).
    const auto extent = [](const auto &matrix) {
        return checked_add(checked_product(matrix.cols() - 1, matrix.outerStride()), matrix.rows());
    };
    const std::less<const Scalar *> less;
    return less(rhs.data(), dest.data() + extent(dest)) && less(dest.data(), rhs.data() + extent(rhs));
}

} // namespace toeplitz_new_detail

/** RAII owner of a block Toeplitz representation, independent of arithmetic.
 * Owns all dimensions, block metadata, and numerical storage. Construction
 * checks the shape and consolidates dense leaves; destruction releases them.
 * Copying shares immutable coefficient buffers, moving leaves an empty owner.
 * No matvec/matmul or multiplication workspace belongs to this class.
 */
template <ToeplitzMatrixBlock Block> class ToeplitzStorage {
  public:
    using block_type = Block;
    using scalar_type = typename ToeplitzBlockTraits<Block>::scalar_type;
    using index_type = Eigen::Index;
    static constexpr bool has_packed_coefficients = toeplitz_new_detail::PackedTree<Block>::value;

  private:
    template <typename> friend class toeplitz_new_detail::Storage;
    using payload_type = toeplitz_new_detail::Storage<Block>;
    index_type block_rows_ = 0, block_cols_ = 0;
    index_type rows_in_block_ = 0, cols_in_block_ = 0;
    payload_type payload_;

    template <typename Factory>
    static std::vector<Block> generate(index_type rows, index_type cols, Factory &&factory) {
        const auto count = toeplitz_new_detail::unique_count(rows, cols);
        std::vector<Block> blocks;
        blocks.reserve(static_cast<std::size_t>(count));
        for (index_type j = 0; j < cols; ++j)
            blocks.emplace_back(std::invoke(factory, 0, j));
        for (index_type i = 1; i < rows; ++i)
            blocks.emplace_back(std::invoke(factory, i, 0));
        return blocks;
    }

    // Only a parent storage can rebind a subtree during construction.
    [[nodiscard]] index_type coefficient_count() const noexcept { return payload_.coefficient_count(); }
    [[nodiscard]] std::span<const scalar_type> coefficients() const requires(has_packed_coefficients) {
        return payload_.coefficients();
    }
    void repack(const std::shared_ptr<toeplitz_new_detail::Buffer<scalar_type>> &arena, index_type &offset)
        requires(has_packed_coefficients) {
        payload_.repack(arena, offset);
    }

  public:
    ToeplitzStorage() = default;

    ToeplitzStorage(index_type block_rows, index_type block_cols, std::vector<Block> blocks)
        : block_rows_(block_rows), block_cols_(block_cols) {
        const auto count = toeplitz_new_detail::unique_count(block_rows, block_cols);
#ifndef NDEBUG
        if (std::cmp_not_equal(blocks.size(), count))
            throw std::invalid_argument("Toeplitz needs block_rows + block_cols - 1 blocks");
#endif
        if (count == 0)
            return;
        rows_in_block_ = toeplitz_new_detail::positive_extent(blocks.front().rows());
        cols_in_block_ = toeplitz_new_detail::positive_extent(blocks.front().cols());
#ifndef NDEBUG
        for (const auto &block : blocks)
            if (std::cmp_not_equal(block.rows(), rows_in_block_) ||
                std::cmp_not_equal(block.cols(), cols_in_block_))
                throw std::invalid_argument("Toeplitz blocks must have identical dimensions");
        (void)toeplitz_new_detail::checked_product(block_rows_, rows_in_block_);
        (void)toeplitz_new_detail::checked_product(block_cols_, cols_in_block_);
#endif
        payload_ = payload_type(std::move(blocks));
    }

    template <typename Factory>
        requires std::invocable<Factory &, index_type, index_type> &&
                 std::convertible_to<std::invoke_result_t<Factory &, index_type, index_type>, Block>
    ToeplitzStorage(index_type block_rows, index_type block_cols, Factory &&factory)
        : ToeplitzStorage(block_rows, block_cols, generate(block_rows, block_cols, std::forward<Factory>(factory))) {}

    ToeplitzStorage(const ToeplitzStorage &) = default;
    ToeplitzStorage(ToeplitzStorage &&other) noexcept
        : block_rows_(std::exchange(other.block_rows_, 0)), block_cols_(std::exchange(other.block_cols_, 0)),
          rows_in_block_(std::exchange(other.rows_in_block_, 0)),
          cols_in_block_(std::exchange(other.cols_in_block_, 0)), payload_(std::move(other.payload_)) {}

    ToeplitzStorage &operator=(ToeplitzStorage other) noexcept {
        swap(other);
        return *this;
    }

    void swap(ToeplitzStorage &other) noexcept {
        using std::swap;
        swap(block_rows_, other.block_rows_);
        swap(block_cols_, other.block_cols_);
        swap(rows_in_block_, other.rows_in_block_);
        swap(cols_in_block_, other.cols_in_block_);
        swap(payload_, other.payload_);
    }

    [[nodiscard]] index_type rows() const noexcept { return block_rows_ * rows_in_block_; }
    [[nodiscard]] index_type cols() const noexcept { return block_cols_ * cols_in_block_; }
    [[nodiscard]] index_type block_rows() const noexcept { return block_rows_; }
    [[nodiscard]] index_type block_cols() const noexcept { return block_cols_; }
    [[nodiscard]] index_type rows_in_block() const noexcept { return rows_in_block_; }
    [[nodiscard]] index_type cols_in_block() const noexcept { return cols_in_block_; }
    [[nodiscard]] index_type unique_blocks() const noexcept {
        return block_rows_ == 0 ? 0 : block_rows_ + (block_cols_ - 1);
    }

    [[nodiscard]] index_type get_linear_index(index_type row, index_type col) const {
#ifndef NDEBUG
        if (row < 0 || row >= block_rows_ || col < 0 || col >= block_cols_)
            throw std::out_of_range("Toeplitz block index out of range");
#endif
        return col >= row ? col - row : block_cols_ + (row - col - 1);
    }
    [[nodiscard]] decltype(auto) get_block(index_type row, index_type col) const {
        return payload_.get(get_linear_index(row, col));
    }
    [[nodiscard]] decltype(auto) unique_block(index_type i) const {
#ifndef NDEBUG
        if (i < 0 || i >= unique_blocks())
            throw std::out_of_range("Toeplitz stored block index out of range");
#endif
        return payload_.get(i);
    }
    [[nodiscard]] std::span<const scalar_type> packed_coefficients() const requires(has_packed_coefficients) {
        return payload_.coefficients();
    }
};

/** Caller-owned scratch: one instance per simultaneous multiplication.
 * Buffers grow as needed and never shrink during multiplication. They are
 * reused across unique blocks AND across calls; child scratch is shared by
 * sequential applications of all blocks on the next level. Eigen/BLAS may
 * still allocate its own GEMM packing buffers. Do not use these buffers as RHS
 * or destination. No mutable scratch is stored inside the matrix.
 */
template <typename Scalar, typename ChildWorkspace> struct ToeplitzWorkspace {
    toeplitz_new_detail::Buffer<Scalar> rhs_panel;
    toeplitz_new_detail::Buffer<Scalar> product_panel;
    toeplitz_new_detail::Buffer<Scalar> prepared_alias_copy;
    ToeplitzDenseMatrix<Scalar> alias_copy;
    [[no_unique_address]] ChildWorkspace child{};
};

/** Persistent gathered representation of one vector for a square, uniform
 * nested Toeplitz layout. Every logical leaf segment is duplicated into the
 * panel of every Toeplitz leaf block that consumes it. That redundancy removes
 * RHS gathering from repeated matvecs. A prepared result has the same layout,
 * so two instances can be alternated in a chain x <- A*x without converting
 * through the ordinary vector layout.
 *
 * Instances are created by ToeplitzTemplateNew::prepare_vector() or
 * make_prepared_vector(). The numerical storage is intentionally opaque: its
 * layout belongs to the matrix which created it.
 */
template <typename Scalar> class ToeplitzPreparedVector {
    template <ToeplitzMatrixBlock> friend class ToeplitzTemplateNew;

    toeplitz_new_detail::Buffer<Scalar> values_;
    std::vector<toeplitz_new_detail::index> layout_signature_;
    toeplitz_new_detail::index logical_size_ = 0;

  public:
    using scalar_type = Scalar;
    using index_type = Eigen::Index;

    ToeplitzPreparedVector() = default;

    [[nodiscard]] index_type size() const noexcept { return logical_size_; }
    [[nodiscard]] std::size_t storage_size() const noexcept { return values_.size(); }
    [[nodiscard]] const scalar_type *data() const noexcept { return values_.data(); }
};

/** Block Toeplitz matrix, with scalar type inferred from Block.
 * Arithmetic facade composed with a ToeplitzStorage<Block> RAII owner. All
 * shape/ownership/copy/move logic belongs to that owner. The existing factory
 * and vector constructors delegate to it; an independently constructed owner
 * can also be passed by value. get_storage() provides read-only inspection.
 *
 * Stored order matches ToeplitzContainer: B(0,0), B(0,1), ..., B(0,C-1),
 * B(1,0), ..., B(R-1,0). No symmetry, reciprocity or square blocks assumed.
 * For a nonempty grid, a factory is invoked exactly R+C-1 times in that order.
 *
 * For one diagonal starting at (i0,j0) with length L, form
 *   X = [x_j0 ... x_(j0+L-1)], then add B*X to [y_i0 ... y_(i0+L-1)].
 * For a vector these panels are zero-copy Maps. For multiple RHS, column
 * q*L+t of the packed panel is RHS.block((j0+t)*cols_in_block,q,cols_in_block,1).
 * Each unique block is applied ONCE per call, including at nested levels.
 * This is GEMM batching, not an FFT or a reduction in arithmetic complexity:
 * a single level still costs O(R*C*m*n*Q) for m-by-n dense blocks and Q RHS.
 * The largest gathered panels use (m+n)*min(R,C)*Q scalars per level; Q can
 * grow by min(R,C) at each recursive level. There is no implicit tiling.
 * Uniform nested trees over dense Eigen leaves use a flattened hot path: the
 * hierarchy is reduced at construction to linear arrays containing every
 * compressed gather/scatter run, its source/destination offsets, panel offset,
 * and occurrence count. Dense leaves are addressed by one linear index. Thus
 * matvec performs no hierarchy traversal or mixed-radix index conversion, and
 * a leaf value is gathered/scattered once rather than once per nesting level.
 * Heterogeneous trees retain the recursive fallback.
 *
 * Dense leaves of a nested tree share one contiguous immutable coefficient
 * arena, in depth-first stored-block order. Copies share immutable numerical
 * data; metadata has value semantics. Construction consolidates factory
 * results, so temporary allocations/copies are possible during construction.
 *
 * For a square uniform tree, prepare_vector() can build a persistent gathered
 * RHS. matvec_prepared() consumes and produces that same redundant layout. Its
 * GEMM inputs are contiguous and need no per-call gather. The result reduction
 * and duplication into the panels required by the next multiplication remain
 * O(number of logical block interactions): no permutation alone can avoid
 * them, because different Toeplitz diagonals overlap in their sources and
 * destinations. Prepared storage contains one leaf segment per logical block
 * interaction, versus one segment per logical column in an ordinary vector.
 * Generic user-defined blocks have contiguous descriptors only: no allocator
 * can force an opaque block's internally allocated payload to be contiguous.
 *
 * get_block() returns a const reference for generic/nested blocks and a const
 * Eigen Map BY VALUE for dense blocks. Views must not outlive their owner.
 * There is intentionally no mutable block accessor that could invalidate
 * dimensions or break shared coefficient ownership.
 *
 * Example (all declarations are in this namespace):
 *   using T1 = ToeplitzTemplateNew<Eigen::MatrixXcd>;
 *   using T2 = ToeplitzTemplateNew<T1>;
 *   using T3 = ToeplitzTemplateNew<T2>;
 *   T1 a(3, 4, [](Eigen::Index i, Eigen::Index j) -> Eigen::MatrixXcd {
 *       return Eigen::MatrixXcd::Constant(2, 2, {double(i+j+1), double(i-j)});
 *   });
 *   T2 b(2, 3, [&](auto, auto) { return a; });
 *   T2::workspace_type scratch;
 *   Eigen::VectorXcd x = Eigen::VectorXcd::Ones(b.cols()), y(b.rows());
 *   b.matvec(x, y, scratch);      // overwrite y
 *   b.matvec_add(x, y, scratch);  // y += b*x
 *
 * The implementation launches no threads and changes no global Eigen/BLAS
 * settings. Backend threading, if any, remains under the caller's control.
 */
template <ToeplitzMatrixBlock Block> class ToeplitzTemplateNew {
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
    using prepared_vector_type = ToeplitzPreparedVector<scalar_type>;
    using storage_type = ToeplitzStorage<Block>;
    static constexpr bool has_packed_coefficients = toeplitz_new_detail::PackedTree<Block>::value;

  private:
    using packed_matrices_type = toeplitz_new_detail::PackedArithmeticMatrices<Block>;
    using packed_leaf_matrix_t = typename packed_matrices_type::leaf_type;
    using packed_rhs_panel_t = typename packed_matrices_type::rhs_panel_type;
    using packed_product_panel_t = typename packed_matrices_type::product_panel_type;

    template <typename> friend class toeplitz_new_detail::Storage;
    storage_type storage_;
    toeplitz_new_detail::FlatToeplitzIndex flat_index_;

    template <typename Node>
    static bool build_flat_index(const Node &node, std::size_t depth,
                                 toeplitz_new_detail::FlatToeplitzIndex &index,
                                 const scalar_type *coefficient_base, index_type &leaf_number) {
        if constexpr (requires {
                          node.block_rows();
                          node.block_cols();
                          node.rows_in_block();
                          node.cols_in_block();
                          node.unique_blocks();
                          node.unique_block(index_type{});
                      }) {
            const toeplitz_new_detail::FlatToeplitzLevel level{
                node.block_rows(), node.block_cols(), node.rows_in_block(), node.cols_in_block(),
                node.unique_blocks()};
            if (index.levels.size() == depth) {
                index.levels.push_back(level);
            } else {
                const auto &expected = index.levels[depth];
                if (expected.block_rows != level.block_rows || expected.block_cols != level.block_cols ||
                    expected.rows_in_block != level.rows_in_block || expected.cols_in_block != level.cols_in_block ||
                    expected.unique_blocks != level.unique_blocks)
                    return false;
            }
            for (index_type block = 0; block < level.unique_blocks; ++block)
                if (!build_flat_index(node.unique_block(block), depth + 1, index, coefficient_base, leaf_number))
                    return false;
            return true;
        } else if constexpr (requires { node.rows(); node.cols(); node.data(); }) {
            const index_type rows = toeplitz_new_detail::positive_extent(node.rows());
            const index_type cols = toeplitz_new_detail::positive_extent(node.cols());
            if (leaf_number == 0) {
                index.leaf_rows = rows;
                index.leaf_cols = cols;
            } else if (index.leaf_rows != rows || index.leaf_cols != cols) {
                return false;
            }
            const index_type stride = toeplitz_new_detail::checked_product(rows, cols);
            if (node.data() != coefficient_base + toeplitz_new_detail::checked_product(leaf_number, stride))
                return false;
            ++leaf_number;
            return true;
        } else {
            return false;
        }
    }

    void initialize_flat_index() {
        flat_index_ = {};
        if constexpr (has_packed_coefficients) {
            const auto coefficients = storage_.packed_coefficients();
            if (coefficients.empty())
                return;
            index_type leaves = 0;
            if (!build_flat_index(*this, 0, flat_index_, coefficients.data(), leaves)) {
                flat_index_ = {};
                return;
            }
            flat_index_.leaf_blocks = leaves;
            // A single level already has a zero-copy vector path. Flattening is
            // useful when it removes at least one intermediate recursive panel.
            flat_index_.enabled = flat_index_.levels.size() > 1 && leaves > 0;
            flat_index_.recurrent_layout = leaves > 0 && flat_index_.leaf_rows == flat_index_.leaf_cols &&
                                           std::ranges::all_of(flat_index_.levels, [](const auto &level) {
                                               return level.block_rows == level.block_cols &&
                                                      level.rows_in_block == level.cols_in_block;
                                           });
            if (flat_index_.enabled || flat_index_.recurrent_layout) {
                toeplitz_new_detail::build_flat_runs(flat_index_);
                if (flat_index_.recurrent_layout) {
                    // Leaf zero is the main diagonal at every level. Its panel
                    // is the canonical vector in ordinary contiguous order.
                    const index_type begin = flat_index_.leaf_run_begin.front();
                    const index_type end = flat_index_.leaf_run_begin[1];
                    for (index_type run = begin; run < end; ++run) {
                        const auto position = static_cast<std::size_t>(run);
                        const index_type canonical = toeplitz_new_detail::checked_product(
                            flat_index_.run_panel_offset[position], flat_index_.leaf_cols);
                        if (flat_index_.run_source_offset[position] != canonical ||
                            flat_index_.run_destination_offset[position] != canonical) {
                            flat_index_.recurrent_layout = false;
                            break;
                        }
                    }
                    flat_index_.recurrent_layout =
                        flat_index_.recurrent_layout &&
                        flat_index_.leaf_panel_offset[1] == cols();
                }
            }
        }
    }

    [[nodiscard]] std::vector<index_type> prepared_layout_signature() const {
        std::vector<index_type> result;
        result.reserve(4 + 5 * flat_index_.levels.size());
        result.push_back(static_cast<index_type>(flat_index_.levels.size()));
        result.push_back(flat_index_.leaf_rows);
        result.push_back(flat_index_.leaf_cols);
        result.push_back(flat_index_.leaf_blocks);
        for (const auto &level : flat_index_.levels) {
            result.push_back(level.block_rows);
            result.push_back(level.block_cols);
            result.push_back(level.rows_in_block);
            result.push_back(level.cols_in_block);
            result.push_back(level.unique_blocks);
        }
        return result;
    }

#ifndef NDEBUG
    void check_prepared_vector(const prepared_vector_type &vector) const {
        if (!flat_index_.recurrent_layout)
            throw std::logic_error("Prepared Toeplitz vectors require a square uniform dense-leaf tree");
        if (vector.logical_size_ != cols() || vector.layout_signature_ != prepared_layout_signature() ||
            std::cmp_not_equal(vector.values_.size(), flat_index_.leaf_panel_offset.back()))
            throw std::invalid_argument("Prepared vector has an incompatible Toeplitz layout");
    }
#endif

    void repack_prepared_from_canonical(prepared_vector_type &vector) const {
        const scalar_type *canonical = vector.values_.data();
        // Leaf zero already is the canonical main-diagonal panel.
        for (index_type leaf = 1; leaf < flat_index_.leaf_blocks; ++leaf) {
            const auto leaf_position = static_cast<std::size_t>(leaf);
            scalar_type *panel = vector.values_.data() + flat_index_.leaf_panel_offset[leaf_position];
            const index_type run_begin = flat_index_.leaf_run_begin[leaf_position];
            const index_type run_end = flat_index_.leaf_run_begin[leaf_position + 1];
            const index_type run_length = flat_index_.leaf_run_length[leaf_position];
            for (index_type run = run_begin; run < run_end; ++run) {
                const auto run_position = static_cast<std::size_t>(run);
                Eigen::Map<const packed_rhs_panel_t> source(
                    canonical + flat_index_.run_source_offset[run_position],
                    flat_index_.leaf_cols, run_length);
                Eigen::Map<packed_rhs_panel_t> target(
                    panel + toeplitz_new_detail::checked_product(
                                flat_index_.run_panel_offset[run_position], flat_index_.leaf_cols),
                    flat_index_.leaf_cols, run_length);
                target = source;
            }
        }
    }

    void multiply_add_prepared_noalias(const scalar_type *rhs, prepared_vector_type &dest,
                                       workspace_type &workspace) const
        requires(has_packed_coefficients)
    {
        scalar_type *canonical_dest = dest.values_.data();
        const auto coefficients = storage_.packed_coefficients();
        const index_type leaf_stride =
            toeplitz_new_detail::checked_product(flat_index_.leaf_rows, flat_index_.leaf_cols);
        const scalar_type *coefficient = coefficients.data();

        for (index_type leaf = 0; leaf < flat_index_.leaf_blocks; ++leaf) {
            const auto leaf_position = static_cast<std::size_t>(leaf);
            const index_type occurrences = flat_index_.leaf_occurrences[leaf_position];
            const index_type run_length = flat_index_.leaf_run_length[leaf_position];
            const index_type run_begin = flat_index_.leaf_run_begin[leaf_position];
            const index_type run_end = flat_index_.leaf_run_begin[leaf_position + 1];

            Eigen::Map<const packed_leaf_matrix_t> block(
                coefficient, flat_index_.leaf_rows, flat_index_.leaf_cols);
            coefficient += leaf_stride;
            Eigen::Map<const packed_rhs_panel_t> x(
                rhs + flat_index_.leaf_panel_offset[leaf_position],
                flat_index_.leaf_cols, occurrences);
            if (run_end == run_begin + 1) {
                Eigen::Map<packed_product_panel_t> target(
                    canonical_dest + flat_index_.run_destination_offset[static_cast<std::size_t>(run_begin)],
                    flat_index_.leaf_rows, occurrences);
                target.noalias() += block * x;
                continue;
            }
            toeplitz_new_detail::grow(
                workspace.product_panel,
                toeplitz_new_detail::checked_product(flat_index_.leaf_rows, occurrences));
            Eigen::Map<packed_product_panel_t> y(
                workspace.product_panel.data(), flat_index_.leaf_rows, occurrences);
            y.noalias() = block * x;

            for (index_type run = run_begin; run < run_end; ++run) {
                const auto run_position = static_cast<std::size_t>(run);
                Eigen::Map<packed_product_panel_t> target(
                    canonical_dest + flat_index_.run_destination_offset[run_position],
                    flat_index_.leaf_rows, run_length);
                target += y.middleCols(flat_index_.run_panel_offset[run_position], run_length);
            }
        }
        repack_prepared_from_canonical(dest);
    }

    void apply_prepared(const prepared_vector_type &rhs, prepared_vector_type &dest,
                        workspace_type &workspace, bool add) const
        requires(has_packed_coefficients)
    {
#ifndef NDEBUG
        check_prepared_vector(rhs);
        check_prepared_vector(dest);
#endif
        const scalar_type *input = rhs.values_.data();
        if (&rhs == &dest) {
            workspace.prepared_alias_copy.assign(rhs.values_.begin(), rhs.values_.end());
            input = workspace.prepared_alias_copy.data();
        }
        if (!add)
            std::fill_n(dest.values_.data(), static_cast<std::size_t>(rows()), scalar_type{});
        multiply_add_prepared_noalias(input, dest, workspace);
    }

    void multiply_add_flat(const_matrix_ref rhs, matrix_ref dest, workspace_type &workspace) const
        requires(has_packed_coefficients)
    {
        const auto coefficients = storage_.packed_coefficients();
        const index_type leaf_stride =
            toeplitz_new_detail::checked_product(flat_index_.leaf_rows, flat_index_.leaf_cols);
        const index_type q = rhs.cols();
        const scalar_type *coefficient = coefficients.data();

        for (index_type leaf = 0; leaf < flat_index_.leaf_blocks; ++leaf) {
            const auto leaf_position = static_cast<std::size_t>(leaf);
            const index_type occurrences = flat_index_.leaf_occurrences[leaf_position];
            const index_type run_length = flat_index_.leaf_run_length[leaf_position];
            const index_type run_begin = flat_index_.leaf_run_begin[leaf_position];
            const index_type run_end = flat_index_.leaf_run_begin[leaf_position + 1];

            Eigen::Map<const packed_leaf_matrix_t> block(
                coefficient, flat_index_.leaf_rows, flat_index_.leaf_cols);
            coefficient += leaf_stride;
            if (q == 1 && run_end == run_begin + 1) {
                Eigen::Map<const packed_rhs_panel_t> x(
                    rhs.data() + flat_index_.run_source_offset[static_cast<std::size_t>(run_begin)],
                    flat_index_.leaf_cols, occurrences);
                Eigen::Map<packed_product_panel_t> y(
                    dest.data() + flat_index_.run_destination_offset[static_cast<std::size_t>(run_begin)],
                    flat_index_.leaf_rows, occurrences);
                y.noalias() += block * x;
                continue;
            }

            const index_type width = toeplitz_new_detail::checked_product(occurrences, q);
            toeplitz_new_detail::grow(
                workspace.rhs_panel, toeplitz_new_detail::checked_product(flat_index_.leaf_cols, width));
            toeplitz_new_detail::grow(
                workspace.product_panel, toeplitz_new_detail::checked_product(flat_index_.leaf_rows, width));
            Eigen::Map<packed_rhs_panel_t> x(
                workspace.rhs_panel.data(), flat_index_.leaf_cols, width);
            Eigen::Map<packed_product_panel_t> y(
                workspace.product_panel.data(), flat_index_.leaf_rows, width);
            for (index_type column = 0; column < q; ++column)
                for (index_type run = run_begin; run < run_end; ++run) {
                    const index_type panel_column =
                        flat_index_.run_panel_offset[static_cast<std::size_t>(run)];
                    Eigen::Map<const packed_rhs_panel_t> source(
                        rhs.data() + column * rhs.outerStride() +
                            flat_index_.run_source_offset[static_cast<std::size_t>(run)],
                        flat_index_.leaf_cols, run_length);
                    x.middleCols(column * occurrences + panel_column, run_length) = source;
                }
            y.noalias() = block * x;
            for (index_type column = 0; column < q; ++column)
                for (index_type run = run_begin; run < run_end; ++run) {
                    const index_type panel_column =
                        flat_index_.run_panel_offset[static_cast<std::size_t>(run)];
                    Eigen::Map<packed_product_panel_t> target(
                        dest.data() + column * dest.outerStride() +
                            flat_index_.run_destination_offset[static_cast<std::size_t>(run)],
                        flat_index_.leaf_rows, run_length);
                    target += y.middleCols(column * occurrences + panel_column, run_length);
                }
        }
    }

#ifndef NDEBUG
    void check_product_shapes(const_matrix_ref rhs, matrix_ref dest) const {
        if (rhs.rows() != cols() || dest.rows() != rows() || rhs.cols() != dest.cols())
            throw std::invalid_argument("Toeplitz multiplication dimension mismatch");
        if (rhs.outerStride() < 0 || dest.outerStride() < 0)
            throw std::invalid_argument("Toeplitz multiplication requires nonnegative column strides");
    }
#endif

    void apply(const_matrix_ref rhs, matrix_ref dest, workspace_type &workspace, bool add) const {
#ifndef NDEBUG
        check_product_shapes(rhs, dest);
#endif
        // Copy BEFORE clearing destination, including partial overlap of views.
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
    ToeplitzTemplateNew() = default; // canonical empty 0-by-0 matrix

    explicit ToeplitzTemplateNew(storage_type storage) : storage_(std::move(storage)) { initialize_flat_index(); }

    ToeplitzTemplateNew(index_type block_rows, index_type block_cols, std::vector<Block> blocks)
        : storage_(block_rows, block_cols, std::move(blocks)) {
        initialize_flat_index();
    }

    template <typename Factory>
        requires std::invocable<Factory &, index_type, index_type> &&
                 std::convertible_to<std::invoke_result_t<Factory &, index_type, index_type>, Block>
    ToeplitzTemplateNew(index_type block_rows, index_type block_cols, Factory &&factory)
        : storage_(block_rows, block_cols, std::forward<Factory>(factory)) {
        initialize_flat_index();
    }

    ToeplitzTemplateNew(const ToeplitzTemplateNew &) = default;
    ToeplitzTemplateNew(ToeplitzTemplateNew &&other) noexcept
        : storage_(std::move(other.storage_)), flat_index_(std::move(other.flat_index_)) {
        other.flat_index_.reset();
    }

    ToeplitzTemplateNew &operator=(ToeplitzTemplateNew other) noexcept {
        swap(other);
        return *this;
    }

    void swap(ToeplitzTemplateNew &other) noexcept {
        storage_.swap(other.storage_);
        flat_index_.swap(other.flat_index_);
    }

    [[nodiscard]] const storage_type &get_storage() const noexcept { return storage_; }
    [[nodiscard]] index_type rows() const noexcept { return storage_.rows(); }
    [[nodiscard]] index_type cols() const noexcept { return storage_.cols(); }
    [[nodiscard]] index_type block_rows() const noexcept { return storage_.block_rows(); }
    [[nodiscard]] index_type block_cols() const noexcept { return storage_.block_cols(); }
    [[nodiscard]] index_type rows_in_block() const noexcept { return storage_.rows_in_block(); }
    [[nodiscard]] index_type cols_in_block() const noexcept { return storage_.cols_in_block(); }
    [[nodiscard]] index_type unique_blocks() const noexcept { return storage_.unique_blocks(); }
    [[nodiscard]] bool uses_flat_indexing() const noexcept { return flat_index_.enabled; }
    [[nodiscard]] std::size_t flattened_levels() const noexcept { return flat_index_.levels.size(); }
    [[nodiscard]] index_type flattened_leaf_blocks() const noexcept { return flat_index_.leaf_blocks; }
    [[nodiscard]] bool supports_prepared_matvec() const noexcept {
        return flat_index_.recurrent_layout;
    }

    [[nodiscard]] index_type get_linear_index(index_type row, index_type col) const {
        return storage_.get_linear_index(row, col);
    }
    [[nodiscard]] decltype(auto) get_block(index_type row, index_type col) const {
        return storage_.get_block(row, col);
    }
    [[nodiscard]] decltype(auto) unique_block(index_type i) const {
        return storage_.unique_block(i);
    }
    [[nodiscard]] std::span<const scalar_type> packed_coefficients() const requires(has_packed_coefficients) {
        return storage_.packed_coefficients();
    }

    [[nodiscard]] scalar_type operator()(index_type row, index_type col) const {
#ifndef NDEBUG
        if (row < 0 || row >= rows() || col < 0 || col >= cols())
            throw std::out_of_range("Toeplitz coefficient index out of range");
#endif
        return get_block(row / rows_in_block(), col / cols_in_block())(row % rows_in_block(), col % cols_in_block());
    }

    [[nodiscard]] dense_matrix_t to_dense() const {
        (void)toeplitz_new_detail::checked_product(rows(), cols());
        dense_matrix_t result(rows(), cols());
        for (index_type j = 0; j < cols(); ++j)
            for (index_type i = 0; i < rows(); ++i)
                result(i, j) = (*this)(i, j);
        return result;
    }

    [[nodiscard]] vector_t diagonal() const {
        vector_t result(std::min(rows(), cols()));
        for (index_type i = 0; i < result.size(); ++i)
            result(i) = (*this)(i, i);
        return result;
    }

    /** Low-level static block interface: dest += this*rhs.
     * RHS/destination/workspace MUST NOT overlap. Public matmul/matvec wrappers
     * below handle RHS/destination aliasing. Shapes are checked in all builds.
     */
    void multiply_add_noalias(const_matrix_ref rhs, matrix_ref dest, workspace_type &workspace) const {
#ifndef NDEBUG
        check_product_shapes(rhs, dest);
#endif
        if (unique_blocks() == 0 || rhs.cols() == 0)
            return;
        if constexpr (has_packed_coefficients) {
            if (flat_index_.enabled) {
                multiply_add_flat(rhs, dest, workspace);
                return;
            }
        }
        const index_type block_rows_ = storage_.block_rows(), block_cols_ = storage_.block_cols();
        const index_type rows_in_block_ = storage_.rows_in_block(), cols_in_block_ = storage_.cols_in_block();
        const index_type q = rhs.cols();
        if (q > 1) {
            const auto width = toeplitz_new_detail::checked_product(std::min(block_rows_, block_cols_), q);
            toeplitz_new_detail::grow(workspace.rhs_panel,
                                     toeplitz_new_detail::checked_product(cols_in_block_, width));
            toeplitz_new_detail::grow(workspace.product_panel,
                                     toeplitz_new_detail::checked_product(rows_in_block_, width));
        }

        for (index_type d = 0; d < unique_blocks(); ++d) {
            const index_type i0 = d < block_cols_ ? 0 : d - block_cols_ + 1;
            const index_type j0 = d < block_cols_ ? d : 0;
            const index_type repeats = std::min(block_rows_ - i0, block_cols_ - j0);
            const auto &block = storage_.unique_block(d);
            if (q == 1) {
                // Consecutive vector segments are already a column-major panel.
                Eigen::Map<const dense_matrix_t> x(rhs.data() + j0 * cols_in_block_, cols_in_block_, repeats);
                Eigen::Map<dense_matrix_t> y(dest.data() + i0 * rows_in_block_, rows_in_block_, repeats);
                traits_type::multiply_add(block, x, y, workspace.child);
            } else {
                const index_type width = repeats * q; // checked above using maximum repeats
                Eigen::Map<dense_matrix_t> x(workspace.rhs_panel.data(), cols_in_block_, width);
                Eigen::Map<dense_matrix_t> y(workspace.product_panel.data(), rows_in_block_, width);
                for (index_type column = 0; column < q; ++column) {
                    // One contiguous run per RHS, including for strided outer views.
                    Eigen::Map<const dense_matrix_t> source(
                        rhs.data() + column * rhs.outerStride() + j0 * cols_in_block_, cols_in_block_, repeats);
                    x.middleCols(column * repeats, repeats) = source;
                }
                y.setZero();
                traits_type::multiply_add(block, x, y, workspace.child);
                for (index_type column = 0; column < q; ++column) {
                    Eigen::Map<dense_matrix_t> target(
                        dest.data() + column * dest.outerStride() + i0 * rows_in_block_, rows_in_block_, repeats);
                    target += y.middleCols(column * repeats, repeats);
                }
            }
        }
    }

    // Overwrite and additive APIs. No output resizing; dimensions must match.
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
        matmul(rhs, result);
        return result;
    }

    void matvec(const_vector_ref rhs, vector_ref dest, workspace_type &workspace) const {
        matmul(rhs, dest, workspace);
    }
    void matvec_add(const_vector_ref rhs, vector_ref dest, workspace_type &workspace) const {
        matmul_add(rhs, dest, workspace);
    }
    void matvec(const_vector_ref rhs, vector_ref dest) const { matmul(rhs, dest); }
    void matvec_add(const_vector_ref rhs, vector_ref dest) const { matmul_add(rhs, dest); }
    [[nodiscard]] vector_t matvec(const_vector_ref rhs) const {
        vector_t result(rows());
        matvec(rhs, result);
        return result;
    }
    [[nodiscard]] vector_t operator*(const_vector_ref rhs) const { return matvec(rhs); }

    /** Allocate a zero prepared vector compatible with this matrix. The
     * precondition is supports_prepared_matvec()==true.
     */
    [[nodiscard]] prepared_vector_type make_prepared_vector() const
        requires(has_packed_coefficients)
    {
#ifndef NDEBUG
        if (!flat_index_.recurrent_layout)
            throw std::logic_error("Prepared Toeplitz vectors require a square uniform dense-leaf tree");
#endif
        prepared_vector_type result;
        result.logical_size_ = cols();
        result.layout_signature_ = prepared_layout_signature();
        result.values_.resize(
            static_cast<std::size_t>(flat_index_.leaf_panel_offset.back()), scalar_type{});
        return result;
    }

    /** Gather an ordinary vector once into all persistent leaf panels. */
    [[nodiscard]] prepared_vector_type prepare_vector(const_vector_ref rhs) const
        requires(has_packed_coefficients)
    {
#ifndef NDEBUG
        if (rhs.size() != cols())
            throw std::invalid_argument("Toeplitz prepared-vector dimension mismatch");
#endif
        prepared_vector_type result = make_prepared_vector();
        std::copy_n(rhs.data(), static_cast<std::size_t>(rhs.size()), result.values_.data());
        repack_prepared_from_canonical(result);
        return result;
    }

    /** Recover the ordinary contiguous vector. This conversion is not needed
     * between prepared matvecs; it is intended for consumers using the usual
     * vector API and for validation/output.
     */
    [[nodiscard]] vector_t unpack_vector(const prepared_vector_type &vector) const
        requires(has_packed_coefficients)
    {
#ifndef NDEBUG
        check_prepared_vector(vector);
#endif
        return Eigen::Map<const vector_t>(vector.values_.data(), vector.logical_size_);
    }

    /** dest = this*rhs, preserving the prepared layout for another call. */
    void matvec_prepared(const prepared_vector_type &rhs, prepared_vector_type &dest,
                         workspace_type &workspace) const
        requires(has_packed_coefficients)
    {
        apply_prepared(rhs, dest, workspace, false);
    }

    /** dest += this*rhs, preserving the prepared layout for another call. */
    void matvec_prepared_add(const prepared_vector_type &rhs, prepared_vector_type &dest,
                             workspace_type &workspace) const
        requires(has_packed_coefficients)
    {
        apply_prepared(rhs, dest, workspace, true);
    }

    void matvec_prepared(const prepared_vector_type &rhs, prepared_vector_type &dest) const
        requires(has_packed_coefficients)
    {
        workspace_type workspace;
        matvec_prepared(rhs, dest, workspace);
    }

    void matvec_prepared_add(const prepared_vector_type &rhs, prepared_vector_type &dest) const
        requires(has_packed_coefficients)
    {
        workspace_type workspace;
        matvec_prepared_add(rhs, dest, workspace);
    }
};

} // namespace EMW::Math::LinAgl::Matrix

#endif // EMW_MATH_MATRIX_TEOPLITZ_TEMPLATE_NEW_HPP
