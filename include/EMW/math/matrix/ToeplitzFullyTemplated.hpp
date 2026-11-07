//
// Created by evgen on 12.02.2025.
//

#ifndef TOEPLITZFULLYTEMPLATED_HPP
#define TOEPLITZFULLYTEMPLATED_HPP

#include "ToeplitzContainer.hpp"
#include "types/Types.hpp"

#include <cassert>
#include <algorithm>
#include <exception>
#include <functional>
#include <iostream>
#include <limits>
#include <vector>

namespace EMW::Math::LinAgl::Matrix {
/**
 * Матрица со структурой тёплиц-общий_вид
 */
template <typename scalar_t, typename block_t> class ToeplitzStructure {
    // Проверка на дефолт контруирование
    static_assert(std::is_default_constructible_v<block_t>);

  public:
    using vector_t = Types::VectorX<scalar_t>;
    using dense_matrix_t = Types::MatrixX<scalar_t>;
    using block_type = block_t;
    using scalar_type = scalar_t;

  private:
    struct MatvecWorkspace {
        std::vector<vector_t> rhs_panels;
        std::vector<vector_t> product_panels;
        vector_t input_alias_copy;
        vector_t factor_first;
        vector_t factor_second;
    };

    ToeplitzContainer<block_t> blocks;
    // Характеристика одного блока
    // Эти поля нужны для удобства
    Types::index rows_in_block_ = 0;
    Types::index cols_in_block_ = 0;
    // matvec_wise_block() is allocation-free after construction. The scratch
    // storage is mutable because multiplication does not change matrix values.
    // Concurrent calls on the same matrix instance require separate matrices
    // (or an external synchronization point).
    mutable MatvecWorkspace matvec_workspace_;

    void initialize_matvec_workspace();

  public:
    ToeplitzStructure() = default;

    /**
     * Конструктор сейчас работает при условии, что в типе block_t есть методы rows() и cols()
     *
     * @param block_rows -- количество строк матрицы, считая в блоках
     * @param block_cols -- количество столбцов матрицы, считая в блоках
     * @param get_block -- функция, которая возвращает квадратную матрицу, размеры одинаковы для любых пар (i, j)
     *
     * Консистентность состояния поддерживается, если функция get_block возвращает блоки
     * одного и того же размера
     */
    ToeplitzStructure(Types::index block_rows, Types::index block_cols,
                      const std::function<block_t(Types::index i, Types::index j)> &get_block);

    ToeplitzStructure(Types::index block_rows, Types::index block_cols, Containers::vector<block_t> &&blocks_);

    /*
     * Построить нулевую матрицу с заданном формате.
     *
     * Это почти как дефолтный конструктор, только не совсем
     *
     */
    ToeplitzStructure(Types::index block_rows, Types::index block_cols, Types::index rows_in_block,
                      Types::index cols_in_block);

    /** Работа с индексами
     */
    [[nodiscard]] inline decltype(auto) get_linear_index(size_t row, size_t col) const noexcept {
        return blocks.get_linear_index(row, col);
    }

    /** Умножение матрицы на вектор */
    [[nodiscard]] vector_t matvec(const vector_t &vec) const noexcept;
    /**
     * dest += matvec(vec);
     */
    void matvec(Eigen::Ref<vector_t> vec, Eigen::Ref<vector_t> dest) const noexcept;

    void matvec(scalar_t *vec, size_t vec_size, Eigen::Ref<vector_t> dest) const noexcept;

    void matvec(scalar_t *vec, size_t vec_size, scalar_t *dest, size_t dest_size) const noexcept;
    // Apply one matrix-matrix operation per unique Toeplitz block.
    void matvec_wise_block(scalar_t *vec, size_t vec_size, scalar_t *dest, size_t dest_size) const noexcept;
    void matvec_wise(scalar_t *vec, size_t vec_size, scalar_t *dest, size_t dest_size) const noexcept;

    /** Умножение матрицы на число с возвращением копии */
    [[nodiscard]] ToeplitzStructure mull(scalar_t value) const noexcept;
    /** Умножение себя на число */
    void mull_inplace(scalar_t value) noexcept;
    /** Взятие диагонали */
    [[nodiscard]] vector_t diagonal() const;

    // --- Selectors --- //
    /**
     * Вернуть блок в матрице под передаваемой нумерацией
     * @param row номер строки (в контексте блочной структуры)
     * @param col номер столбца (в контексте блочной структуры)
     * @return const ref на соотвествующий блок в тёплицевом контейнере
     */
    [[nodiscard]] const block_t &get_block(Types::index row, Types::index col) const noexcept {
        return blocks(row, col);
    }

    [[nodiscard]] block_t &get_block(Types::index row, Types::index col) noexcept { return blocks(row, col); }

    // Возвращают значения строк и столбцов в каждом блоке
    [[nodiscard]] Types::index rows_in_block() const noexcept { return rows_in_block_; }
    [[nodiscard]] Types::index cols_in_block() const noexcept { return cols_in_block_; }
    // Возвращают количество строк и столбцов во всей матрице (то есть в большой матрице, которая тут удобно хранится)
    // имеется ввиду outerSize в терминологии Eigen
    [[nodiscard]] Types::index rows() const noexcept { return rows_in_block_ * blocks.rows(); }
    [[nodiscard]] Types::index cols() const noexcept { return cols_in_block_ * blocks.cols(); }

    [[nodiscard]] Types::index rows_in_toeplitrz() const { return blocks.rows(); };
    [[nodiscard]] Types::index cols_in_toeplitrz() const { return blocks.cols(); };

    // --- Доступ к элементам на чтение --- //

    /** Доступ к настоящим элементам в матрице, а не к блокам! */
    [[nodiscard]] const scalar_t &operator()(Types::index i, Types::index j) const noexcept;
    /** Доступ к блокам в матрице */
    const ToeplitzContainer<block_t> &get_blocks() const noexcept { return blocks; }

    // ---- Static methods --- //
    inline static Types::index get_size_of_container(Types::index rows, Types::index cols) noexcept
        __attribute__((always_inline)) {
        return rows + cols - 1;
    };

    // Приведение к плотной матрице
    Types::MatrixX<scalar_t> to_dense() const noexcept;
};

template <typename scalar_t, typename block_t>
typename ToeplitzStructure<scalar_t, block_t>::vector_t ToeplitzStructure<scalar_t, block_t>::diagonal() const {
    if (rows_in_block_ != cols_in_block_) {
        throw std::runtime_error(
            "To correct use of diagonal method internal blocks of toeplitz structure must be square");
    }

    const Types::index diagonal_size = std::min(rows_in_block_, cols_in_block_);
    const Types::index how_many_diagonal_blocks = std::min(blocks.rows(), blocks.cols());
    const vector_t subdiag = blocks(0, 0).diagonal();

    vector_t result = vector_t::Zero(diagonal_size * how_many_diagonal_blocks);
    for (Types::index i = 0; i < how_many_diagonal_blocks; ++i) {
        result.block(i * diagonal_size, 0, diagonal_size, 1) = subdiag;
    }
    return result;
}

template <typename scalar_t, typename block_t>
[[nodiscard]] const scalar_t &ToeplitzStructure<scalar_t, block_t>::operator()(Types::index i,
                                                                               Types::index j) const noexcept {
    // сначала поймем в каком из "верхних блоков" мы находимся
    const Types::index row_on_current_level = i / rows_in_block_;
    const Types::index col_on_current_level = j / cols_in_block_;
    // далее смотрим индексы на следующем уровне вложенности
    const Types::index i_new = i % rows_in_block_;
    const Types::index j_new = j % cols_in_block_;
    // и идём стучаться в него
    return blocks(row_on_current_level, col_on_current_level)(i_new, j_new);
}

template <typename scalar_t, typename block_t>
ToeplitzStructure<scalar_t, block_t>::ToeplitzStructure(
    Types::index block_rows, Types::index block_cols,
    const std::function<block_t(Types::index i, Types::index j)> &get_block)
    : blocks(block_rows, block_cols, get_block), rows_in_block_(blocks(0, 0).rows()),
      cols_in_block_(blocks(0, 0).cols()) {
    initialize_matvec_workspace();
}

template <typename scalar_t, typename block_t>
ToeplitzStructure<scalar_t, block_t>::ToeplitzStructure(Types::index block_rows, Types::index block_cols,
                                                        Containers::vector<block_t> &&blocks_)
    : blocks(block_rows, block_cols, std::move(blocks_)), rows_in_block_(blocks(0, 0).rows()),
      cols_in_block_(blocks(0, 0).cols()) {
    initialize_matvec_workspace();
}

template <typename scalar_t, typename block_t>
void ToeplitzStructure<scalar_t, block_t>::initialize_matvec_workspace() {
    const auto checked_product = [](size_t lhs, size_t rhs) -> size_t {
        if (rhs != 0 && lhs > std::numeric_limits<size_t>::max() / rhs)
            std::terminate();
        return lhs * rhs;
    };

    std::vector<size_t> rhs_sizes;
    std::vector<size_t> product_sizes;
    size_t factor_workspace_size = 0;
    const auto collect = [&](auto &&self, const auto &matrix, size_t rhs_columns, size_t depth) -> void {
        if constexpr (requires {
                          matrix.rows_in_toeplitrz();
                          matrix.cols_in_toeplitrz();
                          matrix.rows_in_block();
                          matrix.cols_in_block();
                          matrix.get_blocks();
                      }) {
            if (rhs_sizes.size() <= depth) {
                rhs_sizes.resize(depth + 1, 0);
                product_sizes.resize(depth + 1, 0);
            }
            const size_t block_rows = static_cast<size_t>(matrix.rows_in_toeplitrz());
            const size_t block_cols = static_cast<size_t>(matrix.cols_in_toeplitrz());
            const size_t repetitions_max = std::min(block_rows, block_cols);
            if (rhs_columns > 1) {
                const size_t panel_columns = checked_product(repetitions_max, rhs_columns);
                rhs_sizes[depth] = std::max(rhs_sizes[depth], checked_product(
                    static_cast<size_t>(matrix.cols_in_block()), panel_columns));
                product_sizes[depth] = std::max(product_sizes[depth], checked_product(
                    static_cast<size_t>(matrix.rows_in_block()), panel_columns));
            }

            const size_t count = static_cast<size_t>(matrix.get_blocks().get_actual_size());
            for (size_t diagonal = 0; diagonal < count; ++diagonal) {
                const size_t i0 = diagonal < block_cols ? 0 : diagonal - block_cols + 1;
                const size_t j0 = diagonal < block_cols ? diagonal : 0;
                const size_t repetitions = std::min(block_rows - i0, block_cols - j0);
                self(self, matrix.get_blocks()(static_cast<Types::index>(diagonal)),
                     checked_product(repetitions, rhs_columns), depth + 1);
            }
        } else if constexpr (requires { matrix.matmul_wise_workspace_size(rhs_columns); }) {
            factor_workspace_size = std::max(
                factor_workspace_size, matrix.matmul_wise_workspace_size(rhs_columns));
        }
    };
    collect(collect, *this, 1, 0);

    matvec_workspace_.rhs_panels.resize(rhs_sizes.size());
    matvec_workspace_.product_panels.resize(product_sizes.size());
    for (size_t depth = 0; depth < rhs_sizes.size(); ++depth) {
        matvec_workspace_.rhs_panels[depth].resize(static_cast<Eigen::Index>(rhs_sizes[depth]));
        matvec_workspace_.product_panels[depth].resize(static_cast<Eigen::Index>(product_sizes[depth]));
    }
    matvec_workspace_.input_alias_copy.resize(cols());
    matvec_workspace_.factor_first.resize(static_cast<Eigen::Index>(factor_workspace_size));
    matvec_workspace_.factor_second.resize(static_cast<Eigen::Index>(factor_workspace_size));
}

template <typename scalar_t, typename block_t>
typename ToeplitzStructure<scalar_t, block_t>::vector_t
ToeplitzStructure<scalar_t, block_t>::matvec(const vector_t &vec) const noexcept {
    assert(vec.size() == cols());
    // создаем нулевой вектор результата, в который будем записывать ответ
    vector_t result = vector_t::Zero(rows());
    // Далее итерируемся по всем блокам (потому что обычное умножение, а не потому что бесструктурная матрица!)
#pragma omp parallel for num_threads(14) shared(result, blocks, vec) schedule(static)
    for (Types::index i = 0; i < blocks.rows(); ++i) {
        for (Types::index j = 0; j < blocks.cols(); ++j) {
            // Достаем ссылку на текущий блок (тут как раз проявляется тёплицевость)
            const block_t &current_block = blocks(i, j);
            // Теперь умножаем на соответствующий подвектор
            const vector_t &sub_vector = vec.block(j * cols_in_block_, 0, cols_in_block_, 1);
            // тут пришлось скопировать, потому что block -- это не вектор, а block-expression внутри Eigen
            // тут просто получилось несоответствие типов для вызова матвека
            vector_t local_res = current_block * sub_vector;
            // Складываем результат
            { result.block(i * rows_in_block_, 0, rows_in_block_, 1) += local_res; }
        }
    }
    return result;
}

template <typename scalar_t, typename block_t>
void ToeplitzStructure<scalar_t, block_t>::matvec(Eigen::Ref<vector_t> vec, Eigen::Ref<vector_t> dest) const noexcept {
    assert(vec.size() == cols());
    for (Types::index i = 0; i < blocks.rows(); ++i) {
        for (Types::index j = 0; j < blocks.cols(); ++j) {
            // Достаем ссылку на текущий блок (тут как раз проявляется тёплицевость)
            auto &&current_block = blocks(i, j);
            // Теперь умножаем на соответствующий подвектор
            auto sub_vector = vec.block(j * cols_in_block_, 0, cols_in_block_, 1);
            if constexpr (std::is_same_v<block_t, Types::MatrixX<scalar_t>>) {
                dest.block(i * rows_in_block_, 0, rows_in_block_, 1).noalias() += current_block * sub_vector;
            } else {
                current_block.matvec(sub_vector, dest.block(i * rows_in_block_, 0, rows_in_block_, 1));
            }
        }
    }
}

/** dest += this * vec, with one matrix-matrix application per unique block.
 * Consecutive vector segments already form a column-major panel, so the first
 * level uses zero-copy Maps. Nested Toeplitz levels gather multiple RHS into
 * reusable panels and continue recursively; dense leaves use GEMM. Blocks
 * exposing matmul_wise retain their own matrix-matrix implementation.
 * Input/output overlap is supported by snapshotting the input before writing.
 * All scratch buffers are sized by the constructor, outside the timed matvec.
 * The existing noexcept API is preserved: invalid buffer sizes and size
 * overflow terminate.
 */
template <typename scalar_t, typename block_t>
void ToeplitzStructure<scalar_t, block_t>::matvec_wise_block(scalar_t *vec, size_t vec_size,
                                                           scalar_t *dest, size_t dest_size) const noexcept {
    using index_t = Eigen::Index;
    using const_matrix_ref = Eigen::Ref<const dense_matrix_t>;
    using matrix_ref = Eigen::Ref<dense_matrix_t>;
    const auto eigen_index = [](size_t size) -> index_t {
        if (size > static_cast<size_t>(std::numeric_limits<index_t>::max()))
            std::terminate();
        return static_cast<index_t>(size);
    };
    const auto product = [](index_t a, index_t b) -> index_t {
        if (a < 0 || b < 0 || (b != 0 && a > std::numeric_limits<index_t>::max() / b))
            std::terminate();
        return a * b;
    };

    if (vec_size != cols() || dest_size != rows() || (vec_size != 0 && vec == nullptr) ||
        (dest_size != 0 && dest == nullptr))
        std::terminate();
    const index_t input_size = eigen_index(vec_size), output_size = eigen_index(dest_size);
    if (input_size == 0 || output_size == 0)
        return;

    const scalar_t *input = vec;
    const std::less<const scalar_t *> less;
    if (less(vec, dest + output_size) && less(dest, vec + input_size)) {
        Eigen::Map<vector_t>(matvec_workspace_.input_alias_copy.data(), input_size) =
            Eigen::Map<const vector_t>(vec, input_size);
        input = matvec_workspace_.input_alias_copy.data();
    }
    const auto multiply_add = [&](auto &&self, const auto &matrix, const_matrix_ref rhs,
                                  matrix_ref output, size_t depth) -> void {
        if (rhs.rows() != eigen_index(matrix.cols()) || output.rows() != eigen_index(matrix.rows()) ||
            rhs.cols() != output.cols())
            std::terminate();
        if constexpr (requires { matrix.rows_in_toeplitrz(); matrix.cols_in_toeplitrz(); matrix.get_blocks(); }) {
            const index_t nr = eigen_index(matrix.rows_in_toeplitrz());
            const index_t nc = eigen_index(matrix.cols_in_toeplitrz());
            const index_t m = eigen_index(matrix.rows_in_block());
            const index_t n = eigen_index(matrix.cols_in_block());
            const index_t q = rhs.cols();
            const index_t count = eigen_index(matrix.get_blocks().get_actual_size());
            for (index_t d = 0; d < count; ++d) {
                const index_t i0 = d < nc ? 0 : d - nc + 1;
                const index_t j0 = d < nc ? d : 0;
                const index_t repeats = std::min(nr - i0, nc - j0);
                const auto &block = matrix.get_blocks()(static_cast<Types::index>(d));
                if (q == 1) {
                    Eigen::Map<const dense_matrix_t> x(rhs.data() + j0 * n, n, repeats);
                    Eigen::Map<dense_matrix_t> y(output.data() + i0 * m, m, repeats);
                    self(self, block, x, y, depth + 1);
                } else {
                    const index_t width = repeats * q; // bounded by the checked maximum above
                    assert(depth < matvec_workspace_.rhs_panels.size());
                    assert(depth < matvec_workspace_.product_panels.size());
                    assert(matvec_workspace_.rhs_panels[depth].size() >= product(n, width));
                    assert(matvec_workspace_.product_panels[depth].size() >= product(m, width));
                    Eigen::Map<dense_matrix_t> x(matvec_workspace_.rhs_panels[depth].data(), n, width);
                    Eigen::Map<dense_matrix_t> y(matvec_workspace_.product_panels[depth].data(), m, width);
                    for (index_t column = 0; column < q; ++column) {
                        Eigen::Map<const dense_matrix_t> source(rhs.data() + column * rhs.outerStride() + j0 * n,
                                                                n, repeats);
                        x.middleCols(column * repeats, repeats) = source;
                    }
                    y.setZero();
                    self(self, block, x, y, depth + 1);
                    for (index_t column = 0; column < q; ++column) {
                        Eigen::Map<dense_matrix_t> target(output.data() + column * output.outerStride() + i0 * m,
                                                         m, repeats);
                        target += y.middleCols(column * repeats, repeats);
                    }
                }
            }
        } else if constexpr (requires {
                                 matrix.matmul_wise_noalloc(
                                     rhs, output, static_cast<scalar_t *>(nullptr),
                                     static_cast<scalar_t *>(nullptr), size_t{});
                             }) {
            matrix.matmul_wise_noalloc(
                rhs, output, matvec_workspace_.factor_first.data(),
                matvec_workspace_.factor_second.data(),
                static_cast<size_t>(matvec_workspace_.factor_first.size()));
        } else if constexpr (requires { matrix.matmul_wise(rhs, output); }) {
            matrix.matmul_wise(rhs, output);
        } else {
            output.noalias() += matrix * rhs;
        }
    };

    Eigen::Map<const dense_matrix_t> rhs(input, input_size, 1);
    Eigen::Map<dense_matrix_t> output(dest, output_size, 1);
    multiply_add(multiply_add, *this, rhs, output, 0);
}

template <typename scalar_t, typename block_t>
void ToeplitzStructure<scalar_t, block_t>::matvec_wise(
    scalar_t *vec, size_t vec_size, scalar_t *dest, size_t dest_size) const noexcept {
    matvec_wise_block(vec, vec_size, dest, dest_size);
}

template <typename scalar_t, typename block_t>
ToeplitzStructure<scalar_t, block_t> ToeplitzStructure<scalar_t, block_t>::mull(scalar_t value) const noexcept {
    std::cout << "Mul with one copy" << std::endl;
    return ToeplitzStructure(*this).mull_inplace(value);
}

template <typename scalar_t, typename block_t>
void ToeplitzStructure<scalar_t, block_t>::mull_inplace(scalar_t value) noexcept {
    std::cout << "Mull inplace" << std::endl;
    for (Types::index index = 0, total = blocks.get_actual_size(); index != total; ++index) {
        if constexpr (std::is_same_v<block_t, Types::MatrixX<scalar_t>>) {
            blocks(index) *= value;
        } else {
            blocks(index).mull_inplace(value);
        }
    }
}

template <typename scalar_t, typename block_t>
Types::MatrixX<scalar_t> ToeplitzStructure<scalar_t, block_t>::to_dense() const noexcept {
    Types::MatrixX<scalar_t> result(rows(), cols());
    for (Types::index i = 0; i < blocks.rows(); ++i) {
        for (Types::index j = 0; j < blocks.cols(); ++j) {
            const auto &block = blocks(i, j);
            if constexpr (std::is_same_v<block_t, Types::MatrixX<scalar_t>>) {
                result.block(i * rows_in_block_, j * cols_in_block_, rows_in_block_, cols_in_block_) = block;
            } else {
                // std::cout << block.rows() << " " << rows_in_block_ << std::endl;
                result.block(i * rows_in_block_, j * cols_in_block_, rows_in_block_, cols_in_block_) = block.to_dense();
            }
        }
    }
    return result;
}

// --- Defined binary operators --- //
#if 0
template <typename scalar_t, typename block_t>
ToeplitzStructure<scalar_t, block_t> operator*(const ToeplitzStructure<scalar_t, block_t> &matrix, scalar_t value) {
    return matrix.mull(value);
}

template<typename scalar_t, typename block_t>
ToeplitzStructure<scalar_t, block_t> operator*(scalar_t value, const ToeplitzStructure<scalar_t, block_t> &matrix) {
    return matrix.mull(value);
}

template<typename scalar_t, typename block_t>
ToeplitzStructure<scalar_t, block_t> & operator*=(ToeplitzStructure<scalar_t, block_t>& matrix, scalar_t value) {
    matrix.mull_inplace(value);
    return matrix;
}
#endif

template <typename scalar_t, typename block_t>
Types::VectorX<scalar_t> operator*(const ToeplitzStructure<scalar_t, block_t> &matrix,
                                   const Types::VectorX<scalar_t> &vector) noexcept {
    Types::VectorXc result = Types::VectorXc::Zero(matrix.cols());
    // TODO: убрать const_cast.
    // matrix.matvec_wise(const_cast<scalar_t*>(vector.data()), vector.size(), result.data(), result.size());
    result = matrix.matvec(vector);
    return result;
}
} // namespace EMW::Math::LinAgl::Matrix

#if 0
// Code from codex

template <typename scalar_t, typename block_t>
    void ToeplitzStructure<scalar_t, block_t>::matvec_wise_block(scalar_t* vec, size_t vec_size,
                                                                 scalar_t* dest, size_t dest_size) const noexcept
    {
        matmul_wise(vec, vec_size, 1, dest, dest_size, 1);
    }

    template <typename scalar_t, typename block_t>
    void ToeplitzStructure<scalar_t, block_t>::matmul_wise(scalar_t* rhs, size_t rhs_rows, size_t rhs_cols,
                                                           scalar_t* dest, size_t dest_rows,
                                                           size_t dest_cols) const noexcept
    {
        assert(rhs != nullptr);
        assert(dest != nullptr);
        assert(rhs_rows == static_cast<size_t>(cols()));
        assert(dest_rows == static_cast<size_t>(rows()));
        assert(rhs_cols == dest_cols);

        Eigen::Map<const dense_matrix_t> rhs_matrix(rhs, static_cast<Eigen::Index>(rhs_rows),
                                                    static_cast<Eigen::Index>(rhs_cols));
        Eigen::Map<dense_matrix_t> dest_matrix(dest, static_cast<Eigen::Index>(dest_rows),
                                               static_cast<Eigen::Index>(dest_cols));
        matmul_wise(rhs_matrix, dest_matrix);
    }

    template <typename scalar_t, typename block_t>
    void ToeplitzStructure<scalar_t, block_t>::matmul_wise(Eigen::Ref<const dense_matrix_t> rhs,
                                                           Eigen::Ref<dense_matrix_t> dest) const noexcept
    {
        assert(rhs.rows() == cols());
        assert(dest.rows() == rows());
        assert(rhs.cols() == dest.cols());

        const Types::index block_rows = blocks.rows();
        const Types::index block_cols = blocks.cols();
        const Types::index rhs_cols = static_cast<Types::index>(rhs.cols());

        // 1) Уникальные блоки первой строки.
        for (Types::index d = 0; d < block_cols; ++d)
        {
            const auto& current_block = blocks(0, d);
            const Types::index repeats = std::min(block_rows, block_cols - d);

            for (Types::index i = 0; i < repeats; ++i)
            {
                auto x_block = rhs.block(static_cast<Eigen::Index>((i + d) * cols_in_block_), 0,
                                         static_cast<Eigen::Index>(cols_in_block_),
                                         static_cast<Eigen::Index>(rhs_cols));
                auto y_block = dest.block(static_cast<Eigen::Index>(i * rows_in_block_), 0,
                                          static_cast<Eigen::Index>(rows_in_block_),
                                          static_cast<Eigen::Index>(rhs_cols));

                if constexpr (std::is_same_v<block_t, Types::MatrixX<scalar_t>>)
                {
                    y_block.noalias() += current_block * x_block;
                }
                else
                {
                    current_block.matmul_wise(x_block, y_block);
                }
            }
        }

        // 2) Уникальные блоки первого столбца (кроме [0,0]).
        for (Types::index k = 1; k < block_rows; ++k)
        {
            const auto& current_block = blocks(k, 0);
            const Types::index repeats = std::min(block_cols, block_rows - k);

            for (Types::index j = 0; j < repeats; ++j)
            {
                auto x_block = rhs.block(static_cast<Eigen::Index>(j * cols_in_block_), 0,
                                         static_cast<Eigen::Index>(cols_in_block_),
                                         static_cast<Eigen::Index>(rhs_cols));
                auto y_block = dest.block(static_cast<Eigen::Index>((j + k) * rows_in_block_), 0,
                                          static_cast<Eigen::Index>(rows_in_block_),
                                          static_cast<Eigen::Index>(rhs_cols));

                if constexpr (std::is_same_v<block_t, Types::MatrixX<scalar_t>>)
                {
                    y_block.noalias() += current_block * x_block;
                }
                else
                {
                    current_block.matmul_wise(x_block, y_block);
                }
            }
        }
    }

#endif

#endif //TOEPLITZFULLYTEMPLATED_HPP
