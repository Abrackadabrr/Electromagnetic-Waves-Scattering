//
// Created by evgen on 12.02.2025.
//

#ifndef TOEPLITZFULLYTEMPLATED_HPP
#define TOEPLITZFULLYTEMPLATED_HPP

#include "ToeplitzContainer.hpp"
#include "types/Types.hpp"

#include <algorithm>
#include <cassert>
#include <exception>
#include <functional>
#include <iostream>
#include <limits>
#include <vector>

namespace EMW::Math::LinAgl::Matrix {

template <typename scalar_t, typename block_t> class ToeplitzStructure;

template <typename T> struct ToeplitzDepth : std::integral_constant<std::size_t, 0> {};

template <typename scalar_t, typename block_t>
struct ToeplitzDepth<ToeplitzStructure<scalar_t, block_t>>
    : std::integral_constant<std::size_t, 1 + ToeplitzDepth<block_t>::value> {};

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

    // глубина теплицевой структуры в данной матрице
    constexpr static std::size_t toeplitz_depth = 1 + ToeplitzDepth<block_t>::value;

  private:
    ToeplitzContainer<block_t> blocks;
    // Характеристика одного блока
    // Эти поля нужны для удобства
    Types::index rows_in_block_ = 0;
    Types::index cols_in_block_ = 0;

    // -------- MATRIX OPERATIONS KERNELS ------------ //
    void matvec_blocked_kernel(const scalar_t *vec, size_t vec_size, scalar_t *dest, size_t dest_size) const noexcept;
    void matmul_naive(const scalar_t *rhs, size_t rhs_size, scalar_t *dest, size_t dest_size,
                      size_t n_rhs_cols) const noexcept;

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

    // слегка костыльные методы для инициализации, должны быть рекурсивно заданные
    // по-хорошему, в конец нужно передать аргументы фабрики для самого внутреннего типа block_t или готовое значение
    static ToeplitzStructure Zero(Types::index block_rows, Types::index block_cols, Types::index rows_in_block,
                           Types::index cols_in_block);
    static ToeplitzStructure Zero(Types::index block_rows_l2, Types::index block_cols_l2, Types::index block_rows_l1,
                           Types::index block_cols_l1, Types::index rows_in_block, Types::index cols_in_block);

    // ----------- MATVECS --------------- //
    /** dest += matvec(vec); */
    void matvec(const vector_t &vec, vector_t &dest) const noexcept;
    /** dest += matvec(vec); но по-умному */
    void matvec_blocked(const vector_t &vec, vector_t &dest) const noexcept;

    // ----------- MATMULLS ------------------- //
    /** dest += T * A, где A -- матричный блок */
    void matmul(const dense_matrix_t &rhs, dense_matrix_t &dest) const noexcept;
    void matmul(const scalar_t *rhs, size_t rhs_size, scalar_t *dest, size_t dest_size,
                size_t n_rhs_cols) const noexcept;

    // ---------- SOME SHIT --------------- //
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
    [[nodiscard]] inline decltype(auto) get_linear_index(size_t row, size_t col) const noexcept {
        return blocks.get_linear_index(row, col);
    }
    /** Доступ к настоящим элементам в матрице, а не к блокам! */
    [[nodiscard]] const scalar_t &operator()(Types::index i, Types::index j) const noexcept;
    /** Доступ к блокам в матрице */
    const ToeplitzContainer<block_t> &get_blocks() const noexcept { return blocks; }

    // ---- Static methods --- //
    inline static Types::index get_size_of_container(Types::index rows, Types::index cols) noexcept {
        return std::max(0L, static_cast<Types::integer>(rows) + static_cast<Types::integer>(cols) - 1);
    };

    // Приведение к плотной матрице
    dense_matrix_t to_dense() const noexcept;

    // T = L + D + U, это возвращает L + U
    ToeplitzStructure extract_lower_triangle_plus_diagonal();
};

template <typename scalar_t, typename block_t>
ToeplitzStructure<scalar_t, block_t>
ToeplitzStructure<scalar_t, block_t>::Zero(Types::index block_rows, Types::index block_cols, Types::index rows_in_block,
                                           Types::index cols_in_block) {
    static_assert(toeplitz_depth == 1, "method only for single layer toeplitz structure");
    Containers::vector<block_t> zero_blocks(get_size_of_container(block_rows, block_cols),
                                            block_t::Zero(rows_in_block, cols_in_block));
    return {block_rows, block_cols, std::move(zero_blocks)};
}
template <typename scalar_t, typename block_t>
ToeplitzStructure<scalar_t, block_t>
ToeplitzStructure<scalar_t, block_t>::Zero(Types::index block_rows_l2, Types::index block_cols_l2,
                                           Types::index block_rows_l1, Types::index block_cols_l1,
                                           Types::index rows_in_block, Types::index cols_in_block) {
    static_assert(toeplitz_depth == 2, "method only for doubly toeplitz structure");
    Containers::vector<block_t> zero_blocks(get_size_of_container(block_rows_l2, block_cols_l2),
                                            block_t::Zero(block_rows_l1, block_cols_l1, rows_in_block, cols_in_block));
    return {block_rows_l2, block_cols_l2, std::move(zero_blocks)};
}

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
      cols_in_block_(blocks(0, 0).cols()) {}

template <typename scalar_t, typename block_t>
ToeplitzStructure<scalar_t, block_t>::ToeplitzStructure(Types::index block_rows, Types::index block_cols,
                                                        Containers::vector<block_t> &&blocks_)
    : blocks(block_rows, block_cols, std::move(blocks_)), rows_in_block_(blocks(0, 0).rows()),
      cols_in_block_(blocks(0, 0).cols()) {}

// как и выше, но с уже выделенным вектором, куда складывается результат
template <typename scalar_t, typename block_t>
void ToeplitzStructure<scalar_t, block_t>::matvec(const vector_t &vec, vector_t &dest) const noexcept {
    assert(vec.size() == cols());
    vector_t additional_space = vector_t::Zero(vec.size() / blocks.cols());
    for (Types::index i = 0; i < blocks.rows(); ++i) {
        for (Types::index j = 0; j < blocks.cols(); ++j) {
            // Достаем ссылку на текущий блок (тут как раз проявляется тёплицевость)
            auto &&current_block = blocks(i, j);
            // Теперь умножаем на соответствующий подвектор
            vector_t sub_vector = vec.block(j * cols_in_block_, 0, cols_in_block_, 1);
            if constexpr (std::is_same_v<block_t, Types::MatrixX<scalar_t>>) {
                // вызывается матвек для плотного блока
                dest.block(i * rows_in_block_, 0, rows_in_block_, 1).noalias() += current_block * sub_vector;
            } else {
                // иначе вызывается соотвествующая функция для внутреннего блока
                additional_space = dest.block(i * rows_in_block_, 0, rows_in_block_, 1);
                current_block.matvec(sub_vector, additional_space);
            }
        }
    }
}

template <typename scalar_t, typename block_t>
void ToeplitzStructure<scalar_t, block_t>::matvec_blocked(const vector_t &vec, vector_t &dest) const noexcept {
    matvec_blocked_kernel(const_cast<scalar_t *>(vec.data()), vec.size(), dest.data(), dest.size());
}

template <typename scalar_t, typename block_t>
void ToeplitzStructure<scalar_t, block_t>::matvec_blocked_kernel(const scalar_t *vec, size_t vec_size, scalar_t *dest,
                                                                 size_t dest_size) const noexcept {
    Types::integer cols_in_each_toeplitz_block = blocks(0).cols();
    Types::integer rows_in_each_toeplitz_block = blocks(0).rows();
    assert(vec_size == blocks.cols() * cols_in_each_toeplitz_block); // проверка на то, что я могу умножить это дело

    const scalar_t *rhs_head_pointer = vec;
    const Types::integer cols_in_toeplitz_structure =
        blocks.cols(); // количество блоков в строке (или блочных столбцов)
    const Types::integer rows_in_toeplitz_structure = blocks.rows();
    // опеределяем окна столбцов для решейпа
    const Types::integer rhs_window_len = rows_in_toeplitz_structure;
    std::pair<Types::integer, Types::integer> rhs_window{-rhs_window_len + 1, 0};
    // это окно, которое показывает ОТРЕЗОК столбцов, которые нужно взять для умножения на укальный теплицев блок
    const Types::integer lhs_window_size = cols_in_toeplitz_structure;
    std::pair<Types::integer, Types::integer> lhs_window{rows_in_toeplitz_structure - 1,
                                                         rows_in_toeplitz_structure + lhs_window_size - 2};

    // проходимся по всем блокам в toeplitz_container и делаем уникальные матмуллы
    size_t n_unique_blocks_in_toeplitz_structure = blocks.get_actual_size();
    for (size_t block_idx = 0; block_idx < n_unique_blocks_in_toeplitz_structure; block_idx++) {
        // данные для матмула с уникальным блоков
        const Types::integer rhs_col_start = std::max(0L, rhs_window.first);
        const Types::integer rhs_col_end = std::min(rhs_window.second, cols_in_toeplitz_structure - 1);
        // блок по теплицевой структуре
        const auto &block = blocks(block_idx); // взяли блок по теплицевой структуре,
                                               // тут ВАЖНО, что мы знаем как устроен внутри ToeplitzContainer
        // решейп вектора в матрицу
        const Eigen::Index n_cols_in_rhs_block = rhs_col_end - rhs_col_start + 1;
        assert(n_cols_in_rhs_block <= std::min(rows_in_toeplitz_structure, cols_in_toeplitz_structure));
        // решейп правой части
        Eigen::Map<const Eigen::MatrixXcd> reshaped_vector_block(rhs_head_pointer +
                                                                     rhs_col_start * cols_in_each_toeplitz_block,
                                                                 cols_in_each_toeplitz_block, n_cols_in_rhs_block);
        // решейп левой части
        const Types::integer lhs_col_start = std::max(0L, lhs_window.first);
        const Types::integer lhs_col_end = std::min(lhs_window.second, rows_in_toeplitz_structure - 1);
        assert(lhs_col_end - lhs_col_start == rhs_col_end - rhs_col_start); // одинаковое число столбцов в окнах
        Eigen::Map<Eigen::MatrixXcd> reshaped_dest(dest + lhs_col_start * rows_in_each_toeplitz_block,
                                                   rows_in_each_toeplitz_block, n_cols_in_rhs_block);

        // в зависимости от типа внутреннего блока нужно вызывать разное
        if constexpr (std::is_same_v<block_t, dense_matrix_t>) {
            // плотня матрица -> плотный матмул (пока что через Eigen)
            reshaped_dest.noalias() += block * reshaped_vector_block;
        } else {
            // иначе вызывается внутренний матмул блока
            block.matmul(reshaped_vector_block.data(), reshaped_vector_block.size(), reshaped_dest.data(),
                         reshaped_dest.size(), reshaped_vector_block.cols());
        }
        // двигаем окна на один столбец/строку
        rhs_window.first++;
        rhs_window.second++;
        lhs_window.first--;
        lhs_window.second--;
    }
    // ну и все, умножение закончилось
}

template <typename scalar_t, typename block_t>
void ToeplitzStructure<scalar_t, block_t>::matmul(const ToeplitzStructure::dense_matrix_t &rhs,
                                                  ToeplitzStructure::dense_matrix_t &dest) const noexcept {
    matmul_naive(rhs.data(), rhs.size(), dest.data(), dest.size(), rhs.cols());
}

template <typename scalar_t, typename block_t>
void ToeplitzStructure<scalar_t, block_t>::matmul(const scalar_t *rhs, size_t rhs_size, scalar_t *dest,
                                                  size_t dest_size, size_t n_rhs_cols) const noexcept {
    matmul_naive(rhs, rhs_size, dest, dest_size, n_rhs_cols);
}

template <typename scalar_t, typename block_t>
void ToeplitzStructure<scalar_t, block_t>::matmul_naive(const scalar_t *rhs, size_t rhs_size, scalar_t *dest,
                                                        size_t dest_size, size_t n_rhs_cols) const noexcept {
    /** Берем из правой  */
    const size_t rhs_col_size = blocks.cols() * blocks(0).cols();
    const size_t dest_col_size = blocks.rows() * blocks(0).rows();

    assert(rhs_size == rhs_col_size * n_rhs_cols);
    assert(dest_size == dest_col_size * n_rhs_cols);

    for (size_t col = 0; col < n_rhs_cols; ++col) {
        matvec_blocked_kernel(rhs + col * rhs_col_size, rhs_col_size, dest + col * dest_col_size, dest_col_size);
    }
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

template <typename scalar_t, typename block_t>
Types::VectorX<scalar_t> operator*(const ToeplitzStructure<scalar_t, block_t> &matrix,
                                   const Types::VectorX<scalar_t> &vector) noexcept {
    Types::VectorXc result = Types::VectorXc::Zero(matrix.cols());
    matrix.matvec_blocked(vector, result);
    return result;
}
} // namespace EMW::Math::LinAgl::Matrix

#endif // TOEPLITZFULLYTEMPLATED_HPP
