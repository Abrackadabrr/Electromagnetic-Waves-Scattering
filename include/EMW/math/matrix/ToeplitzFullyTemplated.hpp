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
    ToeplitzContainer<block_t> blocks;
    // Характеристика одного блока
    // Эти поля нужны для удобства
    Types::index rows_in_block_ = 0;
    Types::index cols_in_block_ = 0;

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

    /** Работа с индексами */
    [[nodiscard]] inline decltype(auto) get_linear_index(size_t row, size_t col) const noexcept {
        return blocks.get_linear_index(row, col);
    }

    /** Умножение матрицы на вектор c возвращением копии */
    [[nodiscard]] vector_t matvec(const vector_t &vec) const noexcept;
    /** dest += matvec(vec); */
    void matvec(const vector_t &vec, vector_t &dest) const noexcept;
    /** dest += matvec(vec); но по-умному */
    void matvec_blocked(const vector_t &vec, vector_t &dest) const noexcept;
    void matvec_blocked_kernel(scalar_t *vec, size_t vec_size, scalar_t *dest, size_t dest_size) const noexcept;

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
    inline static Types::index get_size_of_container(Types::index rows, Types::index cols) noexcept {
        return rows + cols - 1;
    };

    // Приведение к плотной матрице
    dense_matrix_t to_dense() const noexcept;
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
      cols_in_block_(blocks(0, 0).cols()) {}

template <typename scalar_t, typename block_t>
ToeplitzStructure<scalar_t, block_t>::ToeplitzStructure(Types::index block_rows, Types::index block_cols,
                                                        Containers::vector<block_t> &&blocks_)
    : blocks(block_rows, block_cols, std::move(blocks_)), rows_in_block_(blocks(0, 0).rows()),
      cols_in_block_(blocks(0, 0).cols()) {}

// обычный вариант блочного матвека, в котором выделяются блоки в матрице и в векторе, а дальше происходит умножение
// по обычным блочным формулами
template <typename scalar_t, typename block_t>
auto ToeplitzStructure<scalar_t, block_t>::matvec(const vector_t &vec) const noexcept -> vector_t {
    assert(vec.size() == cols());
    // создаем нулевой вектор результата, в который будем записывать ответ
    vector_t result = vector_t::Zero(rows());
    // Далее итерируемся по всем блокам (потому что обычное умножение, а не потому что бесструктурная матрица!)
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
            {
                result.block(i * rows_in_block_, 0, rows_in_block_, 1) += local_res;
            }
        }
    }
    return result;
}

// как и выше, но с уже выделенным вектором, куда складывается результат
template <typename scalar_t, typename block_t>
void ToeplitzStructure<scalar_t, block_t>::matvec(const vector_t &vec, vector_t &dest) const noexcept {
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

template <typename scalar_t, typename block_t>
void ToeplitzStructure<scalar_t, block_t>::matvec_blocked(const vector_t &vec, vector_t &dest) const noexcept {
    matvec_blocked_kernel(const_cast<scalar_t *>(vec.data()), vec.size(), dest.data(), dest.size());
}

template <typename scalar_t, typename block_t>
void ToeplitzStructure<scalar_t, block_t>::matvec_blocked_kernel(scalar_t *vec, size_t vec_size, scalar_t *dest,
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
        if constexpr (std::is_same_v<block_t, Types::MatrixX<scalar_t>>) {
            // делается плотный матмул
            const Eigen::Index n_cols_in_rhs_block = rhs_col_end - rhs_col_start + 1;
            assert(n_cols_in_rhs_block <= std::min(rows_in_toeplitz_structure, cols_in_toeplitz_structure));
            // решейп правой части
            Eigen::Map<const Eigen::MatrixXcd> reshaped_vector_block(rhs_head_pointer +
                                                                         rhs_col_start * cols_in_each_toeplitz_block,
                                                                     cols_in_each_toeplitz_block, n_cols_in_rhs_block);
            // решейп левой части
            const Types::integer lhs_col_start = std::max(0L, lhs_window.first);
            const Types::integer lhs_col_end = std::min(lhs_window.second, rows_in_toeplitz_structure - 1);

            assert(lhs_col_end - lhs_col_start ==
                   rhs_col_end - rhs_col_start); // окна содержат одинаковое число столбцов

            Eigen::Map<Eigen::MatrixXcd> reshaped_dest(dest + lhs_col_start * rows_in_each_toeplitz_block,
                                                       rows_in_each_toeplitz_block, n_cols_in_rhs_block);
            reshaped_dest.noalias() += block * reshaped_vector_block;
        } else {
            // иначе вызывается внутренний матмул блока
            throw std::runtime_error("not implemented yet");
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

#endif // TOEPLITZFULLYTEMPLATED_HPP
