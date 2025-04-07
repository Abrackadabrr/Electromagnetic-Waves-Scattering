//
// Created by evgen on 12.02.2025.
//

#ifndef TOEPLITZFULLYTEMPLATED_HPP
#define TOEPLITZFULLYTEMPLATED_HPP

#include "ToeplitzContainer.hpp"
#include "types/Types.hpp"

#include <cassert>
#include <iostream>

namespace EMW::Math::LinAgl::Matrix {

/**
 * Матрица со структурой тёплиц-тёплиц-общий_вид
 */
template <typename scalar_t, typename block_t> class ToeplitzStructure {
    // Проверка на дефолт контруирование
    static_assert(std::is_default_constructible_v<block_t>);
  public:
    using vector_t = Types::VectorX<scalar_t>;
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

    /** Умножение матрицы на вектор */
    [[nodiscard]] vector_t matvec(const vector_t &vec) const noexcept;
    /** Умножение матрицы на число с возвращением копии */
    [[nodiscard]] ToeplitzStructure mull(scalar_t value) const noexcept;
    /** Умножение себя на число */
    const ToeplitzStructure &mull_inplace(scalar_t value) noexcept;
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
    // имеется ввиду что-то типа outerSize в терминологии Eigen
    [[nodiscard]] Types::index rows() const noexcept { return rows_in_block_ * blocks.rows(); }
    [[nodiscard]] Types::index cols() const noexcept { return cols_in_block_ * blocks.cols(); }

    [[nodiscard]] Types::index rows_in_toeplitrz() { return blocks.rows(); };
    [[nodiscard]] Types::index cols_in_toeplitrz() { return blocks.cols(); };

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
    const Types::index i_new = i % rows_in_block_;
    const Types::index j_new = j % cols_in_block_;
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

template <typename scalar_t, typename block_t>
typename ToeplitzStructure<scalar_t, block_t>::vector_t
ToeplitzStructure<scalar_t, block_t>::matvec(const vector_t &vec) const noexcept {

    assert(vec.size() == cols());
    // создаем нулевой вектор результата, в который будем записывать ответ
    vector_t result = vector_t::Zero(rows());

    // Далее итерируемся по всем блокам (потому что обычное умножение, а не потому что бесструктурная матрица!)
#pragma omp parallel for schedule(static) collapse(2) num_threads(14)
    for (Types::index i = 0; i < blocks.rows(); ++i) {
        for (Types::index j = 0; j < blocks.cols(); ++j) {
            // Достаем ссылку на текущий блок (тут как раз проявляется тёплицевость)
            const block_t &current_block = blocks(i, j);
            // Теперь умножаем на соответствующий подвектор
            const vector_t &sub_vector = vec.block(j * cols_in_block_, 0, cols_in_block_, 1);
            // тут пришлось скопировать, потому что block -- это не вектор, а block-expression внутри Eigen
            vector_t local_res = current_block * sub_vector;
            // Складываем результат
#pragma omp critical
            { result.block(i * rows_in_block_, 0, rows_in_block_, 1) += local_res; }
        }
    }
    return result;
}

template <typename scalar_t, typename block_t>
ToeplitzStructure<scalar_t, block_t> ToeplitzStructure<scalar_t, block_t>::mull(scalar_t value) const noexcept {
    std::cout << "Mul with one copy" << std::endl;
    return ToeplitzStructure(*this).mull_inplace(value);
}

template <typename scalar_t, typename block_t>
const ToeplitzStructure<scalar_t, block_t> &
ToeplitzStructure<scalar_t, block_t>::mull_inplace(scalar_t value) noexcept {
    std::cout << "Mull inplace" << std::endl;
    for (Types::index index = 0, total = blocks.get_actual_size(); index != total; ++index) {
        blocks(index) *= value;
    }
    return *this;
}

template <typename scalar_t, typename block_t>
Types::MatrixX<scalar_t> ToeplitzStructure<scalar_t, block_t>::to_dense() const noexcept {
    Types::MatrixX<scalar_t> result(rows(), cols());
    for (Types::index i = 0; i < blocks.rows(); ++i) {
        for (Types::index j = 0; j < blocks.cols(); ++j) {
            const auto& block = blocks(i, j);
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

template <typename scalar_t, typename block_t>
Types::VectorX<scalar_t> operator*(const ToeplitzStructure<scalar_t, block_t> &matrix, const Types::VectorX<scalar_t> &vector) noexcept {
    return matrix.matvec(vector);
}



}


#endif //TOEPLITZFULLYTEMPLATED_HPP
