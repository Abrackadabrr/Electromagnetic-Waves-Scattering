//
// Created by evgen on 12.02.2025.
//

#ifndef TOEPLITZFULLYTEMPLATED_HPP
#define TOEPLITZFULLYTEMPLATED_HPP

#include "types/Types.hpp"

#include "ToeplitzContainer.hpp"
#include <omp.h>

#include <cassert>
#include <iostream>

namespace EMW::Math::LinAgl::Matrix {

/**
 * Матрица со структурой тёплиц-тёплиц-общий_вид
 */
template <typename scalar_t, typename block_t> class ToeplitzStructure {
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

    // --- Selectors --- //
    [[nodiscard]] const block_t &get_block(Types::index row, Types::index col) const noexcept {
        return blocks(row, col);
    }
    [[nodiscard]] block_t &get_block(Types::index row, Types::index col) noexcept { return blocks(row, col); }

    // Возвращают значения строк и столбцов в каждом блоке
    [[nodiscard]] Types::index rows_in_block() const noexcept { return rows_in_block_; }
    [[nodiscard]] Types::index cols_in_block() const noexcept { return cols_in_block_; }
    // Возвращают значение строк и столбцов во всей матрице
    [[nodiscard]] Types::index rows() const noexcept { return rows_in_block_ * blocks.rows(); }
    [[nodiscard]] Types::index cols() const noexcept { return cols_in_block_ * blocks.cols(); }

    // ---- Static methods --- //
    inline static Types::index get_size_of_container(Types::index rows, Types::index cols) noexcept
        __attribute__((always_inline)) {
        return rows + cols - 1;
    };
};

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
#pragma omp parallel for schedule(dynamic) collapse(2)
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
