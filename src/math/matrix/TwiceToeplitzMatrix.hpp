//
// Created by evgen on 18.01.2025.
//

#ifndef TWICETOEPLITZMATRIX_HPP
#define TWICETOEPLITZMATRIX_HPP

#include "types/Types.hpp"

#include "Toeplitz.hpp"
#include "ToeplitzContainer.hpp"

namespace EMW::Math::LinAgl::Matrix {

/**
* Матрица со структурой тёплиц-тёплиц-общий_вид
*/
template<typename T> class TwiceToeplitzBlock {

    ToeplitzContainer<BlockToeplitz<T>> blocks;
    // Характеристика одного блока
    // Эти поля нужны для удобства
    Types::index rows_in_block_;
    Types::index cols_in_block_;

public:
    using block_t = Types::MatrixX<T>;

    TwiceToeplitzBlock() = default;
    // Консистентность состояния поддерживается, если функция возвращает блоки
    // одного и того же размера
    /**
     * @param block_rows -- количество блоков в матрице
     * @param block_cols -- количество блоков в матрице
     * @param get_block -- функция, которая возвращает квадратную матрицу, размеры одинаковы для любых пар (i, j)
     */
    TwiceToeplitzBlock(Types::index block_rows, Types::index block_cols,
                  const std::function<block_t(Types::index i, Types::index j)> &get_block)
        : blocks(block_rows, block_cols, get_block), rows_in_block_(blocks(0, 0).rows()),
          cols_in_block_(blocks(0, 0).cols()) {}

    [[nodiscard]] Types::VectorX<T> matvec(Types::VectorX<T> vec) const;

    // --- Selectors --- //
    [[nodiscard]] const block_t &get_block(Types::index row, Types::index col) const {
        return blocks.get_block(row, col);
    }

    // Возвращают значения строк и столбцов в каждом блоке
    [[nodiscard]] Types::index rows_in_block() const { return rows_in_block_; }
    [[nodiscard]] Types::index cols_in_block() const { return cols_in_block_; }
    // Возвращают значение строк и столбцов во всей матрице
    [[nodiscard]] Types::index rows() const { return rows_in_block_ * blocks.rows(); }
    [[nodiscard]] Types::index cols() const { return cols_in_block_ * blocks.cols(); }
};

template <typename scalar_t>
Types::VectorX<scalar_t> TwiceToeplitzBlock<scalar_t>::matvec(Types::VectorX<scalar_t> vec) const {
    assert(vec.size() == cols());
    // создаем нулевой вектор результата, в который будем записывать ответ
    Types::VectorX<scalar_t> result = Types::VectorX<scalar_t>::Zero(rows());

    // Вектор для результата локального блочного умножения
    Types::VectorX<scalar_t> local_res = Types::VectorX<scalar_t>::Zero(rows_in_block_);

    // Далее итерируемся по всем блокам (потому что обычное умножение, а не потому что бесструктурная матрица!)
    for (Types::index i = 0; i < blocks.rows(); ++i) {
        for (Types::index j = 0; j < blocks.cols(); ++j) {
            // Достаем ссылку на текущий блок
            const block_t &current_block = blocks.get_block(i, j);
            // Теперь умножаем на соответствующий подвектор
            local_res = current_block * vec.block(j * cols_in_block_, 0, cols_in_block_, 1);;
            // Складываем результат
            result.block(i * rows_in_block_, 0, rows_in_block_, 1) += local_res;
        }
    }
    return result;
}

template <typename scalar_t>
Types::VectorX<scalar_t> operator*(const TwiceToeplitzBlock<scalar_t> &matrix, const Types::VectorX<scalar_t> &vector) {
    return matrix.matvec(vector);
}

}

#endif //TWICETOEPLITZMATRIX_HPP
