//
// Created by evgen on 10.02.2025.
//

#ifndef BASETOEPLITZ_HPP
#define BASETOEPLITZ_HPP

#include "types/Types.hpp"

namespace EMW::Math::LinAgl::Matrix {
/**
 * Тёплицева матрица со специальным форматом хранения
 * Элементы матрицы -- это блоки общего вида
 *
 * @tparam data_type -- тип чисел, которые будут храниться в матрице
 */
template <typename data_type> class ToeplitzContainer {
    // Проверяем наличие дефолт-конструктора
    static_assert(std::is_default_constructible_v<data_type>);

    // Структура хранения
    Types::index rows_ = 0;
    Types::index cols_ = 0;
    // Контейнер для элементов в специальном виде. Маппинг нумерации приведен ниже (c 1 начиная почему-то)
    // 1  2  3  4  5
    // 6  1  2  3  4
    // 7  6  1  2  3
    // 8  7  6  1  2
    Containers::vector<data_type> values;

  public:
    ToeplitzContainer() = default;
    ToeplitzContainer(Types::index rows, Types::index cols): rows_(rows), cols_(cols) {};
    ToeplitzContainer(Types::index rows, Types::index cols, data_type value): rows_(rows), cols_(cols), values(rows_ + cols_ - 1, value) {};
    ToeplitzContainer(Types::index rows, Types::index cols, const std::function<data_type(Types::index row, Types::index col)>& function);

    // --- Selectors --- //
    /** Доступ к элементу матрицы по заданным координатам */
    data_type operator[](Types::index row, Types::index col) const noexcept;
    /** Доступ к элементу матрицы по заданным координатам */
    data_type &operator[](Types::index row, Types::index col) noexcept;
    /** Доступ к элементу матрицы по заданным координатам. Определен для обратной совместимости */
    data_type operator()(Types::index row, Types::index col) const noexcept;
    /** Доступ к элементу матрицы по заданным координатам. Определен для обратной совместимости */
    data_type &operator()(Types::index row, Types::index col) noexcept;
};

template<typename data_type>
ToeplitzContainer<data_type>::ToeplitzContainer(Types::index rows, Types::index cols,
    const std::function<data_type(Types::index row, Types::index col)>& function) : ToeplitzContainer<data_type>(rows, cols) {
    values.resize(rows_ + cols_ - 1);
    for (Types::index index = 0; index < cols_; ++index) {
        values.push_back(function(0, index));
    }
    for (Types::index index = 1; index < rows_; ++index) {
        values.push_back(function(index, 0));
    }
    assert(values.size() == rows_ + cols_ - 1);
}

template <typename data_type>
data_type ToeplitzContainer<data_type>::operator[](Types::index row_index, Types::index col_index) const noexcept {
    return col_index >= row_index ? values[col_index - row_index] : values[row_index - col_index + cols_ - 1];
}

template <typename data_type>
data_type &ToeplitzContainer<data_type>::operator[](Types::index row_index, Types::index col_index) noexcept {
    return col_index >= row_index ? values[col_index - row_index] : values[row_index - col_index + cols_ - 1];
}

template <typename data_type>
data_type ToeplitzContainer<data_type>::operator()(Types::index row, Types::index col) const noexcept{
    return operator[](row, col);
}

template <typename data_type>
data_type& ToeplitzContainer<data_type>::operator()(Types::index row, Types::index col) noexcept {
    return operator[](row, col);
}

}

#endif //BASETOEPLITZ_HPP
