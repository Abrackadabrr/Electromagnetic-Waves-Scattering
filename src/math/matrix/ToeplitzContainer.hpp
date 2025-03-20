//
// Created by evgen on 10.02.2025.
//

#ifndef BASETOEPLITZ_HPP
#define BASETOEPLITZ_HPP

#include "types/Types.hpp"

namespace EMW::Math::LinAgl::Matrix {
/**
 * Тёплицев контейнер для данных со специальным форматом хранения
 * Элементы контейнера -- это что угодно, хоть май классы
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
    Containers::vector<data_type> values;

    // Инвариант класса: rows_ + cols_ - 1 == values.size();

  public:
    ToeplitzContainer() = default;
    ToeplitzContainer(Types::index rows, Types::index cols, data_type value = data_type{})
        : rows_(rows), cols_(cols), values(rows_ + cols_ - 1, value){};
    ToeplitzContainer(Types::index rows, Types::index cols,
                      const std::function<data_type(Types::index row, Types::index col)> &function);
    ToeplitzContainer(Types::index rows, Types::index cols,
                      const std::function<data_type(Types::index index)> &row_function,
                      const std::function<data_type(Types::index index)> &col_function);
    ToeplitzContainer(const Containers::vector<data_type> &row, const Containers::vector<data_type> &col);
    ToeplitzContainer(Types::index rows, Types::index cols, Containers::vector<data_type> &&values);

    // --- Selectors --- //
    [[nodiscard]] Types::index rows() const noexcept { return rows_; }
    [[nodiscard]] Types::index cols() const noexcept { return cols_; }

    [[nodiscard]] Types::index get_linear_index(Types::index row_index, Types::index col_index) const noexcept {
        return (col_index >= row_index) ? (col_index - row_index) : (row_index - col_index + cols_ - 1);
    }

    [[nodiscard]] Types::index get_actual_size() const noexcept { return values.size(); }

    /** Доступ к элементу контейнера по заданным координатам */
    const data_type &operator[](Types::index row, Types::index col) const noexcept;
    /** Доступ к элементу контейнера по заданным координатам */
    data_type &operator[](Types::index row, Types::index col) noexcept;
    /** Доступ к элементу контейнера по заданным координатам. Определен для обратной совместимости */
    const data_type &operator()(Types::index row, Types::index col) const noexcept;
    /** Доступ к элементу контейнера по заданным координатам. Определен для обратной совместимости */
    data_type &operator()(Types::index row, Types::index col) noexcept;

    /** Далее идут служебные функции, не для пользовательского использования */
    /** Доступ к элементу контейнера по линейному индексу */
    const data_type &operator()(Types::index index) const noexcept;
    /** Доступ к элементу контейнера по линейному индексу */
    data_type &operator()(Types::index index) noexcept;
    /** Доступ на чтение к внутренней структуре */
    Containers::vector<data_type> get_values() const noexcept { return values; }
};

template <typename data_type>
ToeplitzContainer<data_type>::ToeplitzContainer(
    Types::index rows, Types::index cols, const std::function<data_type(Types::index row, Types::index col)> &function)
    : ToeplitzContainer(rows, cols) {
    for (Types::index index = 0; index < cols_; ++index) {
        values[index] = function(0, index);
    }
    for (Types::index index = 1; index < rows_; ++index) {
        values[index + cols_ - 1] = function(index, 0);
    }
    if (!(values.size() == rows_ + cols_ - 1))
        throw std::invalid_argument("ToeplitzContainer constructor with one function failed");
}

template <typename data_type>
ToeplitzContainer<data_type>::ToeplitzContainer(Types::index rows, Types::index cols,
                                                const std::function<data_type(Types::index index)> &row_function,
                                                const std::function<data_type(Types::index index)> &col_function)
    : ToeplitzContainer(rows, cols) {
    for (Types::index index = 0; index < cols_; ++index) {
        values[index] = row_function(index);
    }
    for (Types::index index = 1; index < rows_; ++index) {
        values[index + cols_ - 1] = col_function(index);
    }
    if (!(values.size() == rows_ + cols_ - 1))
        throw std::invalid_argument("ToeplitzContainer constructor with two functions failed");
}

template <typename data_type>
ToeplitzContainer<data_type>::ToeplitzContainer(const Containers::vector<data_type> &row,
                                                const Containers::vector<data_type> &col)
    : rows_(row.size()), cols_(col.size()) {
    values.insert(values.end(), row.begin(), row.end());
    values.insert(values.end(), col.begin() + 1, col.end());
    if (!(values.size() == rows_ + cols_ - 1))
        throw std::invalid_argument("ToeplitzContainer constructor with &row and &col failed");
}

template <typename data_type>
ToeplitzContainer<data_type>::ToeplitzContainer(Types::index rows, Types::index cols,
                                                Containers::vector<data_type> &&values_)
    : rows_(rows), cols_(cols), values(std::move(values_)) {
    if (!(values.size() == rows_ + cols_ - 1))
        throw std::invalid_argument("ToeplitzContainer constructor with &&vector failed");
}

template <typename data_type>
const data_type &ToeplitzContainer<data_type>::operator[](Types::index row_index,
                                                          Types::index col_index) const noexcept {
    return col_index >= row_index ? values[col_index - row_index] : values[row_index - col_index + cols_ - 1];
}

template <typename data_type>
data_type &ToeplitzContainer<data_type>::operator[](Types::index row_index, Types::index col_index) noexcept {
    return col_index >= row_index ? values[col_index - row_index] : values[row_index - col_index + cols_ - 1];
}

template <typename data_type>
const data_type &ToeplitzContainer<data_type>::operator()(Types::index row, Types::index col) const noexcept{
    return operator[](row, col);
}

template <typename data_type>
data_type& ToeplitzContainer<data_type>::operator()(Types::index row, Types::index col) noexcept {
    return operator[](row, col);
}

template <typename data_type>
const data_type &ToeplitzContainer<data_type>::operator()(Types::index index) const noexcept{
    return values[index];
}

template <typename data_type>
data_type& ToeplitzContainer<data_type>::operator()(Types::index index) noexcept {
    return values[index];
}

}

#endif //BASETOEPLITZ_HPP
