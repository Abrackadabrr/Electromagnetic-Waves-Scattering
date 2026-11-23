//
// Created by evgen on 10.02.2025.
//

#ifndef BASETOEPLITZ_HPP
#define BASETOEPLITZ_HPP

#include "types/Types.hpp"

namespace EMW::Math::LinAgl::Matrix {

/**
 * Тёплицев контейнер для данных со специальным форматом хранения.
 *
 * @tparam data_type -- тип элементов матрицы
 */
template <typename data_type> class ToeplitzContainer {
    static_assert(std::is_default_constructible_v<data_type>);

    // Структура хранения
    Types::index rows_ = 0;
    Types::index cols_ = 0;

    // Маппинг нумерации (с 1):
    // 3  4  5  6  7
    // 2  3  4  5  6
    // 1  2  3  4  5
    Containers::vector<data_type> values;

    // Инвариант:
    // rows_ == cols_ == 0 => values.empty()
    // rows_ > 0 && cols_ > 0 => values.size() == rows_ + cols_ - 1

    static Types::index storage_size(Types::index rows, Types::index cols) {
        if (rows == 0 && cols == 0)
            return 0;

        if (rows <= 0 || cols <= 0) {
            throw std::invalid_argument("ToeplitzContainer: either both dimensions are zero, "
                                        "or both are positive");
        }

        if (rows - 1 > std::numeric_limits<Types::index>::max() - cols) {
            throw std::length_error("ToeplitzContainer: storage size overflow");
        }

        return rows - 1 + cols;
    }

  public:
    ToeplitzContainer() = default;

    ToeplitzContainer(Types::index rows, Types::index cols, data_type value = data_type{})
        : rows_(rows), cols_(cols), values(storage_size(rows, cols), value) {}

    ToeplitzContainer(Types::index rows, Types::index cols,
                      const std::function<data_type(Types::index row, Types::index col)> &function);

    ToeplitzContainer(Types::index rows, Types::index cols, const std::function<data_type(Types::index)> &row_function,
                      const std::function<data_type(Types::index)> &col_function);

    ToeplitzContainer(const Containers::vector<data_type> &row, const Containers::vector<data_type> &col);

    ToeplitzContainer(Types::index rows, Types::index cols, Containers::vector<data_type> &&values);

    // --- Selectors --- //

    [[nodiscard]] Types::index rows() const noexcept { return rows_; }

    [[nodiscard]] Types::index cols() const noexcept { return cols_; }

    [[nodiscard]] Types::index get_linear_index(Types::index row_index, Types::index col_index) const {

        if (std::cmp_less(row_index, 0) || std::cmp_less(col_index, 0) || row_index >= rows_ || col_index >= cols_) {

            throw std::out_of_range("ToeplitzContainer: row or column index out of range");
        }

        return rows_ - 1 - row_index + col_index;
    }

    [[nodiscard]] Types::index get_actual_size() const noexcept { return values.size(); }

    [[nodiscard]] Types::index diagonal_length() const noexcept { return std::min(rows_, cols_);};

    /** Доступ к элементу контейнера по заданным координатам */
    const data_type &operator[](Types::index row, Types::index col) const;

    /** Доступ к элементу контейнера по заданным координатам */
    data_type &operator[](Types::index row, Types::index col);

    /** Доступ к элементу контейнера по заданным координатам */
    const data_type &operator()(Types::index row, Types::index col) const;

    /** Доступ к элементу контейнера по заданным координатам */
    data_type &operator()(Types::index row, Types::index col);

    /** Доступ к элементу контейнера по линейному индексу */
    const data_type &operator()(Types::index index) const;

    /** Доступ к элементу контейнера по линейному индексу */
    data_type &operator()(Types::index index);

    /** Iterator-type assess */
    auto begin() { return values.begin(); };
    auto end() { return values.end(); };
    auto cbegin() { return values.cbegin(); };
    auto cend() { return values.cend(); };

    /** Доступ на чтение к внутренней структуре */
    const Containers::vector<data_type> &get_values() const noexcept { return values; }
};

// --- Constructors --- //

template <typename data_type>
ToeplitzContainer<data_type>::ToeplitzContainer(
    Types::index rows, Types::index cols, const std::function<data_type(Types::index row, Types::index col)> &function)
    : ToeplitzContainer(rows, cols) {

    // Первый столбец снизу вверх
    for (Types::index index = 0; index < rows_; ++index) {
        values[index] = function(rows_ - 1 - index, 0);
    }

    // Первая строка слева направо, без первого элемента
    for (Types::index index = 1; index < cols_; ++index) {
        values[rows_ - 1 + index] = function(0, index);
    }
}

template <typename data_type>
ToeplitzContainer<data_type>::ToeplitzContainer(Types::index rows, Types::index cols,
                                                const std::function<data_type(Types::index)> &row_function,
                                                const std::function<data_type(Types::index)> &col_function)
    : ToeplitzContainer(rows, cols) {

    // Первый столбец снизу вверх
    for (Types::index index = 0; index < rows_; ++index) {
        values[index] = col_function(rows_ - 1 - index);
    }

    // Первая строка слева направо, без первого элемента
    for (Types::index index = 1; index < cols_; ++index) {
        values[rows_ - 1 + index] = row_function(index);
    }
}

template <typename data_type>
ToeplitzContainer<data_type>::ToeplitzContainer(const Containers::vector<data_type> &row,
                                                const Containers::vector<data_type> &col)
    : rows_(col.size()), cols_(row.size()) {

    values.reserve(storage_size(rows_, cols_));

    // Первый столбец снизу вверх
    for (Types::index index = 0; index < rows_; ++index) {
        values.push_back(col[rows_ - 1 - index]);
    }

    // Первая строка слева направо, без первого элемента
    for (Types::index index = 1; index < cols_; ++index) {
        values.push_back(row[index]);
    }
}

template <typename data_type>
ToeplitzContainer<data_type>::ToeplitzContainer(Types::index rows, Types::index cols,
                                                Containers::vector<data_type> &&values_)
    : rows_(rows), cols_(cols), values(std::move(values_)) {

    if (values.size() != storage_size(rows_, cols_)) {
        throw std::invalid_argument("ToeplitzContainer: incorrect number of stored values");
    }
}

// --- Element access --- //

template <typename data_type>
const data_type &ToeplitzContainer<data_type>::operator[](Types::index row, Types::index col) const {

    return values[get_linear_index(row, col)];
}

template <typename data_type> data_type &ToeplitzContainer<data_type>::operator[](Types::index row, Types::index col) {

    return values[get_linear_index(row, col)];
}

template <typename data_type>
const data_type &ToeplitzContainer<data_type>::operator()(Types::index row, Types::index col) const {

    return operator[](row, col);
}

template <typename data_type> data_type &ToeplitzContainer<data_type>::operator()(Types::index row, Types::index col) {

    return operator[](row, col);
}

template <typename data_type> const data_type &ToeplitzContainer<data_type>::operator()(Types::index index) const {

    if (std::cmp_less(index, 0) || std::cmp_greater_equal(index, values.size())) {

        throw std::out_of_range("ToeplitzContainer: linear index out of range");
    }

    return values[index];
}

template <typename data_type> data_type &ToeplitzContainer<data_type>::operator()(Types::index index) {

    if (std::cmp_less(index, 0) || std::cmp_greater_equal(index, values.size())) {

        throw std::out_of_range("ToeplitzContainer: linear index out of range");
    }

    return values[index];
}

} // namespace EMW::Math::LinAgl::Matrix

#endif // BASETOEPLITZ_HPP
