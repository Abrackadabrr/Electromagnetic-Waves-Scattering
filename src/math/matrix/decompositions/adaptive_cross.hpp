//
// Created by evgen on 15.05.2025.
//

#ifndef ADAPTIVE_CROSS_HPP
#define ADAPTIVE_CROSS_HPP

#include "types/Types.hpp"
#include "math/matrix/DynamicFactoredMatrix.hpp"


namespace EMW::Math::LinAgl::Decompositions {
template <typename matrix_t, typename vector_t, typename value_t> struct ACA {

    struct element_indexing {
        Types::index i, j;
    };

    /*
     * Критерий остановки для метода адаптивной крестовой аппроксимации
     */
    static inline bool stop_criterion(Types::scalar uv_prev_norm, Types::scalar e, Types::index n, Types::index m,
                               Types::index rank, value_t ij_element) {
        return e * uv_prev_norm >= std::abs(ij_element) * std::sqrt((m - rank) * (n - rank));
    }

    /*
     * Вычислить строку в матрице
     */
    template <typename MatrixElementFunction>
    static value_t get_row(Types::index i, Types::index n, Types::index m, const MatrixElementFunction &&element_function) {
        auto result = vector_t::Zero(m);
        for (Types::index j = 0; j < m; ++j) {
            result[j] = std::forward<MatrixElementFunction>(element_function)(i, j);
        }
        return result;
    }

    /*
     * Вычислить столбец в матрице
     */
    template <typename MatrixElementFunction>
    static value_t get_col(Types::index j, Types::index n, Types::index m, const MatrixElementFunction &&element_function) {
        auto result = vector_t::Zero(n);
        for (Types::index i = 0; i < m; ++i) {
            result[i] = std::forward<MatrixElementFunction>(element_function)(i, j);
        }
        return result;
    }

    /*
     * Вычислить приближение к максимальному элементу в матрице
     */
    template <typename MatrixElementFunction>
    static element_indexing get_argmax_of_matrix_elements(const MatrixElementFunction &&element_function, Types::index n,
                                                   Types::index m, Types::index start_column) {
        int iterations = 2;
        int current_column = start_column;
        int current_row = 0;
        for (int i = 0; i < iterations; ++i) {
            // Проходка вдоль столбца
            const auto &&col = get_col(current_column, n, m, element_function);
            // Нашли строку максимального элемента в столбце
            current_row = std::max_element(col.begin(), col.end()) - col.begin();

            // Проходка вдоль строки
            const auto &&row = get_row(current_row, n, m, element_function);
            // Нашли столбец максимального элемента в строке
            current_column = std::max_element(row.begin(), row.end()) - row.begin();

            // TODO: вставить проверку на то, что колонка не изменилась, а значит, что пора останавливать функцию
        }
        return {current_row, current_column};
    }

    template <typename MatrixElementFunction>
    static element_indexing get_argmax_of_residual_matrix(const MatrixElementFunction &&element_function, const matrix_t &U,
                                                   const matrix_t &V, Types::index start_column) {
        return {0, 0};
    }

    static Types::scalar norm_update(Types::scalar past_norm, const matrix_t &U, const matrix_t &V, const vector_t &u_new, const vector_t &v_new) {
        // Обновляем норму по аналитической формуле
        return std::sqrt(past_norm * past_norm +
                         2 * std::real((U.transpose() * u_new).dot(V * v_new)) +
                         u_new.squared_norm() * v_new.squared_norm());
    }

    template <typename MatrixElementFunction>
    static LinAgl::Matrix::DynamicFactoredMatrix<matrix_t> compute(const MatrixElementFunction&& element,
                                                                   Types::index n, Types::index m,
                                                                   Types::scalar error_control_parameter) {

    };

}

#endif //ADAPTIVE_CROSS_HPP
