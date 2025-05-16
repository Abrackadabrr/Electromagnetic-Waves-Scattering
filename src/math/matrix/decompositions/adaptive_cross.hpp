//
// Created by evgen on 15.05.2025.
//

#ifndef ADAPTIVE_CROSS_HPP
#define ADAPTIVE_CROSS_HPP

#include "math/matrix/DynamicFactoredMatrix.hpp"
#include "types/Types.hpp"

#include <numeric>

namespace EMW::Math::LinAgl::Decompositions {
template <typename matrix_t, typename vector_t, typename value_t> struct ACA {

    struct element_indexing {
        Types::index i, j;
        value_t value;
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
    static vector_t get_row(Types::index i, Types::index n, Types::index m,
                            const MatrixElementFunction &&element_function) {
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
    static vector_t get_col(Types::index j, Types::index n, Types::index m,
                            const MatrixElementFunction &&element_function) {
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
    static element_indexing get_argmax_of_matrix_elements(const MatrixElementFunction &&element_function,
                                                          Types::index n, Types::index m, Types::index start_column) {
        Types::index iterations = 2; // это подкручиваемый параметр
        Types::index current_column = start_column;
        Types::index current_row = 0;
        value_t max_element_value = 0;
        for (int i = 0; i < iterations; ++i) {
            // Проходка вдоль столбца
            const auto &&col = get_col(current_column, n, m, element_function);
            const auto&& col_abs = col.cwiseAbs();
            // Нашли строку максимального элемента в столбце
            current_row = std::max_element(col.begin(), col.end()) - col.begin();

            // Проходка вдоль строки
            const auto &&row = get_row(current_row, n, m, element_function);
            const auto&& row_abs = row.cwiseAbs();
            // Нашли столбец максимального элемента в строке
            current_column = std::max_element(row.begin(), row.end()) - row.begin();

            // TODO: вставить проверку на то, что колонка не изменилась, а значит, что пора останавливать функцию
        }
        return {current_row, current_column, };
    }

    /**
     *  Поиск максимального элемента в новой матрице, то есть в матрице,
     *  скорректированной на матрицы U и V
     *
     *  Пока что функция не написана
     */
    template <typename MatrixElementFunction>
    static element_indexing get_argmax_of_residual_matrix(const MatrixElementFunction &&element_function,
                                                          Types::index n, Types::index m, const matrix_t &U,
                                                          const matrix_t &V, Types::index start_column) {
        int iterations = 2; // это подкручиваемый параметр
        int current_column = start_column;
        int current_row = 0;
        for (int i = 0; i < iterations; ++i) {
            // Проходка вдоль столбца
            const auto &&col = get_col(current_column, n, m, element_function) - get_col_UV(U, V, current_column);
            const auto&& col_abs = col.cwiseAbs();
            // Нашли строку максимального элемента в столбце
            current_row = std::max_element(col_abs.begin(), col_abs.end()) - col_abs.begin();

            // Проходка вдоль строки
            const auto &&row = get_row(current_row, n, m, element_function) - get_row_UV(U, V, current_row);
            const auto&& row_abs = row.cwiseAbs();
            // Нашли столбец максимального элемента в строке
            current_column = std::max_element(row_abs.begin(), row_abs.end()) - row_abs.begin();

            // TODO: вставить проверку на то, что колонка не изменилась, а значит, что пора останавливать функцию
        }
        return {current_row, current_column};
    }

    static Types::scalar norm_update(Types::scalar past_norm, const matrix_t &U, const matrix_t &V,
                                     const vector_t &u_new, const vector_t &v_new) {
        // Обновляем норму по аналитической формуле
        return std::sqrt(past_norm * past_norm + 2 * std::real((U.transpose() * u_new).dot(V * v_new)) +
                         u_new.squared_norm() * v_new.squared_norm());
    }

    static vector_t get_col_UV(const matrix_t &U, const matrix_t &V, Types::index j) {
        return U * V.row(j);
    }

    static vector_t get_row_UV(const matrix_t &U, const matrix_t &V, Types::index i) {
        return V * U.row(i);
    }

    template <typename MatrixElementFunction>
    static LinAgl::Matrix::DynamicFactoredMatrix<matrix_t> compute(const MatrixElementFunction &&element,
                                                                   Types::index n, Types::index m,
                                                                   Types::scalar error_control_parameter) {
        Containers::vector<Types::index> I;
        I.reserve(n);
        Containers::vector<Types::index> J;
        J.reserve(m);
        Containers::vector<Types::index> I_z;
        I_z.reserve(n);
        Containers::vector<Types::index> J_z;
        J_z.reserve(m);
        matrix_t U(n, 1), V(m, 1);

        std::iota(I.begin(), I.end(), 0);
        std::iota(J.begin(), J.end(), 0);

        // Первая итерация с поиском максимального элемента (TODO: семплинг нескольких начальных колонок и выбор
        // наибольшего элемента из них)
        auto [i, j] = get_argmax_of_matrix_elements<MatrixElementFunction>(element, m, n, m, 0);
        J_z.push_back(j);
        I_z.push_back(i);
        I.erase(std::remove(I.begin(), I.end(), i));
        J.erase(std::remove(I.begin(), I.end(), j));
        const value_t ij_element = element(i, j);
        // Добавление столбцов и строк в матрицы U и V (TODO: решить чето с памятью)
        U.col(0) = get_col(j, n, m, element) / std::abs(ij_element);
        V.col(0) = get_row(i, n, m, element) * std::abs(ij_element) / ij_element;

        // Основной цикл работы
        Types::index rank = 1;
        Types::scalar norm = U.col(0).norm() * V.col(0).norm();
        while (!stop_criterion(norm, error_control_parameter, n, m, rank, ij_element)) {
            rank += 1;
            // Поиск элемента "максимального объема" (maxvol для 1d)
            std::tie(i, j) = get_argmax_of_residual_matrix(element, n, m, U, V, J.front());
            J_z.push_back(j);
            I_z.push_back(i);
            I.erase(std::remove(I.begin(), I.end(), i));
            J.erase(std::remove(I.begin(), I.end(), j));
            // Вычисляем новые столбцы для матриц U и V
            const auto&& u = (get_col(j, n, m, element) - get_col_UV(U, V, j)) / std::abs(ij_element);
            const auto&& v = (get_row(i, n, m, element) - get_row_UV(U, V, i)) * std::abs(ij_element) / ij_element;
            // Апдейтим норму
            norm = norm_update(norm, U, V, u, v);
            // И дальше цикл заново, но сначала обновим матрицы
            U.conservativeResize(Eigen::NoChange_t::NoChange, U.cols() + 1);
            V.conservativeResize(Eigen::NoChange_t::NoChange, V.cols() + 1);
            U.col(U.cols() - 1) = u;
            V.col(V.cols() - 1) = v;
        }

        return {std::move(U), std::move(V.transpose())};
    }
};

}
#endif //ADAPTIVE_CROSS_HPP
