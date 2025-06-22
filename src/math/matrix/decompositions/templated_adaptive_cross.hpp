//
// Created by evgen on 22.06.2025.
//

#ifndef TEMPLATED_ADAPTIVE_CROSS_HPP
#define TEMPLATED_ADAPTIVE_CROSS_HPP


#include "math/matrix/DynamicFactoredMatrix.hpp"
#include "types/Types.hpp"

#include <bits/random.h>
#include <numeric>

namespace EMW::Math::LinAgl::Decompositions {
template <typename Callable, typename matrix_t, typename vector_t, typename value_t> struct T_ACA {
#ifndef NDEBUG
    /*
     * Вычислить полную матрицy
     */
    template <typename MatrixElementFunction>
    static matrix_t get_matrix(Types::index n, Types::index m, const MatrixElementFunction &element_function) {
        matrix_t result = matrix_t::Zero(n, m);
        for (int row = 0; row < n; ++row)
            for (int j = 0; j < m; ++j) {
                result(row, j) = element_function(row, j);
            }
        return result;
    }
#endif

    struct maxvol_solution {
        Types::index i, j;
        value_t value;
        vector_t row;
        vector_t col;

        bool operator<(const maxvol_solution &other) const { return std::abs(value) < std::abs(other.value); }
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
                            MatrixElementFunction &&element_function) {
        vector_t result = vector_t::Zero(m);
        for (Types::index j = 0; j < m; ++j) {
            result[j] = element_function(i, j);
        }
        return result;
    }

    /*
     * Вычислить столбец в матрице
     */
    template <typename MatrixElementFunction>
    static vector_t get_col(Types::index j, Types::index n, Types::index m,
                            MatrixElementFunction &&element_function) {
        vector_t result = vector_t::Zero(n);
        for (Types::index i = 0; i < n; ++i) {
            result[i] = element_function(i, j);
        }
        return result;
    }

    /*
     * Вычислить приближение к максимальному элементу в матрице
     */
    template <typename MatrixElementFunction>
    static maxvol_solution get_argmax_of_matrix_elements(MatrixElementFunction &&element_function, Types::index n,
                                                         Types::index m, Types::index start_column) {
        Types::index iterations = 4; // это подкручиваемый параметр
        Types::index current_column = start_column;
        Types::index current_row = 0;
        value_t max_element_value = 0;
        vector_t row;
        vector_t col;
        for (int i = 0; i < iterations; ++i) {
            // Проходка вдоль столбца
            col = get_col(current_column, n, m, element_function);
            const Types::VectorXd col_abs = col.cwiseAbs();
            // Нашли строку максимального элемента в столбце
            current_row = std::max_element(col_abs.begin(), col_abs.end()) - col_abs.begin();

            // Проходка вдоль строки
            row = get_row(current_row, n, m, element_function);
            const Types::VectorXd row_abs = row.cwiseAbs();
            // Нашли столбец максимального элемента в строке
            current_column = std::max_element(row_abs.begin(), row_abs.end()) - row_abs.begin();

            // TODO: вставить проверку на то, что колонка не изменилась, а значит, что пора останавливать функцию
        }

        return {current_row, current_column, row(current_column), row, col};
    }

    static vector_t get_col_UV(const matrix_t &U, const matrix_t &V, Types::index j) {
        return U * (V.row(j).transpose());
    }

    static vector_t get_row_UV(const matrix_t &U, const matrix_t &V, Types::index i) {
        return V * (U.row(i).transpose());
    }

    /**
     *  Поиск максимального элемента в новой матрице, то есть в матрице,
     *  скорректированной на матрицы U и V
     *
     *  Пока что функция не написана
     */
    template <typename MatrixElementFunction>
    static maxvol_solution get_argmax_of_residual_matrix(MatrixElementFunction &&element_function, Types::index n,
                                                         Types::index m, const matrix_t &U, const matrix_t &V,
                                                         Types::index start_column) {
        int iterations = 4; // это подкручиваемый параметр
        Types::index current_column = start_column;
        Types::index current_row = 0;
        value_t max_element_value;
        vector_t row;
        vector_t col;
        for (int i = 0; i < iterations; ++i) {
            // Проходка вдоль столбца
            col = get_col(current_column, n, m, element_function) - get_col_UV(U, V, current_column);
            const Types::VectorXd col_abs = col.cwiseAbs();
            // Нашли строку максимального элемента в столбце
            current_row = std::max_element(col_abs.begin(), col_abs.end()) - col_abs.begin();

            // Проходка вдоль строки
            row = get_row(current_row, n, m, element_function) - get_row_UV(U, V, current_row);
            const Types::VectorXd row_abs = row.cwiseAbs();
            // Нашли столбец максимального элемента в строке
            current_column = std::max_element(row_abs.begin(), row_abs.end()) - row_abs.begin();

            // TODO: вставить проверку на то, что колонка не изменилась, а значит, что пора останавливать функцию
        }
        return {current_row, current_column, row(current_column), row, col};
    }

    static Types::scalar squared_norm_update(Types::scalar past_norm_square, const matrix_t &U, const matrix_t &V,
                                             const vector_t &u_new, const vector_t &v_new) {
        const vector_t vec1 = (U.transpose() * u_new);
        // Обновляем норму по аналитической формуле
        return (past_norm_square + 2 * std::real(vec1.dot(V.transpose() * v_new)) +
                u_new.squaredNorm() * v_new.squaredNorm());
    }

    template <typename MatrixElementFunction>
    static Matrix::DynamicFactoredMatrix<matrix_t> compute(MatrixElementFunction &&element, Types::index n,
                                                           Types::index m, Types::scalar error_control_parameter) {
        Containers::vector<Types::index> I;
        I.resize(n);
        Containers::vector<Types::index> J;
        J.resize(m);
        Containers::vector<Types::index> I_z;
        I_z.reserve(n);
        Containers::vector<Types::index> J_z;
        J_z.reserve(m);
        matrix_t U(n, 1), V(m, 1);

        std::iota(I.begin(), I.end(), 0);
        std::iota(J.begin(), J.end(), 0);

        // Первая итерация с поиском максимального элемента путем семплинга нескольких колонок
        // со случайной нумерацией элементов, так будет более вероятно найти максимальный элемент
        // в матрице (решение задачи максвола с k = 1)

        // 1) Семплинг некоторого случайного количества колонок
        const int N_SAMPLES = 6; // как раз эта константа
        std::array<Types::index, N_SAMPLES> initial_cols;
        initial_cols[0] = 0;
        std::random_device rd;  // a seed source for the random number engine
        std::mt19937 gen(rd()); // mersenne_twister_engine seeded with rd()
        std::uniform_int_distribution<> distrib(0, m);
        std::generate(initial_cols.begin() + 1, initial_cols.end(), [&distrib, &gen]() {return distrib(gen);});
        // Расчет максвола для каждой их этих колонок
        Containers::set<maxvol_solution> initial_samples;
        for (auto &&start_col : initial_cols)
            initial_samples.emplace(get_argmax_of_matrix_elements<MatrixElementFunction>(element, n, m, start_col));

        // 2) Выбор наилучшего приближения для максимального элемента
        auto [i, j, ij_element, row, col] = *initial_samples.rbegin();
        // 3) Работа уже с этим приближением
        J_z.push_back(j);
        I_z.push_back(i);
        I.erase(std::remove(I.begin(), I.end(), i));
        J.erase(std::remove(J.begin(), J.end(), j));
        // Добавление столбцов и строк в матрицы U и V (TODO: решить чето с памятью)
        U.col(0) = col / std::abs(ij_element);
        V.col(0) = row * std::abs(ij_element) / ij_element;

#ifndef NDEBUG
        std::cout << "Rank = " << 1 << "; residual matrix is \n"
                  << get_matrix(n, m, element) - U * V.transpose() << std::endl;
        std::cout << "U\n" << U << std::endl;
        std::cout << "V^T\n" << V.transpose() << std::endl;
        std::cout << "Element was " << i << ' ' << j << " with value:" << ij_element << std::endl;
        std::cout << "// ------- //" << std::endl;
#endif

        // Основной цикл работы
        Types::index rank = 1;
        Types::scalar sq_norm = U.col(0).squaredNorm() * V.col(0).squaredNorm();
        while (!stop_criterion(std::sqrt(sq_norm), error_control_parameter, n, m, rank, ij_element)) {
            // Поиск элемента "максимального объема" (maxvol для 1d)
            const auto maxvol_sol = get_argmax_of_residual_matrix(element, n, m, U, V, J.front());
            i = maxvol_sol.i;
            j = maxvol_sol.j;
            ij_element = maxvol_sol.value;
            if (std::abs(ij_element) < std::numeric_limits<Types::scalar>::epsilon()) {
                break;
            } // TODO: перестроить цикл так, чтобы убрать этот лишний if
            J_z.push_back(j);
            I_z.push_back(i);
            I.erase(std::remove(I.begin(), I.end(), i));
            J.erase(std::remove(J.begin(), J.end(), j));
            // Вычисляем новые столбцы для матриц U и V
            const vector_t u = maxvol_sol.col / std::abs(ij_element);
            const vector_t v = maxvol_sol.row * std::abs(ij_element) / ij_element;
            // Апдейтим норму
            sq_norm = squared_norm_update(sq_norm, U, V, u, v);
            // И дальше цикл заново, но сначала обновим матрицы
            U.conservativeResize(Eigen::NoChange_t::NoChange, U.cols() + 1);
            V.conservativeResize(Eigen::NoChange_t::NoChange, V.cols() + 1);
            U.col(U.cols() - 1) = u;
            V.col(V.cols() - 1) = v;
            rank += 1;

#ifndef NDEBUG
            std::cout << "Rank = " << rank << "; residual matrix is \n" << get_matrix(n, m, element) - U * V.transpose() << std::endl;
            std::cout << "U\n" << U << std::endl;
            std::cout << "V^T\n" << V.transpose() << std::endl;
            std::cout << "Element was " << i << ' ' << j << " with value:" << ij_element << std::endl;
            std::cout << "// ------- //" << std::endl;
#endif
        }

        return {Containers::vector{std::move(U), std::move(V)}, Containers::vector{false, true}};
    }
};

}

#endif //TEMPLATED_ADAPTIVE_CROSS_HPP
