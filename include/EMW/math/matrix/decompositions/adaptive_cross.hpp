//
// Created by evgen on 15.05.2025.
//

#ifndef ADAPTIVE_CROSS_HPP
#define ADAPTIVE_CROSS_HPP

#include "math/matrix/DynamicFactoredMatrix.hpp"
#include "types/Types.hpp"

#include <Eigen/QR>
#include <Eigen/SVD>

#include <numeric>
#include <random>

namespace EMW::Math::LinAgl::Decompositions {
template <typename matrix_t, typename vector_t, typename value_t> struct ACA {

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
        // std::cout << uv_prev_norm << std::endl;
        // std::cout << e << std::endl;
        // std::cout << n << std::endl;
        // std::cout << m << std::endl;
        // std::cout << rank << std::endl;
        // std::cout << ij_element << std::endl;
        return e * uv_prev_norm >= std::abs(ij_element) * std::sqrt((m - rank) * (n - rank));
    }

    /*
     * Вычислить строку в матрице
     */
    template <typename MatrixElementFunction>
    static vector_t get_row(Types::index i, Types::index n, Types::index m, MatrixElementFunction &&element_function) {
        vector_t result = vector_t::Zero(m);
#pragma omp parallel for num_threads(14)
        for (Types::index j = 0; j < m; ++j) {
            result[j] = std::forward<MatrixElementFunction>(element_function)(i, j);
        }
        return result;
    }

    /*
     * Вычислить столбец в матрице
     */
    template <typename MatrixElementFunction>
    static vector_t get_col(Types::index j, Types::index n, Types::index m, MatrixElementFunction &&element_function) {
        vector_t result = vector_t::Zero(n);
#pragma omp parallel for num_threads(14)
        for (Types::index i = 0; i < n; ++i) {
            result[i] = std::forward<MatrixElementFunction>(element_function)(i, j);
        }
        return result;
    }

    /*
     * Вычислить приближение к максимальному элементу в исходной матрице
     * по методу адаптивной крестовой аппроксимации
     *
     * compute_row/compute_col принимает в качестве аргумента только индекс строки/столбца
     */
    template <typename MatrixRowFunction, typename MatrixColumnFunсtion>
    static maxvol_solution iterate_pure_matrix(MatrixRowFunction &&compute_row, MatrixColumnFunсtion &&compute_col,
                                               Types::index start_column) {
        constexpr Types::index iterations = 2; // это подкручиваемый параметр
        // (на самом деле есть теорема о том, что 2 раза почти всегда достаточно)
        Types::index current_column = start_column;
        Types::index current_row = 0;
        value_t max_element_value = 0;
        vector_t row;
        vector_t col;
        for (int i = 0; i < iterations; ++i) {
            // Проходка вдоль столбца
            col = std::forward<MatrixColumnFunсtion>(compute_col)(current_column);
            const Types::VectorXd col_abs = col.cwiseAbs();
            // Нашли строку максимального элемента в столбце
            current_row = std::ranges::max_element(col_abs) - col_abs.begin();

            // Проходка вдоль строки
            row = std::forward<MatrixRowFunction>(compute_row)(current_row);
            const Types::VectorXd row_abs = row.cwiseAbs();
            // Нашли столбец максимального элемента в строке
            const Types::index new_column = std::ranges::max_element(row_abs) - row_abs.begin();

            // Проверка на то, что крест остановился раньше
            if (new_column == current_column) {
                return {current_row, current_column, row(current_column), std::move(row), std::move(col)};
            }
            // А если не остановился, то продолжаем итерацию с новой колонки
            current_column = new_column;
        }

        return {current_row, current_column, row(current_column), std::move(row), std::move(col)};
    }

    // --------- Функции для поиска хорошего креста в остаточной матрице А - U_r V_r^T
    template <typename MatrixU, typename MatrixV>
    static vector_t get_col_UV(const MatrixU &U, const MatrixV &V, Types::index j) {
        return U * V.row(j).transpose();
    }

    template <typename MatrixU, typename MatrixV>
    static vector_t get_row_UV(const MatrixU &U, const MatrixV &V, Types::index i) {
        return V * U.row(i).transpose();
    }

    /**
     *  Поиск максимального элемента в новой матрице, то есть в матрице,
     *  скорректированной на матрицы U и V
     *
     *  Пока что неясно нужна она или нет
     */
    template <typename MatrixRowFunction, typename MatrixColumnFunсtion, typename MatrixU, typename MatrixV>
    static maxvol_solution iterate_residual_matrix(MatrixRowFunction &&compute_row, MatrixColumnFunсtion &&compute_col,
                                                   const MatrixU &U, const MatrixV &V, Types::index start_column) {
        constexpr Types::index iterations = 2; // это подкручиваемый параметр
        Types::index current_column = start_column;
        Types::index current_row = 0;
        value_t max_element_value;
        vector_t row;
        vector_t col;
        for (int i = 0; i < iterations; ++i) {
            // Проходка вдоль столбца
            col = std::forward<MatrixColumnFunсtion>(compute_col)(current_column) - get_col_UV(U, V, current_column);
            const Types::VectorXd col_abs = col.cwiseAbs();
            // Нашли строку максимального элемента в столбце
            current_row = std::ranges::max_element(col_abs) - col_abs.begin();

            // Проходка вдоль строки
            row = std::forward<MatrixRowFunction>(compute_row)(current_row) - get_row_UV(U, V, current_row);
            const Types::VectorXd row_abs = row.cwiseAbs();
            // Нашли столбец максимального элемента в строке
            const Types::index new_column = std::ranges::max_element(row_abs) - row_abs.begin();

            // Проверка на то, что крест остановился раньше
            if (new_column == current_column) {
                return {current_row, current_column, row(current_column), std::move(row), std::move(col)};
            }
            // А если не остановился, то продолжаем итерацию с новой колонки
            current_column = new_column;
        }
        return {current_row, current_column, row(current_column), std::move(row), std::move(col)};
    }

    template <typename MatrixU, typename MatrixV>
    static Types::scalar squared_norm_update(Types::scalar past_norm_square, const MatrixU &U, const MatrixV &V,
                                             const vector_t &u_new, const vector_t &v_new) {
        const Types::scalar u_new_sq_norm = u_new.squaredNorm();
        const Types::scalar v_new_sq_norm = v_new.squaredNorm();
        // Обновляем норму по аналитической формуле
        return (past_norm_square + 2 * std::real((u_new.transpose() * U).dot((v_new.transpose() * V))) + u_new_sq_norm * v_new_sq_norm);
    }

    template <typename MatrixRowFunction, typename MatrixColumnFunсtion>
    static Matrix::DynamicFactoredMatrix<matrix_t> compute(MatrixRowFunction &&compute_row,
                                                           MatrixColumnFunсtion &&compute_col, Types::index n,
                                                           Types::index m, Types::scalar error_control_parameter) {
        // странный параметер, который получается из соображения, что если матрица размером больше,
        // то посчитать её полностью будет дешевле, чем строить крест
        // А ещё это помогает не релоцировать память для матриц евери тайм
        Types::index max_rank = std::min(n, m) / 8;

        // Начало алгоритма
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
        constexpr int N_SAMPLES = 2; // как раз эта константа
        std::array<Types::index, N_SAMPLES> initial_cols{};

        std::random_device rd;  // a seed source for the random number engine
        std::mt19937 gen(rd()); // mersenne_twister_engine seeded with rd()
        std::uniform_int_distribution<> distrib(0, m - 1);
        initial_cols[0] = distrib(gen);
        initial_cols[1] = distrib(gen);

        // Расчет максвола для каждой их этих колонок
        Containers::set<maxvol_solution> initial_samples;
        for (auto &&start_col : initial_cols)
            initial_samples.emplace(iterate_pure_matrix(compute_row, compute_col, start_col));

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

        // Основной цикл работы
        Types::index rank = 1;
        Types::scalar sq_norm = U.col(0).squaredNorm() * V.col(0).squaredNorm();
        while (!stop_criterion(std::sqrt(sq_norm), error_control_parameter, n, m, rank, ij_element)) {
            // Поиск элемента "максимального объема" (maxvol для 1d)
            const auto maxvol_sol = iterate_residual_matrix(compute_row, compute_col, U,
                                                            V, J.front());
            i = maxvol_sol.i;
            j = maxvol_sol.j;
            ij_element = maxvol_sol.value;
            if (std::abs(ij_element) < std::numeric_limits<Types::scalar>::epsilon())
                break; // TODO: перестроить цикл так, чтобы убрать этот лишний if
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
            U.conservativeResize(Eigen::NoChange, U.cols() + 1);
            V.conservativeResize(Eigen::NoChange, V.cols() + 1);
            U.col(rank) = u;
            V.col(rank) = v;
            rank += 1;
        }
        // вот и кончился алгоритм
        return {Containers::vector<matrix_t>{std::move(U), matrix_t{V.transpose()}}, Containers::vector{false, false}};
    }

    template<typename MatrixElementFunction>
    static Matrix::DynamicFactoredMatrix<matrix_t> compute(MatrixElementFunction &&compute_element,
                                                           Types::index n, Types::index m,
                                                           Types::scalar error_control_parameter) {
        const auto compute_row = [&](Types::index i) { return get_row(i, n, m, std::forward<MatrixElementFunction>(compute_element)); };
        const auto compute_col = [&](Types::index j) { return get_col(j, n, m, std::forward<MatrixElementFunction>(compute_element)); };

        return compute(compute_row, compute_col, n, m, error_control_parameter);
    };

    static Matrix::DynamicFactoredMatrix<matrix_t> svd_postcompression(
        Matrix::DynamicFactoredMatrix<matrix_t>&& matrix, Types::scalar tolerance) {
        if (matrix.factor_number() != 2) {
            return std::move(matrix);
        }
        // tolerance = tolerance / matrix.template get<0>().cols();
        const matrix_t U = matrix.template get<0>();
        const matrix_t V_transposed = matrix.template get<1>(); // stored as V^T

        if (U.cols() != V_transposed.rows()) {
            throw std::invalid_argument("ACA::svd_postcompression expects factors U (n x r) and V^T (r x m)");
        }

        if (U.cols() == 0) {
            return std::move(matrix);
        }

        const matrix_t V = V_transposed.transpose();
        const Types::index rank = U.cols();

        // U = Q1 * R1
        Eigen::HouseholderQR<matrix_t> qr_u(U);
        const matrix_t Q1 = qr_u.householderQ() * matrix_t::Identity(U.rows(), rank);
        const matrix_t R1 = qr_u.matrixQR().topRows(rank).template triangularView<Eigen::Upper>();

        // V = Q2 * R2, where input stores V^T
        Eigen::HouseholderQR<matrix_t> qr_v(V);
        const matrix_t Q2 = qr_v.householderQ() * matrix_t::Identity(V.rows(), rank);
        const matrix_t R2 = qr_v.matrixQR().topRows(rank).template triangularView<Eigen::Upper>();

        // Since A = U * V^T = Q1 * R1 * R2^T * Q2^T
        const matrix_t compressed_core = R1 * R2.transpose();

        Eigen::BDCSVD<matrix_t> svd(compressed_core, Eigen::ComputeThinU | Eigen::ComputeThinV);
        const Types::VectorXd singular_values = svd.singularValues();

        if (singular_values.size() == 0) {
            return std::move(matrix);
        }

        const Types::scalar safe_tolerance = std::max<Types::scalar>(0, tolerance);
        const Types::scalar sv_threshold = safe_tolerance * singular_values(0);

        Types::index truncated_rank = 0;
        for (Types::index i = 0; i < singular_values.size(); ++i) {
            if (singular_values(i) >= sv_threshold) {
                truncated_rank += 1;
            } else {
                break;
            }
        }
        if (truncated_rank == 0) {
            truncated_rank = 1;
        }

        const matrix_t P = svd.matrixU().leftCols(truncated_rank);
        const matrix_t H = svd.matrixV().leftCols(truncated_rank);

        const auto S =
            singular_values.head(truncated_rank).template cast<value_t>().asDiagonal();

        const matrix_t left_factor = Q1 * P * S;
        // A = (Q1 * P * S) * (H* * Q2^T), i.e. factorized as U_new * V_new^T
        const matrix_t right_factor_transposed = H.adjoint() * Q2.transpose();

        Containers::vector<matrix_t> factors;
        factors.reserve(2);
        factors.emplace_back(std::move(left_factor));
        factors.emplace_back(std::move(right_factor_transposed));
        return {std::move(factors), Containers::vector<bool>{false, false}};
    }
};
}
#endif //ADAPTIVE_CROSS_HPP
