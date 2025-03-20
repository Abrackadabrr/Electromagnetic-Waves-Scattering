//
// Created by evgen on 12.02.2025.
//

#include "math/matrix/Matrix.hpp"
#include "toeplitz_matrix_tests.hpp"
#include "types/Types.hpp"

#include <Utils.hpp>

Types::scalar simple_toeplitz_matrix_1(Types::index i, Types::index j) {
    return (static_cast<scalar>(i) - static_cast<scalar>(j) + 10);
};

Types::scalar simple_toeplitz_matrix_2(Types::index i, Types::index j) {
    return (static_cast<scalar>(i) - static_cast<scalar>(j) - 10);
};

/**
 * Функция, которая рассчитывает один БЛОК.
 * Матрица, построенная с помощью таких блоков будет ОДИН РАЗ ТЁПЛИЦЕВА
 * @param rows
 * @param cols
 * @return
 */
Types::MatrixX<scalar> get_toeplitz_matrix(Types::index rows, Types::index cols) {
    Types::MatrixX<scalar> result = Types::MatrixX<scalar>::Zero(rows, cols);
    for (Types::index i = 0; i < rows; i++)
        for (Types::index j = 0; j < cols; j++)
            result(i, j) = simple_toeplitz_matrix_1(i, j);
    return result;
}

// Дальше идут функции, которые описывают псевдо дважды блочную тёплицевость

Types::MatrixX<scalar> get_internal_block(Types::index i, Types::index j) {
    Types::MatrixX<scalar> result = Types::MatrixX<scalar>::Zero(1, 1);
    result(0, 0) = simple_toeplitz_matrix_1(i, j);
    return result;
}

template <typename matrix_t>
void final_check_for_vectors(const Types::MatrixX<scalar>& matrix, const matrix_t &toeplitz,
                             const Types::VectorX<scalar> &vec) {
    auto start = std::chrono::system_clock::now();

    const Types::VectorX<scalar> result_straightforward = matrix * vec;

    auto end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "Прямое произведение: " << elapsed.count() << '\n';

    start = std::chrono::system_clock::now();

    const Types::VectorX<scalar> result_special = toeplitz * vec;

    end = std::chrono::system_clock::now();
    elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "Произведение с новой матрицей: " << elapsed.count() << '\n';

    // Сравниваем результаты умножения
    ASSERT_NEAR((result_straightforward - result_special).norm(), 0, 0);
}

/**
 * Обычная тёплиц-блочная матрица тоже является дважды блочно теплицевой
 * Этот тест проверяет тривиальную корректность
 *
 * Здесь очень плохой пример с точки зрения времени работы: прямое умножение в 80 раз быстрее
 */
TEST_F(TOEPLITZ_MATRIX_TESTS, TWICE_TOEPLITZ_SIMPLE) {
    // Описываем структуру матрицы
    const Types::index internal_block_rows = 1;
    const Types::index internal_block_cols = 1;

    const Types::index first_layer_rows = 20;
    const Types::index first_layer_cols = 20;

    const Types::index second_layer_rows = 1;
    const Types::index second_layer_cols = 1;

    const Types::index total_rows = internal_block_rows * first_layer_rows * second_layer_rows;
    const Types::index total_cols = internal_block_cols * first_layer_cols * second_layer_cols;

    // Строим обычнцую матрицу, без специального хранения элементов
    const auto toeplitz_matrix = get_toeplitz_matrix(total_rows, total_cols, internal_block_cols);
    // Собираем дважды тёплицеву матрицу
    const Math::LinAgl::Matrix::ToeplitzToeplitzBlock<complex_d> test_matrix(
        second_layer_rows, second_layer_cols, [&first_layer_cols, &first_layer_rows](Types::index i, Types::index j) {
            return Math::LinAgl::Matrix::ToeplitzBlock<complex_d>(first_layer_rows, first_layer_cols, get_internal_block);
        });
    // Собираем вектор для тестового умножения
    Types::VectorX<scalar> vec = Types::VectorX<scalar>::Zero(total_cols);
    for (Types::index i = 0; i < total_cols; i++)
        vec(i) = static_cast<scalar>(i);

    final_check_for_vectors(toeplitz_matrix, test_matrix, vec);
}

// А теперь настоящие дважды тёплицевы тесты

// Пусть теперь внутренние блоки будут одинаковые, притом с заданным размером
// Как и положено, блоки нумеруются четырьмя индексами, притом зависят на самом деле только от их разности
Types::MatrixX<complex_d> get_internal_block(Types::index rows, Types::index cols, Types::index i, Types::index j,
                                          Types::index k, Types::index m) {
    Types::MatrixX<scalar> result = Types::MatrixX<scalar>::Zero(rows, cols);
    for (Types::index p = 0; p < cols; p++)
        for (Types::index t = 0; t < rows; t++)
            result(t, p) = simple_toeplitz_matrix_1(i, j) * simple_toeplitz_matrix_2(k, m);
    return result;
}

Types::MatrixX<complex_d> get_full_matrix(Types::index s_r, Types::index s_c, Types::index f_r, Types::index f_c,
                                       Types::index i_r, Types::index i_c) {
    const Types::index primary_block_rows = f_r * i_r;
    const Types::index primary_block_cols = f_c * i_c;

    const Types::index total_rows = s_r * primary_block_rows;
    const Types::index total_cols = s_c * primary_block_cols;

    Types::MatrixX<complex_d> result(total_rows, total_cols);

    // Заполняем матрицу
    for (Types::index i = 0; i < s_r; i++)
        for (Types::index j = 0; j < s_c; j++) {
            for (Types::index k = 0; k < f_r; k++)
                for (Types::index l = 0; l < f_c; l++) {
                    result.block(i * primary_block_rows, j * primary_block_cols, primary_block_rows, primary_block_cols)
                        .block(k * i_r, l * i_c, i_r, i_c) = get_internal_block(i_r, i_c, i, j, k, l);
                }
        }
    return result;
}

TEST_F(TOEPLITZ_MATRIX_TESTS, TWICE_TOEPLITZ_REAL_1_NO_PARALLEL) {
    Eigen::setNbThreads(1);
    // Описываем структуру матрицы
    const Types::index internal_block_rows = 10;
    const Types::index internal_block_cols = 20;

    const Types::index first_layer_rows = 20;
    const Types::index first_layer_cols = 10;

    const Types::index second_layer_rows = 10;
    const Types::index second_layer_cols = 10;

    const Types::index total_rows = internal_block_rows * first_layer_rows * second_layer_rows;
    const Types::index total_cols = internal_block_cols * first_layer_cols * second_layer_cols;

    // Строим обычнцую матрицу, без специального хранения элементов
    const Types::MatrixX<complex_d> full_matrix =
        get_full_matrix(second_layer_rows, second_layer_cols, first_layer_rows, first_layer_cols, internal_block_rows,
                        internal_block_cols);
    // std::cout << full_matrix << std::endl;

    // Собираем дважды тёплицеву матрицу
    const auto get_toeplitz_block_by_indexes = [&internal_block_rows, &internal_block_cols, &first_layer_rows,
                                                &first_layer_cols](Types::index i, Types::index j) {
        const auto get_internal_block_by_indexes = [&internal_block_rows, &internal_block_cols, &i,
                                                    &j](Types::index k, Types::index l) {
            return get_internal_block(internal_block_rows, internal_block_cols, i, j, k, l);
        };
        return Math::LinAgl::Matrix::ToeplitzBlock<complex_d>{first_layer_rows, first_layer_cols,
                                                           get_internal_block_by_indexes};
    };
    const auto test_matrix = Math::LinAgl::Matrix::ToeplitzToeplitzBlock<complex_d>{second_layer_rows, second_layer_cols,
                                                                                 get_toeplitz_block_by_indexes};

    // Собираем вектор для тестового умножения
    Types::VectorX<scalar> vec = Types::VectorX<scalar>::Zero(total_cols);
    for (Types::index i = 0; i < total_cols; i++)
        vec(i) = static_cast<scalar>(i);
    // Сравниваем результаты умножения

    final_check_for_vectors(full_matrix, test_matrix, vec);

    // Проверка совпадения матриц
    ASSERT_NEAR((full_matrix - test_matrix.to_dense()).norm(), 0, 1e-14);
}

TEST_F(TOEPLITZ_MATRIX_TESTS, TWICE_TOEPLITZ_REAL_2_NO_PARALLEL) {
    Eigen::setNbThreads(1);

    // Описываем структуру матрицы
    const Types::index internal_block_rows = 100;
    const Types::index internal_block_cols = 100;

    const Types::index first_layer_rows = 7;
    const Types::index first_layer_cols = 7;

    const Types::index second_layer_rows = 7;
    const Types::index second_layer_cols = 7;

    const Types::index total_rows = internal_block_rows * first_layer_rows * second_layer_rows;
    const Types::index total_cols = internal_block_cols * first_layer_cols * second_layer_cols;

    const Types::scalar element_in_gb = 16. / (1024 * 1024 * 1024);
    const Types::scalar memory_for_full_matrix = total_cols * total_rows * element_in_gb;
    const Types::scalar memory_for_toeplitz_matrix = internal_block_cols * internal_block_rows *
                                                     (2 * second_layer_rows - 1) * (2 * second_layer_cols - 1) *
                                                     element_in_gb;

    std::cout << "Количество памяти для полной матрицы: " << memory_for_full_matrix << " Gb" << std::endl;
    std::cout << "Количество памяти для новой матрицы: " << memory_for_toeplitz_matrix << " Gb" << std::endl;

    // Строим обычнцую матрицу, без специального хранения элементов
    const Types::MatrixX<complex_d> full_matrix =
        get_full_matrix(second_layer_rows, second_layer_cols, first_layer_rows, first_layer_cols, internal_block_rows,
                        internal_block_cols);

    // Собираем дважды тёплицеву матрицу
    const auto get_toeplitz_block_by_indexes = [&internal_block_rows, &internal_block_cols, &first_layer_rows,
                                                &first_layer_cols](Types::index i, Types::index j) {
        const auto get_internal_block_by_indexes = [&internal_block_rows, &internal_block_cols, &i,
                                                    &j](Types::index k, Types::index l) {
            return get_internal_block(internal_block_rows, internal_block_cols, i, j, k, l);
        };
        return Math::LinAgl::Matrix::ToeplitzBlock<complex_d>{first_layer_rows, first_layer_cols,
                                                           get_internal_block_by_indexes};
    };
    const auto test_matrix = Math::LinAgl::Matrix::ToeplitzToeplitzBlock<complex_d>{second_layer_rows, second_layer_cols,
                                                                                 get_toeplitz_block_by_indexes};

    // Собираем вектор для тестового умножения
    Types::VectorX<scalar> vec = Types::VectorX<scalar>::Random(total_cols);

    // Проверка умножения с измерением времени
    final_check_for_vectors(full_matrix, test_matrix, vec);

    // Проверка совпадения матриц
    ASSERT_NEAR((full_matrix - test_matrix.to_dense()).norm(), 0, 1e-14);
}

/**
 * Попробуем в этих тестах собрать матрицу по-другому
 * то есть через конструктор, который принимает вектор значений сразу
 */
using ToeplitzBlock = Math::LinAgl::Matrix::ToeplitzBlock<complex_d>;
using block_t = ToeplitzBlock::block_type;

TEST_F(TOEPLITZ_MATRIX_TESTS, TWICE_TOEPLITZ_NEW_CONSTRUCT) {
    // Eigen::setNbThreads(14);
    // Описываем структуру матрицы
    const Types::index internal_block_rows = 2400;
    const Types::index internal_block_cols = 2400;

    const Types::index first_layer_rows = 3;
    const Types::index first_layer_cols = first_layer_rows;

    const Types::index second_layer_rows = 3;
    const Types::index second_layer_cols = second_layer_rows;

    const Types::index total_rows = internal_block_rows * first_layer_rows * second_layer_rows;
    const Types::index total_cols = internal_block_cols * first_layer_cols * second_layer_cols;

    // Строим обычнцую матрицу, без специального хранения элементов
    const Types::MatrixX<complex_d> full_matrix =
        get_full_matrix(second_layer_rows, second_layer_cols, first_layer_rows, first_layer_cols, internal_block_rows,
                        internal_block_cols);

    // Собираем дважды тёплицеву матрицу
    Containers::vector<ToeplitzBlock> internal_blocks;
    internal_blocks.reserve(ToeplitzBlock::get_size_of_container(second_layer_rows, second_layer_cols));
    // Собираем большой вектор из тёплицевых блоков
    // Этот сбор опирается на то, как именно хранится вектор,
    // который задает тёплицеву матрицу (сначала хранится --, затем |)
    for (int i = 0; i < second_layer_cols; i++) {
        Containers::vector<block_t> blocks;
        Types::index first_layer_size = ToeplitzBlock::get_size_of_container(first_layer_cols, first_layer_cols);
        blocks.reserve(first_layer_size);
        for (int k = 0; k < first_layer_cols; k++) {
            // а тут обычный расчет
            blocks.emplace_back(get_internal_block(internal_block_rows, internal_block_cols, 0, i, 0, k));
        }
        for (int k = 1; k < first_layer_rows; k++) {
            // тут снова обычный расчет
            blocks.emplace_back(get_internal_block(internal_block_rows, internal_block_cols, 0, i, k, 0));
        }
        internal_blocks.emplace_back(first_layer_rows, first_layer_cols, std::move(blocks));
        ASSERT_EQ(blocks.size(), 0);
    }
    for (int i = 1; i < second_layer_rows; i++) {
        Containers::vector<block_t> blocks;
        Types::index first_layer_size = ToeplitzBlock::get_size_of_container(first_layer_cols, first_layer_cols);
        blocks.reserve(first_layer_size);
        for (int k = 0; k < first_layer_cols; k++) {
            blocks.emplace_back(get_internal_block(internal_block_rows, internal_block_cols, i, 0, 0, k));
        }
        for (int k = 1; k < first_layer_rows; k++) {
            blocks.emplace_back(get_internal_block(internal_block_rows, internal_block_cols, i, 0, k, 0));
        }
        internal_blocks.emplace_back(first_layer_rows, first_layer_cols, std::move(blocks));
        ASSERT_EQ(blocks.size(), 0);
    }

    const Math::LinAgl::Matrix::ToeplitzToeplitzBlock test_matrix
                            {second_layer_rows, second_layer_cols, std::move(internal_blocks)};

    ASSERT_EQ(internal_blocks.size(), 0);

    // Собираем вектор для тестового умножения
    Types::VectorX<scalar> vec = Types::VectorX<scalar>::Zero(total_cols);
    for (Types::index i = 0; i < total_cols; i++)
        vec(i) = static_cast<scalar>(i);

    // Проверка умножения с измерением времени
    final_check_for_vectors(full_matrix, test_matrix, vec);

    // Проверка методов rows_in_block и rows()
    ASSERT_EQ(test_matrix.cols(), total_cols);
    ASSERT_EQ(test_matrix.rows(), total_rows);
    ASSERT_EQ(test_matrix.cols_in_block(), internal_block_rows * first_layer_rows);
    ASSERT_EQ(test_matrix.rows_in_block(), internal_block_cols * first_layer_cols);

    std::cout << "Количество памяти для полной матрицы: " << Utils::get_memory_usage(test_matrix).full_matrix << " Gb" << std::endl;
    std::cout << "Количество памяти для новой матрицы: " << Utils::get_memory_usage(test_matrix).toeplitz_matrix << " Gb" << std::endl;
}

/** Тест на доступ к элементам матрицы */
TEST_F(TOEPLITZ_MATRIX_TESTS, TWICE_TOEPLITZ_ELEMENTS_ACCESS) {
    // Описываем структуру матрицы
    const Types::index internal_block_rows = 400;
    const Types::index internal_block_cols = 400;

    const Types::index first_layer_rows = 3;
    const Types::index first_layer_cols = first_layer_rows;

    const Types::index second_layer_rows = 4;
    const Types::index second_layer_cols = second_layer_rows;

    const Types::index total_rows = internal_block_rows * first_layer_rows * second_layer_rows;
    const Types::index total_cols = internal_block_cols * first_layer_cols * second_layer_cols;

    // Строим обычнцую матрицу, без специального хранения элементов
    const Types::MatrixX<complex_d> full_matrix =
        get_full_matrix(second_layer_rows, second_layer_cols, first_layer_rows, first_layer_cols, internal_block_rows,
                        internal_block_cols);

    // Собираем дважды тёплицеву матрицу
    Containers::vector<ToeplitzBlock> internal_blocks;
    internal_blocks.reserve(ToeplitzBlock::get_size_of_container(second_layer_rows, second_layer_cols));
    // Собираем большой вектор из тёплицевых блоков
    // Этот сбор опирается на то, как именно хранится вектор,
    // который задает тёплицеву матрицу (сначала хранится --, затем |)
    for (int i = 0; i < second_layer_cols; i++) {
        Containers::vector<block_t> blocks;
        Types::index first_layer_size = ToeplitzBlock::get_size_of_container(first_layer_cols, first_layer_cols);
        blocks.reserve(first_layer_size);
        for (int k = 0; k < first_layer_cols; k++) {
            // а тут обычный расчет
            blocks.emplace_back(get_internal_block(internal_block_rows, internal_block_cols, 0, i, 0, k));
        }
        for (int k = 1; k < first_layer_rows; k++) {
            // тут снова обычный расчет
            blocks.emplace_back(get_internal_block(internal_block_rows, internal_block_cols, 0, i, k, 0));
        }
        internal_blocks.emplace_back(first_layer_rows, first_layer_cols, std::move(blocks));
        ASSERT_EQ(blocks.size(), 0);
    }
    for (int i = 1; i < second_layer_rows; i++) {
        Containers::vector<block_t> blocks;
        Types::index first_layer_size = ToeplitzBlock::get_size_of_container(first_layer_cols, first_layer_cols);
        blocks.reserve(first_layer_size);
        for (int k = 0; k < first_layer_cols; k++) {
            blocks.emplace_back(get_internal_block(internal_block_rows, internal_block_cols, i, 0, 0, k));
        }
        for (int k = 1; k < first_layer_rows; k++) {
            blocks.emplace_back(get_internal_block(internal_block_rows, internal_block_cols, i, 0, k, 0));
        }
        internal_blocks.emplace_back(first_layer_rows, first_layer_cols, std::move(blocks));
        ASSERT_EQ(blocks.size(), 0);
    }

    const Math::LinAgl::Matrix::ToeplitzToeplitzBlock test_matrix
                            {second_layer_rows, second_layer_cols, std::move(internal_blocks)};

    ASSERT_EQ(internal_blocks.size(), 0);

    // Проверка методов rows_in_block и rows()
    ASSERT_EQ(test_matrix.cols(), total_cols);
    ASSERT_EQ(test_matrix.rows(), total_rows);
    ASSERT_EQ(test_matrix.cols_in_block(), internal_block_rows * first_layer_rows);
    ASSERT_EQ(test_matrix.rows_in_block(), internal_block_cols * first_layer_cols);

    // Проверка оператора доступа к элементам
    for (Types::index i = 0; i < total_rows; i++)
        for (Types::index k = 0; k < total_cols; k++)
            ASSERT_EQ(test_matrix(i, k), full_matrix(i, k)) << i << " " << k << std::endl;
}
