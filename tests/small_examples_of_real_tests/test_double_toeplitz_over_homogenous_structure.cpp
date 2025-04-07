//
// Created by evgen on 25.02.2025.
//

#include <gtest/gtest.h>

#include "types/Types.hpp"

#include "mesh/Parser.hpp"
#include "mesh/SurfaceMesh.hpp"

#include "geometry/PeriodicStructure.hpp"

#include "math/MathConstants.hpp"

#include "experiment/PhysicalCondition.hpp"

#include "research/lattice/GeneralizedEquations.hpp"

#include "Utils.hpp"

#include "equations/EquationsOverGeometry.hpp"

using namespace EMW;

class DOUBLE_TOEPLITZ_OVER_HOMOGENOUS_STRUCTURE : public ::testing::Test {
  protected:
    Mesh::SurfaceMesh mesh_base;
    const Types::scalar a = 0.07;
    const Types::scalar freq = Math::Constants::c / 1e8;
    const Types::complex_d k{Physics::get_k_on_frquency(freq), 0};

    // Блоки на первом уровне
    using ToeplitzBlock = Math::LinAgl::Matrix::ToeplitzBlock<Types::complex_d>;
    // Блоки на внутреннем уровне
    using block_t = ToeplitzBlock::block_type;

  public:
    void SetUp() override {
        const std::string cellsFile = "/home/evgen/Education/MasterDegree/thesis/Electromagnetic-Waves-Scattering/"
                                      "tests/meshes_for_tests/lattice_redused/80_cells.csv";
        const std::string nodesFile = "/home/evgen/Education/MasterDegree/thesis/Electromagnetic-Waves-Scattering/"
                                      "tests/meshes_for_tests/lattice_redused/320_nodes.csv";
        const EMW::Types::index nNodes = 320;
        const EMW::Types::index nCells = 80;

        // собираем сетки
        const auto parser_out = EMW::Parser::parseMesh(nodesFile, cellsFile, nNodes, nCells);
        mesh_base = EMW::Mesh::SurfaceMesh{parser_out.first, parser_out.second};
    };

    template <typename matrix_t>
    void final_check_for_vectors(const Types::MatrixX<Types::complex_d> &matrix, const matrix_t &toeplitz,
                                 const Types::VectorX<Types::complex_d> &vec) {
        auto start = std::chrono::system_clock::now();

        std::cout << "Норма вектора: " << vec.norm() << std::endl;

        const Types::VectorX<Types::complex_d> result_straightforward = matrix * vec;

        auto end = std::chrono::system_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        std::cout << "Прямое произведение: " << elapsed.count() << '\n';

        start = std::chrono::system_clock::now();

        const Types::VectorX<Types::complex_d> result_special = toeplitz * vec;

        end = std::chrono::system_clock::now();
        elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        std::cout << "Произведение с новой матрицей: " << elapsed.count() << '\n';

        // Сравниваем результаты умножения
        ASSERT_NEAR((result_straightforward - result_special).norm() / result_special.norm(), 0, 1e-14);
    }

    Types::MatrixXc get_inside_block(const Math::LinAgl::Matrix::ToeplitzToeplitzBlock<Types::complex_d> &mat,
                                     Types::index i, Types::index j, Types::index k, Types::index m) {
        return mat.get_block(i, j).get_block(k, m);
    }
};

template <Types::index N1, Types::index N2> using Scene = Geometry::PeriodicStructure<N1, N2>;

TEST_F(DOUBLE_TOEPLITZ_OVER_HOMOGENOUS_STRUCTURE, ASSEMBLING) {
    // Для сборки матрицы необходимо знать геометрию решетки, которую мы собрались считать
    constexpr Types::index N1 = 6;
    constexpr Types::index N2 = 7;
    const Geometry::PeriodicStructure<N1, N2> real_geometry{0.1, 0.1, mesh_base};

    // Теперь соберем ту структуру, по которой мы будем производить расчет
    constexpr auto n1 = 2 * N1 - 1;
    constexpr auto n2 = 2 * N2 - 1;
    // тут геометрия содержит сетки с индексами [0, n1 - 1 == 2 N1 - 2] x [0, n2 - 1 == 2 N2 - 2]
    const Geometry::PeriodicStructure<n1, n2> expanded_geometry{0.1, 0.1, mesh_base};

    // Тут обязательно есть сетка с origin == {0, 0, 0}, вот относительно неё и будем считать
    // Теперь нам нужно собрать все необходимые блоки в вектор и std::move-нуть его в нужное место
    // в структуре теплицевой матрицы

    constexpr auto first_layer_rows = N2;  // количество строк и столбцов
    constexpr auto first_layer_cols = N2;  // количество строк и столбцов
    constexpr auto second_layer_rows = N1; // количество строк и столбцов
    constexpr auto second_layer_cols = N1; // количество строк и столбцов

    // Собираем дважды тёплицеву матрицу
    Containers::vector<ToeplitzBlock> internal_blocks;
    // Собираем большой вектор из тёплицевых блоков
    internal_blocks.resize(ToeplitzBlock::get_size_of_container(second_layer_rows, second_layer_cols));
    // Находим сетку в середине, относительно которой и будем все считать
    constexpr int reference_row = N1 - 1;
    constexpr int reference_col = N2 - 1;
    const auto &reference_mesh = expanded_geometry.get(reference_row, reference_col);

    // заполняем горизонтальную часть вектора, отвечающего за внешнюю теплицевость
    for (int i = 0; i < second_layer_cols; i++) {
        Containers::vector<block_t> blocks;
        Types::index first_layer_size = ToeplitzBlock::get_size_of_container(first_layer_rows, first_layer_cols);
        blocks.resize(first_layer_size);
        int row_current_line = reference_row + i;
        // заполняем горизонтальную часть вектора, отвечающего за внутреннюю теплицевость
        for (int m = 0; m < first_layer_cols; m++) {
            if ((i == 0) && (m == 0)) {
                blocks[0] = std::move(WaveGuideWithActiveSection::diagonal(reference_mesh, a, k));
                std::cout << row_current_line << ' ' <<  reference_col + m << std::endl;
            } else {
                // а тут обычный расчет
                // надо найти сетку, по которой считаем
                int col = reference_col + m;
                const auto &another_mesh = expanded_geometry.get(row_current_line, col);
                std::cout << row_current_line << ' ' <<  col << std::endl;
                blocks[m] = std::move(WaveGuideWithActiveSection::submatrix(another_mesh, reference_mesh, a, k));
            }
        }
        // заполняем вертикальную часть вектора, отвечающего за внутреннюю теплицевость
        for (int m = 1; m < first_layer_rows; m++) {
            // тут снова обычный расчет
            // надо найти сетку, по которой считаем
            int col = reference_col - m;
            const int index_to_insert = first_layer_cols + m - 1;

            std::cout << row_current_line << ' ' <<  col << std::endl;
            const auto &another_mesh = expanded_geometry.get(row_current_line, col);
            blocks[index_to_insert] =
                std::move(WaveGuideWithActiveSection::submatrix(another_mesh, reference_mesh, a, k));
        }
        internal_blocks[i] = ToeplitzBlock{first_layer_rows, first_layer_cols, std::move(blocks)};
        // Проверка на то, что все корректно мувнулось куда нужно, а не скопировалось
        ASSERT_EQ(blocks.size(), 0);
    }
    // заполняем вертикальную часть вектора, отвечающего за внешнюю теплицевость
    for (int i = 0; i < second_layer_rows - 1; i++) {
        Containers::vector<block_t> blocks;
        Types::index first_layer_size = ToeplitzBlock::get_size_of_container(first_layer_rows, first_layer_cols);
        blocks.resize(first_layer_size);
        int row = reference_row - (i + 1);
        // заполняем горизонтальную часть вектора, отвечающего за внутреннюю теплицевость
        for (int m = 0; m < first_layer_cols; m++) {
            // а тут обычный расчет
            // надо найти сетку, по которой считаем
            int col = reference_col + m;
            const auto &another_mesh = expanded_geometry.get(row, col);
            std::cout << row << ' ' <<  col << std::endl;
            blocks[m] = std::move(WaveGuideWithActiveSection::submatrix(another_mesh, reference_mesh, a, k));
        }
        // заполняем вертикальную часть вектора, отвечающего за внутреннюю теплицевость
        for (int m = 1; m < first_layer_rows; m++) {
            // тут снова обычный расчет
            // надо найти сетку, по которой считаем
            int col = reference_col - m;
            const int index_to_insert = first_layer_cols + m - 1;
            std::cout << row << ' ' <<  col << std::endl;
            const auto &another_mesh = expanded_geometry.get(row, col);
            blocks[index_to_insert] =
                std::move(WaveGuideWithActiveSection::submatrix(another_mesh, reference_mesh, a, k));
        }
        internal_blocks[second_layer_cols + i] = ToeplitzBlock{first_layer_rows, first_layer_cols, std::move(blocks)};
        // Проверка на то, что все корректно мувнулось куда нужно, а не скопировалось
        ASSERT_EQ(blocks.size(), 0);
    }

    const Math::LinAgl::Matrix::ToeplitzToeplitzBlock test_matrix{second_layer_rows, second_layer_cols,
                                                                  std::move(internal_blocks)};

    ASSERT_NEAR((test_matrix.get_block(0, 0).get_block(0, 0) -
                            test_matrix.get_block(1,1).get_block(0, 0)).norm(), 0, 1e-14);

    // Проверка на то, что все корректно мувнулось куда нужно, а не скопировалось
    ASSERT_EQ(internal_blocks.size(), 0);

    // Проверка количества памяти для матрицы
    std::cout << Utils::get_memory_usage(test_matrix) << std::endl;

    // быстренько собираем матрицу, в которую свято верим
    const auto &full_matrix = Equations::MatrixForStructure::compute(
        real_geometry, WaveGuideWithActiveSection::diagonal, WaveGuideWithActiveSection::submatrix, {a, k});
    const auto& dense_test_matrix = test_matrix.to_dense();
    std::cout << full_matrix.rows()  << ' ' << full_matrix.cols() << std::endl;
    std::cout << dense_test_matrix.rows()  << ' ' << dense_test_matrix.cols() << std::endl;

#if 0    // Поблочная проверка матриц
    const Types::index internal_block_rows = get_inside_block(test_matrix, 0, 0, 0, 0).rows();
    const Types::index internal_block_cols = get_inside_block(test_matrix, 0, 0, 0, 0).cols();
    const Types::index first_layer_block_rows = internal_block_rows * first_layer_rows;
    const Types::index first_layer_block_cols = internal_block_cols * first_layer_cols;
    for (int i = 0; i < second_layer_rows; i++)
        for (int j = 0; j < second_layer_cols; j++)
            for (int k = 0; k < first_layer_rows; k++)
                for (int m = 0; m < first_layer_cols; m++) {
                    const auto &ref = full_matrix.block(i * first_layer_rows + k * internal_block_rows,
                                                        j * first_layer_cols + m * internal_block_cols,
                                                        internal_block_rows, internal_block_cols);
                    EXPECT_NEAR((ref - get_inside_block(test_matrix, i, j, k, m)).norm(), 0, 1e-12)
                        << i << " " << j << " " << k << " " << m;
                }
#endif
    // Проверка совпадения матриц
    const Types::scalar error = (full_matrix - dense_test_matrix).norm();
    ASSERT_NEAR(error, 0, 1e-10) << "Норма разницы матриц: " << error << std::endl;

#if 0  // Проверка умножения
    std::array<Types::complex_d, N1 * N2> phases{};
    phases.fill(Types::complex_d{1., 0.});
    const auto& vec = WaveGuideWithActiveSection::getRhs(phases, mesh_base, a, k);
    auto vec1 = Types::VectorXc{test_matrix.cols()};
    for (int i = 0; i < vec1.rows(); i++) vec1[i] = 1;
    Types::VectorXc vec2{test_matrix.cols()};
    vec2.setRandom();

    ASSERT_EQ(vec.size(), test_matrix.cols());

    final_check_for_vectors(full_matrix, test_matrix, vec);
    final_check_for_vectors(full_matrix, test_matrix, vec1);
    final_check_for_vectors(full_matrix, test_matrix, vec2);
#endif // Проверка умножения
}

#include "math/matrix/iterative_solvers_coverage//MatrixTraits.hpp"
#include "math/matrix/iterative_solvers_coverage/MatrixReplacement.hpp"

#include <unsupported/Eigen/IterativeSolvers>

template <typename MatrixType>
Types::VectorXc solve(const MatrixType &A, const Types::VectorXc &b, Types::scalar tolerance) {

    Eigen::GMRES<MatrixType, Eigen::IdentityPreconditioner> method;

    Types::index max_iterations = 10000;

    method.setMaxIterations(max_iterations);
    method.setTolerance(tolerance);
    method.set_restart(max_iterations);

    Types::VectorXc x{};

    auto start = std::chrono::high_resolution_clock::now();

    method.compute(A);
    x = method.solve(b);

    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed_seconds = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "Время решения GMRES: " << elapsed_seconds.count() << std::endl;
    std::cout << "total iterations: " << method.iterations() << std::endl;
    std::cout << "tolerance: " << method.tolerance() << std::endl;
    std::cout << "total error: " << method.error() << std::endl;
    std::cout << "Info: " << static_cast<int>(method.info()) << std::endl;
    return x;
}

TEST_F(DOUBLE_TOEPLITZ_OVER_HOMOGENOUS_STRUCTURE, GMRES_TESTING) {
    // Для сборки матрицы необходимо знать геометрию решетки, которую мы собрались считать
    constexpr Types::index N1 = 2;
    constexpr Types::index N2 = 2;
    const Geometry::PeriodicStructure<N1, N2> real_geometry{0.1, 0.1, mesh_base};

    // Теперь соберем ту структуру, по которой мы будем производить расчет
    constexpr auto n1 = 2 * N1 - 1;
    constexpr auto n2 = 2 * N2 - 1;
    // тут геометрия содержит сетки с индексами [0, n1 - 1 == 2 N1 - 2] x [0, n2 - 1 == 2 N2 - 2]
    const Geometry::PeriodicStructure<n1, n2> expanded_geometry{0.1, 0.1, mesh_base};

    // Тут обязательно есть сетка с origin == {0, 0, 0}, вот относительно неё и будем считать
    // Теперь нам нужно собрать все необходимые блоки в вектор и std::move-нуть его в нужное место
    // в структуре теплицевой матрицы

    constexpr auto first_layer_rows = N2;  // количество строк и столбцов
    constexpr auto first_layer_cols = N2;  // количество строк и столбцов
    constexpr auto second_layer_rows = N1; // количество строк и столбцов
    constexpr auto second_layer_cols = N1; // количество строк и столбцов

    // Собираем дважды тёплицеву матрицу
    Containers::vector<ToeplitzBlock> internal_blocks;
    internal_blocks.reserve(ToeplitzBlock::get_size_of_container(second_layer_rows, second_layer_cols));
    // Собираем большой вектор из тёплицевых блоков

    // Находим сетку в середине, относительно которой и будем все считать
    constexpr auto reference_row = N1 - 1;
    constexpr auto reference_col = N2 - 1;
    const auto &reference_mesh = expanded_geometry.get(reference_row, reference_col);

    // заполняем горизонтальную часть вектора, отвечающего за внешнюю теплицевость

    std::vector<std::pair<int, int>> indexes;
    indexes.reserve(n1 * n2);

    for (int i = 0; i < second_layer_cols; i++) {
        Containers::vector<block_t> blocks;
        Types::index first_layer_size = ToeplitzBlock::get_size_of_container(first_layer_cols, first_layer_cols);
        blocks.reserve(first_layer_size);
        // заполняем горизонтальную часть вектора, отвечающего за внутреннюю теплицевость
        for (int m = 0; m < first_layer_cols; m++) {
            if ((i == 0) && (m == 0)) {
                int row = reference_row;
                int col = reference_col;
                indexes.emplace_back(row, col);

                blocks.emplace_back(WaveGuideWithActiveSection::diagonal(reference_mesh, a, k));
            } else {
                // а тут обычный расчет
                // надо найти сетку, по которой считаем
                int row = reference_row + i;
                int col = reference_col + m;
                indexes.emplace_back(row, col);

                const auto &another_mesh = expanded_geometry.get(row, col);
                blocks.emplace_back(WaveGuideWithActiveSection::submatrix(another_mesh, reference_mesh, a, k));
            }
        }
        // заполняем вертикальную часть вектора, отвечающего за внутреннюю теплицевость
        for (int m = first_layer_rows - 1; m >= 1; m--) {
            // тут снова обычный расчет
            // надо найти сетку, по которой считаем
            int row = reference_row + i;
            int col = reference_col - m;
            indexes.emplace_back(row, col);

            const auto &another_mesh = expanded_geometry.get(row, col);
            blocks.emplace_back(WaveGuideWithActiveSection::submatrix(another_mesh, reference_mesh, a, k));
        }
        internal_blocks.emplace_back(first_layer_rows, first_layer_cols, std::move(blocks));
        // Проверка на то, что все корректно мувнулось куда нужно, а не скопировалось
        ASSERT_EQ(blocks.size(), 0);
    }
    // заполняем вертикальную часть вектора, отвечающего за внешнюю теплицевость
    for (int i = 0; i < second_layer_rows - 1; i++) {
        Containers::vector<block_t> blocks;
        Types::index first_layer_size = ToeplitzBlock::get_size_of_container(first_layer_cols, first_layer_cols);
        blocks.reserve(first_layer_size);
        // заполняем горизонтальную часть вектора, отвечающего за внутреннюю теплицевость
        for (int m = 0; m < first_layer_cols; m++) {
            // а тут обычный расчет
            // надо найти сетку, по которой считаем
            int row = i;
            int col = reference_col + m;
            indexes.emplace_back(row, col);

            const auto &another_mesh = expanded_geometry.get(row, col);
            blocks.emplace_back(WaveGuideWithActiveSection::submatrix(another_mesh, reference_mesh, a, k));
        }
        // заполняем вертикальную часть вектора, отвечающего за внутреннюю теплицевость
        for (int m = first_layer_rows - 1; m >= 1; m--) {
            // тут снова обычный расчет
            // надо найти сетку, по которой считаем
            int row = i;
            int col = reference_col - m;
            indexes.emplace_back(row, col);

            const auto &another_mesh = expanded_geometry.get(row, col);
            blocks.emplace_back(WaveGuideWithActiveSection::submatrix(another_mesh, reference_mesh, a, k));
        }
        internal_blocks.emplace_back(first_layer_rows, first_layer_cols, std::move(blocks));
        // Проверка на то, что все корректно мувнулось куда нужно, а не скопировалось
        ASSERT_EQ(blocks.size(), 0);
    }

    const Math::LinAgl::Matrix::ToeplitzToeplitzBlock test_matrix
                            {second_layer_rows, second_layer_cols, std::move(internal_blocks)};

    // Проверка на то, что все корректно мувнулось куда нужно, а не скопировалось
    ASSERT_EQ(internal_blocks.size(), 0);

    // Проверка правильности прохода по индексам
    ASSERT_EQ(indexes.size(), n1 * n2);
    // Проверка того, что каждый индекс встречается один раз
    for (int i = 0; i < n1 * n2; i++) {
        ASSERT_EQ(indexes.begin() + i, std::find_if(indexes.begin(), indexes.end(),
                    [value = indexes[i]](const auto& v) {return v == value;}));
    }

    // Проверка количества памяти для матрицы
    std::cout << Utils::get_memory_usage(test_matrix) << std::endl;

    std::array<Types::complex_d, N1 * N2> phases{};
    phases.fill(Types::complex_d{1., 0.});
    const auto& vec = WaveGuideWithActiveSection::getRhs(phases, mesh_base, a, k);

    ASSERT_EQ(vec.size(), test_matrix.cols());
    // быстренько собираем матрицу, в которую свято верим
    const auto& full_matrix = Equations::MatrixForStructure::compute(
        real_geometry, WaveGuideWithActiveSection::diagonal, WaveGuideWithActiveSection::submatrix, {a, k});

    // Дальше подключаем работу с GMRES
    Types::scalar tolerance = 1e-2;
    // Решаем с полной матрицей
    const auto res_full = solve(full_matrix, vec, tolerance);
    // Решаем с тёплицевой структурой
    const auto res_toeplitz = solve(Math::LinAgl::Matrix::Wrappers::MatrixReplacement{test_matrix}, vec, tolerance);
    // Сравниваем решения между собой
    ASSERT_NEAR((res_full - res_toeplitz).norm() / res_full.norm(), 0, 1e-13);
}
