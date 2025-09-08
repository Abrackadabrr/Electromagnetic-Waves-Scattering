//
// Created by evgen on 15.05.2025.
//

#include <gtest/gtest.h>

#include "math/matrix/decompositions/Decompositions.hpp"

#include "mat_decomp.hpp"

TEST_F(MATRIX_DECOMPOSITIONS_TESTS, ACA_AUX_TEST1) {
    const Types::index N = 1000;
    const Types::index M = 1200;
    const auto mat = get_sinus_mat(N, M);
    const auto function_for_element = [&mat](const Types::index i, const Types::index j) { return mat(i, j); };

    // Совпадение стоблцов
    for (int i = 0; i < M; i++) {
        const auto &&col = Math::LinAgl::Decompositions::RealACA::get_col(i, N, M, function_for_element);
        ASSERT_NEAR((col - mat.col(i)).norm(), 0, 1e-14);
    }

    // Совпадение строк
    for (int i = 0; i < N; i++) {
        const auto &&row = Math::LinAgl::Decompositions::RealACA::get_row(i, N, M, function_for_element);
        ASSERT_NEAR((row - mat.row(i).transpose()).norm(), 0, 1e-14);
    }

    // 1. Проверка расчета строчек и столбцов от произведения матриц
    const Types::MatrixXd UVT = mat * (mat.transpose()); // матрицы N x N

    // Совпадение стоблцов (всего в итоговой матрице N стоблцов)
    for (int i = N; i < N; i++) {
        const auto &&col = Math::LinAgl::Decompositions::RealACA::get_col_UV(mat, mat, i);
        ASSERT_NEAR((col - UVT.col(i)).norm() / col.norm(), 0, 1e-14);
    }

    // Совпадение строк
    for (int i = 0; i < N; i++) {
        const auto &&row = Math::LinAgl::Decompositions::RealACA::get_row_UV(mat, mat, i);
        ASSERT_NEAR((row - UVT.row(i).transpose()).norm() / row.norm(), 0, 1e-14);
    }

    // 2. Вычисление нормы матрицы
}

TEST_F(MATRIX_DECOMPOSITIONS_TESTS, ACA_SINUS) {
    const Types::index N = 4400;
    const Types::MatrixXd mat = get_sinus_mat(N);

    const auto function_for_element = [&mat](const Types::index i, const Types::index j) { return mat(i, j); };

    const Types::scalar error = 1;
    const auto factored_matrix = Math::LinAgl::Decompositions::RealACA::compute(function_for_element, N, N, error);

    std::cout << "Rank = " << factored_matrix.get<0>().cols() << std::endl;

    const Types::scalar err = (factored_matrix.compute() - mat).norm() / mat.norm();
    std::cout << err << std::endl;
    std::cout << "Full mem usage: " << N * N * 8. / (1024 * 1024) << " Mb" << std::endl;
    std::cout << "Compressed mem usage: " << factored_matrix.memory_usage() * 8 / (1024 * 1024) << " Mb" << std::endl;
}

TEST_F(MATRIX_DECOMPOSITIONS_TESTS, ACA_HILBERT) {
    const Types::index N = 4400;
    const Types::MatrixXd mat = get_hilbert_mat(N);

    const auto function_for_element = [&mat](const Types::index i, const Types::index j) { return mat(i, j); };

    const Types::scalar error = 1e-5;
    const auto factored_matrix = Math::LinAgl::Decompositions::RealACA::compute(function_for_element, N, N, error);

    std::cout << "Rank = " << factored_matrix.get<0>().cols() << std::endl;

    const Types::scalar err = (factored_matrix.compute() - mat).norm() / mat.norm();
    std::cout << "Relative error: " << err << std::endl;
    std::cout << "Full mem usage: " << N * N * 8. / (1024 * 1024) << " Mb" << std::endl;
    std::cout << "Compressed mem usage: " << factored_matrix.memory_usage() * 8 / (1024 * 1024) << " Mb" << std::endl;
}

#include "mesh/Parser.hpp"
#include "mesh/SurfaceMesh.hpp"

#include "math/matrix/iterative_solvers_coverage/MatrixReplacement.hpp"

#include "geometry/PeriodicStructure.hpp"

#include "VTKFunctions.hpp"

#include "research/lattice/FieldCalculation.hpp"
#include "research/lattice/SpecificLatticeEquations.hpp"

#include <unsupported/Eigen/IterativeSolvers>

#include "experiment/PhysicalCondition.hpp"

#include <chrono>

TEST_F(MATRIX_DECOMPOSITIONS_TESTS, ACA_REAL_CASE_TEST) {
    const std::string nodesFile = "/home/evgen/Education/MasterDegree/thesis/Electromagnetic-Waves-Scattering/meshes/"
                                  "lattice/8000_nodes.csv";
    const std::string cellsFile = "/home/evgen/Education/MasterDegree/thesis/Electromagnetic-Waves-Scattering/meshes/"
                                  "lattice/2000_cells.csv";

    // собираем сетки
    const auto parser_out = EMW::Parser::parse_mesh_without_tag(nodesFile, cellsFile);
    auto mesh_base = Mesh::SurfaceMesh{parser_out.nodes, parser_out.cells};

    // сделаем ещё одну сетку для анализа внедиагонлаьных элементов
    const auto &mesh_1 = mesh_base;

    const Types::scalar between = 4 * 0.07;
    const Types::Vector3d origin{0.07 + between, 0, 0.0};
    const auto mesh_2 = Mesh::Utils::move_by_vector(mesh_base, origin);

    // Геометрические параметры антенн
    // Короткая сторона волновода
    const Types::scalar a = 0.07;
    // Физика волны в пространстве
    // частота в гигагерцах
    const Types::scalar freq = 2 * Math::Constants::c / 1e8;
    const Types::complex_d k{Physics::get_k_on_frquency(freq), 0};

    // Собираем матрицу и смотрим время
    auto start = std::chrono::system_clock::now();

    const auto out_diagonal_block = WaveGuideWithActiveSection::submatrix(mesh_2, mesh_1, a, k);
    // const auto diagonal_block = WaveGuideWithActiveSection::diagonal(mesh_1, a, k);

    auto end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

    std::cout << "Matrix assembling: " << elapsed.count() << " ms" << std::endl;

    // А теперь делаем крестовую аппроксимацию для этого блока

    // Вспомогательная функция для эмулирвоания функционального задания матрциы
    const auto element = [&out_diagonal_block](Types::index i, Types::index j) { return out_diagonal_block(i, j); };

    const auto compute_col = [&out_diagonal_block](Types::index i) -> Types::VectorXc {
        return out_diagonal_block.col(i);
    };
    const auto compute_row = [&out_diagonal_block](Types::index i) -> Types::VectorXc {
        return out_diagonal_block.row(i);
    };

    Types::scalar error_control = 0.1;

    start = std::chrono::system_clock::now();
#if 1
    const auto out_factored_matrix = Math::LinAgl::Decompositions::ComplexACA::compute(
        compute_row, compute_col, out_diagonal_block.rows(), out_diagonal_block.cols(), error_control);
#else
    const auto out_factored_matrix = Math::LinAgl::Decompositions::ComplexACA::compute(
        element, out_diagonal_block.rows(), out_diagonal_block.cols(), error_control);
#endif
    end = std::chrono::system_clock::now();
    elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

    std::cout << "Elapsed time: " << elapsed.count() << " ms" << std::endl;

    const Types::scalar err = (out_diagonal_block - out_factored_matrix.compute()).norm() / out_diagonal_block.norm();
    std::cout << "Relative error: " << err << std::endl;
    std::cout << "Rank = " << out_factored_matrix.get<0>().cols() << std::endl;

    std::cout << "Память для хранения общего блока: " << out_diagonal_block.size() * 16. / (1024 * 1024 * 1024)
              << std::endl;
    std::cout << "Память для хранения фактора: " << out_factored_matrix.memory_usage() * 16. / (1024 * 1024 * 1024)
              << std::endl;

    // Проверка всего что можно проверить у матрицы с факторами
    ASSERT_EQ(out_factored_matrix.rows(), out_diagonal_block.rows());
    ASSERT_EQ(out_factored_matrix.cols(), out_diagonal_block.cols());
    ASSERT_EQ(out_factored_matrix.factor_number(), 2);

    // В том числе умножение
    const Types::VectorXc vec = Types::VectorXc::Random(out_diagonal_block.rows());
    const Types::VectorXc res1 = out_diagonal_block * vec;
    const Types::VectorXc res2 = out_factored_matrix * vec;

    std::cout << "Relative error of matvec: " << (res1 - res2).norm() / vec.norm() << "; error_control = " << error_control << std::endl;
}
