//
// Created by evgen on 18.03.2025.
//

#include "math/matrix/decompositions/Decompositions.hpp"
#include "types/Types.hpp"

#include <gtest/gtest.h>

using namespace EMW;

class RSVD_TESTS: public ::testing::Test {
protected:
    Types::scalar sin_mat_el(Types::index i, Types::index j) {
        return std::sin(i + j);
    }

    Types::MatrixXd get_sinus_mat(Types::index N) {
        Types::MatrixXd sin_mat = Types::MatrixXd::Zero(N, N);
        for (Types::index i = 0; i < N; i++) {
            for (Types::index j = 0; j < N; j++) {
                sin_mat(i, j) = sin_mat_el(i, j);
            }
        }
        return sin_mat;
    }
};

TEST_F(RSVD_TESTS, TEST_SINUS_MATRIX) {
    const Types::index N = 1000;
    const auto mat = get_sinus_mat(N);

    const auto rsvd = Math::LinAgl::Decompositions::RealRSVD::compute(mat, 2, 2);

#if 0
    std::cout << rsvd.matrixQ() * rsvd.matrixU() << std::endl;
    std::cout << Types::DiagonalMatrixXd(rsvd.singularValues())  << std::endl;
    std::cout << rsvd.matrixVh() << std::endl;

    std::cout << rsvd.matrixQ() * rsvd.matrixU() * Types::DiagonalMatrixXd(rsvd.singularValues()) * rsvd.matrixVh() << std::endl;
#endif

    const auto& factored_matrix = rsvd;
    const Types::scalar err =  (factored_matrix.compute() - mat).norm() / mat.norm();
    ASSERT_NEAR(err, 0, 1e-13);
    std::cout << err << std::endl;
}

#include "mesh/Parser.hpp"
#include "mesh/SurfaceMesh.hpp"

#include "math/matrix/iterative_solvers_coverage/MatrixReplacement.hpp"

#include "geometry/PeriodicStructure.hpp"

#include "VTKFunctions.hpp"

#include "research/lattice/FieldCalculation.hpp"
#include "research/lattice/GeneralizedEquations.hpp"

#include <unsupported/Eigen/IterativeSolvers>

#include "experiment/PhysicalCondition.hpp"

#include <chrono>

TEST_F(RSVD_TESTS, REAL_CASE_TEST) {
    const std::string nodesFile = "/home/evgen/Education/MasterDegree/thesis/Electromagnetic-Waves-Scattering/meshes/"
                                  "lattice/8000_nodes.csv";
    const std::string cellsFile = "/home/evgen/Education/MasterDegree/thesis/Electromagnetic-Waves-Scattering/meshes/"
                                  "lattice/2000_cells.csv";
    constexpr EMW::Types::index nNodes = 8000;
    constexpr EMW::Types::index nCells = 2000;

    // собираем сетки
    const auto parser_out = EMW::Parser::parseMesh(nodesFile, cellsFile, nNodes, nCells);
    auto mesh_base = Mesh::SurfaceMesh{parser_out.first, parser_out.second};

    // сделаем ещё одну сетку для анализа внедиагонлаьных элементов
    const auto &mesh_1 = mesh_base;

    const Types::scalar between = 0.07;
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

    // А теперь делаем рандомизированное svd для внедиагонального блока
    start = std::chrono::system_clock::now();

    const Types::index rank = 80 * (a / between);
    std::cout << "RSVD with rank = " << rank << std::endl;
    const auto out_factored_matrix = Math::LinAgl::Decompositions::ComplexRSVD::compute(out_diagonal_block, rank, rank, 2);

    end = std::chrono::system_clock::now();
    elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

    std::cout << "Elapsed time: " << elapsed.count() << " ms" << std::endl;

    const Types::scalar err = (out_diagonal_block - out_factored_matrix.compute()).norm();
    std::cout << err << std::endl;
}
