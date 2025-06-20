//
// Created by evgen on 18.05.2025.
//

#include "mesh/Parser.hpp"
#include "mesh/SurfaceMesh.hpp"

#include "math/matrix/iterative_solvers_coverage/MatrixReplacement.hpp"
#include "math/matrix/decompositions/Decompositions.hpp"

#include "geometry/PeriodicStructure.hpp"

#include "VTKFunctions.hpp"

#include "research/lattice/FieldCalculation.hpp"
#include "research/lattice/GeneralizedEquations.hpp"

#include <unsupported/Eigen/IterativeSolvers>

#include "experiment/PhysicalCondition.hpp"

#include <chrono>
#include <gtest/gtest.h>

#include <math/matrix/decompositions/Decompositions.hpp>

TEST(REAL_CASE, COMPARISON) {
    const std::string nodesFile = "/home/evgen/Education/MasterDegree/thesis/Electromagnetic-Waves-Scattering/meshes/"
                                  "lattice/8000_nodes.csv";
    const std::string cellsFile = "/home/evgen/Education/MasterDegree/thesis/Electromagnetic-Waves-Scattering/meshes/"
                                  "lattice/2000_cells.csv";

    // собираем сетки
    const auto parser_out = EMW::Parser::parseMesh(nodesFile, cellsFile);
    auto mesh_base = Mesh::SurfaceMesh{parser_out.first, parser_out.second};

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

    //const auto block = WaveGuideWithActiveSection::submatrix(mesh_2, mesh_1, a, k);
    const auto block = WaveGuideWithActiveSection::diagonal(mesh_1, a, k);

    auto end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

    std::cout << "Matrix assembling: " << elapsed.count() << " ms" << std::endl;

    // А теперь делаем крестовую аппроксимацию для этого блока

    // Вспомогательная функция для эмулирвоания функционального задания матрциы
    const auto element = [&block](Types::index i, Types::index j) { return block(i, j); };
    Types::scalar error_control = 0.1;

    start = std::chrono::system_clock::now();

    const auto out_factored_matrix = Math::LinAgl::Decompositions::ComplexACA::compute(element,
                                                                                       block.rows(),
                                                                                         block.cols(),
                                                                                         error_control);

    end = std::chrono::system_clock::now();
    elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

    std::cout << "Elapsed time: " << elapsed.count() << " ms" << std::endl;

    const Types::scalar err = (block - out_factored_matrix.compute()).norm();
    std::cout << "Ошибка приближения: " << err << std::endl;
    std::cout << "Rank = " << out_factored_matrix.get<0>().cols() << std::endl;

    std::cout << "Память для хранения общего блока: " << block.size() * 16. / (1024 * 1024 * 1024)
              << std::endl;
    std::cout << "Память для хранения фактора: " << out_factored_matrix.memory_usage() * 16. / (1024 * 1024 * 1024) << std::endl;

    // Проверка всего что можно проверить у матрицы с факторами
    ASSERT_EQ(out_factored_matrix.rows(), block.rows());
    ASSERT_EQ(out_factored_matrix.cols(), block.cols());
    ASSERT_EQ(out_factored_matrix.factor_number(), 2);

    // В том числе умножение
    const Types::VectorXc vec = Types::VectorXc::Random(block.rows());
    const Types::VectorXc res1 = block * vec;
    const Types::VectorXc res2 = out_factored_matrix * vec;

    std::cout << "Relative multiplication error: " << (res1 - res2).norm() / vec.norm() << "; error_control = " << error_control << std::endl;
}
