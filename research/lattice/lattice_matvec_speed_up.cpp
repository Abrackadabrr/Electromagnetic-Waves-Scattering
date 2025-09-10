//
// Created by evgen on 01.04.2025.
//

#include <string>

#include "mesh/MeshTypes.hpp"
#include "mesh/Parser.hpp"
#include "mesh/SurfaceMesh.hpp"

#include "math/fields/SurfaceVectorField.hpp"
#include "math/matrix/iterative_solvers_coverage/DiagonalPreconditioner.hpp"
#include "math/matrix/iterative_solvers_coverage/MatrixReplacement.hpp"
#include "math/matrix/iterative_solvers_coverage/MatrixTraits.hpp"

#include "research/Solve.hpp"

#include "geometry/PeriodicStructure.hpp"

#include "VTKFunctions.hpp"

#include "FieldCalculation.hpp"
#include "FieldOverGeometry.hpp"
#include "GeneralEquation.hpp"

#include "math/matrix/Matrix.hpp"

#include <unsupported/Eigen/IterativeSolvers>

#include "experiment/PhysicalCondition.hpp"

#include <Utils.hpp>
#include <chrono>
#include <experiment/ESA.hpp>
#include <iostream>

namespace eq = WaveGuideWithActiveSection;

template <int N1, int N2> using Scene = Geometry::PeriodicStructure<N1, N2>;

namespace LAMatrix = Math::LinAgl::Matrix;
// каким методом расчитываем матрицу
constexpr Research::Lattice::CalculationMethod calc_method = Research::Lattice::CalculationMethod::ACA;
using TTBMatrix = Research::Lattice::CalcTraits<calc_method>::ReturnType;
// Предобуславливатели для матрицы
using DiagonalPrec = LAMatrix::Preconditioning::DiagonalPreconditioner<Types::complex_d, TTBMatrix>;
using BlockDiagPrec = LAMatrix::Preconditioning::BlockDiagonalPreconditioner<Types::complex_d, TTBMatrix>;
using NoPrec = LAMatrix::Preconditioning::IdentityPreconditioner<Types::complex_d, TTBMatrix>;
using MatrixWrapper = LAMatrix::Wrappers::MatrixReplacement<TTBMatrix, BlockDiagPrec>;

template <typename matrix_t, typename matrix_full_t>
void final_check_for_vectors(const matrix_full_t &matrix, const matrix_t &toeplitz,
                             const Types::VectorX<Types::complex_d> &vec) {

    Types::index n_repeat = 5;
    Types::scalar time_1 = 0;
    Types::scalar time_2 = 0;

    for (int i = 0; i < n_repeat; i++) {
        auto start = std::chrono::system_clock::now();

        const Types::VectorX<Types::complex_d> result_straightforward = matrix * vec;

        auto end = std::chrono::system_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        time_1 += elapsed.count();
    }
    std::cout << "Прямое произведение: " << time_1 / n_repeat << '\n';

    for (int i = 0; i < n_repeat; i++) {
        auto start = std::chrono::system_clock::now();

        const Types::VectorX<Types::complex_d> result_special = toeplitz * vec;

        auto end = std::chrono::system_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        time_2 += elapsed.count();
    }
    std::cout << "Произведение с новой матрицей: " << time_2 / n_repeat << '\n';

    const Types::VectorX<Types::complex_d> result_straightforward = matrix * vec;
    const Types::VectorX<Types::complex_d> result_special = toeplitz * vec;

    // Сравниваем результаты умножения
    std::cout << "Относительная норма ошибки: "
              << (result_straightforward - result_special).norm() / result_special.norm();
}

int main() {
    // считываем сетку на антенне
    const std::string nodesFile = "/home/evgen/Education/MasterDegree/thesis/Electromagnetic-Waves-Scattering/meshes/"
                                  "lattice/8000_nodes.csv";
    const std::string cellsFile = "/home/evgen/Education/MasterDegree/thesis/Electromagnetic-Waves-Scattering/meshes/"
                                  "lattice/2000_cells.csv";

    // собираем сетки
    const auto parser_out = EMW::Parser::parse_mesh_without_tag(nodesFile, cellsFile);
    auto mesh_base = Mesh::SurfaceMesh{parser_out.nodes, parser_out.cells};

    constexpr Types::index N1 = 3; // количество строк в решетке
    constexpr Types::index N2 = 4; // количество столбцов в решетке
    constexpr Types::index N1_x_N2 = N1 * N2;
    constexpr Types::scalar d1 = 0.07; // расстояние между строками
    constexpr Types::scalar d2 = 0.14; // расстояние между столбцами
    const Types::Vector3d dir1 = Types::Vector3d{0, -1, 0}.normalized();
    const Types::Vector3d dir2 = Types::Vector3d{1, 0, 0}.normalized();
    const Geometry::PeriodicStructure<N1, N2> geometry{dir1, dir2, d1, d2, mesh_base};

    // Геометрические параметры антенн
    // Короткая сторона волновода
    const Types::scalar a = 0.07;
    // Физика волны в пространстве
    // частота в гигагерцах
    const Types::scalar freq = Math::Constants::c / 1e8;
    const Types::complex_d k{Physics::get_k_on_frquency(freq), 0};
    // расчет коэффициента импеданса
    const Types::complex_d beta = std::sqrt(k * k - (EMW::Math::Constants::PI_square<Types::scalar>() / (a * a)));
    std::cout << "Волновое число в волноводе: " << beta.real()
              << "; Длина волны в волноводе: " << 2 * Math::Constants::PI<Types::scalar>() / beta.real() << std::endl;
    std::cout << "Волновое число в свободном пространстве: " << k.real()
              << "; Длина волны в свободном пространстве: " << 2 * Math::Constants::PI<Types::scalar>() / k.real()
              << std::endl;

    // замеряем время на сбор аппроксимации
    Types::scalar approximation_calculation_time = 0;
    Types::index n_repeats = 20;

    for (int i = 0; i < n_repeats; i++) {
        auto start = std::chrono::steady_clock::now();
        const auto matrix = Research::Lattice::getMatrix<calc_method>(geometry, a, k);

        auto end = std::chrono::steady_clock::now();
        approximation_calculation_time += std::chrono::duration_cast<std::chrono::seconds>(end - start).count();
        std::cout << i << " iteration, " << matrix.rows() << std::endl;
    }

    const auto matrix = Research::Lattice::getMatrix<calc_method>(geometry, a, k);

    std::cout << Utils::get_memory_usage(matrix) << std::endl;
    std::cout << "Matrix assembled, size: " << matrix.rows()
              << ";  average calculation time: " << approximation_calculation_time / n_repeats << std::endl;

    auto start = std::chrono::steady_clock::now();

    const auto matrix_full = Research::Lattice::getMatrix<Research::Lattice::CalculationMethod::Full>(geometry, a, k);

    auto end = std::chrono::steady_clock::now();
    auto elapsed_full = std::chrono::duration_cast<std::chrono::seconds>(end - start);
    std::cout << "Full Matrix assembled, size: " << matrix_full.rows() << "; time elapsed: " << elapsed_full
              << std::endl;

    Types::scalar err_mat_norm = 0;
    Types::scalar mat_norm = 0;

    for (int i = 0; i < matrix.rows_in_toeplitrz(); i++)
        for (int j = 0; j < matrix.cols_in_toeplitrz(); j++) {
            const auto block_f = matrix_full.get_block(i, j).to_dense();
            const auto block_comp = matrix.get_block(i, j).to_dense();
            err_mat_norm += (block_f - block_comp).norm();
            mat_norm += block_f.norm();
        }

    std::cout << "Relative error is: " << err_mat_norm / mat_norm << std::endl;

    final_check_for_vectors(matrix_full, matrix, Types::VectorX<Types::complex_d>::Random(matrix_full.cols()));

}
