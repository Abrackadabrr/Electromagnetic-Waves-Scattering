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

#include "EquationsWithToeplitzStructure.hpp"
#include "FieldCalculation.hpp"
#include "FieldOverGeometry.hpp"
#include "GeneralizedEquations.hpp"

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
using TTBMatrix = LAMatrix::ToeplitzToeplitzDynFactoredBlock<Types::complex_d>;
using DiagonalPrec = LAMatrix::Preconditioning::DiagonalPreconditioner<Types::complex_d, TTBMatrix>;
using BlockDiagPrec = LAMatrix::Preconditioning::BlockDiagonalPreconditioner<Types::complex_d, TTBMatrix>;
using NoPrec = LAMatrix::Preconditioning::IdentityPreconditioner<Types::complex_d, TTBMatrix>;
using MatrixWrapper = LAMatrix::Wrappers::MatrixReplacement<TTBMatrix, BlockDiagPrec>;

template <typename matrix_t, typename matrix_full_t>
void final_check_for_vectors(const matrix_full_t &matrix, const matrix_t &toeplitz,
                             const Types::VectorX<Types::complex_d> &vec) {
    auto start = std::chrono::system_clock::now();

    std::cout << "Норма вектора: " << vec.norm() << std::endl;

    const Types::VectorX<Types::complex_d>
        result_straightforward = matrix * vec;

    auto end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "Прямое произведение: " << elapsed.count() << '\n';

    start = std::chrono::system_clock::now();

    const Types::VectorX<Types::complex_d> result_special = toeplitz * vec;

    end = std::chrono::system_clock::now();
    elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "Произведение с новой матрицей: " << elapsed.count() << '\n';

    // Сравниваем результаты умножения
    std::cout << "Относительная норма ошибки: " << (result_straightforward - result_special).norm() / result_special.norm();
}

int main() {
    // считываем сетку на антенне
    const std::string nodesFile = "/home/evgen/Education/MasterDegree/thesis/Electromagnetic-Waves-Scattering/meshes/"
                                  "lattice/8000_nodes.csv";
    const std::string cellsFile = "/home/evgen/Education/MasterDegree/thesis/Electromagnetic-Waves-Scattering/meshes/"
                                  "lattice/2000_cells.csv";
    constexpr EMW::Types::index nNodes = 8000;
    constexpr EMW::Types::index nCells = 2000;

    // собираем сетки
    const auto parser_out = EMW::Parser::parseMesh(nodesFile, cellsFile, nNodes, nCells);
    auto mesh_base = Mesh::SurfaceMesh{parser_out.first, parser_out.second};

    constexpr Types::index N1 = 1;
    constexpr Types::index N2 = 10;
    constexpr Types::index N1_x_N2 = N1 * N2;
    const Scene<N1, N2> geometry{0.1, 0.05, mesh_base};

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

    // собираем общую маленькую тёплицеву матрицу, притом сжатую

    auto start = std::chrono::high_resolution_clock::now();

    const auto matrix = Research::Lattice::getMatrixCompressed(geometry, a, k);

    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(end - start);

    std::cout << Utils::get_memory_usage(matrix) << std::endl;
    std::cout << "Matrix assembled, size: " << matrix.rows() << "; time elapsed: " << elapsed << std::endl;

    const auto matrix_full = Research::Lattice::getMatrix(geometry, a, k);

    final_check_for_vectors(matrix_full, matrix, Types::VectorX<Types::complex_d>::Random(matrix_full.cols()));

}
