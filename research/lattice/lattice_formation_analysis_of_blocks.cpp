//
// Created by evgen on 06.03.2025.
//
#include <string>

#include "mesh/Parser.hpp"
#include "mesh/SurfaceMesh.hpp"

#include "math/matrix/iterative_solvers_coverage/MatrixReplacement.hpp"

#include "geometry/PeriodicStructure.hpp"

#include "VTKFunctions.hpp"

#include "FieldCalculation.hpp"
#include "GeneralizedEquations.hpp"

#include <unsupported/Eigen/IterativeSolvers>

#include "experiment/PhysicalCondition.hpp"

#include <Utils.hpp>
#include <chrono>
#include <iostream>

template <typename Matrix> void dump(const Matrix &matrix, const std::string &filename) {
    std::ofstream file(filename);
    file.precision(20);
    for (Types::index row = 0; row < matrix.rows(); ++row) {
        for (Types::index col = 0; col < matrix.cols(); ++col) {
            file << matrix(row, col).real() << "+" << matrix(row, col).imag() << "j";
            if (col !=  matrix.cols() - 1) file << " ";
        }
        file << std::endl;
    }
    file.close();
}

int main() {
    // путь до папки, куда писать результаты
    const std::string path = "/home/evgen/Education/MasterDegree/thesis/results/lattice/";
    // считываем сетку на антенне

#define REAL_CASE 1
#define SVD 0

#if REAL_CASE
    const std::string nodesFile = "/home/evgen/Education/MasterDegree/thesis/Electromagnetic-Waves-Scattering/meshes/"
                                  "lattice/8000_nodes.csv";
    const std::string cellsFile = "/home/evgen/Education/MasterDegree/thesis/Electromagnetic-Waves-Scattering/meshes/"
                                  "lattice/2000_cells.csv";
    constexpr EMW::Types::index nNodes = 8000;
    constexpr EMW::Types::index nCells = 2000;
#else
    const std::string cellsFile = "/home/evgen/Education/MasterDegree/thesis/Electromagnetic-Waves-Scattering/tests/"
                                  "meshes_for_tests/lattice_redused/80_cells.csv";
    const std::string nodesFile = "/home/evgen/Education/MasterDegree/thesis/Electromagnetic-Waves-Scattering/tests/"
                                  "meshes_for_tests/lattice_redused/320_nodes.csv";
    constexpr EMW::Types::index nNodes = 320;
    constexpr EMW::Types::index nCells = 80;
#endif

    // собираем сетки
    const auto parser_out = EMW::Parser::parseMesh(nodesFile, cellsFile, nNodes, nCells);
    auto mesh_base = Mesh::SurfaceMesh{parser_out.first, parser_out.second};

    // сделаем ещё одну сетку для анализа внедиагонлаьных элементов
    const auto &mesh_1 = mesh_base;

    const Types::Vector3d origin{4 * 0.14, 0, 0.0};
    const auto mesh_2 = Mesh::Utils::move_by_vector(mesh_base, origin);

    // Геометрические параметры антенн
    // Короткая сторона волновода
    const Types::scalar a = 0.07;
    // Физика волны в пространстве
    // частота в гигагерцах
    const Types::scalar freq = 2 * Math::Constants::c / 1e8;
    const Types::complex_d k{Physics::get_k_on_frquency(freq), 0};

    const auto out_diagonal_block = WaveGuideWithActiveSection::submatrix(mesh_2, mesh_1, a, k);
    const auto diagonal_block = WaveGuideWithActiveSection::diagonal(mesh_1, a, k);

#if SVD
    const auto svd = out_diagonal_block.bdcSvd();
    const auto SVD_diag = diagonal_block.bdcSvd();

    const auto sing_vals = svd.singularValues();
    const auto sing_vals_diag = SVD_diag.singularValues();
    Containers::vector<int> indexes{};
    indexes.resize(sing_vals.size());
    std::iota(indexes.begin(), indexes.end(), 0);
    // Запись сингулярных чисел
    std::ofstream file(path + "sv_out.csv");
    Utils::to_csv(indexes, sing_vals, "ind", "sv", file);

    std::ofstream file_diag(path + "sv_diagonal.csv");
    Utils::to_csv(indexes, sing_vals_diag, "ind", "sv_out_diag", file);
#endif

    // Запись матриц в файлик для обработке в други месте
    dump(diagonal_block, path + "matrix/diagonal.txt");
    dump(out_diagonal_block, path + "matrix/out_diagonal.txt");
}
