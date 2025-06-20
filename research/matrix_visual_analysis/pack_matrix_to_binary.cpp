//
// Created by evgen on 02.06.2025.
//
#include <string>

#include "mesh/Parser.hpp"
#include "mesh/SurfaceMesh.hpp"

#include "math/matrix/iterative_solvers_coverage/MatrixReplacement.hpp"

#include "geometry/PeriodicStructure.hpp"

#include "VTKFunctions.hpp"

#include "../lattice/FieldCalculation.hpp"
#include "../lattice/GeneralizedEquations.hpp"

#include "experiment/PhysicalCondition.hpp"

#include "Utils.hpp"

#include <chrono>
#include <iostream>
#include <numeric>

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
    const std::string nodesFile_c = "/home/evgen/Education/MasterDegree/thesis/Electromagnetic-Waves-Scattering/meshes/waveguide/coarse_mesh/129_nodes.csv";
    const std::string cellsFile_c = "/home/evgen/Education/MasterDegree/thesis/Electromagnetic-Waves-Scattering/meshes/waveguide/coarse_mesh/127_cells.csv";
    constexpr EMW::Types::index nNodes_c = 8000;
    constexpr EMW::Types::index nCells_c = 2000;

    // путь до папки, куда писать результаты
    const std::string path = "/home/evgen/Education/MasterDegree/thesis/results/blocks_analysis/";
    // собираем сетки
    const auto parser_out = EMW::Parser::parseMesh(nodesFile_c, cellsFile_c);
    auto mesh_base = Mesh::SurfaceMesh{parser_out.first, parser_out.second};

    // сделаем ещё одну сетку для анализа внедиагонлаьных элементов
    const auto &mesh_1 = mesh_base;

    const Types::Vector3d origin{0.075, 0, 0.0};
    const auto mesh_2 = Mesh::Utils::move_by_vector(mesh_base, origin);

    // Геометрические параметры антенн
    // Короткая сторона волновода
    const Types::scalar a = 0.07;
    // Физика волны в пространстве
    // частота в гигагерцах
    const Types::scalar freq = 2 * Math::Constants::c / 1e8;
    const Types::complex_d k{Physics::get_k_on_frquency(freq), 0};

    const auto out_diagonal_block = WaveGuideWithActiveSection::submatrix(mesh_2, mesh_1, a, k);
    const auto diag_block = WaveGuideWithActiveSection::diagonal(mesh_2, a, k);

    dump(diag_block, "/home/evgen/Education/MasterDegree/thesis/results/blocks_analysis/diagonal.txt");
    dump(out_diagonal_block, "/home/evgen/Education/MasterDegree/thesis/results/blocks_analysis/out_diag.txt");
}
