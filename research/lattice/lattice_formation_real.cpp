//
// Created by evgen on 05.02.2025.
//

#include <string>

#include "mesh/MeshTypes.hpp"
#include "mesh/Parser.hpp"
#include "mesh/SurfaceMesh.hpp"

#include "math/fields/SurfaceVectorField.hpp"

#include "GeneralizedEquations.hpp"

#include "equations/EquationsOverGeometry.hpp"


#include "geometry/PeriodicStructure.hpp"

#include "VTKFunctions.hpp"

#include "meshes/plate/PlateGrid.hpp"

#include "FieldCalculation.hpp"

#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include <unsupported/Eigen/IterativeSolvers>

#include "experiment/PhysicalCondition.hpp"

#include "Utils.hpp"

#include <chrono>
#include <iostream>

namespace eq = WaveGuideWithActiveSection;

template <int N1, int N2> using Scene = Geometry::PeriodicStructure<N1, N2>;

Types::VectorXc solve(const Types::MatrixXc &A, const Types::VectorXc &b, Types::scalar tolerance) {
    auto method = Eigen::GMRES<Types::MatrixXc>{};

    Types::index max_iterations = 1000;

    method.setMaxIterations(max_iterations);
    std::cout << method.maxIterations() << std::endl;
    method.setTolerance(tolerance);
    method.set_restart(max_iterations);

    std::chrono::time_point<std::chrono::system_clock> start = std::chrono::system_clock::now();

    method.compute(A);
    auto j_vec = Types::VectorXc{method.solve(b)};

    std::chrono::time_point<std::chrono::system_clock> end = std::chrono::system_clock::now();
    std::chrono::duration<double, std::ratio<1, 1>> elapsed_seconds = end - start;
    std::cout << "Время решения GMRES: " << elapsed_seconds.count() << std::endl;
    std::cout << "total iterations: " << method.iterations() << std::endl;
    std::cout << "tolerance: " << method.tolerance() << std::endl;
    std::cout << "total error: " << method.error() << std::endl;
    std::cout << "Info: " << static_cast<int>(method.info()) << std::endl;
    return j_vec;
}

#if 0
Containers::vector<Types::Vector3c>
calculateField(const Math::SurfaceVectorField &e_c_1, const Math::SurfaceVectorField &e_c_2,
               const Math::SurfaceVectorField &m_c_1, const Math::SurfaceVectorField &m_c_2,
               const Containers::vector<Mesh::point_t> &points, Types::complex_d k) {
    Containers::vector<Types::Vector3c> field;
    field.resize(points.size());
#pragma omp parallel for num_threads(14) schedule(static)
    for (int i = 0; i < points.size(); ++i) {
        field[i] =
            Lattice::getE_in_point(e_c_1, m_c_1, k, points[i]) + Lattice::getE_in_point(e_c_2, m_c_2, k, points[i]);
    }
    return field;
}


void getSigmaValuesXZ(const Types::complex_d k, const Math::SurfaceVectorField &j_e_1, const Math::SurfaceVectorField &j_e_2,
                    const Math::SurfaceVectorField &j_m_1, const Math::SurfaceVectorField &j_m_2) {
    int samples = 360;
    Containers::vector<Types::scalar> esas;
    esas.reserve(samples);
    Containers::vector_d angles;
    angles.reserve(samples);

    for (int i = 0; i < samples; i++) {
        Types::scalar angle = i * Math::Constants::PI<Types::scalar>() * 2 / samples;
        Types::Vector3d tau = {std::sin(angle), 0, std::cos(angle)};
        esas.push_back(ESA::calculateESA(tau, k, {j_e_1, j_e_2}, {j_m_1, j_m_2}));
        angles.push_back(angle);
    }
    std::ofstream sigma(
        "/home/evgen/Education/MasterDegree/thesis/results/lattice/sigmaXZ.csv");

    Utils::to_csv(esas, angles, "sigma", "angle", sigma);
}
#endif

int main() {
    // считываем сетку на антенне
    const std::string nodesFile = "/home/evgen/Education/MasterDegree/thesis/Electromagnetic-Waves-Scattering/meshes/"
                                  "lattice/8000_nodes.csv";
    const std::string cellsFile = "/home/evgen/Education/MasterDegree/thesis/Electromagnetic-Waves-Scattering/meshes/"
                                  "lattice/2000_cells.csv";
    const EMW::Types::index nNodes = 8000;
    const EMW::Types::index nCells = 2000;

    // собираем сетки
    const auto parser_out = EMW::Parser::parseMesh(nodesFile, cellsFile, nNodes, nCells);
    auto mesh_base = Mesh::SurfaceMesh{parser_out.first, parser_out.second};

    constexpr Types::index N1 = 2;
    constexpr Types::index N2 = 2;
    constexpr Types::index N1_x_N2 = N1 * N2;
    const Scene<N1, N2> geometry{0.1, 0.05, mesh_base};

#if 0
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

    // собираем общую матрицу
    const auto matrix = EMW::Equations::MatrixFor<Scene<N1, N2>, eq::diagonal_t, eq::submatrix_t>::compute(
        geometry, eq::diagonal, eq::submatrix, {a, k});

    std::cout << "Matrix assembled, size: " << matrix.rows() << std::endl;

    // собираем правую часть шаманским спобосом (очень шаманским)
    // решаем какой будет фазовый фактор на волноводах
    const Types::scalar phi = 0 * Math::Constants::PI<Types::scalar>();
    const Types::complex_d phase_factor = std::exp(Math::Constants::i * phi);
    const Containers::array<Types::complex_d, N1_x_N2> phase_factors{phase_factor, phase_factor, phase_factor,
                                                                     phase_factor};
    const auto rhs = eq::getRhs(phase_factors, mesh_base, a, k);

    std::cout << "RHS assembled, size: " << rhs.rows() << std::endl;

    const auto result = solve(matrix, rhs, 1e-2);

    // Разбиваем на токи и рисуем на разных многообразиях
    const Types::VectorXc electric_current_1 = result.block(0, 0, 2 * mesh_all.getCells().size(), 1);
    const Types::VectorXc magnetic_current_1 =
        result.block(2 * mesh_all.getCells().size(), 0, 2 * mesh_zero_1.getCells().size(), 1);
    const Types::VectorXc electric_current_2 = result.block(rhs_1.rows(), 0, 2 * mesh_all.getCells().size(), 1);
    const Types::VectorXc magnetic_current_2 =
        result.block(rhs_1.rows() + 2 * mesh_all.getCells().size(), 0, 2 * mesh_zero_1.getCells().size(), 1);
#endif
    const std::string path = "/home/evgen/Education/MasterDegree/thesis/results/lattice/arbitrary_sizes/";

    VTK::geometry_snapshot(geometry, path + "mesh.vtu");

#if 0

    // Рисуем картину поля в плоскости y = 0
    int N1 = 100;
    Types::scalar h1 = 1. / (N1 - 1);
    int N2 = 200;
    Types::scalar h2 = 2. / (N2 - 1);

    std::vector<Mesh::point_t> points;
    points.reserve(N1 * N2);
    Mesh::Utils::cartesian_product_unevenXZ(std::ranges::views::iota(0, N1), std::ranges::views::iota(0, N2),
                                            std::back_inserter(points), N1, N2, h1, h2);

    const auto calculated_field = calculateField(e_c_1, e_c_2, m_c_1, m_c_2, points, k);

    VTK::field_in_points_snapshot({calculated_field}, {"E"}, points, "surrounding_mesh", path);

    // считаем ДН системы из двух антенн
    getSigmaValuesXZ(k, e_c_1, e_c_2, m_c_1, m_c_2);
#endif
}
