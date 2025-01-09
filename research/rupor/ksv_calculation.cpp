//
// Created by evgen on 10.01.2025.
//
//
// Created by evgen on 24.12.2024.
//

#include <string>

#include "mesh/MeshTypes.hpp"
#include "mesh/Parser.hpp"
#include "mesh/SurfaceMesh.hpp"

#include "math/fields/SurfaceVectorField.hpp"

#include "experiment/ESA.hpp"
#include "experiment/SWC.hpp"

#include "Equations.hpp"
#include "VTKFunctions.hpp"

#include "meshes/plate/PlateGrid.hpp"

#include "operators/OperatorK.hpp"
#include "operators/OperatorR.hpp"

#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include <unsupported/Eigen/IterativeSolvers>

#include "experiment/PhysicalCondition.hpp"

#include "Utils.hpp"

#include <math/integration/gauss_quadrature/GaussLegenderPoints.hpp>

#include <chrono>
#include <iostream>

using namespace EMW;

Types::VectorXc solve(const Types::MatrixXc &A, const Types::VectorXc &b, Types::scalar tolerance) {
    auto method = Eigen::GMRES<Types::MatrixXc>{};

    Types::index max_iterations = 500;

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
#if 1

Types::Vector3c operatorK_in_point(const Math::SurfaceVectorField &field, const Types::complex_d k,
                                   const Mesh::point_t &point) {
    return EMW::OperatorK::K1<DefiniteIntegrals::GaussLegendre::Quadrature<4, 4>>(point, k, field) +
           EMW::OperatorK::K0<DefiniteIntegrals::GaussLegendre::Quadrature<4>>(point, k, field);
}

Math::SurfaceVectorField operatorK(const Math::SurfaceVectorField &field, const Types::complex_d k,
                                   const Mesh::SurfaceMesh &targetMesh) {

    const auto analytical = [field, k](const Types::Vector3d &point) -> Types::Vector3c {
        return operatorK_in_point(field, k, point);
    };
    return Math::SurfaceVectorField(targetMesh, analytical);
}

Types::Vector3c operatorR_in_point(const Math::SurfaceVectorField &field, const Types::complex_d k,
                                   const Mesh::point_t &point) {
    return EMW::OperatorR::detail::R<DefiniteIntegrals::GaussLegendre::Quadrature<4, 4>>(point, k, field);
}

Math::SurfaceVectorField operatorR(const Math::SurfaceVectorField &field, const Types::complex_d k,
                                   const Mesh::SurfaceMesh &targetMesh) {
    const auto analytical = [field, k](const Types::Vector3d &point) -> Types::Vector3c {
        return operatorR_in_point(field, k, point);
    };
    return Math::SurfaceVectorField(targetMesh, analytical);
}

Types::Vector3c getE_in_point(const Math::SurfaceVectorField &j_e, const Math::SurfaceVectorField &j_m,
                              const Types::complex_d k, const Mesh::point_t &point) {
    const Types::complex_d mul = Math::Constants::i / k;
    return mul * operatorK_in_point(j_e, k, point) - operatorR_in_point(j_m, k, point);
}
Math::SurfaceVectorField getE(const Math::SurfaceVectorField &j_e, const Math::SurfaceVectorField &j_m,
                              const Types::complex_d k, const Mesh::SurfaceMesh &targetMesh) {
    const Types::complex_d mul = Math::Constants::i / k;
#ifdef OLD
    return mul * operatorK(j_e, k, targetMesh) - operatorR(j_m, k, targetMesh);
#else
    const Types::MatrixXc R_matrix = Matrix::getMatrixR(k, j_m.getManifold(), targetMesh);
    std::cout << "matrix R computation done" << std::endl;
    const Types::MatrixXc K_matrix = Matrix::getMatrixK(k, j_e.getManifold(), targetMesh);
    std::cout << "matrix K computation done" << std::endl;

    std::cout << K_matrix.rows() << std::endl;
    std::cout << K_matrix.cols() << std::endl;
    std::cout << R_matrix.rows() << std::endl;
    std::cout << R_matrix.cols() << std::endl;
    std::cout << j_e.asSLAERHS().rows() << std::endl;
    std::cout << j_m.asSLAERHS().rows() << std::endl;
    const Types::VectorXc e_vector = mul * K_matrix * j_e.asSLAERHS() - R_matrix * j_m.asSLAERHS();
    std::cout << "multiplication done" << std::endl;
    return Math::SurfaceVectorField::TangentField(targetMesh, e_vector);
#endif
}
#endif

template <typename Callable>
void getSWConAxis(const Math::SurfaceVectorField &j_e, const Math::SurfaceVectorField &j_m, const Types::complex_d k,
                  const Callable &direct_field) {
    // собрать точки на оси
    int N = 5 * 41;
    Types::scalar h = 0.17 / (N - 1);
    const auto points_view = std::views::iota(0, N) | std::views::transform([&](auto p) { return Types::Vector3d{0, 0, p * h}; });
    const auto z_values_view = points_view | std::views::transform([&](auto p) { return p.z(); });
    Containers::vector<Mesh::point_t> points{std::begin(points_view), std::end(points_view)};
    const Containers::vector<Types::scalar> z_coordinate{std::begin(z_values_view), std::end(z_values_view)};

    const auto inverse_field = [&](Mesh::point_t &point) {
        return getE_in_point(j_e, j_m, k, point) - direct_field(point);
    };

    const auto result = EngineeringCoefficients::SWC(direct_field, inverse_field, points);

    std::ofstream sigma(
        "/home/evgen/Education/MasterDegree/thesis/Electromagnetic-Waves-Scattering/vtk_files/studies/rupor/ksv.csv");
    Utils::to_csv(result, z_coordinate, "z", "ksv", sigma);
}

int main() {
    // считываем сетку на антенне
    const std::string nodesFile = "/home/evgen/Education/MasterDegree/thesis/Electromagnetic-Waves-Scattering/meshes/"
                                  "rupor/15200_nodes.csv";
    const std::string cellsFile = "/home/evgen/Education/MasterDegree/thesis/Electromagnetic-Waves-Scattering/meshes/"
                                  "rupor/3800_cells.csv";
    const EMW::Types::index nNodes = 15200;
    const EMW::Types::index nCells = 3800;

    const auto parser_out = EMW::Parser::parseMesh(nodesFile, cellsFile, nNodes, nCells);
    auto mesh_all = Mesh::SurfaceMesh{parser_out.first, parser_out.second};
    mesh_all.setName("total_mesh");
    auto mesh_sigma = mesh_all.getSubmesh(Mesh::IndexedCell::Tag::SIGMA);
    mesh_sigma.setName("sigma_mesh");
    auto mesh_zero = mesh_all.getSubmesh(Mesh::IndexedCell::Tag::WAVEGUIDE_CROSS_SECTION);
    mesh_zero.setName("zero_mesh");

    // Короткая сторона волновода
    const Types::scalar a = 0.07;
    // Физика волны в пространстве
    // частота в гигагерцах
    const Types::scalar freq = 3;
    const Types::complex_d k{Physics::get_k_on_frquency(freq), 0};
    // расчет коэффициента импеданса
    const Types::complex_d beta = std::sqrt(k * k - (EMW::Math::Constants::PI_square<Types::scalar>() / (a * a)));
    std::cout << beta << std::endl;
    std::cout << k.real() << std::endl;

    // тут сделать векторное поле (direct wave в волноводе)
    const auto get_e_h10_mode = [&k, &a, &beta](const Mesh::point_t &x) {
        const auto pi = Math::Constants::PI<Types::scalar>();
        const Types::complex_d mult = Math::Constants::i * (a / pi) * k / Math::Constants::e_0_c;
        const Types::complex_d exp = std::exp(Math::Constants::i * beta * x.z());
        const Types::scalar sin = std::cos(pi * x.x() / a);
        return Types::Vector3c{Types::complex_d{0, 0}, mult * sin * exp, Types::complex_d{0, 0}};
    };

    // Считаем матрицу и правую часть и решаем СЛАУ
    const auto rhs = Rupor::getRhs(mesh_all, mesh_zero, get_e_h10_mode);

    const auto matrix = Rupor::getMatrix(mesh_all, mesh_sigma, mesh_zero, a, k);

    std::cout << "Matrix assembled" << std::endl;

    const auto result = solve(matrix, rhs, 1e-2);

    // Разбиваем на два тока и рисуем на разных многообразиях
    const Types::VectorXc electric_current = result.block(0, 0, 2 * mesh_all.getCells().size(), 1);
    const Types::VectorXc magnetic_current =
        result.block(2 * mesh_all.getCells().size(), 0, 2 * mesh_zero.getCells().size(), 1);

    const std::string path =
        "/home/evgen/Education/MasterDegree/thesis/Electromagnetic-Waves-Scattering/vtk_files/studies/rupor/";

    Math::SurfaceVectorField e_c = Math::SurfaceVectorField::TangentField(mesh_all, electric_current);
    e_c.setName("e_c");
    Math::SurfaceVectorField m_c = Math::SurfaceVectorField::TangentField(mesh_zero, magnetic_current);
    m_c.setName("m_c");

    VTK::united_snapshot({e_c}, {}, mesh_all, path);
    VTK::united_snapshot({m_c}, {}, mesh_zero, path);

    getSWConAxis(e_c, m_c, k, get_e_h10_mode);
}
