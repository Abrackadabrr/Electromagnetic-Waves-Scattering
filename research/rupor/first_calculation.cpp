//
// Created by evgen on 24.12.2024.
//

#include <string>

#include "mesh/MeshTypes.hpp"
#include "mesh/Parser.hpp"
#include "mesh/SurfaceMesh.hpp"

#include "math/fields/SurfaceVectorField.hpp"

#include "experiment/ESA.hpp"

#include "Equations.hpp"
#include "VTKFunctions.hpp"

#include "meshes/plate/PlateGrid.hpp"

#include "operators/OperatorK.hpp"
#include "operators/OperatorR.hpp"

#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include <unsupported/Eigen/IterativeSolvers>

#include <chrono>
#include <iostream>
#include <math/integration/gauss_quadrature/GaussLegenderPoints.hpp>

using namespace EMW;

Types::VectorXc solve(const Types::MatrixXc &A, const Types::VectorXc &b, Types::scalar tolerance) {
    auto method = Eigen::GMRES<Types::MatrixXc>{};

    Types::index max_iterations = 200;

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
Math::SurfaceVectorField operatorK(const Math::SurfaceVectorField &field, const Types::complex_d k,
                                   const Mesh::SurfaceMesh &targetMesh) {
    const auto analytical = [field, k](const Types::Vector3d &point) -> Types::Vector3c {
        return EMW::OperatorK::K1<DefiniteIntegrals::GaussLegendre::Quadrature<4, 4>>(point, k, field) +
               EMW::OperatorK::K0<DefiniteIntegrals::GaussLegendre::Quadrature<4>>(point, k, field);
    };

    return Math::SurfaceVectorField(targetMesh, analytical);
}

Math::SurfaceVectorField operatorR(const Math::SurfaceVectorField &field, const Types::complex_d k,
                                   const Mesh::SurfaceMesh &targetMesh) {
    const auto analytical = [field, k](const Types::Vector3d &point) -> Types::Vector3c {
        return EMW::OperatorR::detail::R<DefiniteIntegrals::GaussLegendre::Quadrature<4, 4>>(point, k, field);
    };

    return Math::SurfaceVectorField(targetMesh, analytical);
}

Math::SurfaceVectorField getE(const Math::SurfaceVectorField &j_e, const Math::SurfaceVectorField &j_m,
                              const Types::complex_d k, const Mesh::SurfaceMesh &targetMesh) {
    const Types::complex_d mul = Math::Constants::i * Math::Constants::mu_0_c / k;

    return mul * operatorK(j_e, k, targetMesh) - operatorR(j_m, k, targetMesh);
}
#endif

template <typename Container>
void to_csv(const Container &cont1, const Container &cont2, std::ofstream &str) {
    str << "sigma,alpha\n";
    for (int i = 0; i < cont1.size(); i++) {
        str << cont1[i] << ',' << cont2[i] << '\n';
    }
}

void getSigmaValues(const Types::complex_d k, const Math::SurfaceVectorField &j_e,
                                                 const Math::SurfaceVectorField &j_m) {
    int samples = 360;
    Containers::vector<Types::scalar> esas;
    esas.reserve(samples);
    Containers::vector_d angles;
    angles.reserve(samples);

    for (int i = 0; i < samples; i++) {
        Types::scalar angle = i * Math::Constants::PI<Types::scalar>() * 2 / samples;
        Types::Vector3d tau = {std::cos(angle), 0, std::sin(angle)};
        esas.push_back(ESA::calculateESA(tau, k, j_e, j_m));
        angles.push_back(angle);
    }
    std::ofstream sigma(
        "/home/evgen/Education/MasterDegree/thesis/Electromagnetic-Waves-Scattering/vtk_files/studies/rupor/sigma.csv");

    to_csv(esas, angles, sigma);
}

int main() {
    // считываем сетку на антенне
    const std::string nodesFile = "/home/evgen/Education/MasterDegree/thesis/Electromagnetic-Waves-Scattering/meshes/"
                                  "rupor/15200_nodes.csv";
    const std::string cellsFile = "/home/evgen/Education/MasterDegree/thesis/Electromagnetic-Waves-Scattering/meshes/"
                                  "rupor/3800_cells.csv";
    const EMW::Types::index nNodes = 15200;
    const EMW::Types::index nCells = 3800;
q
    const auto parser_out = EMW::Parser::parseMesh(nodesFile, cellsFile, nNodes, nCells);
    auto mesh_all = Mesh::SurfaceMesh{parser_out.first, parser_out.second};
    mesh_all.setName("total_mesh");
    auto mesh_sigma = mesh_all.getSubmesh(Mesh::IndexedCell::Tag::SIGMA);
    mesh_sigma.setName("sigma_mesh");
    auto mesh_zero = mesh_all.getSubmesh(Mesh::IndexedCell::Tag::WAVEGUIDE_CROSS_SECTION);
    mesh_zero.setName("zero_mesh");

    // Физика волны в волноводе
    const Types::scalar a = 0.07;
    const Types::scalar lambda_z = 0.10;
    //const Types::complex_d k{
    //    EMW::Math::Constants::PI<Types::scalar>() * std::sqrt((1 / (a * a)) + (4 / (lambda_z * lambda_z))), 0};
    const Types::complex_d k{2 * Math::Constants::PI<Types::scalar>() / lambda_z,0};
    // расчет коэффициента импеданса
    const Types::complex_d beta = std::sqrt(k * k - (EMW::Math::Constants::PI_square<Types::scalar>() / (a * a)));
    std::cout << beta << std::endl;
    std::cout << 2 * Math::Constants::PI<Types::scalar>() / lambda_z << std::endl;
    const Types::complex_d z = -k * Math::Constants::one_div_c / beta;
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

    const auto result = solve(matrix, rhs, 1e-5);

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

    // Рисуем картину поля в плоскости y = 0
    int N1 = 41;
    Types::scalar h1 = a / (N1 - 1);
    int N2 = 5 * 41;
    Types::scalar h2 = 2 * 0.17 / (N2 - 1);
    auto set_of_points_for_E = Examples::Plate::generateRectangularMesh(N1, N2, h1, h2);
    set_of_points_for_E.setName("E_mesh");
    auto e_field =
        getE(e_c, m_c, k, set_of_points_for_E);
    //  + Math::SurfaceVectorField{set_of_points_for_E, get_e_h10_mode}
    e_field.setName("e_field");
    VTK::united_snapshot({e_field}, {}, set_of_points_for_E, path);

    getSigmaValues(k, e_c, m_c);
}
