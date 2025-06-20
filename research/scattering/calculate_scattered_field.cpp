//
// Created by evgen on 11.06.2025.
//

#include "mesh/Parser.hpp"
#include "types/Types.hpp"
#include "mesh/SurfaceMesh.hpp"
#include "mesh/MeshTypes.hpp"
#include "visualisation/VTKFunctions.hpp"
#include "mesh/Parser.hpp"
#include "slae_generation/MatrixGeneration.hpp"
#include "math/MathConstants.hpp"
#include "examples/pathes.hpp"
#include "experiment/PhysicalCondition.hpp"
#include "math/integration/gauss_quadrature/GaussLegenderPoints.hpp"
#include "mesh/Utils.hpp"

#include "operators/OperatorK.hpp"
#include "slae_generation/MatrixGeneration.hpp"

#include "../Solve.hpp"

#include <Eigen/Core>
#include <Eigen/IterativeLinearSolvers>
#include <fstream>
#include <iostream>
#include <ranges>
#include <unsupported/Eigen/IterativeSolvers>


using namespace EMW;
using namespace EMW::Types;

inline Types::Vector3c operatorK_in_point(const Math::SurfaceVectorField &field, const Types::complex_d k,
                                          const Mesh::point_t &point) {
    return EMW::OperatorK::K1<DefiniteIntegrals::GaussLegendre::Quadrature<4, 4>>(point, k, field) +
           EMW::OperatorK::K0<DefiniteIntegrals::GaussLegendre::Quadrature<4>>(point, k, field);
}

inline Types::Vector3c getE_in_point(const Math::SurfaceVectorField &j_e, const Types::complex_d k,
                                     const Mesh::point_t &point) {
    const Types::complex_d mul = Math::Constants::i / k;
    const auto value = operatorK_in_point(j_e, k, point);
    // std::cout << value << std::endl;
    return value;
}


int main() {
    const std::string nodesFile = "/home/evgen/Education/MasterDegree/thesis/objects/mobuis/2000_nodes.csv";
    const std::string cellsFile = "/home/evgen/Education/MasterDegree/thesis/objects/mobuis/1881_cells.csv";

    const auto parser_out = EMW::Parser::parseMesh(nodesFile, cellsFile);
    auto mesh_base = Mesh::SurfaceMesh{parser_out.first, parser_out.second};
    mesh_base.setName("Mobuis");

    // задаем падающую волну в виде плоской волны
    const Types::scalar lamdba = 0.5;
    const Types::complex_d k{2 * Math::Constants::PI<scalar>() / lamdba, 0};
    EMW::Physics::planeWaveCase physics(Vector3d{0, 0, 1}.normalized(),     // polarization
                                        k,  // wave figure
                                        Vector3d{1, 0, 0}.normalized());   // wave unit vector
    // смотрим след поля на поверхности
    const auto E_0_field = Math::SurfaceVectorField{mesh_base, [&physics](const Mesh::point_t& p) {return physics.value(p);}};
    const Types::VectorXc b = -1 * E_0_field.asVector();

    // расчет
    const MatrixXc A = EMW::Matrix::getMatrixK(k, mesh_base);
    const auto j_vec = Research::solve<Eigen::GMRES>(A, b, 2000, 1e-2);

    // закидываем ток как поле на поверхности
    const Math::SurfaceVectorField j_e = Math::SurfaceVectorField::TangentField(mesh_base, j_vec, "j_e");

    // расчитываем поле вокруг заданной геометрии
#define FIELD_CALCULATION 1
#if FIELD_CALCULATION
    // Рисуем картину поля в плоскости y = 0
    int k1 = 199;
    Types::scalar h1 = 4. / (k1 - 1);
    int k2 = 199;
    Types::scalar h2 = 4. / (k2 - 1);

    std::vector<Mesh::point_t> points;
    points.reserve(k1 * k2);
    Mesh::Utils::cartesian_product_unevenXY(std::ranges::views::iota(0, k2), std::ranges::views::iota(0, k1),
                                            std::back_inserter(points), k2, k1, h2, h1);

    std::cout << "Surrounding Mesh Constructed" << std::endl;

    Containers::vector<Types::Vector3c> calculated_field; calculated_field.resize(points.size());

#pragma omp parallel for num_threads(14) schedule(dynamic)
    for (int i = 0; i < points.size(); ++i) {
        calculated_field[i] = getE_in_point(j_e, k, points[i]) + physics.value(points[i]);
    }

    // отрисовка полученного поля
    std::transform(calculated_field.begin(), calculated_field.end(), calculated_field.begin(),
        [](auto&& v) {
            if (v.norm() > 10) return Types::Vector3c{{0,0},{0,0},{0,0}};
            return v;
        });

    const std::string path = "/home/evgen/Education/MasterDegree/thesis/results/mobius/";
    VTK::field_in_points_snapshot({calculated_field}, {"E"}, points, "surrounding_mesh", path);
#endif

    // рисуем поле токов
    VTK::united_snapshot({j_e}, {}, mesh_base, path);
}
