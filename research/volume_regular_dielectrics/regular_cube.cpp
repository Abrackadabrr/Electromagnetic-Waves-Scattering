//
// Created by evgen on 10.02.2026.
//

#include "types/Types.hpp"

#include "mesh/volume_mesh/CubeMeshWithData.hpp"

#include "operators/volume/OperatorK.hpp"
#include "operators/volume/ProjectorOnMesh.hpp"

#include "experiment/PhysicalCondition.hpp"

#include "../Solve.hpp"

#include "VTKFunctions.hpp"

#include "math/fields/Utils.hpp"

#include <iostream>

using namespace EMW;

int main() {
    // Параметры куба
    constexpr Types::scalar cube_length = 0.025;
    constexpr Types::index Nx = 3;
    constexpr Types::index Ny = static_cast<size_t>(1 / cube_length) + 1;
    constexpr Types::index Nz = Ny;
    constexpr Types::scalar cube_measure = cube_length * cube_length * cube_length;
    // L2 норма базисных функций
    std::cout << "L2 норма базисных функций = " << cube_measure << std::endl;
    // Сетка
    Mesh::VolumeMesh::CubeMeshWithData mesh{Types::point_t{0, -0.5, -0.5}, (Nx - 1) * cube_length,
                                            (Ny - 1) * cube_length, (Nz - 1) * cube_length, Nx, Ny, Nz};

    // Параметры падающей волны
    constexpr double freq = 0.3; // GGz
    constexpr Types::complex_d k{Physics::get_k_on_frquency(freq), 0.};
    Physics::planeWaveCase incident_field{Types::Vector3d{0, 1, 0}, k, Types::Vector3d{1, 0, 0}};

    std::cout << "Длина волны в свободном пространстве = " << 2 * M_PI / k.real() << std::endl;

    // Галеркинская проекция правой части на сетку (правая часть линейной системы)
    Operators::Volume::ProjectorOnMesh proj{mesh};
    auto rhs = proj([incident_field](Types::point_t p) { return incident_field.value(p); });
    std::cout << "RHS size: " << rhs.rows() << std::endl;

    // Матрицу собираем
    Operators::Volume::operator_K_over_cube_mesh operator_K{k, mesh};
    auto mat = operator_K.get_galerkin_matrix();
    std::cout << "Matrix (" << mat.rows() << ", " << mat.cols() << ')' << std::endl;
    std::cout << "Matrix for norm = " << mat.norm() << std::endl;
    // Собираем матрицу системы
    mat = Types::MatrixXc::Identity(mat.rows(), mat.cols()) - mat / (cube_measure * cube_measure);
    // Поправляем правую часть
    rhs = rhs / (cube_measure * cube_measure);

    // Решаем линейную систему
    const auto solution = Research::solve<Eigen::GMRES>(mat, rhs, 100, 1e-5);

    // Преобразовываем в векторное поле на ячейках
    std::vector<Types::Vector3c> field_on_mesh;
    field_on_mesh.reserve(mesh.getCells().size());
    for (size_t idx = 0; idx < mesh.getCells().size(); ++idx) {
        field_on_mesh.emplace_back(solution[3 * idx + 0], solution[3 * idx + 1], solution[3 * idx + 2]);
    }

    // Добавляем поле
    mesh.setVectorData("solution", std::move(field_on_mesh));

    // Отрисовываем результат
    VTK::volume_mesh_withdata_snapshot(
        mesh,
        "/home/evgen/Education/MasterDegree/thesis/Electromagnetic-Waves-Scattering/research/volume_regular_dielectrics/");

#define CALCULATE_FIELD 0

#if CALCULATE_FIELD
    // Делаем множество точек вокруг куба в плоскости Оху
    int k1 = 199;
    Types::scalar h1 = 4. / (k1 - 1);
    int k2 = 199;
    Types::scalar h2 = 4. / (k2 - 1);

    std::vector<Mesh::point_t> points;
    points.reserve(k1 * k2);
    Mesh::Utils::cartesian_product_unevenYZ(std::ranges::views::iota(0, k2), std::ranges::views::iota(0, k1),
                                            std::back_inserter(points), k2, k1, h2, h1);

    std::cout << "Surrounding Mesh Constructed" << std::endl;

    Containers::vector<Types::Vector3c> calculated_field;
    calculated_field.resize(points.size());

#pragma omp parallel for num_threads(14) schedule(dynamic)
    for (int i = 0; i < points.size(); ++i) {
        calculated_field[i] = getE_in_point(j_e, k, points[i]) + physics.value(points[i]);
    }

    // отрисовка полученного поля
    std::transform(calculated_field.begin(), calculated_field.end(), calculated_field.begin(), [](auto &&v) {
        if (v.norm() > 10)
            return Types::Vector3c{{0, 0}, {0, 0}, {0, 0}};
        return v;
    });
    VTK::field_in_points_snapshot({calculated_field}, {}, {"E"}, {}, points, "surrounding_mesh", path);
#endif
}
