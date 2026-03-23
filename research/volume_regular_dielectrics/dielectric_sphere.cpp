//
// Created by evgen on 22.03.2026.
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

constexpr Types::scalar SPHERE_RADUIS = 0.25;

Types::scalar permittivity_distribution(Types::point_t x) {
    return x.norm() < SPHERE_RADUIS ? 2 : 1;
}

int main() {
    // 1. Рисуем сетку
    constexpr Types::scalar cube_length = 2 * SPHERE_RADUIS;
    constexpr Types::index Nx = 16;
    constexpr Types::index Ny = 16;
    constexpr Types::index Nz = 16;
    constexpr Types::scalar mesh_one_axis_size = SPHERE_RADUIS / (Nx - 1);
    Mesh::VolumeMesh::CubeMeshWithData mesh{Types::point_t{0, 0, 0}, (Nx - 1) * mesh_one_axis_size,
                                            (Ny - 1) * mesh_one_axis_size, (Nz - 1) * mesh_one_axis_size, Nx, Ny, Nz};
    mesh.setName("sphere");
    // Настраиваем диэлектрическую проницаемость
    mesh.invokeScalarData("eps", permittivity_distribution);

    // 2. Параметры падающей волны
    constexpr double freq = 0.3; // GGz
    constexpr Types::complex_d k{freq, 0.};
    Physics::planeWaveCase incident_field{Types::Vector3d{0, 1, 0}, k, Types::Vector3d{1, 0, 0}};
    std::cout << "Длина волны в свободном пространстве = " << 2 * M_PI / k.real() << std::endl;

    // 3. Галеркинская проекция
    Operators::Volume::ProjectorOnMesh proj{mesh};
    auto rhs = proj([incident_field](Types::point_t p) { return incident_field.value(p); });
    std::cout << "Rhs rows = " << rhs.rows() << std::endl;

    // Матрицу собираем
    Operators::Volume::operator_K_over_cube_mesh operator_K{k, mesh};
    auto mat = operator_K.get_galerkin_matrix();
    const Types::scalar cube_measure = mesh.dx() * mesh.dy() * mesh.dz();

    // Изменяем матрицу
    Types::VectorXc diag_eps = Types::VectorXc::Zero(3 * mesh.getCells().size());
    const auto eps_data = mesh.getScalarData("eps");
    for (size_t idx = 0; idx < mesh.getCells().size(); ++idx) {
        diag_eps[3 * idx] = eps_data[idx] - 1.;
        diag_eps[3 * idx + 1] = eps_data[idx] - 1.;
        diag_eps[3 * idx + 2] = eps_data[idx] - 1.;
    }
    const Types::MatrixXc A = Types::MatrixXc::Identity(mat.rows(), mat.cols()) - (mat / cube_measure) * Eigen::Diagonal
                              {diag_eps};
    // Поправляем правую часть
    const Types::VectorXc b = rhs / std::sqrt(cube_measure);

    // 4. Решаем систему
    const auto solution = Research::solve<Eigen::GMRES>(A, b, 1000, 1e-8);

    // 5. Преобразовываем в векторное поле на ячейках и пишем в данные сетки
    std::vector<Types::Vector3c> field_on_mesh;
    field_on_mesh.reserve(mesh.getCells().size());
    for (size_t idx = 0; idx < mesh.getCells().size(); ++idx) {
        field_on_mesh.emplace_back(solution[3 * idx + 0] / std::sqrt(cube_measure),
                                   solution[3 * idx + 1] / std::sqrt(cube_measure),
                                   solution[3 * idx + 2] / std::sqrt(cube_measure));
    }

    // Добавляем поле
    mesh.setVectorData("solution", std::move(field_on_mesh));

    // Отрисовываем результат
    VTK::volume_mesh_withdata_snapshot(
        mesh,
        "/home/evgen/Education/MasterDegree/thesis/Electromagnetic-Waves-Scattering/research/volume_regular_dielectrics/");

}
