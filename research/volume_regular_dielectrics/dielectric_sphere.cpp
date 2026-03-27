//
// Created by evgen on 22.03.2026.
//

#include "types/Types.hpp"

#include "mesh/volume_mesh/CubeMeshWithData.hpp"

#include "operators/volume/OperatorK.hpp"
#include "operators/volume/ProjectorOnMesh.hpp"

#include "experiment/PhysicalCondition.hpp"
#include "experiment/ESA.hpp"

#include "../Solve.hpp"

#include "VTKFunctions.hpp"

#include "MatrixTraits.hpp"
#include "MatrixReplacement.hpp"

#include "math/fields/Utils.hpp"
#include "math/integration/newton_cotess/Rectangular.hpp"
#include "math/integration/gauss_quadrature/GaussLegenderPoints.hpp"

#include "Utils.hpp"

#include <iostream>

using namespace EMW;

constexpr Types::scalar SPHERE_RADUIS = 0.5;
constexpr Types::scalar CUBE_LENGTH = 0.5;

Types::scalar permittivity_distribution(const Types::point_t &x) {
    return x.norm() < SPHERE_RADUIS ? 2.56 : 1;
}

Types::scalar permittivity_distribution_cube(const Types::point_t &x) {
    return 3 * ((std::abs(x.x()) < CUBE_LENGTH / 2) &&
    (std::abs(x.y()) < CUBE_LENGTH / 2) &&
    (std::abs(x.z()) < CUBE_LENGTH / 2)) + 1;
}

int main() {
    // 1. Рисуем сетку
    constexpr Types::scalar cube_length = 2 * CUBE_LENGTH;
    constexpr Types::index Nx = 23;
    constexpr Types::index Ny = 23;
    constexpr Types::index Nz = 23;
    constexpr Types::scalar mesh_one_axis_size = cube_length / (Nx - 1);
    Mesh::VolumeMesh::CubeMeshWithData mesh{Types::point_t{-cube_length / 2, -cube_length / 2, -cube_length / 2},
                                            (Nx - 1) * mesh_one_axis_size,
                                            (Ny - 1) * mesh_one_axis_size, (Nz - 1) * mesh_one_axis_size, Nx, Ny, Nz};
    mesh.setName("cube");
    // Настраиваем диэлектрическую проницаемость
    mesh.smoothScalarData<DefiniteIntegrals::NewtonCotess::Quadrature<4, 4, 4>>("eps", permittivity_distribution_cube);

    // 2. Параметры падающей волны
    constexpr double freq = 1; // GHz
    constexpr Types::complex_d k{2 * M_PI / CUBE_LENGTH, 0.};
    Physics::planeWaveCase incident_field{Types::Vector3d{0, 1, 0}, k, Types::Vector3d{1, 0, 0}};
    std::cout << "Длина волны в свободном пространстве = " << 2 * M_PI / k.real() << std::endl;

    // 3. Галеркинская проекция
    Operators::Volume::ProjectorOnMesh proj{mesh};
    const Types::scalar cube_measure = mesh.dx() * mesh.dy() * mesh.dz();
    std::cout << "Объем элемента сетки: " << cube_measure << std::endl;
    auto rhs = proj([incident_field](Types::point_t p) { return incident_field.value(p); });
    std::cout << "Rhs rows = " << rhs.rows() << std::endl;
    // Поправляем правую часть
    Types::VectorXc b = rhs / std::sqrt(cube_measure);

    // Матрицу собираем
    Operators::Volume::operator_K_over_cube_mesh operator_K{k, mesh};
    auto mat = operator_K.compute_galerkin_matrix(std::sqrt(cube_measure));
    // Собираем матрицу diag(eps - 1) как вектор из значений
    Types::VectorXc diag_eps = Types::VectorXc::Zero(3 * mesh.getCells().size());
    const auto eps_data = mesh.getScalarData("eps");
    for (size_t idx = 0; idx < mesh.getCells().size(); ++idx) {
        diag_eps[3 * idx] = eps_data[idx] - 1.;
        diag_eps[3 * idx + 1] = eps_data[idx] - 1.;
        diag_eps[3 * idx + 2] = eps_data[idx] - 1.;
    }
    Math::LinAgl::Matrix::Wrappers::VolumeOperatorMatrixReplacement A{mat, diag_eps};
    // Поправляем правую часть по маске из фиктивных элементов
    b = A.modify_rhs_according_to_mask(b);
    // 4. Решаем систему
    auto solution = Research::solve<Eigen::GMRES>(A, b, 1000, 2e-3);
    // ВАЖНО, что тут дополнительно делим на sqrt(epsilon),
    // потому что в системе была замена переменных
    solution = A.modify_vec_according_to_inverse_epsilon_sqrt(solution);

    // 5. Преобразовываем в векторное поле на ячейках и пишем в данные сетки

    std::vector<Types::Vector3c> field_on_mesh;
    field_on_mesh.reserve(mesh.getCells().size());
    for (size_t idx = 0; idx < mesh.getCells().size(); ++idx) {
        const auto eps_m_one = std::abs(eps_data[idx]) > 1e-14 ? eps_data[idx] : 1;
        field_on_mesh.emplace_back(solution[3 * idx + 0] / (std::sqrt(cube_measure)),
                                   solution[3 * idx + 1] / (std::sqrt(cube_measure)),
                                   solution[3 * idx + 2] / (std::sqrt(cube_measure)));
    }

    // Добавляем поле
    mesh.setVectorData("solution", std::move(field_on_mesh));

    // Отрисовываем результат
    const std::string path =
        "/home/evgen/Education/MasterDegree/thesis/Electromagnetic-Waves-Scattering/research/volume_regular_dielectrics/";
    VTK::volume_mesh_withdata_snapshot(mesh, path);
    // VTK::volume_mesh_withdata_snapshot(sphere_mesh, path);

    // Расчситываем диаграмму направленности
    int N = 720;
    const auto get_tau_hh = [](Types::scalar phi) { return Types::Vector3d{std::cos(phi), 0, std::sin(phi)}; };
    const auto get_tau_vv = [](Types::scalar phi) { return Types::Vector3d{std::cos(phi), std::sin(phi), 0}; };

    auto view = std::views::iota(0, N) | std::views::transform([N](int i) { return i * M_PI / N; });
    std::vector<Types::scalar> phis{view.begin(), view.end()};

    std::vector<Types::scalar> rsp_hh{};
    rsp_hh.reserve(N);
    std::vector<Types::scalar> rsp_vv{};
    rsp_vv.reserve(N);

    for (auto &&p : phis) {
        rsp_vv.push_back(std::log(ESA::calculateRSP(get_tau_vv(p), k, "solution", mesh)));
        rsp_hh.push_back(std::log(ESA::calculateRSP(get_tau_hh(p), k, "solution", mesh)));
    }
    auto degree_view = view | std::views::transform([&](Types::scalar phi) { return phi * 180 / M_PI; });
    std::vector<Types::scalar> phis_degree{degree_view.begin(), degree_view.end()};
    std::ofstream rsp_vv_file{path + "sigma_vv.csv"};
    std::ofstream rsp_hh_file{path + "sigma_hh.csv"};
    Utils::to_csv(phis_degree, rsp_vv, "angle", "rsp", rsp_vv_file);
    Utils::to_csv(phis_degree, rsp_hh, "angle", "rsp", rsp_hh_file);
};

// Надо сделать:
// 1) Уточнить, что такое VV и HH в Гибсоне
// 2) Попробовать убрать единицу из множителя (eps-1) для хорошей обусловленности.
//    Ну с корнем идея сработала, но кажется, что здесь справится диагональный предобуславливатель.
// 3) Сделать блочно-тёплицев формат полностью, оценить степень сжатия матрицы
// 4) Написать умножалку через FFT и сравнить производительность моей умножалки и умножалки через FFT.

// Наблюдения:
// 1. если честно считать eps на кубах, то количество итераций GMRES учаличивается в 2-2.5 раза
//    при этом если делать через домножение на корень, то количество итераций расчет не сильно.
//    Пример: eps = 2, r = 0.25, freq = 1 GHz.
//
// 2.
