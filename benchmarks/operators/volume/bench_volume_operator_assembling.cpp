//
// Created by evgen on 31.03.2026.
//

#include <benchmark/benchmark.h>

#include "mesh/volume_mesh/CubeMeshWithData.hpp"

#include "operators/volume/OperatorK.hpp"

#include "experiment/PhysicalCondition.hpp"

#include "Utils.hpp"

#include <omp.h>

using namespace EMW;

class OperatorKAssemblingBench : public benchmark::Fixture {
  public:
    void SetUp(const benchmark::State &state) override {
        // 0. Настраиваемые параметры
        constexpr Types::index Nx = 9;
        constexpr Types::scalar freq = 0.3; // GHz
        constexpr Types::scalar rTol = 1e-3;
        constexpr Types::scalar aTol = 1e-21;
        constexpr Types::index lev_2d = 1;
        constexpr Types::index lev_3d = 1;
        constexpr Types::index lev_4d = 1;
        constexpr Types::index lev_6d = 1;
        // 1. Сбор сетки
        Eigen::setNbThreads(1);
        const Types::index Ny = Nx;
        const Types::index Nz = Nx;
        const Types::point_t min_corner = Types::point_t{-cube_length / 2, -cube_length / 2, -cube_length / 2};
        mesh_ptr = new Mesh::VolumeMesh::CubeMeshWithData{min_corner, cube_length, cube_length, cube_length, Nx, Ny, Nz};
        const Types::scalar cube_measure = mesh_ptr->dx() * mesh_ptr->dy() * mesh_ptr->dz();
        basis_fn_module = 1. / sqrt(cube_measure);

        // 2. Параметры падающей волны
        constexpr Types::complex_d k{Physics::get_k_on_frquency(freq), 0.0};

        // 3. Галеркинская проекция оператора
        operator_k = new Operators::Volume::operator_K_over_cube_mesh{k, *mesh_ptr};
        operator_k->set_tolerances(rTol, aTol);
        operator_k->set_adaptive_integration_max_levels({lev_2d, lev_3d, lev_4d, lev_6d});
    }

  protected:
    constexpr static Types::scalar SPHERE_RADUIS = 0.5;
    constexpr static Types::scalar CUBE_LENGTH = 0.5;
    constexpr static Types::scalar SPHERE_EPSILON = 2.56;
    constexpr static Types::scalar cube_length = 2 * SPHERE_RADUIS;

    Types::scalar basis_fn_module;

    Mesh::VolumeMesh::CubeMeshWithData* mesh_ptr = nullptr;
    Operators::Volume::operator_K_over_cube_mesh *operator_k = nullptr;
};

/**
 * Измеряем ускорение расчета матрицы с увеличением числа потоков omp
 */
BENCHMARK_DEFINE_F(OperatorKAssemblingBench, AssemblingAcceleration)(benchmark::State &state) {
    omp_set_num_threads(state.range(0));
    auto warming_result = operator_k->compute_galerkin_matrix(basis_fn_module);
    benchmark::DoNotOptimize(warming_result);
    for (auto _ : state) {
        auto result = operator_k->compute_galerkin_matrix(basis_fn_module);
        benchmark::DoNotOptimize(result);
        benchmark::ClobberMemory();
    }
}

// Полная сборка матрицы оператора К
BENCHMARK_DEFINE_F(OperatorKAssemblingBench, SimpleAssembling)(benchmark::State &state) {
    omp_set_num_threads(1);
    auto warming_result = operator_k->compute_galerkin_matrix(basis_fn_module);
    benchmark::DoNotOptimize(warming_result);
    for (auto _ : state) {
        auto result = operator_k->compute_galerkin_matrix(basis_fn_module);
        benchmark::DoNotOptimize(result);
        benchmark::ClobberMemory();
    }
}

BENCHMARK_REGISTER_F(OperatorKAssemblingBench, SimpleAssembling)->Iterations(5)->Unit(benchmark::kMillisecond)->UseRealTime();

BENCHMARK_REGISTER_F(OperatorKAssemblingBench, AssemblingAcceleration)->Arg(1)->Arg(2)->Arg(4)->Arg(8)->Arg(12)
->Unit(benchmark::kMillisecond)->UseRealTime()->Iterations(3);

BENCHMARK_MAIN();
