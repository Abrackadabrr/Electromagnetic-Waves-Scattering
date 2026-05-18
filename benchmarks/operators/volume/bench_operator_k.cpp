//
// Created by evgen on 31.03.2026.
//

#include <benchmark/benchmark.h>

#include "mesh/volume_mesh/CubeMeshWithData.hpp"

#include "operators/volume/OperatorK.hpp"

#include <Utils.hpp>

using namespace EMW;

class OperatorKBench : public benchmark::Fixture {
public:
    void SetUp(const benchmark::State &state) override {
        Eigen::setNbThreads(1);

        constexpr Types::index nodes_per_axis = 81;
        constexpr Types::scalar cube_length = 2;
        const Types::scalar mesh_one_axis_size = cube_length / (nodes_per_axis - 1);
        const Types::scalar basis_fn_norm = 1. / (mesh_one_axis_size * std::sqrt(mesh_one_axis_size));

        mesh_ = Mesh::VolumeMesh::CubeMeshWithData(
            Types::point_t{-cube_length / 2, -cube_length / 2, -cube_length / 2},
            (nodes_per_axis - 1) * mesh_one_axis_size,
            (nodes_per_axis - 1) * mesh_one_axis_size,
            (nodes_per_axis - 1) * mesh_one_axis_size,
            nodes_per_axis, nodes_per_axis, nodes_per_axis);

        constexpr Types::complex_d k{1, 0};

        operator_k = new Operators::Volume::operator_K_over_cube_mesh(k, mesh_);

        my_idx = 0;
        near_idx = mesh_.cube_idx(1, 1, 1);
        middle_idx = mesh_.cube_idx(mesh_.nCubesX() / 2, 0, 0);
        far_idx = mesh_.cube_idx(mesh_.nCubesX() - 1, mesh_.nCubesY() - 1, mesh_.nCubesZ() - 1);
    }

protected:
    Mesh::VolumeMesh::CubeMeshWithData mesh_;
    Operators::Volume::operator_K_over_cube_mesh* operator_k = nullptr;

    size_t my_idx, near_idx, middle_idx, far_idx;
};

/**
 * Измеряем время, которое затрачивается на расчет матричного элемента для влияния куба самого на себя
 */
BENCHMARK_DEFINE_F(OperatorKBench, SelfInteraction)(benchmark::State &state) {
    for (auto _ : state) {
        auto result = operator_k->galerkin_block_for_cubes(my_idx, my_idx);
        benchmark::DoNotOptimize(result.data());
        benchmark::ClobberMemory();
    }
}

/**
 * Измеряем время, которое затрачивается на расчет матричного элемента для влияния куба соседний по ребру
 */
BENCHMARK_DEFINE_F(OperatorKBench, NearInteraction)(benchmark::State &state) {
    for (auto _ : state) {
        auto result = operator_k->galerkin_block_for_cubes(my_idx, near_idx);
        benchmark::DoNotOptimize(result.data());
        benchmark::ClobberMemory();
    }
}


/**
 * Измеряем время, которое затрачивается на расчет матричного элемента для влияния недостаточно далёкий куб, но и
 * недостаточно близкий куб (то есть через преобразование сингулярного интеграла, но без выделения особенности)
 */
BENCHMARK_DEFINE_F(OperatorKBench, MiddleInteraction)(benchmark::State &state) {
    for (auto _ : state) {
        auto result = operator_k->galerkin_block_for_cubes(my_idx, middle_idx);
        benchmark::DoNotOptimize(result.data());
        benchmark::ClobberMemory();
    }
}

/**
 * Измеряем время, которое затрачивается на расчет матричного элемента для влияния достаточно далёкий куб,
 * то есть интегрирование гладкого ядра в интегральном операторе
 */
BENCHMARK_DEFINE_F(OperatorKBench, FatInteraction)(benchmark::State &state) {
    for (auto _ : state) {
        auto result = operator_k->galerkin_block_for_cubes(my_idx, far_idx);
        benchmark::DoNotOptimize(result.data());
        benchmark::ClobberMemory();
    }
}


BENCHMARK_DEFINE_F(OperatorKBench, SelfMatrix2)(benchmark::State &state) {
    for (auto _ : state) {
        auto result = operator_k->matrix_2_coef(my_idx, my_idx);
        benchmark::DoNotOptimize(result.data());
        benchmark::ClobberMemory();
    }
}


BENCHMARK_DEFINE_F(OperatorKBench, NearMatrix2)(benchmark::State &state) {
    for (auto _ : state) {
        auto result = operator_k->matrix_2_coef(my_idx, near_idx);
        benchmark::DoNotOptimize(result.data());
        benchmark::ClobberMemory();
    }
}


BENCHMARK_DEFINE_F(OperatorKBench, MiddleMatrix2)(benchmark::State &state) {
    for (auto _ : state) {
        auto result = operator_k->matrix_2_coef(my_idx, middle_idx);
        benchmark::DoNotOptimize(result.data());
        benchmark::ClobberMemory();
    }
}


BENCHMARK_DEFINE_F(OperatorKBench, FarMatrix2)(benchmark::State &state) {
    for (auto _ : state) {
        auto result = operator_k->matrix_2_coef(my_idx, far_idx);
        benchmark::DoNotOptimize(result.data());
        benchmark::ClobberMemory();
    }
}

BENCHMARK_DEFINE_F(OperatorKBench, SelfMatrix3)(benchmark::State &state) {
    for (auto _ : state) {
        auto result = operator_k->matrix_3_coef(my_idx, my_idx);
        benchmark::DoNotOptimize(result);
        benchmark::ClobberMemory();
    }
}


BENCHMARK_DEFINE_F(OperatorKBench, NearMatrix3)(benchmark::State &state) {
    for (auto _ : state) {
        auto result = operator_k->matrix_3_coef(my_idx, near_idx);
        benchmark::DoNotOptimize(result);
        benchmark::ClobberMemory();
    }
}


BENCHMARK_DEFINE_F(OperatorKBench, MiddleMatrix3)(benchmark::State &state) {
    for (auto _ : state) {
        auto result = operator_k->matrix_3_coef(my_idx, middle_idx);
        benchmark::DoNotOptimize(result);
        benchmark::ClobberMemory();
    }
}


BENCHMARK_DEFINE_F(OperatorKBench, FarMatrix3)(benchmark::State &state) {
    for (auto _ : state) {
        auto result = operator_k->matrix_3_coef(my_idx, far_idx);
        benchmark::DoNotOptimize(result);
        benchmark::ClobberMemory();
    }
}

BENCHMARK_REGISTER_F(OperatorKBench, SelfInteraction)->Unit(benchmark::kMillisecond);
BENCHMARK_REGISTER_F(OperatorKBench, NearInteraction)->Unit(benchmark::kMillisecond);
BENCHMARK_REGISTER_F(OperatorKBench, MiddleInteraction)->Unit(benchmark::kMillisecond);
BENCHMARK_REGISTER_F(OperatorKBench, FatInteraction)->Unit(benchmark::kMillisecond);

BENCHMARK_REGISTER_F(OperatorKBench, SelfMatrix2)->Unit(benchmark::kMillisecond);
BENCHMARK_REGISTER_F(OperatorKBench, NearMatrix2)->Unit(benchmark::kMillisecond);
BENCHMARK_REGISTER_F(OperatorKBench, MiddleMatrix2)->Unit(benchmark::kMillisecond);
BENCHMARK_REGISTER_F(OperatorKBench, FarMatrix2)->Unit(benchmark::kMillisecond);

BENCHMARK_REGISTER_F(OperatorKBench, SelfMatrix3)->Unit(benchmark::kMillisecond);
BENCHMARK_REGISTER_F(OperatorKBench, NearMatrix3)->Unit(benchmark::kMillisecond);
BENCHMARK_REGISTER_F(OperatorKBench, MiddleMatrix3)->Unit(benchmark::kMillisecond);
BENCHMARK_REGISTER_F(OperatorKBench, FarMatrix3)->Unit(benchmark::kMillisecond);

BENCHMARK_MAIN();
