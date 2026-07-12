//
// Created by evgen on 10.07.2026.
//

#include <benchmark/benchmark.h>

#include "math/fourier/TripleToeplitz3x3Fourier.hpp"

#include <omp.h>

#include <cmath>

using namespace EMW;

namespace {

using fourier_t = Math::Fourier::TripleToeplitz3x3FourierParallel<Types::complex_d>;
using tensor_t = fourier_t::tensor_type;

constexpr Types::index CELLS_PER_AXIS = 100;

tensor_t make_test_levels(Types::index cells_per_axis) {
    const Types::index levels_per_axis = 2 * cells_per_axis - 1;
    const Types::integer center = static_cast<Types::integer>(cells_per_axis - 1);

    tensor_t levels(levels_per_axis, levels_per_axis, levels_per_axis, Types::complex_d{0.0, 0.0});

    for (Types::index lz = 0; lz < levels_per_axis; ++lz) {
        const auto dz = static_cast<Types::integer>(lz) - center;
        for (Types::index ly = 0; ly < levels_per_axis; ++ly) {
            const auto dy = static_cast<Types::integer>(ly) - center;
            for (Types::index lx = 0; lx < levels_per_axis; ++lx) {
                const auto dx = static_cast<Types::integer>(lx) - center;
                const auto distance_sqr = static_cast<Types::scalar>(1 + dx * dx + dy * dy + dz * dz);
                const auto kernel_scale = 1.0 / distance_sqr;

                for (Types::index row = 0; row < 3; ++row) {
                    for (Types::index col = 0; col < 3; ++col) {
                        const auto component_scale = static_cast<Types::scalar>((row + 1) * (col + 2));
                        const auto phase_like =
                            static_cast<Types::scalar>((row + 1) * dx - (col + 1) * dy + (row + col + 1) * dz);
                        levels(lx, ly, lz, row, col) =
                            Types::complex_d{component_scale * kernel_scale, 0.01 * phase_like * kernel_scale};
                    }
                }
            }
        }
    }

    return levels;
}

class TripleToeplitz3x3FourierParallelAccelerationBench : public benchmark::Fixture {
  public:
    void SetUp(const benchmark::State & /*state*/) override {
        Eigen::setNbThreads(1);

        fourier_ = fourier_t(make_test_levels(CELLS_PER_AXIS));
        x_ = Types::VectorXc::Random(fourier_.cols());
        y_ = Types::VectorXc::Zero(fourier_.rows());
    }

  protected:
    fourier_t fourier_;
    Types::VectorXc x_;
    Types::VectorXc y_;
};

BENCHMARK_DEFINE_F(TripleToeplitz3x3FourierParallelAccelerationBench, FourierParallelMatVecAcceleration)
(benchmark::State &state) {
    int n_th = static_cast<int>(state.range(0));

    y_ = fourier_.matvec_fftw_threads(x_, n_th);
    benchmark::DoNotOptimize(y_.data());

    for (auto _ : state) {
        y_ = fourier_.matvec_fftw_threads(x_, n_th);
        benchmark::DoNotOptimize(y_.data());
        benchmark::ClobberMemory();
    }
}

BENCHMARK_REGISTER_F(TripleToeplitz3x3FourierParallelAccelerationBench, FourierParallelMatVecAcceleration)
    ->Arg(1)
    ->Arg(2)
    ->Arg(3)
    ->Arg(4)
    ->Arg(5)
    ->Arg(6)
    ->Arg(7)
    ->Arg(8)
    ->Arg(9)
    ->Arg(10)
    ->Arg(11)
    ->Arg(12)
    ->Unit(benchmark::kMillisecond)
    ->UseRealTime();

} // namespace

BENCHMARK_MAIN();
