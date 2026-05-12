//
// Created by evgen on 11.05.2026.
//

#include "benchmark/benchmark.h"

#include <complex>

#include <random>

namespace {
class BasicMath : public benchmark::Fixture {
};

BENCHMARK_DEFINE_F(BasicMath, DoubleSum)(benchmark::State &state) {
    double v1 = rand() % 42;
    double v2 = rand() % 42;
    for (auto _ : state) {
        double d = v1 + v2;
        benchmark::DoNotOptimize(d);
    }
}
BENCHMARK_REGISTER_F(BasicMath, DoubleSum);


BENCHMARK_DEFINE_F(BasicMath, FloatSum)(benchmark::State &state) {
    float v1 = rand() % 42;
    float v2 = rand() % 42;
    for (auto _ : state) {
        float d = v1 + v2;
        benchmark::DoNotOptimize(d);
    }
}
BENCHMARK_REGISTER_F(BasicMath, FloatSum);

BENCHMARK_DEFINE_F(BasicMath, ComplexDoubleSum)(benchmark::State &state) {
     auto v1 = std::complex<double>{1, static_cast<double>(rand() % 42)};
     auto v2 = std::complex<double>{1, static_cast<double>(rand() % 42)};
    for (auto _ : state) {
        auto d = v1 + v2;
        benchmark::DoNotOptimize(d);
    }
}
BENCHMARK_REGISTER_F(BasicMath, ComplexDoubleSum);

BENCHMARK_DEFINE_F(BasicMath, ComplexFloatSum)(benchmark::State &state) {
    auto v1 = std::complex<float>{1, static_cast<float>(rand() % 42)};
    auto v2 = std::complex<float>{1, static_cast<float>(rand() % 42)};
    for (auto _ : state) {
        auto d = v1 + v2;
        benchmark::DoNotOptimize(d);
    }
}
BENCHMARK_REGISTER_F(BasicMath, ComplexFloatSum);

BENCHMARK_DEFINE_F(BasicMath, DoubleMul)(benchmark::State &state) {
    double v1 = rand() % 42;
    double v2 = rand() % 42;
    for (auto _ : state) {
        double d = v1 * v2;
        benchmark::DoNotOptimize(d);
    }
}
BENCHMARK_REGISTER_F(BasicMath, DoubleMul);

BENCHMARK_DEFINE_F(BasicMath, FloatMul)(benchmark::State &state) {
    float v1 = rand() % 42;
    float v2 = rand() % 42;
    for (auto _ : state) {
        float d = v1 * v2;
        benchmark::DoNotOptimize(d);
    }
}
BENCHMARK_REGISTER_F(BasicMath, FloatMul);

BENCHMARK_DEFINE_F(BasicMath, ComplexDoubleMul)(benchmark::State &state) {
    auto v1 = std::complex<double>{1, static_cast<double>(rand() % 42)};
    auto v2 = std::complex<double>{1, static_cast<double>(rand() % 42)};
    for (auto _ : state) {
        auto d = v1 * v2;
        benchmark::DoNotOptimize(d);
    }
}
BENCHMARK_REGISTER_F(BasicMath, ComplexDoubleMul);

BENCHMARK_DEFINE_F(BasicMath, ComplexFloatMul)(benchmark::State &state) {
    auto v1 = std::complex<float>{1, static_cast<float>(rand() % 42)};
    auto v2 = std::complex<float>{1, static_cast<float>(rand() % 42)};
    for (auto _ : state) {
        auto d = v1 * v2;
        benchmark::DoNotOptimize(d);
    }
}
BENCHMARK_REGISTER_F(BasicMath, ComplexFloatMul);

BENCHMARK_DEFINE_F(BasicMath, ComplexExponent)(benchmark::State &state) {
    auto v1 = std::complex<double>{1, static_cast<double>(rand() % 42)};
    for (auto _ : state) {
        auto d = std::exp(v1);
        benchmark::DoNotOptimize(d);
    }
}
BENCHMARK_REGISTER_F(BasicMath, ComplexExponent);

BENCHMARK_DEFINE_F(BasicMath, RealExponent)(benchmark::State &state) {
    auto v1 = static_cast<double>(rand() % 42);
    for (auto _ : state) {
        auto d = std::exp(v1);
        benchmark::DoNotOptimize(d);
    }
}
BENCHMARK_REGISTER_F(BasicMath, RealExponent);

BENCHMARK_DEFINE_F(BasicMath, ComplexExponentViaReal)(benchmark::State &state) {
    auto v1 = std::complex<double>{1, static_cast<double>(rand() % 42)};
    for (auto _ : state) {
        auto d = std::exp(v1.real()) * std::complex<double>{std::cos(v1.imag()), std::sin(v1.imag())};
        benchmark::DoNotOptimize(d);
    }
}
BENCHMARK_REGISTER_F(BasicMath, ComplexExponentViaReal);

}