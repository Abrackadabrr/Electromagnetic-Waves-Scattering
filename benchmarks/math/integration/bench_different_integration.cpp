//
// Created by evgen on 28.04.2026.
//

#include <benchmark/benchmark.h>

#include "EMW/Utils.hpp"
#include "EMW/mesh/MeshTypes.hpp"

#include "math/integration/NumericalIntegration.hpp"

#include <complex>

using namespace EMW;

namespace {
class IntegrationBenches : public benchmark::Fixture {
};

using namespace std::complex_literals;

Types::complex_d f(Types::point_t p) {
    const double r = p.norm();
    const double k = 10;
    return std::exp(1i * r * k) / r;
}

Types::Vector3c f_vect(const Types::point_t &p) {
    const double r = p.norm();
    const double k = 10;
    Types::Vector3c res = p * std::exp(1i * r * k) / r;
    return res;
}


BENCHMARK_DEFINE_F(IntegrationBenches, OldIntegration)(benchmark::State &state) {
    Containers::vector<Types::point_t> points;
    points.resize(4);
    const double h = 0.1;
    points[0] = Types::point_t{0, 0, h};
    points[1] = Types::point_t{h, 0, h};
    points[2] = Types::point_t{h, h, h};
    points[3] = Types::point_t{0, h, h};
    auto reference_cell = Mesh::IndexedCell{{0, 1, 2, 3}, points};

    const auto integrand = [&reference_cell](double x, double y) {
        return f(reference_cell.parametrization(x, y)) * reference_cell.multiplier(x, y);
    };
    Types::complex_d result{};
    for (auto _ : state) {
        result = DecartIntegration::integrate<DecartIntegration::GaussLegendre::Quadrature<4, 4>>(
            integrand, {0, 0}, {1., 1.});
        benchmark::DoNotOptimize(result);
        benchmark::ClobberMemory();
    }
}

BENCHMARK_DEFINE_F(IntegrationBenches, OldIntegrationVect)(benchmark::State &state) {
    Containers::vector<Types::point_t> points;
    points.resize(4);
    const double h = 0.1;
    points[0] = Types::point_t{0, 0, h};
    points[1] = Types::point_t{h, 0, h};
    points[2] = Types::point_t{h, h, h};
    points[3] = Types::point_t{0, h, h};
    auto reference_cell = Mesh::IndexedCell{{0, 1, 2, 3}, points};

    const auto integrand = [&reference_cell](double x, double y) -> Types::Vector3c {
        const auto point = reference_cell.parametrization(x, y);
        const auto multiplier = reference_cell.multiplier(x, y);
        return f_vect(point) * multiplier;
    };

    Types::Vector3c result{};
    for (auto _ : state) {
        result = DecartIntegration::integrate<DecartIntegration::GaussLegendre::Quadrature<3, 3>>(
            integrand, {0, 0}, {1., 1.});
        benchmark::DoNotOptimize(result);
        benchmark::ClobberMemory();
    }
}

BENCHMARK_DEFINE_F(IntegrationBenches, NewIntegration)(benchmark::State &state) {
    Containers::vector<Types::point_t> points;
    points.resize(4);
    const double h = 0.1;
    points[0] = Types::point_t{0, 0, h};
    points[1] = Types::point_t{h, 0, h};
    points[2] = Types::point_t{h, h, h};
    points[3] = Types::point_t{0, h, h};
    Types::complex_d result{};
    double mes = h * h / 2.;
    for (auto _ : state) {
        result = mes * SimplicialIntegration::quadrature_sum<SimplicialIntegration::GaussianPoints<2, 5>>(f, points[0], points[1], points[2]);
        result += mes * SimplicialIntegration::quadrature_sum<SimplicialIntegration::GaussianPoints<2, 5>>(f, points[2], points[3], points[0]);
        benchmark::DoNotOptimize(result);
        benchmark::ClobberMemory();
    }
}

BENCHMARK_DEFINE_F(IntegrationBenches, NewIntegrationVect)(benchmark::State &state) {
    Containers::vector<Types::point_t> points;
    points.resize(4);
    const double h = 0.1;
    points[0] = Types::point_t{0, 0, h};
    points[1] = Types::point_t{h, 0, h};
    points[2] = Types::point_t{h, h, h};
    points[3] = Types::point_t{0, h, h};
    Types::Vector3c result{};
    double mes = h * h / 2.;
    for (auto _ : state) {
        result = mes * SimplicialIntegration::quadrature_sum<SimplicialIntegration::GaussianPoints<2, 5>>(f_vect, points[0], points[1], points[2]);
        result += mes * SimplicialIntegration::quadrature_sum<SimplicialIntegration::GaussianPoints<2, 5>>(f_vect, points[2], points[3], points[0]);
        benchmark::DoNotOptimize(result);
        benchmark::ClobberMemory();
    }
}

}

BENCHMARK_REGISTER_F(IntegrationBenches, NewIntegration);
BENCHMARK_REGISTER_F(IntegrationBenches, OldIntegration);
BENCHMARK_REGISTER_F(IntegrationBenches, NewIntegrationVect);
BENCHMARK_REGISTER_F(IntegrationBenches, OldIntegrationVect);

BENCHMARK_MAIN();
