//
// Created by evgen on 12.05.2026.
//

#include <benchmark/benchmark.h>

#include "EMW/Utils.hpp"

#include "EMW/mesh/MeshTypes.hpp"
#include "EMW/mesh/Parser.hpp"
#include "EMW/mesh/SurfaceMesh.hpp"
#include "EMW/mesh/Utils.hpp"

#include "math/integration/NumericalIntegration.hpp"
#include "math/integration/analytical/SingularIntegration.hpp"

#include <complex>
#include <ranges>

using namespace EMW;

namespace {

inline EMW::Mesh::SurfaceMesh generateRectangularMesh(int N1, int N2, Types::scalar h1, Types::scalar h2) {
    std::vector<Mesh::point_t> meshgrid;
    meshgrid.reserve(N1 * N2);
    Mesh::Utils::cartesian_product_unevenYZ(std::ranges::views::iota(0, N1), std::ranges::views::iota(0, N2),
                                            std::back_inserter(meshgrid), N1, N2, h1, h2);

    const auto cellsView = std::views::iota(0, (N1 - 1) * (N2 - 1)) | std::views::transform([N1](int index) {
                               Types::index i =
                                   index + index / (N1 - 1); // левый нижний индекс зависит от номера СТРОКИ
                               // в которой располагается ячейка, поэтому
                               // здесь стоит N1
                               const auto point = Mesh::IndexedCell::nodes_t{i, i + 1, i + 1 + N1, i + N1};
                               return point;
                           });

    const auto cells =
        Containers::vector<Mesh::IndexedCell::nodes_t>{std::ranges::begin(cellsView), std::ranges::end(cellsView)};
    return Mesh::SurfaceMesh{meshgrid, cells};
}

class MeshStorageBenchmarks : public benchmark::Fixture {
  protected:
    Containers::vector<Mesh::IndexedCell> indexed_cells;
    Containers::vector<Types::point_t> points;
    Containers::vector<Mesh::Cell> cells;

    void SetUp(benchmark::State &st) override {
        auto mesh = generateRectangularMesh(1000, 1000, 0.001, 0.0001);
        indexed_cells = mesh.getCells();
        points = mesh.getNodes();
        cells.resize(indexed_cells.size());
        std::ranges::transform(indexed_cells, cells.begin(),
                               [](const Mesh::IndexedCell &cell) { return cell.getVertex(); });
    }
};

// интегрирование через формальное сведение поверхностного интеграла к интегралу Римана
BENCHMARK_DEFINE_F(MeshStorageBenchmarks, VeryOldIntegration)(benchmark::State &state) {
    // будем по сетке считать интегралы подряд
    Containers::vector<Types::scalar> result;
    result.resize(indexed_cells.size());
    // интегранд
    const auto integrand = [](const Types::point_t &point) { return std::exp(point.norm()); };
    // замер времени
    for (auto _ : state) {
        for (size_t i = 0; i < indexed_cells.size(); i++) {
            const auto local_integrand = [cell = indexed_cells[i], f = integrand](Types::scalar x, Types::scalar y) {
                return f(cell.parametrization(x, y)) * cell.multiplier(x, y);
            };
            result[i] = Math::Integration::Numerical::Decart::integrate<
                Math::Integration::Numerical::Decart::GaussLegendre::Quadrature<4, 4>>(local_integrand, {0., 0.}, {1., 1.});
        }
        benchmark::DoNotOptimize(result.data());
    }
}

// интегрирование через правильную квадратуру на треугольнике,
// с доставанием точек из общего массива
BENCHMARK_DEFINE_F(MeshStorageBenchmarks, NewIntegrationOldCell)(benchmark::State &state) {
    // будем по сетке считать интегралы подряд
    Containers::vector<Types::scalar> result;
    result.resize(indexed_cells.size());
    // интегранд
    const auto integrand = [](const Types::point_t &point) { return std::exp(point.norm()); };
    // замер времени
    for (auto _ : state) {
        for (size_t i = 0; i < indexed_cells.size(); i++) {

            result[i] = Math::Integration::Numerical::Simplicial::quadrature_sum<
                Math::Integration::Numerical::Simplicial::GaussianPoints<2, 5>>(
                integrand, points[indexed_cells[i].points_[0]], points[indexed_cells[i].points_[1]],
                points[indexed_cells[i].points_[2]]);
            result[i] = Math::Integration::Numerical::Simplicial::quadrature_sum<
                Math::Integration::Numerical::Simplicial::GaussianPoints<2, 5>>(
                integrand, points[indexed_cells[i].points_[2]], points[indexed_cells[i].points_[3]],
                points[indexed_cells[i].points_[0]]);

            result[i] += 0.5 * integrand({0, 0, 0});
            result[i] += 0.5 * integrand({0, 0, 1.});
        }
        benchmark::DoNotOptimize(result.data());
    }
}

// интегрирование через правильную квадратуру на треугольнике,
// с расчетом точек на ходу средствами IndexedCell.
BENCHMARK_DEFINE_F(MeshStorageBenchmarks, NewIntegrationOldCellWithPointsCalculation)(benchmark::State &state) {
    // будем по сетке считать интегралы подряд
    Containers::vector<Types::scalar> result;
    result.resize(indexed_cells.size());
    // интегранд
    const auto integrand = [](const Types::point_t &point) { return std::exp(point.norm()); };
    // замер времени
    for (auto _ : state) {
        for (size_t i = 0; i < indexed_cells.size(); i++) {
            auto local_points = indexed_cells[i].getVertexAsArray();
            result[i] = Math::Integration::Numerical::Simplicial::quadrature_sum<
                Math::Integration::Numerical::Simplicial::GaussianPoints<2, 5>>(integrand, local_points[0],
                                                                                local_points[1], local_points[2]);
            result[i] = Math::Integration::Numerical::Simplicial::quadrature_sum<
                Math::Integration::Numerical::Simplicial::GaussianPoints<2, 5>>(integrand, local_points[2],
                                                                                local_points[3], local_points[0]);

            result[i] += 0.5 * integrand({0, 0, 0});
            result[i] += 0.5 * integrand({0, 0, 1.});
        }
        benchmark::DoNotOptimize(result.data());
    }
}

// Новое интегрирование по хранщему Cell
// с правильной квадратурой на треугольнике
BENCHMARK_DEFINE_F(MeshStorageBenchmarks, NewIntegration)(benchmark::State &state) {
    // будем по сетке считать интегралы подряд
    Containers::vector<Types::scalar> result;
    result.resize(indexed_cells.size());
    // интегранд
    const auto integrand = [](const Types::point_t &point) { return std::exp(point.norm()); };
    // замер времени
    for (auto _ : state) {
        for (size_t i = 0; i < indexed_cells.size(); i++) {

            result[i] = Math::Integration::Numerical::Simplicial::quadrature_sum<
                Math::Integration::Numerical::Simplicial::GaussianPoints<2, 5>>(integrand, cells[i].a, cells[i].b,
                                                                                cells[i].c);
            result[i] = Math::Integration::Numerical::Simplicial::quadrature_sum<
                Math::Integration::Numerical::Simplicial::GaussianPoints<2, 5>>(integrand, cells[i].c, cells[i].d,
                                                                                cells[i].a);

            result[i] += 0.5 * integrand({0, 0, 0});
            result[i] += 0.5 * integrand({0, 0, 1.});
        }
        benchmark::DoNotOptimize(result.data());
    }
}

BENCHMARK_DEFINE_F(MeshStorageBenchmarks, OldIntegrationAn)(benchmark::State &state) {
    // будем по сетке считать интегралы подряд
    Containers::vector<Types::scalar> result;
    result.resize(indexed_cells.size());
    Types::point_t ref_p{1., 0., -10.};
    // замер времени
    for (auto _ : state) {
        for (size_t i = 0; i < indexed_cells.size(); i++) {

            result[i] = Math::Integration::Analytical::integrate_1_div_r(ref_p, indexed_cells[i]);
        }
        benchmark::DoNotOptimize(result.data());
    }
}

BENCHMARK_DEFINE_F(MeshStorageBenchmarks, NewIntegrationAn)(benchmark::State &state) {
    // будем по сетке считать интегралы подряд
    Containers::vector<Types::scalar> result;
    result.resize(indexed_cells.size());
    Types::point_t ref_p{1., 0., -10.};
    // замер времени
    for (auto _ : state) {
        for (size_t i = 0; i < indexed_cells.size(); i++) {

            result[i] = Math::Integration::Analytical::integrate_1_div_r(
                ref_p, Containers::array{cells[i].a, cells[i].b, cells[i].c, cells[i].d});
        }
        benchmark::DoNotOptimize(result.data());
    }
}

} // namespace

BENCHMARK_REGISTER_F(MeshStorageBenchmarks, VeryOldIntegration)->Unit(benchmark::kMillisecond);
BENCHMARK_REGISTER_F(MeshStorageBenchmarks, NewIntegrationOldCell)->Unit(benchmark::kMillisecond);
BENCHMARK_REGISTER_F(MeshStorageBenchmarks, NewIntegrationOldCellWithPointsCalculation)->Unit(benchmark::kMillisecond);
BENCHMARK_REGISTER_F(MeshStorageBenchmarks, NewIntegration)->Unit(benchmark::kMillisecond);
BENCHMARK_REGISTER_F(MeshStorageBenchmarks, OldIntegrationAn)->Unit(benchmark::kMillisecond);
BENCHMARK_REGISTER_F(MeshStorageBenchmarks, NewIntegrationAn)->Unit(benchmark::kMillisecond);
