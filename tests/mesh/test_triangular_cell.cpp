//
// Created by evgen on 02.07.24.
//
#include "gtest/gtest.h"

#include "operators/Operators.hpp"
#include "operators/Functions.hpp"
#include "types/Types.hpp"
#include "mesh/MeshTypes.hpp"
#include "integration/gauss_quadrature/GaussLegenderPoints.hpp"
#include "integration/newton_cotess/Rectangular.hpp"

using namespace EMW;

class TRIANGULAR_CELL_TESTS : public testing::Test {
};

TEST_F(TRIANGULAR_CELL_TESTS, METADATA) {
    // объявляем ячейку, по которой происходит интегрирование
    const Containers::vector<Mesh::point_t> points = {
            Mesh::point_t{0, 0, 0},
            Mesh::point_t{1, 0, 0},
            Mesh::point_t{1. / 2, std::sqrt(3) / 2, 0},
            Mesh::point_t{1. / 2, std::sqrt(3) / 2, 0}
    };
    Mesh::IndexedCell cell{{0, 1, 2, 3}, points};

    // расчетные параметры для треугольной ячейки
    const Types::scalar area = std::sqrt(3) / 4;
    const Types::Vector3d normal{0, 0, 1};
    const Types::Vector3d ort1 = points[1] - points[0];
    const Types::Vector3d ort2 = points[3] - points[0];
    const Types::Vector3d diff = -ort1;
    const Types::Vector3d a = ort1.cross(ort2);
    const Types::Vector3d c = -a;

    ASSERT_NEAR(cell.area_, area, 1e-16);
    ASSERT_NEAR((normal - cell.normal).norm(), 0, 1e-16);

    // ортонормированность локального базиса
    ASSERT_NEAR(cell.tau[0].dot(cell.tau[1]), 0, 1e-16);
    ASSERT_NEAR(cell.tau[0].dot(cell.normal), 0, 1e-16);
    ASSERT_NEAR(cell.normal.dot(cell.tau[1]), 0, 1e-16);
    ASSERT_NEAR(cell.tau[0].norm(), 1, 1e-16);
    ASSERT_NEAR(cell.tau[1].norm(), 1, 1e-16);

    // структура ячейки для параметризации (Cell Struncture)
    ASSERT_NEAR((cell.cellStructure.ort1 - ort1).norm(), 0, 1e-16);
    ASSERT_NEAR((cell.cellStructure.ort2 - ort2).norm(), 0, 1e-16);
    ASSERT_NEAR((cell.cellStructure.diff - diff).norm(), 0, 1e-16);

    // параметры, используемые для интегрирования (Integration Parameters)
    ASSERT_NEAR(cell.integrationParameters.b.norm(), 0, 1e-16);
    ASSERT_NEAR((cell.integrationParameters.a - a).norm(), 0, 1e-16);
    ASSERT_NEAR((cell.integrationParameters.c - c).norm(), 0, 1e-16);
}

TEST_F(TRIANGULAR_CELL_TESTS, SURFACE_INTEGRATION_CONSTANT) {
    // объявляем ячейку, по которой происходит интегрирование
    const Containers::vector<Mesh::point_t> points = {
            Mesh::point_t{0, 0, 0},
            Mesh::point_t{1, 0, 0},
            Mesh::point_t{1, 1, 0},
            Mesh::point_t{1, 1, 1},
            Mesh::point_t{3, 1, 5},
            Mesh::point_t{3, 1.2, std::sqrt(12)}
    };
    Mesh::IndexedCell cell[4] = {{{0, 1, 2, 2}, points},
                                 {{0, 3, 5, 5}, points},
                                 {{0, 2, 4, 4}, points},
                                 {{0, 4, 1, 1}, points}};

    // интегрируем константу
    const auto constant = [](Types::scalar x, Types::scalar y) -> Types::scalar { return 1; };
    // парамтеризация функции
    const auto phi1 = [&](Types::scalar p, Types::scalar q) -> Types::scalar {
        const Types::Vector3d y = cell[0].parametrization(p, q);
        const Types::scalar mul = cell[0].multiplier(p, q);
        return constant(y.x(), y.y()) * mul;
    };
    const auto phi2 = [&](Types::scalar p, Types::scalar q) -> Types::scalar {
        const Types::Vector3d y = cell[1].parametrization(p, q);
        const Types::scalar mul = cell[1].multiplier(p, q);
        return constant(y.x(), y.y()) * mul;
    };
    const auto phi3 = [&](Types::scalar p, Types::scalar q) -> Types::scalar {
        const Types::Vector3d y = cell[2].parametrization(p, q);
        const Types::scalar mul = cell[2].multiplier(p, q);
        return constant(y.x(), y.y()) * mul;
    };
    const auto phi4 = [&](Types::scalar p, Types::scalar q) -> Types::scalar {
        const Types::Vector3d y = cell[3].parametrization(p, q);
        const Types::scalar mul = cell[3].multiplier(p, q);
        return constant(y.x(), y.y()) * mul;
    };
    const Types::scalar result1 =
            DefiniteIntegrals::integrate<DefiniteIntegrals::GaussLegendre::Quadrature<8, 8>>(phi1, {0, 0}, {1., 1.});
    const Types::scalar result2 =
            DefiniteIntegrals::integrate<DefiniteIntegrals::GaussLegendre::Quadrature<8, 8>>(phi2, {0, 0}, {1., 1.});
    const Types::scalar result3 =
            DefiniteIntegrals::integrate<DefiniteIntegrals::GaussLegendre::Quadrature<8, 8>>(phi3, {0, 0}, {1., 1.});
    const Types::scalar result4 =
            DefiniteIntegrals::integrate<DefiniteIntegrals::GaussLegendre::Quadrature<8, 8>>(phi4, {0, 0}, {1., 1.});
    ASSERT_NEAR(result1, cell[0].area_, 1e-15);
    ASSERT_NEAR(result2, cell[1].area_, 1e-15);
    ASSERT_NEAR(result3, cell[2].area_, 1e-15);
    ASSERT_NEAR(result4, cell[3].area_, 1e-15);
}

TEST_F(TRIANGULAR_CELL_TESTS, SURFACE_INTEGRATION) {
    // объявляем ячейку, по которой происходит интегрирование
    const Containers::vector<Mesh::point_t> points = {
            Mesh::point_t{0, 0, 0},
            Mesh::point_t{1, 0, 0},
            Mesh::point_t{1, 1, 0},
    };
    Mesh::IndexedCell cell{{0, 1, 2, 2}, points};

    // интегрируем функцию
    const auto constant = [](Types::scalar x, Types::scalar y) -> Types::scalar {
        return x * x * x * y * std::sin(y * y / 2) * std::sin(y * y / 2);
    };
    // парамтеризация функции
    const auto phi = [&](Types::scalar p, Types::scalar q) -> Types::scalar {
        const Types::Vector3d y = cell.parametrization(p, q);
        const Types::scalar mul = cell.multiplier(p, q);
        return constant(y.x(), y.y()) * mul;
    };
    const Types::scalar result =
            DefiniteIntegrals::integrate<DefiniteIntegrals::GaussLegendre::Quadrature<12, 12>>(phi, {0, 0}, {1., 1.});
    ASSERT_NEAR(result, 1./24 - 1./8 * (std::sin(1) - std::cos(1)), 1e-16);
}

TEST_F(TRIANGULAR_CELL_TESTS, CONTOUR_PARAMETRIZATION) {
    const Containers::vector<Mesh::point_t> points = {
            Mesh::point_t{0, 0, 0},
            Mesh::point_t{1, 0, 0},
            Mesh::point_t{1, 1, 0},
            Mesh::point_t{1, 1, 0}
    };
    Mesh::IndexedCell cell{{0, 1, 2, 3}, points};

    // параметризация контура CD не зависит от параметра и возвращает точку C (== D)
    const auto CD = [&](Types::scalar t) -> Types::Vector3d {
        const Types::Vector3d y = cell.parametrization(1 - t, 1);
        return y;
    };
    int N = 100;
    for (int i = 0; i != N; i++) {
        ASSERT_NEAR((CD(i * 1./N) - points[2]).norm(), 0, 2e-16);
    }

    // проверка мультипликаторов
    const Types::Vector3d mul[3] = {
            {0, -1, 0},
            {1, 0, 0},
            {-1, 1, 0}
    };
    // мультипликатор[2] контурного интеграла равен нулю
    ASSERT_NEAR((cell.integrationParameters.mul[2]).norm(), 0, 1e-16);
    // проверка других мультипликаторов
    ASSERT_NEAR((cell.integrationParameters.mul[0] - mul[0]).norm(), 0, 1e-16);
    ASSERT_NEAR((cell.integrationParameters.mul[1] - mul[1]).norm(), 0, 1e-16);
    ASSERT_NEAR((cell.integrationParameters.mul[3] - mul[2]).norm(), 0, 1e-16);
}

TEST_F(TRIANGULAR_CELL_TESTS, CONTOUR_INTEGRATION) {
    // объявляем ячейку, по которой происходит интегрирование
    const Containers::vector<Mesh::point_t> points = {
            Mesh::point_t{0, 0, 0},
            Mesh::point_t{1, 0, 0},
            Mesh::point_t{1, 1, 0},
            Mesh::point_t{1, 1, 0}
    };
    Mesh::IndexedCell cell{{0, 1, 2, 3}, points};

    const auto AB = [&](Types::scalar t) -> Types::Vector3d {
        const Types::Vector3d y = cell.parametrization(t, 0);
        return cell.integrationParameters.mul[0].normalized();
    };
    const auto BC = [&](Types::scalar t) -> Types::Vector3d {
        const Types::Vector3d y = cell.parametrization(1, t);
        return cell.integrationParameters.mul[1].normalized();
    };
    const auto CD = [&](Types::scalar t) -> Types::Vector3d {
        const Types::Vector3d y = cell.parametrization(1 - t, 1);
        return Types::Vector3d{1, 0, 0};
    };
    const auto DA = [&](Types::scalar t) -> Types::Vector3d {
        const Types::Vector3d y = cell.parametrization(0, 1 - t);
        return cell.integrationParameters.mul[3].normalized();
    };

    // сравнение интеграла по трем и четырем отрезкам


}
