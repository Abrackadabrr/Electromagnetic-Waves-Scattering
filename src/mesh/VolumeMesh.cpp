//
// Created by evgen on 08.02.24.
//

#include "VolumeMesh.hpp"
#include "slae_generation/Functions.hpp"
#include "integration/Quadrature.hpp"
#include "integration/gauss_quadrature/GaussLegenderPoints.hpp"
#include "math/MathConstants.hpp"

namespace EMW::Mesh {
    void EMW::Mesh::VolumeMesh::calculateAll(const Types::Vector3d &polarization, const Types::Vector3d &k_vec,
                                             Types::complex_d k) {
        for (auto &node: nodes_) {
            node.E_ = (getZeroPartIntegral(node.point_, k) + getFirstPartIntegral(node.point_, k) +
                       polarization * std::exp(Math::Constants::i * k_vec.dot(node.point_))).real();
        }
    }

    void calculateESS(const Types::Vector3d &polarization, const Types::Vector3d &k_vec,
                      Types::complex_d k) {

    }

    Types::Vector3c VolumeMesh::getFirstPartIntegral(const Point &point, Types::complex_d k) {
        Types::Vector3c result{0, 0, 0};
        // суммирование по всем ячейкам
        for (const auto &cell: surfaceMesh_.getCells()) {
            // расчет на одной ячейке сетки
            const auto phi = [&](Types::scalar p, Types::scalar q) -> Types::complex_d {
                const Types::Vector3d y = cell.parametrization(p, q);
                const Types::scalar mul = cell.multiplier(p, q);
                return Helmholtz::F(k, point, y) * mul;
            };
            result += cell.collPoint_.J_ *
                      DefiniteIntegrals::integrate<DefiniteIntegrals::GaussLegendre::Quadrature<12, 12>>(phi, {0, 0}, {1., 1.});
        }
        return k * k * result;
    }

    Types::Vector3c VolumeMesh::getZeroPartIntegral(const Point &point, Types::complex_d k) {
        Types::Vector3c result{0, 0, 0};
        // суммирование по всем ячейкам
        for (const auto &cell: surfaceMesh_.getCells()) {
            const auto AB = [&](Types::scalar t) -> Types::Vector3c {
                const Types::Vector3d y = cell.parametrization(t, 0);
                return Helmholtz::V(k, point, y);
            };
            const auto BC = [&](Types::scalar t) -> Types::Vector3c {
                const Types::Vector3d y = cell.parametrization(1, t);
                return Helmholtz::V(k, point, y);
            };
            const auto CD = [&](Types::scalar t) -> Types::Vector3c {
                const Types::Vector3d y = cell.parametrization(1 - t, 1);
                return Helmholtz::V(k, point, y);
            };
            const auto DA = [&](Types::scalar t) -> Types::Vector3c {
                const Types::Vector3d y = cell.parametrization(0, 1 - t);
                return Helmholtz::V(k, point, y);
            };
            const Types::Vector3c current_step =
                    (cell.collPoint_.J_.dot(cell.integrationParameters.mul[0])) *
                    DefiniteIntegrals::integrate<DefiniteIntegrals::GaussLegendre::Quadrature<4>>(AB, {0}, {1}) +
                    (cell.collPoint_.J_.dot(cell.integrationParameters.mul[1])) *
                    DefiniteIntegrals::integrate<DefiniteIntegrals::GaussLegendre::Quadrature<4>>(BC, {0}, {1}) +
                    (cell.collPoint_.J_.dot(cell.integrationParameters.mul[2])) *
                    DefiniteIntegrals::integrate<DefiniteIntegrals::GaussLegendre::Quadrature<4>>(CD, {0}, {1}) +
                    (cell.collPoint_.J_.dot(cell.integrationParameters.mul[3])) *
                    DefiniteIntegrals::integrate<DefiniteIntegrals::GaussLegendre::Quadrature<4>>(DA, {0}, {1});

            result += current_step;
        }
        return result;
    }

    Types::Vector3c VolumeMesh::getRIntegral(const Point &point, Types::complex_d k) {
        Types::Vector3c result{0, 0, 0};
        // суммирование по всем ячейкам
        for (const auto &cell: surfaceMesh_.getCells()) {
            // расчет на одной ячейке сетки
            const auto phi = [&](Types::scalar p, Types::scalar q) -> Types::Vector3c {
                const Types::Vector3d y = cell.parametrization(p, q);
                const Types::scalar mul = cell.multiplier(p, q);
                return (-1) * Helmholtz::V(k, point, y) * mul;
            };
            result += cell.collPoint_.J_.cross(
                    DefiniteIntegrals::integrate<DefiniteIntegrals::GaussLegendre::Quadrature<8, 8 >>(phi, {0, 0}, {1., 1.}));
        }
        return result;
    };
}

