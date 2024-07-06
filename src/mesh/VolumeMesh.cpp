//
// Created by evgen on 08.02.24.
//

#include "VolumeMesh.hpp"
#include "operators/Functions.hpp"
#include "operators/Operators.hpp"
#include "integration/Quadrature.hpp"
#include "integration/gauss_quadrature/GaussLegenderPoints.hpp"
#include "math/MathConstants.hpp"
#include "math/Productions.hpp"

namespace EMW::Mesh {
    void EMW::Mesh::VolumeMesh::calculateAll(const Types::Vector3d &polarization, const Types::Vector3c &k_vec,
                                             Types::scalar k) {
        if (surfaceMesh_.jFilled()) {
            for (auto &node: nodes_) {
                node.E_ =
                        Operators::K0<DefiniteIntegrals::GaussLegendre::Quadrature<8>>(node.point_,
                                                                                       surfaceMesh_.getCells(),
                                                                                       k) +
                        Operators::K1<DefiniteIntegrals::GaussLegendre::Quadrature<8, 8>>(node.point_,
                                                                                          surfaceMesh_.getCells(), k) +
                        polarization * std::exp(Math::Constants::i * Math::quasiDot(node.point_, k_vec));
            }
        } else {
            throw std::exception{};
        }
    }

    Types::Vector3c EMW::Mesh::VolumeMesh::sigmaOverCell(Types::complex_d k, const Types::Vector3d &tau,
                                                         const Mesh::IndexedCell &cell) const {
        const auto phi = [&](Types::scalar p, Types::scalar q) -> Types::Vector3c {
            const Mesh::Point y = cell.parametrization(p, q);
            const Types::scalar mul = cell.multiplier(p, q);
            return Helmholtz::sigmaKernel(k, tau, y, cell.collPoint_.J_) * mul;
        };
        return DefiniteIntegrals::integrate<DefiniteIntegrals::GaussLegendre::Quadrature<8, 8>>(phi, {0, 0}, {1, 1});
    }

    Types::scalar EMW::Mesh::VolumeMesh::calculateESS(const Types::Vector3d &tau, Types::complex_d k) const {
        if (surfaceMesh_.jFilled()) {
            Types::Vector3c result = Types::Vector3c::Zero();
            for (auto &cell: surfaceMesh_.getCells()) {
                result += sigmaOverCell(k, tau, cell);
            }
            return Math::Constants::inverse_4PI<Types::scalar>() * result.squaredNorm();
        } else {
            throw std::exception();
        }
    }
};
