//
// Created by evgen on 15.02.24.
//

#include "ESA.hpp"
#include "types/Types.hpp"
#include "mesh/MeshTypes.hpp"
#include "mesh/Mesh.hpp"
#include "integration/Quadrature.hpp"
#include "integration/gauss_quadrature/GaussLegenderPoints.hpp"
#include "operators/Functions.hpp"
#include "math/MathConstants.hpp"

namespace EMW::ESA {
    Types::Vector3c sigmaOverCell(Types::complex_d k, const Types::Vector3d &tau,
                                                         const Mesh::IndexedCell &cell) {
        const auto phi = [&](Types::scalar p, Types::scalar q) -> Types::Vector3c {
            const Mesh::Point y = cell.parametrization(p, q);
            const Types::scalar mul = cell.multiplier(p, q);
            return Helmholtz::sigmaKernel(k, tau, y, cell.collPoint_.J_) * mul;
        };
        return DefiniteIntegrals::integrate < DefiniteIntegrals::GaussLegendre::Quadrature < 8, 8
                >> (phi, {0, 0}, {1, 1});
    }

    Types::scalar calculateESA(const Types::Vector3d &tau, Types::complex_d k,
                                                      const Mesh::SurfaceMesh &mesh) {
        if (mesh.jFilled()) {
            Types::Vector3c result = Types::Vector3c::Zero();
            for (auto &cell: mesh.getCells()) {
                result += sigmaOverCell(k, tau, cell);
            }
            return Math::Constants::inverse_4PI<Types::scalar>() * result.squaredNorm();
        } else {
            throw std::exception();
        }
    }
}