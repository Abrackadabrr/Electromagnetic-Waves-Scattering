//
// Created by evgen on 15.02.24.
//

#include "ESA.hpp"
#include "math/MathConstants.hpp"
#include "math/fields/SurfaceVectorField.hpp"
#include "math/integration/Quadrature.hpp"
#include "math/integration/gauss_quadrature/GaussLegenderPoints.hpp"
#include "mesh/MeshTypes.hpp"
#include "mesh/SurfaceMesh.hpp"
#include "operators/Functions.hpp"
#include "types/Types.hpp"

namespace EMW::ESA {
    Types::Vector3c sigmaOverCell(Types::complex_d k, const Types::Vector3d &tau,
                                  const Mesh::IndexedCell &cell, const Types::Vector3c& j) {
        const auto phi = [&](Types::scalar p, Types::scalar q) -> Types::Vector3c {
            const Mesh::point_t y = cell.parametrization(p, q);
            const Types::scalar mul = cell.multiplier(p, q);
            return Helmholtz::sigmaKernel(k, tau, y, j) * mul;
        };
        return DefiniteIntegrals::integrate<DefiniteIntegrals::GaussLegendre::Quadrature<8, 8
        >>(phi, {0, 0}, {1, 1});
    }

    Types::scalar calculateESA(const Types::Vector3d &tau, Types::complex_d k,
                               const Math::SurfaceVectorField &j_field) {
        const auto & cells = j_field.getManifold().getCells();
        const auto & field = j_field.getField();
        Types::Vector3c result = Types::Vector3c::Zero();
        for (int i = 0; i != j_field.getField().size(); i++) {
            result += sigmaOverCell(k, tau, cells[i], field[i]);
        }
        return Math::Constants::inverse_4PI<Types::scalar>() * result.squaredNorm();
    }
}