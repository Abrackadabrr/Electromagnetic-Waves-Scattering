//
// Created by evgen on 31.01.24.
//

#ifndef ELECTROMAGNETIC_WAVES_SCATTERING_MATRIXGENERATIONFUNCTION_HPP
#define ELECTROMAGNETIC_WAVES_SCATTERING_MATRIXGENERATIONFUNCTION_HPP

#include "types/Types.hpp"
#include "mesh/MeshTypes.hpp"
#include "slae_generation/Functions.hpp"
#include "integration/gauss_quadrature/GaussQuadrature.hpp"
#include "integration/gauss_quadrature/GaussLegenderPoints.hpp"

namespace EMW::Matrix {
    Types::complex_d
    getFirstPart(Types::index i, Types::index j, Types::complex_d k, const Containers::vector<Mesh::IndexedCell> &points) {
        const auto phi = [&](Types::scalar p, Types::scalar q) -> Types::complex_d {
            const Types::Vector3d y = points[j].parametrization(p, q);
            const Types::scalar mul = points[j].multiplier(p, q);
            return Helmholtz::F(k, points[i].collPoint_.point_, y) * mul;
        };
        const Types::scalar tauitauj = points[i].tau1.dot(points[j].tau1);
        return DefiniteIntegrals::integrate<DefiniteIntegrals::Quadrature<6, 6>>(phi, {-1, -1}, {2, 2});
    }
}

#endif //ELECTROMAGNETIC_WAVES_SCATTERING_MATRIXGENERATIONFUNCTION_HPP
