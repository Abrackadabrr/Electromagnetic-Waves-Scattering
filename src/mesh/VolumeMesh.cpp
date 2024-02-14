//
// Created by evgen on 08.02.24.
//

#include "VolumeMesh.hpp"
#include "operators/Functions.hpp"
#include "operators/Operators.hpp"
#include "integration/Quadrature.hpp"
#include "integration/gauss_quadrature/GaussLegenderPoints.hpp"
#include "math/MathConstants.hpp"

namespace EMW::Mesh {
    void EMW::Mesh::VolumeMesh::calculateAll(const Types::Vector3d &polarization, const Types::Vector3c &k_vec,
                                             Types::complex_d k) {
        for (auto &node: nodes_) {
            node.E_ =
                    Operators::K0<DefiniteIntegrals::GaussLegendre::Quadrature<8>>(node.point_, surfaceMesh_.getCells(),
                                                                                   k) +
                    Operators::K1<DefiniteIntegrals::GaussLegendre::Quadrature<8, 8>>(node.point_,
                                                                                      surfaceMesh_.getCells(),
                                                                                      k) +
                    polarization * std::exp(Math::Constants::i * node.point_.dot(k_vec) * k);
        }
    }

    void calculateESS(const Types::Vector3d &polarization, const Types::Vector3d &k_vec,
                      Types::complex_d k) {}
}

