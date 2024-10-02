//
// Created by evgen on 08.02.24.
//

#include "VolumeMesh.hpp"
#include "operators/Functions.hpp"
#include "operators/Operators.hpp"
#include "math/integration/Quadrature.hpp"
#include "math/integration/gauss_quadrature/GaussLegenderPoints.hpp"
#include "math/MathConstants.hpp"
#include "math/Productions.hpp"

namespace EMW::Mesh {
    void EMW::Mesh::VolumeMesh::evaluateOperator(Types::scalar k, const Math::SurfaceField &field) {
        for (auto &node: nodes_) {
            node.E_ = Operators::K0<DefiniteIntegrals::GaussLegendre::Quadrature<8>>(node.point_, k, field) +
                      Operators::K1<DefiniteIntegrals::GaussLegendre::Quadrature<8, 8>>(node.point_, k, field);
        }
    }

    void EMW::Mesh::VolumeMesh::addInitialField(const std::function<Types::Vector3c(Mesh::point_t)> &initial) {
        for (auto &node: nodes_) {
            node.E_ += initial(node.point_);
        }
    }

    void VolumeMesh::calculateFullField(Types::scalar k, const Math::SurfaceField &field,
                                        const std::function<Types::Vector3c(Mesh::point_t)> &initial) {
        evaluateOperator(k, field);
        addInitialField(initial);
    }
};
