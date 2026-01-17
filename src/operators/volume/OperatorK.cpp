//
// Created by evgen on 17.01.2026.
//

#include "OperatorK.hpp"

namespace EMW::Operators::Volume {

Types::Matrix3c operator_K_over_cube_mesh::matrix_3_coef(Types::index k, Types::index p) {

    const auto &cube_k = mesh.getCells()[k];
    const auto &cube_p = mesh.getCells()[p];
    const auto h = mesh.h();

    if ((cube_k.center_ - cube_p.center_).norm() < 2 * h) {
        // интегрирование с вынесением особенности
    } else {
        // обычное интегрирование по объёму дважды
    }

}


}