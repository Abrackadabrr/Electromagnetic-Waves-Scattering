//
// Created by evgen on 23.08.24.
//

#include "SingularIntegration.hpp"
#include <cmath>

namespace EMW::Math::AnalyticalIntegration {

    Types::scalar integrate_1_div_r(const Mesh::point_t &r, const Mesh::IndexedCell &cell) {
        // определяем основные геометрические характеритики
        const Types::Vector3d n = cell.normal;
        const auto vertex = cell.getVertexAsArray();
        const Types::scalar d = std::abs((r - cell.collPoint_.point_).dot(n));
        Types::scalar result = 0;
        for (int i = 0; i != vertex.size(); i++) {
            const auto &rp = vertex[(i + 1) % vertex.size()];
            const auto &rm = vertex[i];
            // расчитываем геометрию, зависящую от сегмента, если этот сегмент не вырожден
            if ((rp - rm).norm() >= 1e-12) {
                const Types::Vector3d l = (rp - rm).normalized();
                // если особая точка лежит не на прямой, содержащей часть границы
                if (std::abs(std::abs((r - rm).dot(l)) - (r - rm).norm()) >= 1e-12) {
                    const Types::Vector3d u = l.cross(n);
                    const Types::scalar l_plus = (rp - r).dot(l);
                    const Types::scalar l_minus = (rm - r).dot(l);
                    const Types::scalar p0 = (rp - r).dot(u);
                    const Types::Vector3d p0_vec = (rp - r - l_plus * l) / p0;
                    const Types::scalar R_plus = (rp - r).norm();
                    const Types::scalar R_minus = (rm - r).norm();
                    const Types::scalar R0_sq = p0 * p0 + d * d;
                    // расчет слагаемого от iй части границы
                    result +=
                            p0_vec.dot(u) * (
                                    p0 * std::log((R_plus + l_plus) / (R_minus + l_minus)) -
                                    d * (
                                            std::atan2(p0 * l_plus, R0_sq + d * R_plus) -
                                            std::atan2(p0 * l_minus, R0_sq + d * R_minus)
                                    )
                            );
                }
            }
        }
        return result;
    }
}