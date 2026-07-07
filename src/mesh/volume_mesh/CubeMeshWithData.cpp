//
// Created by evgen on 10.02.2026.
//

#include "mesh/volume_mesh/CubeMeshWithData.hpp"

#include <ranges>

namespace EMW::Mesh::VolumeMesh {

Types::VectorXc CubeMeshWithData::getScalarDataAsVector(const std::string &name) const {
    const auto& raw_data = getScalarData(name);
    Types::VectorXc result(cellsCount);
    for (auto&& [res_v, v] : std::views::zip(result, raw_data))
        res_v = v;
    return result;
}

Types::VectorXc CubeMeshWithData::getVectorDataAsVector(const std::string &name) const {
    const auto& raw_data = getVectorData(name);
    Types::VectorXc result(3 * cellsCount);
    for (Types::index i = 0; i < cellsCount; ++i)
        result.block(3 * i, 0, 3, 1) = raw_data[i];
    return result;
}

}
