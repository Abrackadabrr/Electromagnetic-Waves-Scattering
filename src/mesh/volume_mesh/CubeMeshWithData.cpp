//
// Created by evgen on 10.02.2026.
//

#include "CubeMeshWithData.hpp"

#include <ranges>

namespace EMW::Mesh::VolumeMesh {

Types::VectorXc CubeMeshWithData::getScalarDataAsVector(const std::string &name) const {
    const auto& raw_data = getScalarData(name);
    Types::VectorXc result(cellsCount);
    for (auto&& [res_v, v] : std::views::zip(result, raw_data.second))
        res_v = v;
    return result;
}

}
