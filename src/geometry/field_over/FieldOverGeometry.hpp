//
// Created by evgen on 06.02.2025.
//

#ifndef FIELDOVERGEOMETRY_HPP
#define FIELDOVERGEOMETRY_HPP
#include <types/Types.hpp>

namespace EMW::Geometry::FieldsColletion {
template<typename TopologicalStructure, typename field_t>
class FieldOverGeometry {
    Containers::vector<field_t> fields;

public:
    FieldOverGeometry(const TopologicalStructure& geometry, const std::function<field_t(Types::index)> &construct);
};

template <typename TopologicalStructure, typename field_t>
FieldOverGeometry<TopologicalStructure, field_t>::FieldOverGeometry(
    const TopologicalStructure &geometry, const std::function<field_t(Types::index)> &construct) {}


}
#endif //FIELDOVERGEOMETRY_HPP
