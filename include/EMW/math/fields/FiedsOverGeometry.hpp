//
// Created by evgen on 05.02.2025.
//

#ifndef FIEDSOVERGEOMETRY_HPP
#define FIEDSOVERGEOMETRY_HPP
#include <tuple>
#include "types/Types.hpp"
namespace EMW::Math {

/**
 *
 * @tparam TopologicalStructure - сборник многообразий, фактически та геометрия, на которой мы считаем
 * @tparam Fields - типы полей, которые мы тут задали
 */
template <typename TopologicalStructure, typename ... Fields>
class FieldsOverGeometry {
    std::tuple<Containers::vector<Fields>...> fields;

public:
    explicit FieldsOverGeometry(const TopologicalStructure& geometry){};

    virtual ~FieldsOverGeometry()= default;
};

}

#endif //FIEDSOVERGEOMETRY_HPP
