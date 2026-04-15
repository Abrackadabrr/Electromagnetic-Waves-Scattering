//
// Created by evgen on 20.03.2025.
//

#ifndef MATRIXTRAITS_HPP
#define MATRIXTRAITS_HPP

#include "types/Types.hpp"

namespace EMW::Math::LinAgl::Matrix {

template<typename matrix_t>
struct MatrixTraits {};

template<typename T>
struct MatrixTraits<Types::MatrixX<T>> {
    using matrix_t = Types::MatrixX<T>;
    using vector_t = Types::VectorX<T>;
    using element_t = T;
    using production_t = matrix_t;
};
}
#endif //MATRIXTRAITS_HPP
