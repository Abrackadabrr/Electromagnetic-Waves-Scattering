//
// Created by evgen on 12.02.2025.
//

#ifndef ASIASES_FOR_MATRIX_HPP
#define ASIASES_FOR_MATRIX_HPP

#include "types/Types.hpp"

#include "ToeplitzFullyTemplated.hpp"
#include "DynamicFactoredMatrix.hpp"

namespace EMW::Math::LinAgl::Matrix {

/** Алиасы для удобной работы с матрицами */

template <typename T> using ToeplitzBlock = ToeplitzStructure<T, Types::MatrixX<T>>;
template <typename T> using ToeplitzToeplitzBlock = ToeplitzStructure<T, ToeplitzBlock<T>>;



} // namespace EMW::Math::LinAgl::Matrix

#endif // ASIASES_FOR_MATRIX_HPP
