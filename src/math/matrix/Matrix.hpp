//
// Created by evgen on 12.02.2025.
//

#ifndef ASIASES_FOR_MATRIX_HPP
#define ASIASES_FOR_MATRIX_HPP

#include "types/Types.hpp"

#include "DynamicFactoredMatrix.hpp"
#include "ToeplitzFullyTemplated.hpp"

namespace EMW::Math::LinAgl::Matrix {

/** Алиасы для удобной работы с матрицами */

template <typename T> using ToeplitzBlock = ToeplitzStructure<T, Types::MatrixX<T>>;
template <typename T> using ToeplitzToeplitzBlock = ToeplitzStructure<T, ToeplitzBlock<T>>;
template <typename T> using TripleToeplitzBlock = ToeplitzStructure<T, ToeplitzToeplitzBlock<T>>;

template <typename T> using ToeplitzDynFactoredBlock = ToeplitzStructure<T, DynamicFactoredMatrix<Types::MatrixX<T>>>;
template <typename T> using ToeplitzToeplitzDynFactoredBlock = ToeplitzStructure<T, ToeplitzDynFactoredBlock<T>>;
template <typename T> using TripleToeplitzFactoredBlock = ToeplitzStructure<T, ToeplitzToeplitzDynFactoredBlock<T>>;

// Zero toeplitz structures with dense internal blocks

template<typename T>
ToeplitzBlock<T> ZeroToeplitzBlock(size_t fl, size_t inner_size) {
    const auto get_zero_inner_block = [&](size_t i, size_t j)->Types::MatrixX<T> {
        return Types::MatrixX<T>::Zero(inner_size, inner_size);
    };
    return ToeplitzBlock<T>(fl, fl, get_zero_inner_block);
}

template<typename T>
ToeplitzToeplitzBlock<T> ZeroDoubleToeplitzBlock(size_t fl, size_t sl, size_t inner_size) {
    const auto get_zero_inner_block = [&](size_t i, size_t j) -> ToeplitzBlock<T> {
        return ZeroToeplitzBlock<T>(fl, inner_size);
    };
    return ToeplitzToeplitzBlock<T>(sl, sl, get_zero_inner_block);
}

template<typename T>
TripleToeplitzBlock<T> ZeroTripleToeplitzBlock(size_t fl, size_t sl, size_t tl, size_t inner_size) {
    const auto get_zero_inner_block = [&](size_t i, size_t j) -> ToeplitzToeplitzBlock<T> {
        return ZeroDoubleToeplitzBlock<T>(fl, sl, inner_size);
    };
    return TripleToeplitzBlock<T>(tl, tl, get_zero_inner_block);
}

// Zero toeplitz structures with factored inside blocks

template<typename T>
decltype(auto) ZeroToeplitzFactoredBlock(size_t fl, size_t inner_size) {
    const auto get_zero_inner_block = [&](size_t i, size_t j) {
        return DynamicFactoredMatrix<Types::MatrixX<T>>{inner_size, inner_size};
    };
    return ToeplitzDynFactoredBlock<T>(fl, fl, get_zero_inner_block);
}

template<typename T>
decltype(auto) ZeroDoubleToeplitzFactoredBlock(size_t fl, size_t sl, size_t inner_size) {
    const auto get_zero_inner_block = [&](size_t i, size_t j) {
        return ZeroToeplitzFactoredBlock<T>(fl, inner_size);
    };
    return ToeplitzToeplitzDynFactoredBlock<T>(sl, sl, get_zero_inner_block);
}

template<typename T>
decltype(auto) ZeroTripleToeplitzFactoredBlock(size_t fl, size_t sl, size_t tl, size_t inner_size) {
    const auto get_zero_inner_block = [&](size_t i, size_t j) {
        return ZeroDoubleToeplitzFactoredBlock<T>(fl, sl, inner_size);
    };
    return TripleToeplitzFactoredBlock<T>(tl, tl, get_zero_inner_block);
}

} // namespace EMW::Math::LinAgl::Matrix

#endif // ASIASES_FOR_MATRIX_HPP
