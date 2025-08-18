//
// Created by evgen on 27.02.2025.
//

#ifndef ITERATIVE_SOLVERS_COVERAGE_MATRIXTRAITS_HPP
#define ITERATIVE_SOLVERS_COVERAGE_MATRIXTRAITS_HPP

#include "IdentityPreconditioner.hpp"
#include "MatrixReplacement.hpp"
#include "types/Types.hpp"
#include <Eigen/Core>

namespace Eigen::internal {
namespace Mat = EMW::Math::LinAgl::Matrix;
namespace My = Mat::Wrappers;
/**
 * Характеритики типа MatrixReplacement<MatrixType>
 * Это необходимо определенить для того, чтобы Eigen верно разрешал типы
 * @tparam MatrixType
 */
template <typename MatrixType, typename PrecondType> struct traits<My::MatrixReplacement<MatrixType, PrecondType>> {
    // Эти три вещи ниже должны зависеть от типа, который мы передаем в качестве MatrixType.
    // Но пока что это просто захардкожено
    typedef EMW::Types::complex_d Scalar;
    typedef EMW::Types::integer StorageIndex;
    static constexpr EMW::Types::integer Options_ = 0;
    // --- //

    typedef Sparse StorageKind;
    typedef MatrixXpr XprKind;
    enum {
        RowsAtCompileTime = Dynamic,
        ColsAtCompileTime = Dynamic,
        MaxRowsAtCompileTime = Dynamic,
        MaxColsAtCompileTime = Dynamic,
        Flags = Options_ | NestByRefBit | LvalueBit | CompressedAccessBit,
    };
};

template <typename MatrixType, typename PrecondType> struct traits<My::MatrixReplacementReal<MatrixType, PrecondType>> {
    // Эти три вещи ниже должны зависеть от типа, который мы передаем в качестве MatrixType.
    // Но пока что это просто захардкожено
    typedef EMW::Types::scalar Scalar;
    typedef EMW::Types::integer StorageIndex;
    static constexpr EMW::Types::integer Options_ = 0;
    // --- //

    typedef Sparse StorageKind;
    typedef MatrixXpr XprKind;
    enum {
        RowsAtCompileTime = Dynamic,
        ColsAtCompileTime = Dynamic,
        MaxRowsAtCompileTime = Dynamic,
        MaxColsAtCompileTime = Dynamic,
        Flags = Options_ | NestByRefBit | LvalueBit | CompressedAccessBit,
    };
};

template <typename MatrixType> struct traits<My::MatrixReplacement<MatrixType>> {
    // Эти три вещи ниже должны зависеть от типа, который мы передаем в качестве MatrixType.
    // Но пока что это просто захардкожено
    typedef EMW::Types::complex_d Scalar;
    typedef EMW::Types::integer StorageIndex;
    static constexpr EMW::Types::integer Options_ = 0;
    // --- //

    typedef Sparse StorageKind;
    typedef MatrixXpr XprKind;
    enum {
        RowsAtCompileTime = Dynamic,
        ColsAtCompileTime = Dynamic,
        MaxRowsAtCompileTime = Dynamic,
        MaxColsAtCompileTime = Dynamic,
        Flags = Options_ | NestByRefBit | LvalueBit | CompressedAccessBit,
    };
};
} // namespace Eigen::internal

#endif // ITERATIVE_SOLVERS_COVERAGE_MATRIXTRAITS_HPP
