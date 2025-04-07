//
// Created by evgen on 20.03.2025.
//

#ifndef DYMANICFACTOREDMATRIX_HPP
#define DYMANICFACTOREDMATRIX_HPP

#include "MatrixTraits.hpp"
#include "types/Types.hpp"

#include <iostream>

namespace EMW::Math::LinAgl::Matrix {

template <typename factor_t> class DynamicFactoredMatrix {

    using factor_traits = MatrixTraits<factor_t>;

    Containers::vector<factor_t> factors_;

  public:
    DynamicFactoredMatrix() = default;

    DynamicFactoredMatrix(Containers::vector<factor_t> &&factors) : factors_(factors){};

    DynamicFactoredMatrix(const Containers::vector<factor_t> &factors) : factors_(factors) {
        std::cout << "FactoredMatrix constructor with copying" << std::endl;
    };

    template <typename vector_t> vector_t matvec(const vector_t &vector) const;

    typename factor_traits::production_t compute() const noexcept;

    template <Types::index I> const factor_t& get() const;

    typename factor_traits::vector_t diagonal() const;

    // ---- Selectors ---- //
    Types::index factor_number() const { return factors_.size(); };
    Types::scalar memory_usage() const;
    Types::index rows() const { return factors_.front().rows(); };
    Types::index cols() const { return factors_.back().cols(); };

    // --- Aux methods --- //
    decltype(auto) to_dense() const noexcept { return compute(); };

    // --- Доступ к элементам матрицы --- //
};

template <typename factor_t>
template <typename vector_t>
vector_t DynamicFactoredMatrix<factor_t>::matvec(const vector_t &vector) const {
    vector_t result = factors_.back() * vector;
    for (Types::integer i = factor_number() - 2; i >= 0; --i) {
        result = factors_[i] * result;
    }
    return result;
}

template <typename factor_t>
typename DynamicFactoredMatrix<factor_t>::factor_traits::vector_t
DynamicFactoredMatrix<factor_t>::diagonal() const {
    if (factor_number() == 1)
        return factors_.front().diagonal();
    std::cout << "diagonal() is potentially long operation while running with more that 1 factor numbers" << std::endl;
    return compute().diagonal();
}

template <typename factor_t>
typename DynamicFactoredMatrix<factor_t>::factor_traits::production_t
DynamicFactoredMatrix<factor_t>::compute() const noexcept {
    factor_t result = factors_[0];
    for (Types::index i = 1; i < factor_number(); ++i) {
        result *= factors_[i];
    }
    return result;
}

template <typename factor_t> template <Types::index I>
const factor_t& DynamicFactoredMatrix<factor_t>::get() const {
    return factors_[I];
}

// --- Selectors --- //
template <typename factor_t> Types::scalar DynamicFactoredMatrix<factor_t>::memory_usage() const {
    Types::index result = 0;
    for (auto &factor : factors_) {
        result += factor.size();
    }
    return result;
}

template <typename factor_t, typename vector_t>
vector_t operator*(const DynamicFactoredMatrix<factor_t> &mat, const vector_t &vector) {
    return mat.matvec(vector);
}
} // namespace EMW::Math::LinAgl::Matrix

#endif // DYMANICFACTOREDMATRIX_HPP
