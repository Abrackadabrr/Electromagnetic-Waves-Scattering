//
// Created by evgen on 18.03.2025.
//

#ifndef FACTOREDMATRIX_HPP
#define FACTOREDMATRIX_HPP

#include "types/TypeTraits.hpp"
#include "types/Types.hpp"

#include <Utils.hpp>
#include <iostream>

namespace EMW::Math::LinAgl::Matrix {

template <typename... factors_t> class FactoredMatrix {
    using factors_pack = Types::TypeTraits::TypePack<factors_t...>;

    Containers::tuple<factors_t...> factors_;

public:
    template <typename vector_t, Types::index... Idx>
    vector_t matvec_impl(const vector_t &vector, const std::index_sequence<Idx...> & /*noname*/) const {
        return (std::get<Idx>(factors_) * ... * vector);
    }

    template <Types::index... Idx> decltype(auto) compute_impl(const std::index_sequence<Idx...> & /*noname*/) const {
        return (std::get<Idx>(factors_) * ...);
    }

  public:
    static constexpr std::size_t number_factors = sizeof...(factors_t);

    /** Warning: it could be universal references with ref collapsing */
    FactoredMatrix(factors_t &&...factors) : factors_(std::forward<factors_t>(factors)...){};

    FactoredMatrix(const factors_t &...factors) : factors_(factors...) {
        std::cout << "FactoredMatrix constructor with copying" << std::endl;
    };

    template<typename vector_t>
    vector_t matvec(const vector_t &vector) const {
        return matvec_impl(vector, std::make_index_sequence<number_factors>{});
    }

    decltype(auto) compute() const {
        return compute_impl(std::make_index_sequence<number_factors>{});
    }

    template <Types::index I> decltype(auto) get() { return std::get<I>(factors_); }
};

template <typename ... factors_t, typename vector_t>
vector_t operator*(const FactoredMatrix<factors_t...> &mat, const vector_t &vector) {
    return mat.matvec(vector);
}

template <typename factor_t>
class DynamicFactoredMatrix {

    Containers::vector<factor_t> factors_;

public:
    DynamicFactoredMatrix(Containers::vector<factor_t> &&factors) : factors_(factors){};

    DynamicFactoredMatrix(const Containers::vector<factor_t> &factors) : factors_(factors) {
        std::cout << "FactoredMatrix constructor with copying" << std::endl;
    };

    template<typename vector_t>
    vector_t matvec(const vector_t &vector) const {
        vector_t result = vector_t::Zero(vector.size());
        for (auto& factor: Utils::reverse(factors_)) {
            result += factor * result;
        }
    }

    decltype(auto) compute() const {
        factor_t result = factors_[0];
        for (Types::index i = 1; i < factor_number(); ++i) {
            result *= factors_[i];
        }
        return result;
    }

    template <Types::index I> decltype(auto) get() { return factors_[I]; }

    // --- Selectors --- //

    Types::index factor_number() const { return factors_.size(); };
    Types::index memory_usage() const {
        Types::index result = 0;
        for (auto& factor: factors_) {
            result += factor.rows() * factor.cols() * 16;
        }
        return result;
    }
};

} // namespace EMW::Math::Matrix

#endif // FACTOREDMATRIX_HPP
