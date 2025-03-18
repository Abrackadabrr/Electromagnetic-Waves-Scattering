//
// Created by evgen on 18.03.2025.
//

#ifndef FACTOREDMATRIX_HPP
#define FACTOREDMATRIX_HPP

#include "types/TypeTraits.hpp"
#include "types/Types.hpp"

namespace EMW::Math::Matrix {

template <typename... factors_t> class FactoredMatrix {
    using factors_pack = Types::TypeTraits::TypePack<factors_t...>;

    Containers::tuple<factors_t...> factors_;

    template <typename vector_t, Types::index... Idx>
    vector_t matvec_impl(const vector_t &vector, const std::index_sequence<Idx...>& /*noname*/) const {
        return (std::get<Idx>(factors_) * ... * vector);
    }

  public:
    /** Warning: it could be universal references with ref collapsing */
    FactoredMatrix(factors_t &&...factors) : factors_(std::forward<factors_t>(factors)...){};

    FactoredMatrix(factors_t &...factors) : factors_(factors...) {
        std::cout << "FactoredMatrix constructor with copying" << std::endl;
    };

    template <typename vector_t> vector_t matvec(const vector_t &vector) const {
        return matvec_impl(vector, std::make_index_sequence<sizeof...(factors_t)>());
    }

    template <Types::index I> decltype(auto) get() { return std::get<I>(factors_); }
};

template<typename ... factors_t, typename vector_t>
vector_t operator*(const FactoredMatrix<factors_t...>& mat, const vector_t &vector) {
    return mat.matvec(vector);
}

} // namespace EMW::Math::Matrix



#endif // FACTOREDMATRIX_HPP
