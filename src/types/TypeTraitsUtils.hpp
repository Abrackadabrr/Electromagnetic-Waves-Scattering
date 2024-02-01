//
// Created by evgen on 01.02.24.
//

#ifndef ELECTROMAGNETIC_WAVES_SCATTERING_TYPETRAITSUTILS_HPP
#define ELECTROMAGNETIC_WAVES_SCATTERING_TYPETRAITSUTILS_HPP

#include <utility>

namespace EMW::Types::TypeTraits::Utils {

    namespace Aux {

        template <typename>
        struct get {};

        template <std::size_t... Is>
        struct get<std::index_sequence<Is...>> {
        template <typename T>
        constexpr static T foo(decltype(Is, reinterpret_cast<void*>(0))..., T*, ...);
    };

}  // namespace Aux

/**
 * Проверяет, содержит ли пачка параметров Ts тип T
 * @tparam T искомый тип
 * @tparam Ts пачка параметров
 *
 * @return true, если Ts содержит T, иначе false
 */
template <typename T, typename... Ts>
constexpr bool contains() noexcept {
    return (std::is_same_v<T, Ts> || ...);
}

/**
 * Возвращает номер первого вхождения типа T в пачку параметров Ts
 * @tparam T искомый тип T
 * @tparam Ts пачка параметров
 *
 * @return индекс типа T в пачке параметров Ts или размер пачки, если она не содержит этот тип
 */
template <typename T, typename... Ts>
constexpr std::size_t find() noexcept {
    constexpr bool is_same[sizeof...(Ts)] = {std::is_same_v<T, Ts>...};
    for (std::size_t i = 0; i < sizeof...(Ts); ++i) {
        if (is_same[i]) {
            return i;
        }
    }
    return sizeof...(Ts);
}

/**
 * Возвращает номер первого типа, удовлетворяющего предикату F в пачке Ts
 * @tparam F предикат
 * @tparam Ts пачка параметров
 *
 * @return индекс типа, удовлетворяющего F, или размер пачки
 */
template <template <typename> typename F, typename... Ts>
constexpr std::size_t find_if() noexcept {
    constexpr bool condition[sizeof...(Ts)] = {F<Ts>::value...};
    for (std::size_t i = 0; i < sizeof...(Ts); ++i) {
        if (condition[i]) {
            return i;
        }
    }
    return sizeof...(Ts);
}

/**
 * Возвращает число вхождений типа T в пачку параметров Ts
 * @tparam T искомый тип T
 * @tparam Ts пачка параметров
 *
 * @return число вхождений
 */
template <typename T, typename... Ts>
constexpr std::size_t count() noexcept {
    constexpr bool is_same[sizeof...(Ts)] = {std::is_same_v<T, Ts>...};
    std::size_t count = 0;
    for (std::size_t i = 0; i < sizeof...(Ts); ++i) {
        count += static_cast<std::size_t>(is_same[i]);
    }
    return count;
}

/**
 * Проверяет, состоит ли пачка параметров из уникальных типов
 * @tparam Ts пачка параметров
 *
 * @return true, если каждый тип входит в пачку только один раз, иначе false
 */
template <typename... Ts>
constexpr bool unique() noexcept {
    return ((count<Ts, Ts...>() == 1) && ...);
}

template <std::size_t I, typename... Ts>
using GetType = decltype(Aux::get<std::make_index_sequence<I>>::foo(reinterpret_cast<Ts*>(0)...));

}

#endif //ELECTROMAGNETIC_WAVES_SCATTERING_TYPETRAITSUTILS_HPP
