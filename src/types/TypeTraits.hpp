//
// Created by evgen on 31.01.24.
//

#ifndef ELECTROMAGNETIC_WAVES_SCATTERING_TYPETRAITS_HPP
#define ELECTROMAGNETIC_WAVES_SCATTERING_TYPETRAITS_HPP

#include "Types.hpp"
#include "TypeTraitsUtils.hpp"

namespace EMW::Types::TypeTraits {
/**
 * Структура для хранения пачки параметров, наследующаяся от них
 * @tparam Ts пачка параметров
 */
    template<typename... Ts>
    struct Pack : Ts ... {
        /** Размер пачки параметров **/
        constexpr static index size = sizeof...(Ts);

        template<template<typename...> typename F>
        using Transformed = Pack<F<Ts>...>;

        template<typename... Us>
        using AddToPack = Pack<Ts..., Us...>;
    };

/**
 * Структура для хранения пачки параметров
 * @tparam Ts пачка параметров
 */
    template<typename... Ts>
    struct TypePack {
    public:
        /** Размер пачки параметров **/
        constexpr static index size = sizeof...(Ts);

        /**
         * Возвращает индекс типа в пачке параметров
         * @tparam T искомый тип
         * @return индекс первого вхождения в пачку или её размер, если он не содержится
         */
        template<typename T>
        [[nodiscard]] constexpr static index getIndex() {
            return Utils::find<T, Ts...>();
        }

        /**
         * Применяет F к пачке параметров
         * @tparam F применяемый тип
         * @return TypePack от F<Ts>...
         */
        template<template<typename...> typename F>
        using Transformed = TypePack<F<Ts>...>;

        /**
         * Подставляет пачку параметров в F
         * @tparam F тип, в который производится подстановка
         * @return тип с подстановкой пачки параметров
         */
        template<template<typename...> typename F>
        using Substitute = F<Ts...>;

        /**
         * Добавляет типы в пачку параметров
         * @tparam Us добавляемая пачка параметров
         * @return тип с объединёнными пачками
         */
        template<typename... Us>
        using AddToPack = TypePack<Ts..., Us...>;

        template<template<typename> typename F>
        using TransformToPack = Pack<F<Ts>...>;
    };

/**
 * Объединяет две пачки параметров
 * @tparam первая пачка параметров
 * @tparam Us вторая пачка параметров
 * @return TypePack с объединёнными пачками
 */
    template<typename... Ts, typename... Us>
    constexpr TypePack<Ts..., Us...> operator+(TypePack<Ts...>, TypePack<Us...>) noexcept {
        return {};
    }

/**
 * Добавляет параметр в пачку, если его там нет (справа)
 * @tparam пачка параметров
 * @tparam U новый параметр
 * @return TypePack с новым типом, если пачка его ещё не содержит
 */
    template<typename... Ts, typename U>
    constexpr auto operator||(TypePack<Ts...>, TypePack<U>) noexcept {
        if constexpr (Utils::contains<U, Ts...>()) {
            return TypePack<Ts...>{};
        } else {
            return TypePack<Ts..., U>{};
        }
    }

/**
 * Добавляет параметр в пачку, если его там нет (слева)
 * @tparam T новый параметр
 * @tparam Ts пачка параметров
 * @return TypePack с новым типом, если он его ещё не содержит
 */
    template<typename T, typename... Ts>
    constexpr auto operator|(TypePack<T>, TypePack<Ts...>) noexcept {
        if constexpr (Utils::contains<T, Ts...>()) {
            return TypePack<Ts...>{};
        } else {
            return TypePack<T, Ts...>{};
        }
    }

    namespace Aux {

        template<typename T, typename... Us>
        constexpr auto cartesianProduct(TypePack<T>, TypePack<Us...>) {
            if constexpr (sizeof...(Us) == 0) {
                return TypePack<TypePack<T>>{};
            } else {
                return TypePack<decltype(TypePack<T>{} + Us{})...>{};
            }
        }
    }  // namespace Aux

/**
 * Выполняет декартово произведение двух пачек параметров
 * @tparam Ts первая пачка параметров
 * @tparam Us вторая пачка параметров (состоит из TypePack'ов)
 * @return TypePack из TypePack'ов с декартовыми произведениями
 */
    template<typename... Ts, typename... Us>
    constexpr auto operator*(TypePack<Ts...>, TypePack<Us...>) {
        return (Aux::cartesianProduct(TypePack<Ts>{}, TypePack<Us...>{}) + ...);
    }

    namespace Aux {

        template<typename... Ts>
        constexpr auto makeUniqueTypePack(TypePack<Ts...>) noexcept {
            return (TypePack<>{} || ... || TypePack<Ts>{});
        }

        template<typename... Ts>
        constexpr auto makeUniqueTypePackReverse(TypePack<Ts...>) noexcept {
            return (TypePack<Ts>{} | ... | TypePack<>{});
        }

        template<typename... Ts, typename... Us>
        constexpr auto notInTypePack(TypePack<Ts...>, TypePack<Us...>) noexcept {
            return (TypePack<>{} + ... + std::conditional_t<Utils::contains<Us, Ts...>(), TypePack<>, TypePack<Us>>{});
        }

        template<typename... Packs>
        constexpr auto mergePacks() noexcept {
            using T = decltype((TypePack<>{} + ... + std::declval<Packs>()));
            return T{};
        }

        template<typename... Packs>
        constexpr auto cartesianProduct() noexcept {
            return (Packs{} * ... * TypePack<>{});
        }

    }  // namespace Aux

/**
 * Объединяет несколько пачек параметров
 * @tparam Packs объединяемые пачки
 * @return TypePack после объединения
 */
    template<typename... Packs>
    using MergeTypePacks = decltype(Aux::mergePacks<Packs...>());

/**
 * Выделяет уникальные типы из TypePack в порядке первого вхождения
 * @tparam Pack пачка параметров, из которой нужно выделить уникальные
 * @return TypePack с уникальными типами
 */
    template<typename Pack>
    using UniqueTypePack = decltype(Aux::makeUniqueTypePack(std::declval<Pack>()));

/**
 * Выделяет уникальные типы из TypePack в порядке последнего вхождения
 * @tparam Pack пачка параметров, из которой нужно выделить уникальные
 * @return TypePack с уникальными типами
 */
    template<typename Pack>
    using UniqueTypePackInverseOrder = decltype(Aux::makeUniqueTypePackReverse(std::declval<Pack>()));

/**
 * Исключает из второго набора типов типы, содержащиеся в первом
 * @tparam Pack1 набор типов для исключения
 * @tparam Pack2 набор типов для фильтрации
 * @return отфильтрованный набор типов
 */
    template<typename Pack1, typename Pack2>
    using NotInTypePack = decltype(Aux::notInTypePack(std::declval<Pack1>(), std::declval<Pack2>()));

/**
 * @tparam Packs TypePack'и для декартова произведения
 * @return TypePack из TypePack'ов с декартовым произведением
 */
    template<typename... Packs>
    using CartesianProduct = decltype(Aux::cartesianProduct<Packs...>());

}

#endif //ELECTROMAGNETIC_WAVES_SCATTERING_TYPETRAITS_HPP
