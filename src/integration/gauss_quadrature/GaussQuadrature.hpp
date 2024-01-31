//
// Created by evgen on 31.01.24.
//

#ifndef ELECTROMAGNETIC_WAVES_SCATTERING_GAUSSQUADRATURE_HPP
#define ELECTROMAGNETIC_WAVES_SCATTERING_GAUSSQUADRATURE_HPP

#include "types/Types.hpp"
#include "types/TypeTraits.hpp"

#include <utility>
#include <tuple>


namespace EMW::DefiniteIntegrals::aux {

    template<typename Arg>
    using DeltaType = decltype(std::declval<Arg>() - std::declval<Arg>());

    // Произвольный объект с оператором ()
    template<typename Callable>
    struct IntegralTypesUnwrap;

    // Специализация под функции
    template<typename FReturnType, typename... Args>
    struct IntegralTypesUnwrap<FReturnType (Args...)> {
        using ArgsTuple = std::tuple<Args...>;
        using DeltasTuple = std::tuple<DeltaType<Args>...>;
        using MeasureType = decltype((std::declval<DeltaType<Args>>() * ...));
        using ResultType = decltype(std::declval<MeasureType>() * std::declval<FReturnType>());
    };

    // Специализация под методы классов
    template <typename CTy, typename RTy, typename... ATy>
    struct IntegralTypesUnwrap<RTy (CTy::*)(ATy...)> {
        using ArgsTuple = std::tuple<ATy...>;
        using DeltasTuple = std::tuple<DeltaType<ATy>...>;
        using MeasureType = decltype((std::declval<DeltaType<ATy>>() * ...));
        using ResultType = decltype(std::declval<MeasureType>() * std::declval<RTy>());
    };

    template<typename Callable>
    struct IntegralTypes;

    template<typename FReturnType, typename... Args>
    struct IntegralTypes<FReturnType (Args...)>: IntegralTypesUnwrap<FReturnType (Args...)> {};

    template <typename CTy, typename RTy, typename... ATy>
    struct IntegralTypes<RTy (CTy::*)(ATy...) const> : IntegralTypesUnwrap<RTy(ATy...)> {};

    template<typename Callable>
    struct IntegralTypes: IntegralTypesUnwrap<decltype(&Callable::operator())> {};

    template<typename ArgsTuple, typename DeltaTuple, typename Points, Types::index... indexes>
    ArgsTuple createPoint(const ArgsTuple &start, const DeltaTuple &deltas, const Points &points,
                          const std::index_sequence<indexes...> & /*sequence*/) {
        return {std::get<indexes>(start) + (std::get<indexes>(deltas) / 2) * (std::get<indexes>(points) + 1) ...};
    }

    template<typename DeltaTuple, Types::index... indexes>
    decltype(auto) calcArea(const DeltaTuple &deltas, const std::index_sequence<indexes...> & /*sequence*/) {
        return (std::get<indexes>(deltas) * ...);
    }

    template<typename Quadrature, typename Callable, Types::index... indexes>
    IntegralTypes<Callable>::ResultType calcQuadratureSum(
            const Callable &f, const typename IntegralTypes<Callable>::ArgsTuple &startArgs,
            const typename IntegralTypes<Callable>::DeltasTuple &deltas,
            const std::index_sequence<indexes...> & /*sumSequence*/) {
        constexpr auto dimSize = std::tuple_size_v<typename IntegralTypes<Callable>::ArgsTuple>;
        static_assert(Quadrature::dim == dimSize);
        constexpr auto pointSequence = std::make_index_sequence<dimSize>();

        return ((Quadrature::nodes[indexes].weight *
                 std::apply(f, createPoint(startArgs, deltas, Quadrature::nodes[indexes].point, pointSequence))) +
                       ...) *
               (calcArea(deltas, pointSequence) / std::pow(2, dimSize));
    }

} // namespace EMW::DefiniteIntegrals::detail

namespace EMW::DefiniteIntegrals {

/** Выполняет интегрирование в многомерной прямоугольной области
 *
 * @tparam Quadrature       структура квадратуры. Содержит поля :
 *                              dim - размерность квадратуры,
 *                              size - количесто точек в квадратуре
 *                              points - точки квадратуры
 *                              weights - веса квадратуры
 *
 * @tparam Callable         Тип функция для интегрирования. Количество аргументов должно совпадать с dim в квадратуре
 *
 * @param f                 Объект функции для интегрирования
 * @param startArgs         начальные точки для интегрирования (координаты нижней левой точки)
 * @param deltas            длины ребер многомерного прямоугольника, в которой происходит интегрирование
 *
 * @return                  Значение интеграла
 */
    template<typename Quadrature, typename Callable>
    typename aux::IntegralTypes<Callable>::ResultType integrate(
            const Callable &f, const typename aux::IntegralTypes<Callable>::ArgsTuple &startArgs,
            const typename aux::IntegralTypes<Callable>::DeltasTuple &deltas) {
        return aux::calcQuadratureSum<Quadrature>(f, startArgs, deltas,
                                                     std::make_index_sequence<Quadrature::size>());
    }
}

#endif //ELECTROMAGNETIC_WAVES_SCATTERING_GAUSSQUADRATURE_HPP
