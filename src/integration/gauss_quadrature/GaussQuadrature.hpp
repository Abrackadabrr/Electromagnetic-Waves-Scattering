//
// Created by evgen on 01.02.24.
//

#ifndef ELECTROMAGNETIC_WAVES_SCATTERING_GAUSSQUADRATURE_HPP
#define ELECTROMAGNETIC_WAVES_SCATTERING_GAUSSQUADRATURE_HPP

#include "types/Types.hpp"
#include "types/TypeTraits.hpp"
#include "types/FunctionExtraction.hpp"

namespace EMW::DefiniteIntegrals {

    namespace detail {

        template <typename Arg>
        using DeltaType = decltype(std::declval<Arg>() - std::declval<Arg>());

        template <typename FReturnType, typename... Args>
        struct IntegralTypes {
            using ArgsTuple = std::tuple<Args...>;
            using DeltasTuple = std::tuple<DeltaType<Args>...>;
            using MeasureType = decltype((std::declval<DeltaType<Args>>() * ...));
            using ResultType = decltype(std::declval<MeasureType>() * std::declval<FReturnType>());
        };

        template <typename FReturnType, typename... Args>
        IntegralTypes<FReturnType, Args...> extractIntegralTypes(const FReturnType& /*fReturnType*/,
                                                                 const std::tuple<Args...>& /*args*/) {
            return {};
        }

        template <typename Callable>
        using ExtructedIntegralTypes =
                decltype(extractIntegralTypes(std::declval<Types::TypeTraits::return_type_of_t<Callable>>(),
        std::declval<Types::TypeTraits::argument_tuple_type_of_t<Callable>>()));

        template <typename ArgsTuple, typename DeltaTuple, typename Points, Types::index... indexes>
        ArgsTuple createPoint(const ArgsTuple& start, const DeltaTuple& deltas, const Points& points,
                              const std::index_sequence<indexes...>& /*sequence*/) {
            return {std::get<indexes>(start)  + (std::get<indexes>(deltas) / 2) * (std::get<indexes>(points) + 1) ...};
        }

        template <typename DeltaTuple, Types::index... indexes>
        decltype(auto) calcArea(const DeltaTuple& deltas, const std::index_sequence<indexes...>& /*sequence*/) {
            return (std::get<indexes>(deltas) * ...);
        }

        template <typename Quadrature, typename Callable, Types::index... indexes>
        typename ExtructedIntegralTypes<Callable>::ResultType calcQuadratureSum(
                const Callable& f, const typename ExtructedIntegralTypes<Callable>::ArgsTuple& startArgs,
                const typename ExtructedIntegralTypes<Callable>::DeltasTuple& deltas,
                const std::index_sequence<indexes...>& /*sumSequence*/) {
            constexpr auto dimSize = std::tuple_size_v<typename ExtructedIntegralTypes<Callable>::ArgsTuple>;
            static_assert(Quadrature::dim == dimSize);
            constexpr auto pointSequence = std::make_index_sequence<dimSize>();

            return ((Quadrature::nodes[indexes].weight *
                     std::apply(f, createPoint(startArgs, deltas, Quadrature::nodes[indexes].point, pointSequence))) +
                           ...) *
                   (calcArea(deltas, pointSequence) / std::pow(2, dimSize));
        }

    }  // namespace detail

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
    template <typename Quadrature, typename Callable>
    typename detail::ExtructedIntegralTypes<Callable>::ResultType integrate(
            const Callable& f, const typename detail::ExtructedIntegralTypes<Callable>::ArgsTuple& startArgs,
            const typename detail::ExtructedIntegralTypes<Callable>::DeltasTuple& deltas) {
        return detail::calcQuadratureSum<Quadrature>(f, startArgs, deltas, std::make_index_sequence<Quadrature::size>());
    }

//template <typename Quadrature, typename Callable>
//typename detail::ExtructedIntegralTypes<Callable>::ResultType integrate(
//    const Callable& f, const typename detail::ExtructedIntegralTypes<Callable>::ArgsTuple& startArgs,
//    const typename detail::ExtructedIntegralTypes<Callable>::ArgsTuple& endArgs,
//    const typename detail::ExtructedIntegralTypes<Callable>::DeltasTuple& steps) {
//
//    return detail::calcQuadratureSum<Quadrature>(f, startArgs, deltas, std::make_index_sequence<Quadrature::size>());
//}

}

#endif //ELECTROMAGNETIC_WAVES_SCATTERING_GAUSSQUADRATURE_HPP
