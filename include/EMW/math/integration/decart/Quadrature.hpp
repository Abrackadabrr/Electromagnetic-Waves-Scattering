//
// Created by evgen on 01.02.24.
//

#ifndef ELECTROMAGNETIC_WAVES_SCATTERING_QUADRATURE_HPP
#define ELECTROMAGNETIC_WAVES_SCATTERING_QUADRATURE_HPP

#include "types/FunctionExtraction.hpp"
#include "types/TypeTraits.hpp"
#include "types/Types.hpp"

namespace EMW::Math::Integration::Numerical::Decart {
namespace detail {
template <typename Arg> using DeltaType = decltype(std::declval<Arg>() - std::declval<Arg>());

template <typename FReturnType, typename... Args> struct IntegralTypes {
    using ArgsTuple = std::tuple<Args...>;
    using DeltasTuple = std::tuple<DeltaType<Args>...>;
    using MeasureType = decltype((std::declval<DeltaType<Args>>() * ...));
    using ResultType = FReturnType;
};

template <typename FReturnType, typename... Args>
IntegralTypes<FReturnType, Args...> extractIntegralTypes(const FReturnType & /*fReturnType*/,
                                                         const std::tuple<Args...> & /*args*/) {
    return {};
}

template <typename Callable>
using ExtructedIntegralTypes =
    decltype(extractIntegralTypes(std::declval<Types::TypeTraits::return_type_of_t<Callable>>(),
                                  std::declval<Types::TypeTraits::argument_tuple_type_of_t<Callable>>()));

template <typename ArgsTuple, typename DeltaTuple, typename Points, Types::index... indexes>
ArgsTuple createPoint(const ArgsTuple &start, const DeltaTuple &deltas, const Points &points,
                      const std::index_sequence<indexes...> & /*sequence*/) {
    return {std::get<indexes>(start) + (std::get<indexes>(deltas) / 2) * (std::get<indexes>(points) + 1)...};
}

template <typename DeltaTuple, Types::index... indexes>
decltype(auto) calcArea(const DeltaTuple &deltas, const std::index_sequence<indexes...> & /*sequence*/) {
    return (std::get<indexes>(deltas) * ...);
}

template <typename Callable, typename ArgTuple, Types::index... index>
Types::TypeTraits::return_type_of_t<Callable> unroll(const Callable &f, const ArgTuple &args,
                                                     std::index_sequence<index...>) {
    return f(std::get<index>(args)...);
}

template <typename Quadrature, typename Callable, Types::index... indexes>
typename ExtructedIntegralTypes<Callable>::ResultType
calcQuadratureSum(const Callable &f, const typename ExtructedIntegralTypes<Callable>::ArgsTuple &startArgs,
                  const typename ExtructedIntegralTypes<Callable>::DeltasTuple &deltas,
                  const std::index_sequence<indexes...> & /*sumSequence*/) {
    constexpr auto dimSize = std::tuple_size_v<typename ExtructedIntegralTypes<Callable>::ArgsTuple>;
    static_assert(Quadrature::dim == dimSize);
    constexpr auto pointSequence = std::make_index_sequence<dimSize>();

    return ((Quadrature::nodes[indexes].weight *
             unroll(f, createPoint(startArgs, deltas, Quadrature::nodes[indexes].point, pointSequence),
                    std::make_index_sequence<Quadrature::dim>{})) +
            ...);
}

/**
 * Переводит линейный индекс в тензорный интекс с правильными размерностями.
 * Сетка по размерностям передается в обратном порядке, начиная с предпоследнего
 */
constexpr std::tuple<size_t> flat_to_volume_index(size_t flat_index) { return flat_index; }

template <typename last_size, typename... sizes>
constexpr decltype(auto) flat_to_volume_index(size_t flat_index, last_size l_s, sizes... sizes_each_dim) {
    static_assert(std::is_same_v<last_size, size_t> && (std::is_same_v<sizes, size_t> && ...),
                  "All size types must be the same size_t");
    size_t dim_size_prod = l_s * (sizes_each_dim * ... * 1);
    size_t last_i = flat_index / dim_size_prod;
    return std::tuple_cat(flat_to_volume_index(flat_index % dim_size_prod, sizes_each_dim...), std::tuple{last_i});
}

/**
 * Аналогичные функции, но для работы с tuples
 *
 */
template <typename tuple_t, size_t... Idx>
constexpr decltype(auto) flat_to_volume_index(size_t flat_index, tuple_t &&index_tuple, std::index_sequence<Idx...>) {
    return flat_to_volume_index(flat_index, std::get<Idx>(index_tuple)...);
}

template <typename... tuple_args>
constexpr decltype(auto) flat_to_volume_index(size_t flat_index, const std::tuple<tuple_args...> &index_tuple) {
    return flat_to_volume_index(flat_index, index_tuple, std::make_index_sequence<sizeof...(tuple_args)>());
}

/**
 * Утилиты для угобного создания tuples
 *
 */
template <typename value_t, size_t... Idx>
constexpr decltype(auto) repeated_tuple_impl(value_t &&v, std::index_sequence<Idx...>) {
    return std::make_tuple([&v](size_t) { return v; }(Idx)...);
}

template <size_t dim, typename value_t> constexpr decltype(auto) repeated_tuple(value_t &&v) {
    return repeated_tuple_impl(std::forward<value_t>(v), std::make_index_sequence<dim>());
}

template <typename tuple_1, typename tuple_2, size_t... Idx>
constexpr decltype(auto) sum_impl(tuple_1 &&t1, tuple_2 &&t2, std::index_sequence<Idx...>) {
    return std::make_tuple((std::get<Idx>(t1) + std::get<Idx>(t2))...);
}

template <typename... args1, typename... args2>
requires requires { sizeof...(args1) == sizeof...(args2); }
constexpr decltype(auto) sum(std::tuple<args1...> t1, std::tuple<args2...> t2) {
    return sum_impl(t1, t2, std::make_index_sequence<sizeof...(args1)>());
}

template <typename tuple_1, typename tuple_2, size_t... Idx>
constexpr decltype(auto) mul_impl(tuple_1 &&t1, tuple_2 &&t2, std::index_sequence<Idx...>) {
    return std::make_tuple((std::get<Idx>(t1) * std::get<Idx>(t2))...);
}

template <typename... args1, typename... args2>
requires requires { sizeof...(args1) == sizeof...(args2); }
constexpr decltype(auto) mul(std::tuple<args1...> t1, std::tuple<args2...> t2) {
    return mul_impl(t1, t2, std::make_index_sequence<sizeof...(args1)>());
}
} // namespace detail

template <typename Quadrature, typename Callable>
constexpr typename detail::ExtructedIntegralTypes<Callable>::ResultType
quadrature_sum(const Callable &f, const typename detail::ExtructedIntegralTypes<Callable>::ArgsTuple &startArgs,
               const typename detail::ExtructedIntegralTypes<Callable>::DeltasTuple &deltas) {
    constexpr auto dimSize = std::tuple_size_v<typename detail::ExtructedIntegralTypes<Callable>::ArgsTuple>;
    static_assert(Quadrature::dim == dimSize);

    // инициализация дефолтным значением (нулем по смыслу)
    typename detail::ExtructedIntegralTypes<Callable>::ResultType result{};
    for (size_t index = 0; index < Quadrature::nodes.size(); ++index) {
        decltype(auto) point =
            detail::createPoint(startArgs, deltas, Quadrature::nodes[index].point, std::make_index_sequence<dimSize>());
        result += Quadrature::nodes[index].weight *
                  detail::unroll(f, point, std::make_index_sequence<dimSize>());
    }
    return result;
}

template <typename Quadrature, typename Callable>
typename detail::ExtructedIntegralTypes<Callable>::ResultType
integrate(const Callable &f, const typename detail::ExtructedIntegralTypes<Callable>::ArgsTuple &startArgs,
          const typename detail::ExtructedIntegralTypes<Callable>::DeltasTuple &deltas) {
    return detail::calcQuadratureSum<Quadrature>(f, startArgs, deltas, std::make_index_sequence<Quadrature::size>()) *
           (detail::calcArea(deltas, std::make_index_sequence<Quadrature::dim>{}) / (1 << Quadrature::dim));
}

/**
 * Квадратура для области == декартово произведение отрезков
 * с заданием уровня подразбиения на 2^level подобластей по каждой из сторон
 */
template <typename Quadrature, typename Callable, typename... ArgsType, typename... DeltasType>
constexpr typename detail::ExtructedIntegralTypes<Callable>::ResultType
quadrature_sum_with_decomposition(const Callable &f, const std::tuple<ArgsType...> &startArgs,
                                  const std::tuple<DeltasType...> &deltas, size_t level) {
    typename detail::ExtructedIntegralTypes<Callable>::ResultType result{};
    // Разбиение по сторонам
    size_t each_side_size = level;
    size_t total_indexes = std::pow(each_side_size, sizeof...(ArgsType));
    // Строим новые дельты
    auto new_deltas = std::apply(
        [each_side_size](auto &&...args) {
            return std::make_tuple(args / static_cast<Types::scalar>(each_side_size)...);
        },
        deltas);
    // Делаем массив из новых агрументов
    for (size_t flat_index = 0; flat_index < total_indexes; ++flat_index) {
        auto tensor_index =
            detail::flat_to_volume_index(flat_index, detail::repeated_tuple<sizeof...(ArgsType) - 1>(each_side_size));
        auto new_start = detail::sum(startArgs, detail::mul(new_deltas, tensor_index));
        result += quadrature_sum<Quadrature>(f, new_start, new_deltas);
    }
    return result / static_cast<Types::scalar>(total_indexes);
}

template <typename Quadrature, typename Callable, typename... ArgsType, typename... DeltasType>
typename detail::ExtructedIntegralTypes<Callable>::ResultType
constexpr integrate_with_decomposition(const Callable &f, const std::tuple<ArgsType...> &startArgs,
                             const std::tuple<DeltasType...> &deltas, size_t level) {
    return quadrature_sum_with_decomposition<Quadrature>(f, startArgs, deltas, level) *
        (detail::calcArea(deltas, std::make_index_sequence<Quadrature::dim>{}) / (1 << Quadrature::dim));
}

template <typename Quadrature, typename Callable, typename StopCriterion>
requires std::invocable<StopCriterion, typename detail::ExtructedIntegralTypes<Callable>::ResultType,
                        typename detail::ExtructedIntegralTypes<Callable>::ResultType>
constexpr std::pair<typename detail::ExtructedIntegralTypes<Callable>::ResultType, size_t>
adaptive_quadrature_sum(const Callable &f, const typename detail::ExtructedIntegralTypes<Callable>::ArgsTuple &startArgs,
                   const typename detail::ExtructedIntegralTypes<Callable>::DeltasTuple &deltas,
                   StopCriterion &&stopCriterion, size_t max_level) {
    decltype(auto) result = quadrature_sum<Quadrature>(f, startArgs, deltas);
    if (max_level == 1)
        return {result, 1};
    size_t level = 2;
    decltype(auto) result_adaptive = quadrature_sum_with_decomposition<Quadrature>(f, startArgs, deltas, level);
    while (!std::invoke(std::forward<StopCriterion>(stopCriterion), result, result_adaptive) && (level < max_level))
    {
        result = result_adaptive;
        result_adaptive = quadrature_sum_with_decomposition<Quadrature>(f, startArgs, deltas, ++level);
    }
    return {result_adaptive, level};
}

template <typename Quadrature, typename Callable, typename StopCriterion>
requires std::invocable<StopCriterion, typename detail::ExtructedIntegralTypes<Callable>::ResultType,
                        typename detail::ExtructedIntegralTypes<Callable>::ResultType>
constexpr std::pair<typename detail::ExtructedIntegralTypes<Callable>::ResultType, size_t>
adaptive_integrate(const Callable &f, const typename detail::ExtructedIntegralTypes<Callable>::ArgsTuple &startArgs,
                   const typename detail::ExtructedIntegralTypes<Callable>::DeltasTuple &deltas,
                   StopCriterion &&stopCriterion, size_t max_level) {
    decltype(auto) result = integrate<Quadrature>(f, startArgs, deltas);
    if (max_level == 1)
        return {result, 1};
    size_t level = 2;
    decltype(auto) result_adaptive = integrate_with_decomposition<Quadrature>(f, startArgs, deltas, level);
    while (!std::invoke(std::forward<StopCriterion>(stopCriterion), result, result_adaptive) && (level < max_level))
        {
            result = result_adaptive;
            result_adaptive = integrate_with_decomposition<Quadrature>(f, startArgs, deltas, ++level);
        }
        return {result_adaptive, level};
    }
}


#endif //ELECTROMAGNETIC_WAVES_SCATTERING_QUADRATURE_HPP
