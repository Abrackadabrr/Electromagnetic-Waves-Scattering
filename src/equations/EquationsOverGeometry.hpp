//
// Created by evgen on 05.02.2025.
//

#ifndef EQUATIONSOVERGEOMETRY_HPP
#define EQUATIONSOVERGEOMETRY_HPP

#include "consepts/Consepts.hpp"
#include "mesh/SurfaceMesh.hpp"
#include "types/FunctionExtraction.hpp"

namespace EMW::Equations::MatrixForStructure {

namespace detail {
/**
 * Структура содержит в себе информацию о типах передаваемых аргументов в функцию, кроме многообразий.
 * Принимаются две функции (callable type): расчет самодействия и расчет взаимодействия двух разных
 */
template <typename SelfAffecting, typename OtherAffecting> struct auxiliary_info {
    // Дальше нужно статически вытащить те аргументы, которые идут после сеток в обеих функциях
    using self_args = Types::TypeTraits::argument_type_of_t<SelfAffecting>;
    using other_args = Types::TypeTraits::argument_type_of_t<OtherAffecting>;
    // Убираем из первой стопки для первых аргумента, из второй -- первый
    using self_additional_args = Types::TypeTraits::RemoveFirst<self_args>;
    using other_additional_args = Types::TypeTraits::RemoveFirstTwo<other_args>;
    // Сравниваем пачки типов дополнительных аргмуентов
    static_assert(std::is_same_v<self_additional_args, other_additional_args>);
    // Собираем tuple из необходимых дополнительных аргументов
    using tuple_additional_args = typename self_additional_args::Substitute<std::tuple>;
    using aux_sequence = std::make_index_sequence<self_additional_args::size>;
};

/**
 * Заготовка на расчет матрицы для обобщенного расчета произвольной сцены
 * Пока что это работает только для PeriodicSrtucture,
 * потому что завязано на один и тот же размер блока, но вскоре это перепишется
 *
 * @tparam TopologicalStructure -- класс, хранящий каким-либо образом множество многообразий
 * @tparam SelfAffecting
 * @tparam OtherAffecting
 */

/**
 * Функция расчета совместной матрицы
 * @param geometry геометрия
 * @param diagonal функиця для расчета самодействия
 * @param submatrix функция для расчета взаимодействия
 * @param args дополнительные аргументы
 * @param unnamed нужен для deducing index_sequence type
 * @return
 */
template <Concepts::GeomtricalStructure TopologicalStructure, typename SelfAffecting, typename OtherAffecting,
          Types::index... IndexSequenceForAdditionalArgumentsExpansion>
Types::MatrixXc compute(const TopologicalStructure &geometry, const SelfAffecting &diagonal,
                        const OtherAffecting &submatrix,
                        const typename auxiliary_info<SelfAffecting, OtherAffecting>::tuple_additional_args &args,
                        const std::index_sequence<IndexSequenceForAdditionalArgumentsExpansion...> &) {
    // Создаем итоговый результат
    // Размер блока в матрице
    const Types::index size_of_block =
        2 * (geometry.get(0).getCells().size() +
             geometry.get(0).getSubmesh(Mesh::IndexedCell::Tag::WAVEGUIDE_CROSS_SECTION).getCells().size());
    // Итоговый размер матрицы
    const Types::index size = size_of_block * geometry.size();
    // std::cout << size << std::endl;
    Types::MatrixXc A(size, size);
    // additional::(Вот тут выше как раз сильная привязка к специфике задачи с расчетом
    // токов на периодической системе волноводов.
    // При обобщении именно тут нужно будет убрать код и писать нормально)

    // 2) Поэтапный расчет внедиагональных блоков без учета того, что какие-то из них одинаковые
    for (int i = 0; i < geometry.size(); i++) {
        for (int j = 0; j < geometry.size(); j++)
            if (i != j) {
                // если у нас недиагональный блок, то считаем матрицу
                // A_ij блок показывает как j сетка влияет на поле в точках коллокации на i сетке
                const auto A_ij = submatrix(geometry.get(j), geometry.get(i),
                                            std::get<IndexSequenceForAdditionalArgumentsExpansion>(args)...);
                // рассчитывем место, где этот блок должен находится
                const Types::index first_row = i * size_of_block;
                const Types::index first_col = j * size_of_block;
                A.block(first_row, first_col, size_of_block, size_of_block) = A_ij;
            }
        A.block(i * size_of_block, i * size_of_block, size_of_block, size_of_block)=
            diagonal(geometry.get(i), std::get<IndexSequenceForAdditionalArgumentsExpansion>(args) ...);
    }
    return A;
}
}

/**
 * Функция, рассчитывающая дискретизованную систему уравнений поверх произвольной топологической структуры
 * @param geometry геометрия
 * @param diagonal функиця для расчета самодействия
 * @param submatrix функция для расчета взаимодействия
 * @param args дополнительные аргументы
 * @return
 */
template <Concepts::GeomtricalStructure TopologicalStructure, typename SelfAffecting, typename OtherAffecting>
Types::MatrixXc compute(const TopologicalStructure &geometry, const SelfAffecting &diagonal, const OtherAffecting &submatrix,
                               const typename detail::auxiliary_info<SelfAffecting, OtherAffecting>::tuple_additional_args &args) {
    return detail::compute(geometry, diagonal, submatrix, args, typename detail::auxiliary_info<SelfAffecting, OtherAffecting>::aux_sequence{});
}

}

#endif //EQUATIONSOVERGEOMETRY_HPP
