//
// Created by evgen on 05.02.2025.
//

#ifndef EQUATIONSOVERGEOMETRY_HPP
#define EQUATIONSOVERGEOMETRY_HPP

#include "mesh/SurfaceMesh.hpp"
#include "types/FunctionExtraction.hpp"
#include "consepts/Consepts.hpp"

namespace EMW::Equations {
/**
 * Заготовка на расчет матрицы для обобщенного расчета произвольной сцены
 * Пока что это работает только для PeriodicSrtucture,
 * потому что завязано на один и тот же размер блока, но вскоре это перепишется
 *
 * @tparam TopologicalStructure -- класс, хранящий каким-либо образом многообразия
 * @tparam SelfAffecting
 * @tparam OtherAffecting
 */
template <Concepts::GeomtricalStructure TopologicalStructure, typename SelfAffecting, typename OtherAffecting> class MatrixFor {
    // Дальше нужно статически вытащить те аргументы, которые идут после сеток в обеих функциях
    using self_args = Types::TypeTraits::argument_type_of_t<SelfAffecting>;

    static_assert(std::is_same_v<
                  self_args, Types::TypeTraits::TypePack<const Mesh::SurfaceMesh &, Types::scalar, Types::complex_d>>);

    using other_args = Types::TypeTraits::argument_type_of_t<OtherAffecting>;
    using matfold_pack = Types::TypeTraits::TypePack<const Mesh::SurfaceMesh &>;
    // Убираем из первой стопки для первых аргумента, из второй -- первый
    using self_additional_args = Types::TypeTraits::RemoveFirst<self_args>;
    using other_additional_args = Types::TypeTraits::RemoveFirstTwo<other_args>;
    // Сравниваем пачки типов дополнительных аргмуентов
    static_assert(std::is_same_v<self_additional_args, other_additional_args>);
    // Собираем tuple из необходимых дополнительных аргументов
    using tuple_additional_args = typename self_additional_args::Substitute<std::tuple>;

    static_assert(std::is_same_v<tuple_additional_args, std::tuple<Types::scalar, Types::complex_d>>);

    using aux_sequence = std::make_index_sequence<self_additional_args::size>;

    // И вот тут уже нормальные вычисления, а не на типах
    /**
     *
     * @param geometry геометрия
     * @param diagonal функиця для расчета самодействия
     * @param submatrix функция для расчета взаимодействия
     * @param args дополнительные аргументы
     * @return
     */
    template<Types::index ... IndexSequenceForAdditionalArgumentsExpansion>
    static Types::MatrixXc compute(const TopologicalStructure &geometry, SelfAffecting diagonal,
                                   OtherAffecting submatrix, const tuple_additional_args &args,
                                   const std::index_sequence<IndexSequenceForAdditionalArgumentsExpansion...>&) {
        // Создаем итоговый результат
        // Размер итоговой матрицы
        const Types::index size_of_block =
            2 * (geometry.get(0).getCells().size() +
                 geometry.get(0).getSubmesh(Mesh::IndexedCell::Tag::WAVEGUIDE_CROSS_SECTION).getCells().size());
        const Types::index size = size_of_block * geometry.size();
        // std::cout << size << std::endl;
        Types::MatrixXc A(size, size);

        // 2) Поэтапный расчет внедиагональных блоков без учета того, что какие-то из них одинаковые
        for (int i = 0; i < geometry.size(); i++) {
            for (int j = 0; j < geometry.size(); j++)
                if (i != j) {
                    // если у нас недиагональный блок, то считаем матрицу
                    // A_ij блок показывает как k сетка влияет на поле в точках коллокации на j сетке
                    const auto A_ij =
                        submatrix(geometry.get(j), geometry.get(i),
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

public:
    /**
     * Функция, рассчитывающая дискретизованную систему уравнений поверх произвольной топологической структуры
     * @param geometry геометрия
     * @param diagonal функиця для расчета самодействия
     * @param submatrix функция для расчета взаимодействия
     * @param args дополнительные аргументы
     * @return
     */
    static Types::MatrixXc compute(const TopologicalStructure &geometry, SelfAffecting diagonal,
                               OtherAffecting submatrix, const tuple_additional_args &args) {
        return compute(geometry, diagonal, submatrix, args, aux_sequence{});
    }
};

}

#endif //EQUATIONSOVERGEOMETRY_HPP
