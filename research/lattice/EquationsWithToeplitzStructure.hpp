//
// Created by evgen on 03.03.2025.
//

#ifndef EQUATIONSWITHTOEPLITZSTRUCTURE_HPP
#define EQUATIONSWITHTOEPLITZSTRUCTURE_HPP

#include "types/Types.hpp"

#include "math/matrix/Matrix.hpp"
#include "math/matrix/DynamicFactoredMatrix.hpp"
#include "math/matrix/decompositions/Decompositions.hpp"

#include "geometry/PeriodicStructure.hpp"

#include "research/lattice/GeneralizedEquations.hpp"

namespace Research::Lattice {

using namespace EMW;

template <Types::index N1, Types::index N2>
Math::LinAgl::Matrix::ToeplitzToeplitzBlock<Types::complex_d>
getMatrix(const Geometry::PeriodicStructure<N1, N2> &geometry, Types::scalar a, Types::complex_d k) {

    // алиасы для удобного пользования
    using ToeplitzBlock = Math::LinAgl::Matrix::ToeplitzBlock<Types::complex_d>;
    using block_t = ToeplitzBlock::block_type;

    // Теперь соберем ту структуру, по которой мы будем производить расчет
    constexpr auto n1 = 2 * N1 - 1;
    constexpr auto n2 = 2 * N2 - 1;
    // тут геометрия содержит сетки с индексами [0, n1 - 1 == 2 N1 - 2] x [0, n2 - 1 == 2 N2 - 2]
    const Geometry::PeriodicStructure<n1, n2> expanded_geometry = geometry.expand();

    // Тут обязательно есть сетка с origin == {0, 0, 0}, вот относительно неё и будем считать
    // Теперь нам нужно собрать все необходимые блоки в вектор и std::move-нуть его в нужное место
    // в структуре теплицевой матрицы

    constexpr auto first_layer_rows = N2;  // количество строк и столбцов
    constexpr auto first_layer_cols = N2;  // количество строк и столбцов
    constexpr auto second_layer_rows = N1; // количество строк и столбцов
    constexpr auto second_layer_cols = N1; // количество строк и столбцов

    // Собираем дважды тёплицеву матрицу
    Containers::vector<ToeplitzBlock> internal_blocks;
    internal_blocks.reserve(ToeplitzBlock::get_size_of_container(second_layer_rows, second_layer_cols));
    // Собираем большой вектор из тёплицевых блоков

    // Находим сетку в середине, относительно которой и будем все считать
    constexpr auto reference_row = N1 - 1;
    constexpr auto reference_col = N2 - 1;
    const auto &reference_mesh = expanded_geometry.get(reference_row, reference_col);

    // заполняем горизонтальную часть вектора, отвечающего за внешнюю теплицевость

    std::vector<std::pair<int, int>> indexes;
    indexes.reserve(n1 * n2);

    for (int i = 0; i < second_layer_cols; i++) {
        Containers::vector<block_t> blocks;
        Types::index first_layer_size = ToeplitzBlock::get_size_of_container(first_layer_cols, first_layer_cols);
        blocks.reserve(first_layer_size);
        // заполняем горизонтальную часть вектора, отвечающего за внутреннюю теплицевость
        for (int m = 0; m < first_layer_cols; ++m) {
            if ((i == 0) && (m == 0)) {
                int row = reference_row;
                int col = reference_col;
                indexes.emplace_back(row, col);

                blocks.emplace_back(WaveGuideWithActiveSection::diagonal(reference_mesh, a, k));
            } else {
                // а тут обычный расчет
                // надо найти сетку, по которой считаем
                int row = reference_row + i;
                int col = reference_col + m;
                indexes.emplace_back(row, col);

                const auto &another_mesh = expanded_geometry.get(row, col);
                blocks.emplace_back(WaveGuideWithActiveSection::submatrix(another_mesh, reference_mesh, a, k));
            }
        }
        // заполняем вертикальную часть вектора, отвечающего за внутреннюю теплицевость
        for (int m = first_layer_rows - 1; m >= 1; --m) {
            // тут снова обычный расчет
            // надо найти сетку, по которой считаем
            int row = reference_row + i;
            int col = reference_col - m;
            indexes.emplace_back(row, col);

            const auto &another_mesh = expanded_geometry.get(row, col);
            blocks.emplace_back(WaveGuideWithActiveSection::submatrix(another_mesh, reference_mesh, a, k));
        }
        internal_blocks.emplace_back(first_layer_rows, first_layer_cols, std::move(blocks));
        // Проверка на то, что все корректно мувнулось куда нужно, а не скопировалось
        if (blocks.size() != 0)
            throw std::runtime_error("Bad memory usage in TTB matrix assembling for lattice while assembling a row");
    }
    // заполняем вертикальную часть вектора, отвечающего за внешнюю теплицевость
    for (int i = 0; i < second_layer_rows - 1; i++) {
        Containers::vector<block_t> blocks;
        Types::index first_layer_size = ToeplitzBlock::get_size_of_container(first_layer_cols, first_layer_cols);
        blocks.reserve(first_layer_size);
        // заполняем горизонтальную часть вектора, отвечающего за внутреннюю теплицевость
        for (int m = 0; m < first_layer_cols; ++m) {
            // а тут обычный расчет
            // надо найти сетку, по которой считаем
            int row = i;
            int col = reference_col + m;
            indexes.emplace_back(row, col);

            const auto &another_mesh = expanded_geometry.get(row, col);
            blocks.emplace_back(WaveGuideWithActiveSection::submatrix(another_mesh, reference_mesh, a, k));
        }
        // заполняем вертикальную часть вектора, отвечающего за внутреннюю теплицевость
        for (int m = first_layer_rows - 1; m >= 1; --m) {
            // тут снова обычный расчет
            // надо найти сетку, по которой считаем
            int row = i;
            int col = reference_col - m;
            indexes.emplace_back(row, col);

            const auto &another_mesh = expanded_geometry.get(row, col);
            blocks.emplace_back(WaveGuideWithActiveSection::submatrix(another_mesh, reference_mesh, a, k));
        }
        internal_blocks.emplace_back(first_layer_rows, first_layer_cols, std::move(blocks));
        // Проверка на то, что все корректно мувнулось куда нужно, а не скопировалось
        if (blocks.size() != 0)
            throw std::runtime_error("Bad memory usage in TTB matrix assembling for lattice while assembling a col");
    }

    return {second_layer_rows, second_layer_cols, std::move(internal_blocks)};
};

template <Types::index N1, Types::index N2>
Math::LinAgl::Matrix::ToeplitzToeplitzDynFactoredBlock<Types::complex_d>
getMatrixCompressed(const Geometry::PeriodicStructure<N1, N2> &geometry, Types::scalar a, Types::complex_d k) {

#define INFO 1

    // алиасы для удобного пользования
    using ToeplitzBlock = Math::LinAgl::Matrix::ToeplitzDynFactoredBlock<Types::complex_d>;
    using block_t = ToeplitzBlock::block_type;
    using compression = EMW::Math::LinAgl::Decompositions::ComplexRSVD;

    // ПАРАМЕТРЫ RSVD
    const Types::index rank = 100;

    // Теперь соберем ту структуру, по которой мы будем производить расчет
    constexpr auto n1 = 2 * N1 - 1;
    constexpr auto n2 = 2 * N2 - 1;
    // тут геометрия содержит сетки с индексами [0, n1 - 1 == 2 N1 - 2] x [0, n2 - 1 == 2 N2 - 2]
    const Geometry::PeriodicStructure<n1, n2> expanded_geometry = geometry.expand();

    // Тут обязательно есть сетка с origin == {0, 0, 0}, вот относительно неё и будем считать
    // Теперь нам нужно собрать все необходимые блоки в вектор и std::move-нуть его в нужное место
    // в структуре теплицевой матрицы

    constexpr auto first_layer_rows = N2;  // количество строк и столбцов
    constexpr auto first_layer_cols = N2;  // количество строк и столбцов
    constexpr auto second_layer_rows = N1; // количество строк и столбцов
    constexpr auto second_layer_cols = N1; // количество строк и столбцов

    // Собираем дважды тёплицеву матрицу
    Containers::vector<ToeplitzBlock> internal_blocks;
    internal_blocks.reserve(ToeplitzBlock::get_size_of_container(second_layer_rows, second_layer_cols));
    // Собираем большой вектор из тёплицевых блоков

    // Находим сетку в середине, относительно которой и будем все считать
    constexpr auto reference_row = N1 - 1;
    constexpr auto reference_col = N2 - 1;
    const auto &reference_mesh = expanded_geometry.get(reference_row, reference_col);

    // заполняем горизонтальную часть вектора, отвечающего за внешнюю теплицевость

    std::vector<std::pair<int, int>> indexes;
    indexes.reserve(n1 * n2);

#if INFO
    Types::index counter = 0;
#endif

    for (int i = 0; i < second_layer_cols; i++) {
        Containers::vector<block_t> blocks;
        Types::index first_layer_size = ToeplitzBlock::get_size_of_container(first_layer_cols, first_layer_cols);
        blocks.reserve(first_layer_size);
        // заполняем горизонтальную часть вектора, отвечающего за внутреннюю теплицевость
        for (int m = 0; m < first_layer_cols; ++m) {
            if ((i == 0) && (m == 0)) {
                // тут считаем диагональный блок
                int row = reference_row;
                int col = reference_col;
                indexes.emplace_back(row, col);
                blocks.emplace_back(block_t({WaveGuideWithActiveSection::diagonal(reference_mesh, a, k)}));

#if INFO
    std::cout << "Block " << counter++ << " out of " << n1 * n2 << std::endl;
#endif
            } else {
                // а тут обычный расчет
                // надо найти сетку, по которой считаем
                int row = reference_row + i;
                int col = reference_col + m;
                indexes.emplace_back(row, col);

                const auto &another_mesh = expanded_geometry.get(row, col);
                blocks.emplace_back(compression::compute(
                    WaveGuideWithActiveSection::submatrix(another_mesh, reference_mesh, a, k), rank, rank));

#if INFO
                std::cout << "Block " << counter++ << " out of " << n1 * n2 << std::endl;
#endif

            }
        }
        // заполняем вертикальную часть вектора, отвечающего за внутреннюю теплицевость
        for (int m = first_layer_rows - 1; m >= 1; --m) {
            // тут снова обычный расчет
            // надо найти сетку, по которой считаем
            int row = reference_row + i;
            int col = reference_col - m;
            indexes.emplace_back(row, col);

            const auto &another_mesh = expanded_geometry.get(row, col);
            blocks.emplace_back(compression::compute(
                WaveGuideWithActiveSection::submatrix(another_mesh, reference_mesh, a, k), rank, rank));
#if INFO
            std::cout << "Block " << counter++ << " out of " << n1 * n2 << std::endl;
#endif

        }
        internal_blocks.emplace_back(first_layer_rows, first_layer_cols, std::move(blocks));
        // Проверка на то, что все корректно мувнулось куда нужно, а не скопировалось
        if (blocks.size() != 0)
            throw std::runtime_error("Bad memory usage in TTB matrix assembling for lattice while assembling a row");
    }
    // заполняем вертикальную часть вектора, отвечающего за внешнюю теплицевость
    for (int i = 0; i < second_layer_rows - 1; i++) {
        Containers::vector<block_t> blocks;
        Types::index first_layer_size = ToeplitzBlock::get_size_of_container(first_layer_cols, first_layer_cols);
        blocks.reserve(first_layer_size);
        // заполняем горизонтальную часть вектора, отвечающего за внутреннюю теплицевость
        for (int m = 0; m < first_layer_cols; ++m) {
            // а тут обычный расчет
            // надо найти сетку, по которой считаем
            int row = i;
            int col = reference_col + m;
            indexes.emplace_back(row, col);

            const auto &another_mesh = expanded_geometry.get(row, col);
            blocks.emplace_back(compression::compute(
                WaveGuideWithActiveSection::submatrix(another_mesh, reference_mesh, a, k), rank, rank));

#if INFO
            std::cout << "Block " << counter++ << " out of " << n1 * n2 << std::endl;
#endif

        }
        // заполняем вертикальную часть вектора, отвечающего за внутреннюю теплицевость
        for (int m = first_layer_rows - 1; m >= 1; --m) {
            // тут снова обычный расчет
            // надо найти сетку, по которой считаем
            int row = i;
            int col = reference_col - m;
            indexes.emplace_back(row, col);

            const auto &another_mesh = expanded_geometry.get(row, col);
            blocks.emplace_back(compression::compute(
                WaveGuideWithActiveSection::submatrix(another_mesh, reference_mesh, a, k), rank, rank));

#if INFO
            std::cout << "Block " << counter++ << " out of " << n1 * n2 << std::endl;
#endif

        }
        internal_blocks.emplace_back(first_layer_rows, first_layer_cols, std::move(blocks));
        // Проверка на то, что все корректно мувнулось куда нужно, а не скопировалось
        if (blocks.size() != 0) throw std::runtime_error("Bad memory usage in TTB matrix assembling for lattice while assembling a col");
    }

    return {second_layer_rows, second_layer_cols, std::move(internal_blocks)};
};
}

#endif //EQUATIONSWITHTOEPLITZSTRUCTURE_HPP
