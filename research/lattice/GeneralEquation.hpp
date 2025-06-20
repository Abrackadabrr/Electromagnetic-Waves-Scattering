//
// Created by evgen on 07.04.2025.
//

#ifndef GENERALEQUATION_HPP
#define GENERALEQUATION_HPP

#include "types/Types.hpp"

#include "math/matrix/DynamicFactoredMatrix.hpp"
#include "math/matrix/Matrix.hpp"
#include "math/matrix/decompositions/Decompositions.hpp"

#include "geometry/PeriodicStructure.hpp"

#include "research/lattice/SpecificLatticeEquations.hpp"

namespace Research::Lattice {

using namespace EMW;

enum class CalculationMethod {
    Full, RSVD, ACA
};

template<CalculationMethod CalculationMethod>
struct CalcTraits {};

template<>
struct CalcTraits<CalculationMethod::Full> {
    using ReturnType = Math::LinAgl::Matrix::ToeplitzToeplitzBlock<Types::complex_d>;
    using ToeplitzBlock = Math::LinAgl::Matrix::ToeplitzBlock<Types::complex_d>;
    using block_t = ToeplitzBlock::block_type;

    inline block_t create_a_block(const block_t& block) {
        return block;
    }
};

template<>
struct CalcTraits<CalculationMethod::RSVD> {
    using ReturnType = Math::LinAgl::Matrix::ToeplitzToeplitzDynFactoredBlock<Types::complex_d>;
    using ToeplitzBlock = Math::LinAgl::Matrix::ToeplitzDynFactoredBlock<Types::complex_d>;
    using block_t = ToeplitzBlock::block_type;
    using compression = EMW::Math::LinAgl::Decompositions::ComplexRSVD;

    // ПАРАМЕТРЫ RSVD
    static Types::index rank;
    static constexpr Types::index calculate_rank(Types::scalar diff, Types::scalar a) { return static_cast<Types::index>(rank * (a / diff)); };
};

Types::index CalcTraits<CalculationMethod::RSVD>::rank = 100;

template<>
struct CalcTraits<CalculationMethod::ACA> {
    using ReturnType = Math::LinAgl::Matrix::ToeplitzToeplitzDynFactoredBlock<Types::complex_d>;
    using ToeplitzBlock = Math::LinAgl::Matrix::ToeplitzDynFactoredBlock<Types::complex_d>;
    using block_t = ToeplitzBlock::block_type;
    using compression = EMW::Math::LinAgl::Decompositions::ComplexACA;

    // ПАРАМЕТРЫ ACA
    static constexpr Types::index error_controller = 1;
};

template <CalculationMethod Calc, Types::index N1, Types::index N2>
typename CalcTraits<Calc>::ReturnType
getMatrix(const Geometry::PeriodicStructure<N1, N2> &geometry, Types::scalar a, Types::complex_d k) {
    // алиасы для удобного пользования
    using CT = CalcTraits<Calc>;
    using ToeplitzBlock = typename CT::ToeplitzBlock;
    using block_t = typename CT::block_t;

    // Теперь соберем ту структуру, по которой мы будем производить расчет
    constexpr auto n1 = 2 * N1 - 1;
    constexpr auto n2 = 2 * N2 - 1;
    // тут геометрия содержит сетки с индексами [0, n1 - 1 == 2 N1 - 2] x [0, n2 - 1 == 2 N2 - 2]
    const Geometry::PeriodicStructure<n1, n2> expanded_geometry = geometry.expand_without_saving_nice_origin();

    // Тут обязательно есть сетка с origin == {0, 0, 0}, вот относительно неё и будем считать
    // Теперь нам нужно собрать все необходимые блоки в вектор и std::move-нуть его в нужное место
    // в структуре теплицевой матрицы

    constexpr auto first_layer_rows = N2;  // количество строк и столбцов
    constexpr auto first_layer_cols = N2;  // количество строк и столбцов
    constexpr auto second_layer_rows = N1; // количество строк и столбцов
    constexpr auto second_layer_cols = N1; // количество строк и столбцов

    // Собираем дважды тёплицеву матрицу
    Containers::vector<ToeplitzBlock> internal_blocks;
    // Собираем большой вектор из тёплицевых блоков
    internal_blocks.resize(ToeplitzBlock::get_size_of_container(second_layer_rows, second_layer_cols));
    // Находим сетку в середине, относительно которой и будем все считать
    constexpr int reference_row = N1 - 1;
    constexpr int reference_col = N2 - 1;
    const auto &reference_mesh = expanded_geometry.get(reference_row, reference_col);
    const Types::Vector3d reference_origin = expanded_geometry.get_origin(reference_row, reference_col);

    // размеры самых внутренних блоков матрицы (нужня для АСА, но могут и в будущем понадобиться)
    const Types::index N_cols = reference_mesh.getCells().size();
    const Types::index K_cols = reference_mesh..getSubmesh(Mesh::IndexedCell::Tag::WAVEGUIDE_CROSS_SECTION).getCells().size();
    const Types::index cols = 2 * (N_cols + K_cols);
    const Types::index rows = cols;

    // заполняем горизонтальную часть вектора, отвечающего за внешнюю теплицевость
    for (int i = 0; i < second_layer_cols; i++) {
        Containers::vector<block_t> blocks;
        Types::index first_layer_size = ToeplitzBlock::get_size_of_container(first_layer_rows, first_layer_cols);
        blocks.resize(first_layer_size);
        int row_current_line = reference_row + i;
        // заполняем горизонтальную часть вектора, отвечающего за внутреннюю теплицевость
        for (int m = 0; m < first_layer_cols; m++) {
            if ((i == 0) && (m == 0)) {
                // вычисление диагонального блока, единственного несжимаемого
                // во всей большой матрице

                const auto block = WaveGuideWithActiveSection::diagonal(reference_mesh, a, k);

                blocks[0] = std::move(block_t({WaveGuideWithActiveSection::diagonal(reference_mesh, a, k)}));

//                std::cout << row_current_line << ' ' <<  reference_col + m << std::endl;

            } else {
                // а тут обычный расчет
                // надо найти сетку, по которой считаем
                int col = reference_col + m;
                const auto &another_mesh = expanded_geometry.get(row_current_line, col);
                std::cout << row_current_line << ' ' <<  col << std::endl;
                auto &where_to_insert = blocks[m];

                // Этот кусок выносим в отдельный код для каждого метода декомпозиции блока и радуемся уму-разуму
                if constexpr (Calc == CalculationMethod::Full)
                    // пихаем полную матрицу
                    where_to_insert = std::move(WaveGuideWithActiveSection::submatrix(another_mesh, reference_mesh, a, k));
                else if constexpr (Calc == CalculationMethod::RSVD) {
                    // считаем ранг, с которым надо приблизить
                    const Types::scalar difference = (expanded_geometry.get_origin(row_current_line, col) - reference_origin).norm();
                    const Types::index r = CT::calculate_rank(difference, a);
                    // считаем полную матрицу
                    const auto full_matrix = WaveGuideWithActiveSection::submatrix(another_mesh, reference_mesh, a, k);
                    // делаем rsvd с заднным рангом
                    where_to_insert = std::move(CT::compression::compute(full_matrix, r, r));
                } else if constexpr (Calc == CalculationMethod::ACA) {
                    // тут дурацкая функция для расчета элемента матрицы
                    const auto ref_mesh_zero & = reference_mesh.getSubmesh(Mesh::IndexedCell::Tag::WAVEGUIDE_CROSS_SECTION);
                    const auto another_mesh_zero & = another_mesh.getSubmesh(Mesh::IndexedCell::Tag::WAVEGUIDE_CROSS_SECTION);
                    const auto element_function = [&reference_mesh, &ref_mesh_zero,
                                                   &another_mesh, &another_mesh_zero,
                                                   k, a] (Types::index row, Types::index col)-> Types::complex_d {
                        return WaveGuideWithActiveSection::element_of_submatrix(row, col, another_mesh, another_mesh_zero,
                                                                                reference_mesh, ref_mesh_zero, a, k);
                    };

                    // тут непосредственно делаем адаптивный крест
                    where_to_insert = std::move(Math::LinAgl::Decompositions::ComplexACA::compute(element_function,
                                                                                   rows,
                                                                                     cols,
                                                                                     CT::error_controller));
            }
        }
        // заполняем вертикальную часть вектора, отвечающего за внутреннюю теплицевость
        for (int m = 1; m < first_layer_rows; m++) {
            // тут снова обычный расчет
            // надо найти сетку, по которой считаем
            int col = reference_col - m;
            const int index_to_insert = first_layer_cols + m - 1;

            std::cout << row_current_line << ' ' <<  col << std::endl;
            const auto &another_mesh = expanded_geometry.get(row_current_line, col);
            auto &where_to_insert = blocks[index_to_insert];
            if constexpr (is_full)
                where_to_insert =
                std::move(WaveGuideWithActiveSection::submatrix(another_mesh, reference_mesh, a, k));
            else {
                const Types::scalar difference = (expanded_geometry.get_origin(row_current_line, col) - reference_origin).norm();
                const Types::index r = CT::calculate_rank(difference, a);
                const auto full_matrix = WaveGuideWithActiveSection::submatrix(another_mesh, reference_mesh, a, k);
                where_to_insert = std::move(CT::compression::compute(full_matrix, r, r));
            }
        }
        internal_blocks[i] = ToeplitzBlock{first_layer_rows, first_layer_cols, std::move(blocks)};
        // Проверка на то, что все корректно мувнулось куда нужно, а не скопировалось
        if (blocks.size() != 0)
            throw std::runtime_error("Bad move in generalized equations");
    }
    // заполняем вертикальную часть вектора, отвечающего за внешнюю теплицевость
    for (int i = 0; i < second_layer_rows - 1; i++) {
        Containers::vector<block_t> blocks;
        Types::index first_layer_size = ToeplitzBlock::get_size_of_container(first_layer_rows, first_layer_cols);
        blocks.resize(first_layer_size);
        int row_current_line = reference_row - (i + 1);
        // заполняем горизонтальную часть вектора, отвечающего за внутреннюю теплицевость
        for (int m = 0; m < first_layer_cols; m++) {
            // а тут обычный расчет
            // надо найти сетку, по которой считаем
            int col = reference_col + m;
            const auto &another_mesh = expanded_geometry.get(row_current_line, col);
            std::cout << row_current_line << ' ' <<  col << std::endl;
            auto &where_to_insert = blocks[m];
            if constexpr (is_full)
                where_to_insert = std::move(WaveGuideWithActiveSection::submatrix(another_mesh, reference_mesh, a, k));
            else {
                const Types::scalar difference = (expanded_geometry.get_origin(row_current_line, col) - reference_origin).norm();
                const Types::index r = CT::calculate_rank(difference, a);
                const auto full_matrix = WaveGuideWithActiveSection::submatrix(another_mesh, reference_mesh, a, k);
                where_to_insert = std::move(CT::compression::compute(full_matrix, r, r));
            }
        }
        // заполняем вертикальную часть вектора, отвечающего за внутреннюю теплицевость
        for (int m = 1; m < first_layer_rows; m++) {
            // тут снова обычный расчет
            // надо найти сетку, по которой считаем
            int col = reference_col - m;
            const int index_to_insert = first_layer_cols + m - 1;
            const auto &another_mesh = expanded_geometry.get(row_current_line, col);

            std::cout << row_current_line << ' ' <<  col << std::endl;

            auto &where_to_insert = blocks[index_to_insert];
            if constexpr (is_full)
                where_to_insert =
                std::move(WaveGuideWithActiveSection::submatrix(another_mesh, reference_mesh, a, k));
            else {
                const Types::scalar difference = (expanded_geometry.get_origin(row_current_line, col) - reference_origin).norm();
                const Types::index r = CT::calculate_rank(difference, a);
                const auto full_matrix = WaveGuideWithActiveSection::submatrix(another_mesh, reference_mesh, a, k);
                where_to_insert = std::move(CT::compression::compute(full_matrix, r, r));
            }
        }
        internal_blocks[second_layer_cols + i] = ToeplitzBlock{first_layer_rows, first_layer_cols, std::move(blocks)};
        // Проверка на то, что все корректно мувнулось куда нужно, а не скопировалось
        if (blocks.size() != 0)
            throw std::runtime_error("Bad move in generalized equations");
    }

    return {second_layer_rows, second_layer_cols, std::move(internal_blocks)};
};

template <CalculationMethod Calc, Types::index N1, Types::index N2>
typename CalcTraits<Calc>::ReturnType
getMatrix(const Geometry::PeriodicStructure<N1, N2> &geometry, Types::scalar a, Types::complex_d k) {
    // алиасы для удобного пользования
    using CT = CalcTraits<Calc>;
    using ToeplitzBlock = typename CT::ToeplitzBlock;
    using block_t = typename CT::block_t;

    std::cout << "Compressed with rank = " << CalcTraits<Calc>::rank << std::endl;

    constexpr bool is_full = Calc == CalculationMethod::Full;

    // Теперь соберем ту структуру, по которой мы будем производить расчет
    constexpr auto n1 = 2 * N1 - 1;
    constexpr auto n2 = 2 * N2 - 1;
    // тут геометрия содержит сетки с индексами [0, n1 - 1 == 2 N1 - 2] x [0, n2 - 1 == 2 N2 - 2]
    const Geometry::PeriodicStructure<n1, n2> expanded_geometry = geometry.expand_without_saving_nice_origin();

    // Тут обязательно есть сетка с origin == {0, 0, 0}, вот относительно неё и будем считать
    // Теперь нам нужно собрать все необходимые блоки в вектор и std::move-нуть его в нужное место
    // в структуре теплицевой матрицы

    constexpr auto first_layer_rows = N2;  // количество строк и столбцов
    constexpr auto first_layer_cols = N2;  // количество строк и столбцов
    constexpr auto second_layer_rows = N1; // количество строк и столбцов
    constexpr auto second_layer_cols = N1; // количество строк и столбцов

    // Собираем дважды тёплицеву матрицу
    Containers::vector<ToeplitzBlock> internal_blocks;
    // Собираем большой вектор из тёплицевых блоков
    internal_blocks.resize(ToeplitzBlock::get_size_of_container(second_layer_rows, second_layer_cols));
    // Находим сетку в середине, относительно которой и будем все считать
    constexpr int reference_row = N1 - 1;
    constexpr int reference_col = N2 - 1;
    const auto &reference_mesh = expanded_geometry.get(reference_row, reference_col);
    const Types::Vector3d reference_origin = expanded_geometry.get_origin(reference_row, reference_col);

    // заполняем горизонтальную часть вектора, отвечающего за внешнюю теплицевость
    for (int i = 0; i < second_layer_cols; i++) {
        Containers::vector<block_t> blocks;
        Types::index first_layer_size = ToeplitzBlock::get_size_of_container(first_layer_rows, first_layer_cols);
        blocks.resize(first_layer_size);
        int row_current_line = reference_row + i;
        // заполняем горизонтальную часть вектора, отвечающего за внутреннюю теплицевость
        for (int m = 0; m < first_layer_cols; m++) {
            if ((i == 0) && (m == 0)) {
                // вычисление матрицы
                const auto block = WaveGuideWithActiveSection::diagonal(reference_mesh, a, k);

                if constexpr (is_full)
                    blocks[0] = std::move(block);
                else
                    blocks[0] = std::move(block_t({WaveGuideWithActiveSection::diagonal(reference_mesh, a, k)}));

                std::cout << row_current_line << ' ' <<  reference_col + m << std::endl;

            } else {
                // а тут обычный расчет
                // надо найти сетку, по которой считаем
                int col = reference_col + m;
                const auto &another_mesh = expanded_geometry.get(row_current_line, col);
                std::cout << row_current_line << ' ' <<  col << std::endl;
                auto &where_to_insert = blocks[m];
                if constexpr (is_full)
                    where_to_insert = std::move(WaveGuideWithActiveSection::submatrix(another_mesh, reference_mesh, a, k));
                else {
                    const Types::scalar difference = (expanded_geometry.get_origin(row_current_line, col) - reference_origin).norm();
                    const Types::index r = CT::calculate_rank(difference, a);
                    const auto full_matrix = WaveGuideWithActiveSection::submatrix(another_mesh, reference_mesh, a, k);
                    where_to_insert = std::move(CT::compression::compute(full_matrix, r, r));
                }
            }
        }
        // заполняем вертикальную часть вектора, отвечающего за внутреннюю теплицевость
        for (int m = 1; m < first_layer_rows; m++) {
            // тут снова обычный расчет
            // надо найти сетку, по которой считаем
            int col = reference_col - m;
            const int index_to_insert = first_layer_cols + m - 1;

            std::cout << row_current_line << ' ' <<  col << std::endl;
            const auto &another_mesh = expanded_geometry.get(row_current_line, col);
            auto &where_to_insert = blocks[index_to_insert];
            if constexpr (is_full)
                where_to_insert =
                std::move(WaveGuideWithActiveSection::submatrix(another_mesh, reference_mesh, a, k));
            else {
                const Types::scalar difference = (expanded_geometry.get_origin(row_current_line, col) - reference_origin).norm();
                const Types::index r = CT::calculate_rank(difference, a);
                const auto full_matrix = WaveGuideWithActiveSection::submatrix(another_mesh, reference_mesh, a, k);
                where_to_insert = std::move(CT::compression::compute(full_matrix, r, r));
            }
        }
        internal_blocks[i] = ToeplitzBlock{first_layer_rows, first_layer_cols, std::move(blocks)};
        // Проверка на то, что все корректно мувнулось куда нужно, а не скопировалось
        if (blocks.size() != 0)
            throw std::runtime_error("Bad move in generalized equations");
    }
    // заполняем вертикальную часть вектора, отвечающего за внешнюю теплицевость
    for (int i = 0; i < second_layer_rows - 1; i++) {
        Containers::vector<block_t> blocks;
        Types::index first_layer_size = ToeplitzBlock::get_size_of_container(first_layer_rows, first_layer_cols);
        blocks.resize(first_layer_size);
        int row_current_line = reference_row - (i + 1);
        // заполняем горизонтальную часть вектора, отвечающего за внутреннюю теплицевость
        for (int m = 0; m < first_layer_cols; m++) {
            // а тут обычный расчет
            // надо найти сетку, по которой считаем
            int col = reference_col + m;
            const auto &another_mesh = expanded_geometry.get(row_current_line, col);
            std::cout << row_current_line << ' ' <<  col << std::endl;
            auto &where_to_insert = blocks[m];
            if constexpr (is_full)
                where_to_insert = std::move(WaveGuideWithActiveSection::submatrix(another_mesh, reference_mesh, a, k));
            else {
                const Types::scalar difference = (expanded_geometry.get_origin(row_current_line, col) - reference_origin).norm();
                const Types::index r = CT::calculate_rank(difference, a);
                const auto full_matrix = WaveGuideWithActiveSection::submatrix(another_mesh, reference_mesh, a, k);
                where_to_insert = std::move(CT::compression::compute(full_matrix, r, r));
            }
        }
        // заполняем вертикальную часть вектора, отвечающего за внутреннюю теплицевость
        for (int m = 1; m < first_layer_rows; m++) {
            // тут снова обычный расчет
            // надо найти сетку, по которой считаем
            int col = reference_col - m;
            const int index_to_insert = first_layer_cols + m - 1;
            const auto &another_mesh = expanded_geometry.get(row_current_line, col);

            std::cout << row_current_line << ' ' <<  col << std::endl;

            auto &where_to_insert = blocks[index_to_insert];
            if constexpr (is_full)
                where_to_insert =
                std::move(WaveGuideWithActiveSection::submatrix(another_mesh, reference_mesh, a, k));
            else {
                const Types::scalar difference = (expanded_geometry.get_origin(row_current_line, col) - reference_origin).norm();
                const Types::index r = CT::calculate_rank(difference, a);
                const auto full_matrix = WaveGuideWithActiveSection::submatrix(another_mesh, reference_mesh, a, k);
                where_to_insert = std::move(CT::compression::compute(full_matrix, r, r));
            }
        }
        internal_blocks[second_layer_cols + i] = ToeplitzBlock{first_layer_rows, first_layer_cols, std::move(blocks)};
        // Проверка на то, что все корректно мувнулось куда нужно, а не скопировалось
        if (blocks.size() != 0)
            throw std::runtime_error("Bad move in generalized equations");
    }

    return {second_layer_rows, second_layer_cols, std::move(internal_blocks)};
};


}

}

#endif //GENERALEQUATION_HPP
