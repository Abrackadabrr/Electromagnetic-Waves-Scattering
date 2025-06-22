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

enum class CalculationMethod { Full, RSVD, ACA };

template <CalculationMethod CalculationMethod> struct CalcTraits {};

template <> struct CalcTraits<CalculationMethod::Full> {
    using ReturnType = Math::LinAgl::Matrix::ToeplitzToeplitzBlock<Types::complex_d>;
    using ToeplitzBlock = Math::LinAgl::Matrix::ToeplitzBlock<Types::complex_d>;
    using block_t = ToeplitzBlock::block_type;

    // Функция расчета
    static inline block_t create_a_block(const Mesh::SurfaceMesh &reference_mesh, const Mesh::SurfaceMesh &another_mesh,
                                         Types::complex_d k, Types::scalar a) {
        return WaveGuideWithActiveSection::submatrix(another_mesh, reference_mesh, a, k);
    }
};

#define COMPUTE_ERRORS_OF_APPROXIMATION 1

template <> struct CalcTraits<CalculationMethod::RSVD> {
    using ReturnType = Math::LinAgl::Matrix::ToeplitzToeplitzDynFactoredBlock<Types::complex_d>;
    using ToeplitzBlock = Math::LinAgl::Matrix::ToeplitzDynFactoredBlock<Types::complex_d>;
    using block_t = ToeplitzBlock::block_type;
    using compression = EMW::Math::LinAgl::Decompositions::ComplexRSVD;

    // ПАРАМЕТРЫ RSVD
    static Types::index rank;
    static constexpr Types::index calculate_rank(Types::scalar diff, Types::scalar a) {
        return static_cast<Types::index>(rank * (a / diff));
    };

    // Функция расчета
    static inline block_t create_a_block(const Mesh::SurfaceMesh &reference_mesh, const Mesh::SurfaceMesh &another_mesh,
                                         Types::complex_d k, Types::scalar a) {
        // считаем ранг, с которым надо приблизить
        const Types::scalar difference =
            (reference_mesh.getCells()[0].cellStructure.A - another_mesh.getCells()[0].cellStructure.A).norm();
        const Types::index r = calculate_rank(difference, a);
        // считаем полную матрицу
        const auto full_matrix = WaveGuideWithActiveSection::submatrix(another_mesh, reference_mesh, a, k);
        // делаем rsvd с заднным рангом
        return compression::compute(full_matrix, r, r);
    }
};

Types::index CalcTraits<CalculationMethod::RSVD>::rank = 100;

template <> struct CalcTraits<CalculationMethod::ACA> {
    using ReturnType = Math::LinAgl::Matrix::ToeplitzToeplitzDynFactoredBlock<Types::complex_d>;
    using ToeplitzBlock = Math::LinAgl::Matrix::ToeplitzDynFactoredBlock<Types::complex_d>;
    using block_t = ToeplitzBlock::block_type;
    using compression = EMW::Math::LinAgl::Decompositions::ComplexACA;

    // ПАРАМЕТРЫ ACA
    static constexpr Types::index error_controller = 0.1;

    // Функция расчета
    static inline block_t create_a_block(const Mesh::SurfaceMesh &reference_mesh, const Mesh::SurfaceMesh &another_mesh,
                                         Types::complex_d k, Types::scalar a) {
        // размеры самых внутренних блоков матрицы (нужня для АСА)
        const Types::index N_cols = reference_mesh.getCells().size();
        const Types::index K_cols =
            reference_mesh.getSubmesh(Mesh::IndexedCell::Tag::WAVEGUIDE_CROSS_SECTION).getCells().size();
        const Types::index cols = 2 * (N_cols + K_cols);
        const Types::index rows = cols;

        // тут дурацкая функция для расчета элемента матрицы
        const auto &ref_mesh_zero = reference_mesh.getSubmesh(Mesh::IndexedCell::Tag::WAVEGUIDE_CROSS_SECTION);
        const auto &another_mesh_zero = another_mesh.getSubmesh(Mesh::IndexedCell::Tag::WAVEGUIDE_CROSS_SECTION);
        // cells extraction
        const auto &ref_cells = reference_mesh.getCells();
        const auto &ref_zero_cells = ref_mesh_zero.getCells();
        const auto &an_cells = another_mesh.getCells();
        const auto &an_zero_cells = another_mesh_zero.getCells();
        const auto element_function = [&ref_cells, &ref_zero_cells, &an_cells, &an_zero_cells, k,
                                       a](Types::index row, Types::index col) -> Types::complex_d {
            return WaveGuideWithActiveSection::element_of_submatrix(row, col, ref_cells, ref_zero_cells, an_cells,
                                                                    an_zero_cells, a, k);
        };

        // тут непосредственно делаем адаптивный крест
        const auto result =
            Math::LinAgl::Decompositions::ComplexACA::compute(element_function, rows, cols, error_controller);

#if COMPUTE_ERRORS_OF_APPROXIMATION
        const auto full_block = CalcTraits<CalculationMethod::Full>::create_a_block(reference_mesh, another_mesh, k, a);
        const auto approximation = result.compute();
        std::cout << result.rows() << " " << result.cols() << std::endl;
        std::cout << full_block.rows() << " " << full_block.cols() << std::endl;
        const Types::MatrixXc error_matrix = (full_block - approximation);
        std::cout << 'm' << std::endl;
        std::cout << "Rank = " << result.get<0>().cols() << std::endl;
        std::cout << "Absolute error of approximation = " << error_matrix.norm() << std::endl;
        std::cout << "Relative error of approximation = " << error_matrix.norm() / full_block.norm() << std::endl;
#endif
        return result;
    }
};

template <CalculationMethod Calc, Types::index N1, Types::index N2>
typename CalcTraits<Calc>::ReturnType getMatrix(const Geometry::PeriodicStructure<N1, N2> &geometry, Types::scalar a,
                                                Types::complex_d k) {
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
                blocks[0] = std::move(block_t(std::vector{WaveGuideWithActiveSection::diagonal(reference_mesh, a, k)}));
            } else {
                // а тут обычный расчет
                // надо найти сетку, по которой считаем
                int col = reference_col + m;
                const auto &another_mesh = expanded_geometry.get(row_current_line, col);
                std::cout << row_current_line << ' ' << col << std::endl;
                // Расчет в зависимости от метода сжатия
                blocks[m] = CT::create_a_block(reference_mesh, another_mesh, k, a);
            }
        }
        // заполняем вертикальную часть вектора, отвечающего за внутреннюю теплицевость
        for (int m = 1; m < first_layer_rows; m++) {
            // тут снова обычный расчет
            // надо найти сетку, по которой считаем
            int col = reference_col - m;
            const int index_to_insert = first_layer_cols + m - 1;
            std::cout << row_current_line << ' ' << col << std::endl;
            const auto &another_mesh = expanded_geometry.get(row_current_line, col);
            // Расчет в зависимости от метода сжатия
            blocks[index_to_insert] = CT::create_a_block(reference_mesh, another_mesh, k, a);
        }
        internal_blocks[i] = std::move(ToeplitzBlock{first_layer_rows, first_layer_cols, std::move(blocks)});
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
            std::cout << row_current_line << ' ' << col << std::endl;
            // Расчет в зависимости от метода сжатия
            blocks[m] = CT::create_a_block(reference_mesh, another_mesh, k, a);
            }
            // заполняем вертикальную часть вектора, отвечающего за внутреннюю теплицевость
            for (int m = 1; m < first_layer_rows; m++) {
                // тут снова обычный расчет
                // надо найти сетку, по которой считаем
                int col = reference_col - m;
                const int index_to_insert = first_layer_cols + m - 1;
                const auto &another_mesh = expanded_geometry.get(row_current_line, col);
                std::cout << row_current_line << ' ' << col << std::endl;
                // Расчет в зависимости от метода сжатия
                blocks[index_to_insert] = CT::create_a_block(reference_mesh, another_mesh, k, a);
            }
            internal_blocks[second_layer_cols + i] =
                ToeplitzBlock{first_layer_rows, first_layer_cols, std::move(blocks)};
            // Проверка на то, что все корректно мувнулось куда нужно, а не скопировалось
            if (blocks.size() != 0)
                throw std::runtime_error("Bad move in generalized equations");
        }
    return {second_layer_rows, second_layer_cols, std::move(internal_blocks)};
};

}

#endif //GENERALEQUATION_HPP
