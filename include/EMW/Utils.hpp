//
// Created by evgen on 10.01.2025.
//

#ifndef ELECTROMAGNETIC_WAVES_SCATTERING_SRC_UTILS_HPP
#define ELECTROMAGNETIC_WAVES_SCATTERING_SRC_UTILS_HPP

#include "math/matrix/Matrix.hpp"
#include "types/Types.hpp"

namespace EMW::Utils
{
    template <typename T>
    class reverse
    {
        T& iterable_;

    public:
        explicit reverse(T& iterable) : iterable_{iterable}
        {
        }

        decltype(auto) begin() const { return std::rbegin(iterable_); }
        decltype(auto) end() const { return std::rend(iterable_); }
    };

    template <typename Container1, typename Container2>
    void to_csv(const Container1& cont1, const Container2& cont2, const std::string& name1, const std::string& name2,
                std::ostream& str, char delimiter = ',')
    {
        str << name1 << delimiter << name2 << "\n";
        for (int i = 0; i < cont1.size(); i++)
        {
            str << cont1[i] << delimiter << cont2[i] << '\n';
        }
    }


    struct MemoryUsage
    {
        Types::scalar full_matrix;
        Types::scalar toeplitz_matrix;
        Types::scalar toeplitz_and_factored_matrix;
    };

    inline std::ostream& operator<<(std::ostream& str, const MemoryUsage& usage)
    {
        str << "Full matrix memory usage: " << usage.full_matrix << " Gb \n";
        str << "Toeplitz matrix memory usage: " << usage.toeplitz_matrix << " Gb \n";
        if (usage.toeplitz_and_factored_matrix > 0)
            str << "Toeplitz and compressed matrix memory usage: " << usage.toeplitz_and_factored_matrix << " Gb \n";
        str << "Profit with toeplitz: " << usage.full_matrix / usage.toeplitz_matrix << "\n";
        if (usage.toeplitz_and_factored_matrix > 0)
            str << "Profit with toeplitz and compressed: " << usage.full_matrix / usage.toeplitz_and_factored_matrix;
        return str;
    }

    template <typename T>
    MemoryUsage get_memory_usage(const Math::LinAgl::Matrix::ToeplitzToeplitzBlock<T>& matrix)
    {
        const Types::index total_rows = matrix.rows();
        const Types::index total_cols = matrix.cols();

        const Types::index rows_in_big_block = matrix.rows_in_block();
        const Types::index cols_in_big_block = matrix.cols_in_block();

        const Types::index second_layer_rows = total_rows / rows_in_big_block;
        const Types::index second_layer_cols = total_cols / cols_in_big_block;

        const Types::index first_layer_rows = matrix.get_block(0, 0).rows() / matrix.get_block(0, 0).rows_in_block();
        const Types::index first_layer_cols = matrix.get_block(0, 0).cols() / matrix.get_block(0, 0).cols_in_block();

        const Types::index internal_block_rows = matrix.get_block(0, 0).rows_in_block();
        const Types::index internal_block_cols = matrix.get_block(0, 0).cols_in_block();

        const Types::scalar element_in_gb = static_cast<Types::scalar>(sizeof(T)) / (1024 * 1024 * 1024);
        const Types::scalar memory_for_full_matrix = total_cols * total_rows * element_in_gb;
        const Types::scalar memory_for_toeplitz_matrix = internal_block_cols * internal_block_rows *
            (2 * second_layer_rows - 1) * (2 * first_layer_rows - 1) *
            element_in_gb;
        return {memory_for_full_matrix, memory_for_toeplitz_matrix, -1};
    }

    template <typename T>
    MemoryUsage get_memory_usage(const Math::LinAgl::Matrix::TripleToeplitzBlock<T>& matrix)
    {
        const Types::index total_rows = matrix.rows();
        const Types::index total_cols = matrix.cols();

        const Types::index rows_in_big_block = matrix.rows_in_block();
        const Types::index third_layer_rows = total_rows / rows_in_big_block;

        const auto& second_layer = matrix.get_block(0, 0);
        const Types::index second_layer_rows = second_layer.rows() / second_layer.rows_in_block();

        const auto& first_layer = second_layer.get_block(0, 0);
        const Types::index first_layer_rows = first_layer.rows() / first_layer.rows_in_block();

        const Types::index internal_block_rows = first_layer.rows_in_block();
        const Types::index internal_block_cols = first_layer.cols_in_block();

        const Types::scalar element_in_gb = static_cast<Types::scalar>(sizeof(T)) / (1024 * 1024 * 1024);
        const Types::scalar memory_for_full_matrix = total_cols * total_rows * element_in_gb;
        const Types::scalar memory_for_toeplitz_matrix = internal_block_cols * internal_block_rows *
            (2 * third_layer_rows - 1) * (2 * second_layer_rows - 1) *
            (2 * first_layer_rows - 1) * element_in_gb;
        return {memory_for_full_matrix, memory_for_toeplitz_matrix, -1};
    }

    template <typename T>
    MemoryUsage get_memory_usage(const Math::LinAgl::Matrix::ToeplitzToeplitzDynFactoredBlock<T>& matrix)
    {
        const Types::index total_rows = matrix.rows();
        const Types::index total_cols = matrix.cols();

        const Types::index rows_in_big_block = matrix.rows_in_block();
        const Types::index cols_in_big_block = matrix.cols_in_block();

        const Types::index second_layer_rows = total_rows / rows_in_big_block;
        const Types::index second_layer_cols = total_cols / cols_in_big_block;

        const Types::index first_layer_rows = matrix.get_block(0, 0).rows() / matrix.get_block(0, 0).rows_in_block();
        const Types::index first_layer_cols = matrix.get_block(0, 0).cols() / matrix.get_block(0, 0).cols_in_block();

        const Types::index internal_block_rows = matrix.get_block(0, 0).rows_in_block();
        const Types::index internal_block_cols = matrix.get_block(0, 0).cols_in_block();

        const Types::scalar element_in_gb = static_cast<Types::scalar>(sizeof(T)) / (1024 * 1024 * 1024);
        const Types::scalar memory_for_full_matrix = total_cols * total_rows * element_in_gb;
        const Types::scalar memory_for_toeplitz_matrix = internal_block_cols * internal_block_rows *
            (2 * second_layer_rows - 1) * (2 * first_layer_rows - 1) *
            element_in_gb;

        Types::scalar toeplitz_factored_matrix = 0;

        for (const auto& second_layer_block : matrix.get_blocks().get_values())
            for (const auto& first_layer_block : second_layer_block.get_blocks().get_values())
                toeplitz_factored_matrix += first_layer_block.memory_usage();

        return {memory_for_full_matrix, memory_for_toeplitz_matrix, toeplitz_factored_matrix * element_in_gb};
    }

    template <typename T>
    MemoryUsage get_memory_usage(const Math::LinAgl::Matrix::TripleToeplitzFactoredBlock<T>& matrix)
    {
        const Types::index total_rows = matrix.rows();
        const Types::index total_cols = matrix.cols();

        const Types::index rows_in_big_block = matrix.rows_in_block();

        const Types::index third_layer_rows = total_rows / rows_in_big_block;

        const auto& second_layer = matrix.get_block(0, 0);
        const Types::index second_layer_rows = second_layer.rows() / second_layer.rows_in_block();

        const auto& first_layer = second_layer.get_block(0, 0);
        const Types::index first_layer_rows = first_layer.rows() / first_layer.rows_in_block();

        const Types::index internal_block_rows = first_layer.rows_in_block();
        const Types::index internal_block_cols = first_layer.cols_in_block();

        const Types::scalar element_in_gb = static_cast<Types::scalar>(sizeof(T)) / (1024 * 1024 * 1024);
        const Types::scalar memory_for_full_matrix = total_cols * total_rows * element_in_gb;
        const Types::scalar memory_for_toeplitz_matrix = internal_block_cols * internal_block_rows *
            (2 * third_layer_rows - 1) * (2 * second_layer_rows - 1) *
            (2 * first_layer_rows - 1) * element_in_gb;

        Types::scalar toeplitz_factored_matrix = 0;
        for (const auto& third_layer_block : matrix.get_blocks().get_values())
            for (const auto& second_layer_block : third_layer_block.get_blocks().get_values())
                for (const auto& first_layer_block : second_layer_block.get_blocks().get_values())
                    toeplitz_factored_matrix += first_layer_block.memory_usage();

        return {memory_for_full_matrix, memory_for_toeplitz_matrix, toeplitz_factored_matrix * element_in_gb};
    }

    template <typename T>
    size_t get_elements_for_parametrization(const Math::LinAgl::Matrix::TripleToeplitzFactoredBlock<T>& matrix)
    {
        const auto& second_layer = matrix.get_block(0, 0);

        const auto& first_layer = second_layer.get_block(0, 0);

        size_t toeplitz_factored_matrix = 0;
        for (const auto& third_layer_block : matrix.get_blocks().get_values())
            for (const auto& second_layer_block : third_layer_block.get_blocks().get_values())
                for (const auto& first_layer_block : second_layer_block.get_blocks().get_values())
                    toeplitz_factored_matrix += first_layer_block.memory_usage();
        return toeplitz_factored_matrix;
    }

    template <typename T>
    Types::scalar relative_frobenius_error(
        const Math::LinAgl::Matrix::TripleToeplitzBlock<T>& reference,
        const Math::LinAgl::Matrix::TripleToeplitzFactoredBlock<T>& approximated)
    {
        Types::scalar ref_frobenius_sq = 0;
        Types::scalar diff_frobenius_sq = 0;

#pragma omp parallel for collapse(2) reduction(+:ref_frobenius_sq) reduction(+:diff_frobenius_sq)
        for (Types::index i3 = 0; i3 < reference.rows_in_toeplitrz(); ++i3)
            for (Types::index j3 = 0; j3 < reference.cols_in_toeplitrz(); ++j3)
            {
                const auto& ref_second_layer = reference.get_block(i3, j3);
                const auto& appr_second_layer = approximated.get_block(i3, j3);

                for (Types::index i2 = 0; i2 < ref_second_layer.rows_in_toeplitrz(); ++i2)
                    for (Types::index j2 = 0; j2 < ref_second_layer.cols_in_toeplitrz(); ++j2)
                    {
                        const auto& ref_first_layer = ref_second_layer.get_block(i2, j2);
                        const auto& appr_first_layer = appr_second_layer.get_block(i2, j2);

                        for (Types::index i1 = 0; i1 < ref_first_layer.rows_in_toeplitrz(); ++i1)
                            for (Types::index j1 = 0; j1 < ref_first_layer.cols_in_toeplitrz(); ++j1)
                            {
                                const auto& ref_block = ref_first_layer.get_block(i1, j1);
                                const Types::MatrixX<T> appr_block = appr_first_layer.get_block(i1, j1).to_dense();

                                ref_frobenius_sq += ref_block.squaredNorm();
                                diff_frobenius_sq += (ref_block - appr_block).squaredNorm();
                            }
                    }
            }

        const Types::scalar ref_frobenius = std::sqrt(ref_frobenius_sq);
        const Types::scalar diff_frobenius = std::sqrt(diff_frobenius_sq);

        if (ref_frobenius == 0)
            return diff_frobenius;
        return diff_frobenius / ref_frobenius;
    }
} // namespace EMW::Utils

#endif // ELECTROMAGNETIC_WAVES_SCATTERING_SRC_UTILS_HPP
