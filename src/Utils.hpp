//
// Created by evgen on 10.01.2025.
//

#ifndef UTILS_HPP
#define UTILS_HPP

#include "math/matrix/Matrix.hpp"
#include "types/Types.hpp"

namespace EMW::Utils {

template <typename T> class reverse {
    T &iterable_;

  public:
    explicit reverse(T &iterable) : iterable_{iterable} {}
    decltype(auto) begin() const { return std::rbegin(iterable_); }
    decltype(auto) end() const { return std::rend(iterable_); }
};

template <typename Container1, typename Container2>
void to_csv(const Container1 &cont1, const Container2 &cont2, const std::string &name1, const std::string &name2,
            std::ostream &str, char delimiter=',') {
    str << name1 << delimiter << name2 << "\n";
    for (int i = 0; i < cont1.size(); i++) {
        str << cont1[i] << delimiter << cont2[i] << '\n';
    }
}

struct MemoryUsage {
    Types::scalar full_matrix;
    Types::scalar toeplitz_matrix;
    Types::scalar toeplitz_and_factored_matrix;
};

inline std::ostream &operator<<(std::ostream &str, const MemoryUsage &usage) {
    str << "Full matrix memory usage: " << usage.full_matrix << " Gb \n";
    str << "Toeplitz matrix memory usage: " << usage.toeplitz_matrix << " Gb \n";
    if (usage.toeplitz_and_factored_matrix > 0)
        str << "Toeplitz and compressed matrix memory usage: " << usage.toeplitz_and_factored_matrix << " Gb \n";
    str << "Profit with toeplitz: " << usage.full_matrix / usage.toeplitz_matrix << "\n";
    if (usage.toeplitz_and_factored_matrix > 0)
        str << "Profit with toeplitz and compressed: " << usage.full_matrix / usage.toeplitz_and_factored_matrix;
    return str;
}

template <typename T> MemoryUsage get_memory_usage(const Math::LinAgl::Matrix::ToeplitzToeplitzBlock<T> &matrix) {
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
MemoryUsage get_memory_usage(const Math::LinAgl::Matrix::ToeplitzToeplitzDynFactoredBlock<T> &matrix) {
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

    for (const auto &second_layer_block : matrix.get_blocks().get_values())
        for (const auto& first_layer_block : second_layer_block.get_blocks().get_values())
            toeplitz_factored_matrix += first_layer_block.memory_usage();

    return {memory_for_full_matrix, memory_for_toeplitz_matrix, toeplitz_factored_matrix * element_in_gb};
}

} // namespace EMW::Utils

#endif // UTILS_HPP
