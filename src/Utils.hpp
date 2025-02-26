//
// Created by evgen on 10.01.2025.
//

#ifndef UTILS_HPP
#define UTILS_HPP

#include "types/Types.hpp"

#include "math/matrix/Matrix.hpp"

namespace EMW::Utils {

template <typename Container>
void to_csv(const Container &cont1, const Container &cont2, const std::string &name1, const std::string &name2,
            std::ostream &str) {
    str << name1 << "," << name2 << "\n";
    for (int i = 0; i < cont1.size(); i++) {
        str << cont1[i] << ',' << cont2[i] << '\n';
    }
}

struct MemoryUsage {
    Types::scalar full_matrix;
    Types::scalar toeplitz_matrix;
};

inline std::ostream& operator<<(std::ostream &str, const MemoryUsage &usage) {
    str << "Full matrix memory usage: " << usage.full_matrix << " Gb \n";
    str << "Toeplitz matrix memory usage: " << usage.toeplitz_matrix << " Gb \n";
    str << "Profit: " << usage.full_matrix / usage.toeplitz_matrix;
    return str;
}

template <typename T> MemoryUsage get_memory_usage(const Math::LinAgl::Matrix::ToeplitzToeplitzBlock<T> &matrix) {
    const Types::index total_rows = matrix.rows();
    const Types::index total_cols = matrix.cols();

    const Types::index rows_in_big_block = matrix.rows_in_block();
    const Types::index cols_in_big_block = matrix.cols_in_block();

    const Types::index second_layer_rows = total_rows / rows_in_big_block;
    const Types::index second_layer_cols = total_cols / cols_in_big_block;

    const Types::index internal_block_rows = matrix.get_block(0, 0).rows_in_block();
    const Types::index internal_block_cols = matrix.get_block(0, 0).cols_in_block();

    const Types::scalar element_in_gb = static_cast<Types::scalar>(sizeof(T)) / (1024 * 1024 * 1024);
    const Types::scalar memory_for_full_matrix = total_cols * total_rows * element_in_gb;
    const Types::scalar memory_for_toeplitz_matrix = internal_block_cols * internal_block_rows *
                                                     (2 * second_layer_rows - 1) * (2 * second_layer_cols - 1) *
                                                     element_in_gb;
    return {memory_for_full_matrix, memory_for_toeplitz_matrix};
}
}


#endif //UTILS_HPP
