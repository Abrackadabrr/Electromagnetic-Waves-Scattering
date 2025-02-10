//
// Created by evgen on 10.02.2025.
//

#ifndef TOEPLITZ_HPP
#define TOEPLITZ_HPP

#include "types/Types.hpp"

namespace EMW::Math::LinAgl::Matrix {



template <typename scalar_t>
const typename ToeplitzMatrix<scalar_t>::block_t &ToeplitzMatrix<scalar_t>::getBlock(Types::index row, Types::index col) const {

}


template <typename scalar_t>
Types::VectorX<scalar_t> operator*(const ToeplitzMatrix<scalar_t> &matrix, const Types::VectorX<scalar_t> &vector) {
    // создаем нулевой вектор результата, в который будем записывать ответ
    typename ToeplitzMatrix<scalar_t>::col_t result;

    DataType *sub_vec = static_cast<DataType *>(std::malloc(msub_col * sizeof(DataType)));
    DataType *local_res = static_cast<DataType *>(std::malloc(msub_row * sizeof(DataType)));

    for (uint64_t i = 0; i < mblock_row; ++i) {
        for (uint64_t j = 0; j < mblock_col; ++j) {

            BLAS::copy(msub_col, x + j * msub_col, 1, sub_vec, 1);
            int block_id = get_block(i, j);

            // std::cout << "res for block " << i << " " << j << " : #" << get_block(i, j) << std::endl;
            // mblock[block_id].print_mat();

            // std::cout << "sub_vector of x :" << std::endl;
            // for (uint64_t k = 0; k < msub_col; ++k) {
            // std::cout << sub_vec[k] << "  ";
            // }
            // std::cout << std::endl;

            mblock[block_id].matvec(sub_vec, local_res, Option::N);

            // std::cout << "local result of matvec : " << std::endl;
            // for (uint64_t k = 0; k < msub_row; ++k) {
            // std::cout << local_res[k] << "  ";
            // }
            // std::cout << std::endl;

            // std::cout << "y before axpy : " << std::endl;
            // for (uint64_t k = 0; k < total_row(); ++k) {
            // std::cout << y[k] << "  ";
            // }
            // std::cout << std::endl;

            BLAS::axpy(msub_row, DataType(1.0), local_res, 1, y + i * msub_row, 1);

            // std::cout << "y after axpy : " << std::endl;
            // for (uint64_t k = 0; k < total_row(); ++k) {
            // std::cout << y[k] << "  ";
            // }
            // std::cout << std::endl;
            // std::cout << "-------------------------------- " << std::endl;
        }
        }

        std::free(sub_vec);
        std::free(local_res);
    }
}
}
#endif //TOEPLITZ_HPP
