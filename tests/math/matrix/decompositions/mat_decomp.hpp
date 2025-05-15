//
// Created by evgen on 15.05.2025.
//

#ifndef MAT_DECOMP_HPP
#define MAT_DECOMP_HPP

#include <gtest/gtest.h>
#include "types/Types.hpp"

using namespace EMW;

class MATRIX_DECOMPOSITIONS_TESTS : public ::testing::Test {
protected:
    Types::scalar sin_mat_el(Types::index i, Types::index j) { return std::sin(i + j); }

    Types::MatrixXd get_sinus_mat(Types::index N) {
        Types::MatrixXd sin_mat = Types::MatrixXd::Zero(N, N);
        for (Types::index i = 0; i < N; i++) {
            for (Types::index j = 0; j < N; j++) {
                sin_mat(i, j) = sin_mat_el(i, j);
            }
        }
        return sin_mat;
    }
};


#endif //MAT_DECOMP_HPP
