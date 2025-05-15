//
// Created by evgen on 15.05.2025.
//

#include <gtest/gtest.h>

#include "math/matrix/decompositions/adaptive_cross.hpp"

#include "mat_decomp.hpp"

TEST_F(MATRIX_DECOMPOSITIONS_TESTS, ACA_SINUS) {
    const Types::index N = 1000;
    const auto mat = get_sinus_mat(N);

    const auto rsvd = Math::LinAgl::Decompositions::RealRSVD::compute(mat, 2, 2);

#if 0
        std::cout << rsvd.matrixQ() * rsvd.matrixU() << std::endl;
        std::cout << Types::DiagonalMatrixXd(rsvd.singularValues())  << std::endl;
        std::cout << rsvd.matrixVh() << std::endl;

        std::cout << rsvd.matrixQ() * rsvd.matrixU() * Types::DiagonalMatrixXd(rsvd.singularValues()) * rsvd.matrixVh() << std::endl;
#endif

    const auto &factored_matrix = rsvd;
    const Types::scalar err = (factored_matrix.compute() - mat).norm() / mat.norm();
    ASSERT_NEAR(err, 0, 1e-13);
    std::cout << err << std::endl;
}
