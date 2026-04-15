//
// Created by evgen on 18.03.2025.
//

#ifndef RSVD_HPP
#define RSVD_HPP

#include "math/matrix/DynamicFactoredMatrix.hpp"

#include <Eigen/Dense>

namespace EMW::Math::LinAgl::Decompositions {

template <typename matrix_t, typename vector_t, typename value_t>
class RSVD {
    static vector_t get_ith_vector(Types::index i, Types::index N) {
        vector_t result = vector_t::Zero(N);
        result(i) = static_cast<value_t>(1);
        return result;
    }

  public:
    /**
     *  Function that returns DENSE matrix
     *  Function that returns householder sequence and eliminator can be used
     *  It's not obvious which one is faster

      Finds a set of orthonormal vectors that approximates the range of A
      Basic idea is that finding orthonormal basis vectors for A*W, where W is set of some
      random vectors w_i, can approximate the range of A
      Most of the time/computation in the randomized SVD is spent here
    */
    static matrix_t FindRandomizedRange(const matrix_t &A, int size, int iter) {
        int nr = A.rows(), nc = A.cols();
        // 1) create sketch matrix
        matrix_t Omega = matrix_t::Random(nc, size);

        // Conduct normalized power iterations
        // Intuition: multiply by A a few times to find a matrix Q that's "more in the range of A"
        //  Simply multiplying by A repeatedly makes alg unstable, so use LU to "normalize"
        // From Facebook implementation: "Please note that even n_iter=1 guarantees superb accuracy, whether or not
        // there is any gap in the singular values of the matrix A being approximated"
        matrix_t L(nr, size);
        Eigen::FullPivLU<matrix_t> lu1(nr, size);
        Eigen::FullPivLU<matrix_t> lu2(nc, nr);

        for (int i = 0; i < iter; ++i) {
            lu1.compute(A * Omega);
            L.setIdentity();
            L.block(0, 0, nr, size).template triangularView<Eigen::StrictlyLower>() = lu1.matrixLU();

            lu2.compute(A.transpose() * L);
            Omega.setIdentity();
            Omega.block(0, 0, nc, size).template triangularView<Eigen::StrictlyLower>() = lu2.matrixLU();
        }

        // 2) Orthogonalize columns
        Eigen::ColPivHouseholderQR<matrix_t> qr(A * Omega);
        const matrix_t r = qr.matrixR().template triangularView<Eigen::Upper>();

        // 3) Recovering skinny Q matrix
        matrix_t result = matrix_t::Zero(nr, size);
        const auto &householder_sequence = qr.householderQ();
        for (int i = 0; i < size; ++i) {
            result.col(i) = householder_sequence * get_ith_vector(i, nr);
        }
        return result;
    }

    static LinAgl::Matrix::DynamicFactoredMatrix<matrix_t> compute(const matrix_t &A, Types::index rank, Types::index oversamples, int iter = 0) {

        using diag_t = Types::DiagonalMatrixX<value_t>;

        Containers::vector<matrix_t> factors; factors.reserve(3);

        // If matrix is too small for desired rank/oversamples
        if ((rank + oversamples) > std::min(A.rows(), A.cols())) {
            rank = std::min(A.rows(), A.cols());
            oversamples = 0;
        }

        factors.emplace_back(FindRandomizedRange(A, rank + oversamples, iter));

        matrix_t B = factors[0].adjoint() * A;

        // 4) Compute the SVD on the thin matrix (much cheaper than SVD on original)
        Eigen::BDCSVD<matrix_t> svd(B, Eigen::ComputeThinU | Eigen::ComputeThinV);

        // 5) Remove oversampled eigenvalues and construct matrixes
        factors.push_back(svd.matrixU().block(0, 0, B.rows(), rank) * diag_t(svd.singularValues().head(rank)));
        factors.push_back(svd.matrixV().block(0, 0, A.cols(), rank).adjoint());

        // 6) Construct factored matrix and return
        return LinAgl::Matrix::DynamicFactoredMatrix(std::move(factors));
    }
};

} // namespace EMW::Math::Matrix::Decompositions

#endif // RSVD_HPP
