//
// Created by evgen on 18.03.2025.
//

#ifndef RSVD_HPP
#define RSVD_HPP

#include <algorithm>
#include <chrono>
#include <iostream>
#include <cmath>
#include "Eigen/Dense"

namespace EMW::Math::Matrix::Decopositions {

template<typename matrix_t, typename vector_t>
class RandomizedSvd {
public:
  RandomizedSvd(const matrix_t& m, int rank, int oversamples = 10, int iter = 2)
      : U_(), V_(), S_() {
    ComputeRandomizedSvd(m, rank, oversamples, iter);
  }

  vector_t singularValues() { return S_; }
  matrix_t matrixU() { return U_; }
  matrix_t matrixV() { return V_; }

private:
  matrix_t U_, V_;
  vector_t S_;

  /*
    Main function for randomized svd
    oversamples: additional samples/rank for accuracy, to account for random sampling
  */
  void ComputeRandomizedSvd(const matrix_t& A, int rank, int oversamples,
                            int iter) {
    using namespace std::chrono;

    // If matrix is too small for desired rank/oversamples
    if((rank + oversamples) > min(A.rows(), A.cols())) {
      rank = min(A.rows(), A.cols());
      oversamples = 0;
    }

    matrix_t Q = FindRandomizedRange(A, rank + oversamples, iter);
    matrix_t B = Q.transpose() * A;

    // Compute the SVD on the thin matrix (much cheaper than SVD on original)
    Eigen::JacobiSVD<matrix_t> svd(B, Eigen::ComputeThinU | Eigen::ComputeThinV);

    U_ = (Q * svd.matrixU()).block(0, 0, A.rows(), rank);
    V_ = svd.matrixV().block(0, 0, A.cols(), rank);
    S_ = svd.singularValues().head(rank);
  }

  /*
    Finds a set of orthonormal vectors that approximates the range of A
    Basic idea is that finding orthonormal basis vectors for A*W, where W is set of some
    random vectors w_i, can approximate the range of A
    Most of the time/computation in the randomized SVD is spent here
  */
  matrix_t FindRandomizedRange(const matrix_t& A, int size, int iter) {
    int nr = A.rows(), nc = A.cols();
    matrix_t L(nr, size);
    Eigen::FullPivLU<matrix_t> lu1(nr, size);
    matrix_t Q = matrix_t::Random(nc, size); // TODO should this be stack or dynamic allocation?
    Eigen::FullPivLU<matrix_t> lu2(nc, nr);

    // Conduct normalized power iterations
    // Intuition: multiply by A a few times to find a matrix Q that's "more in the range of A"
    //  Simply multiplying by A repeatedly makes alg unstable, so use LU to "normalize"
    // From Facebook implementation: "Please note that even n_iter=1 guarantees superb accuracy, whether or not there is any gap in the singular values of the matrix A being approximated"
    for (int i = 0; i < iter; ++i) {
      lu1.compute(A * Q);
      L.setIdentity();
      L.block(0, 0, nr, size).triangularView<Eigen::StrictlyLower>() =
          lu1.matrixLU();

      lu2.compute(A.transpose() * L);
      Q.setIdentity();
      Q.block(0, 0, nc, size).triangularView<Eigen::StrictlyLower>() =
          lu2.matrixLU();
    }

    Eigen::ColPivHouseholderQR<matrix_t> qr(A * Q);
    return qr.householderQ() * matrix_t::Identity(nr, size); // recover skinny Q matrix
  }
};

/*
  Computes spectral norm of error in reconstruction, from SVD matrices.

  Spectral norm = square root of maximum eigenvalue of matrix. Intuitively: the maximum 'scale', by which a matrix can 'stretch' a vector.
  Note: The definition of an eigenvalue is for square matrices. For non square matrices, we can define singular values: Definition: The singular values of a m×n matrix A are the positive square roots of the nonzero eigenvalues of the corresponding matrix A'A. The corresponding eigenvectors are called the singular vectors.
*/
double diff_spectral_norm(matrix_t A, matrix_t U, vector_t s, matrix_t V, int n_iter=20) {
  int nr = A.rows();

  vector_t y = vector_t::Random(nr);
  y.normalize();

  matrix_t B = (A - U*s.asDiagonal()*V.transpose());

  // TODO make this more efficient (don't explicitly calculate B)
  if(B.rows() != B.cols())
     B = B*B.transpose();

 // Run n iterations of the power method
 // TODO implement and compare fbpca's method
  for(int i=0; i<n_iter; ++i) {
    y = B*y;
    y.normalize();
  }
  double eigval = abs((B*y).dot(y) / y.dot(y));
  if(eigval==0) return 0;

  return sqrt(eigval);
}

}

#endif //RSVD_HPP
