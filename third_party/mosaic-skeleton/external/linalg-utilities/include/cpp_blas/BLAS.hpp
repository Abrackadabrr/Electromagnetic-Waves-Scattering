/*
 * Copyright 2023-2026 Bulat Valiakhmetov
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */


#ifndef LINALG_UTILITIES_CPP_BLAS_BLAS_HPP
#define LINALG_UTILITIES_CPP_BLAS_BLAS_HPP


#include <omp.h>

#include <cmath>
#include <cstdint>
#include <complex>
#include <type_traits>


namespace BLAS
{
    //Level 1: scalar and vector operations
    
    // Swaps vectors x and y
    template <class DataType>
        void swap(const uint64_t &n,
                  const DataType *x, const int64_t &incx,
                  DataType *y, const int64_t &incy)
        {
#pragma omp parallel for
            for (uint64_t i = 0; i < n; i++)
            {
                std::swap(x[i * incx], y[i * incy]);
            }
        }

    // Adds x scaled by alpha to y:
    //       y = \alpha x + y
    template <class DataType>
        void axpy(const uint64_t &n,
                  const DataType &alpha,
                  const DataType *x, const int64_t &incx,
                  DataType *y, const int64_t &incy)
        {
#pragma omp parallel for
            for (uint64_t i = 0; i < n; i++)
            {
                y[i * incy] += x[i * incx] * alpha;
            }
        }

    // Copies strided vector x to strided vector y
    template <class DataType>
        void copy(const uint64_t &n,
                  const DataType *x, const int64_t &incx,
                  DataType *y, const int64_t &incy)
        {
#pragma omp parallel for
            for (uint64_t i = 0; i < n; i++)
            {
                y[i * incy] = x[i * incx];
            }
        }

    // Computes dot product of conjugated strided vector x by strided vector y
    template <class DataType>
        DataType dotc(const uint64_t &n,
                  const DataType *x, const int64_t &incx,
                  const DataType *y, const int64_t &incy)
        {
            DataType r{0.0};
#pragma omp parallel
            {
                DataType lr{0.0};
#pragma omp for nowait
                for (uint64_t i = 0; i < n; i++)
                {
                    lr += std::conj(x[i * incx]) * y[i * incy];
                }
#pragma omp critical
                {
                    r += lr;
                }
            }
            return r;
        }

    // Computes dot product of UNconjugated strided vector x by strided vector y
    template <class DataType>
    DataType dotu(const uint64_t &n,
                  const DataType *x, const int64_t &incx,
                  const DataType *y, const int64_t &incy)
    {
        DataType r{0.0};
#pragma omp parallel
        {
            DataType lr{0.0};
#pragma omp for nowait
            for (uint64_t i = 0; i < n; i++)
            {
                lr += x[i * incx] * y[i * incy];
            }
#pragma omp critical
            {
                r += lr;
            }
        }
        return r;
    }

    // Computes Euclidian (second, mean square) norm of strided vector x
    template <class DataType>
        decltype(std::abs(DataType{})) nrm2(const uint64_t &n,
                                            const DataType *x, const int64_t &incx)
        {
            typedef decltype(std::abs(DataType{})) RealType;
            RealType r{0.0};
#pragma omp parallel
            {
                RealType lr{0.0};
#pragma omp for nowait
                for (uint64_t i = 0; i < n; i++)
                {
                    lr += std::norm(x[i * incx]);
                }
#pragma omp critical
                {
                    r += lr;
                }
            }
            return std::sqrt(r);
        }

    // Scales strided vector by real number (for real DataTypes is the same as scal)
    template <class DataType>
        void rscal(const uint64_t &n,
                  const decltype(std::abs(DataType{})) &alpha,
                  DataType *x, const int64_t &incx)
        {
#pragma omp parallel for
            for (uint64_t i = 0; i < n; i++)
            {
                x[i * incx] *= alpha;
            }
        }

    // Scales strided vector
    template <class DataType>
        void scal(const uint64_t &n,
                  const DataType &alpha,
                  DataType *x, const int64_t &incx)
        {
#pragma omp parallel for
            for (uint64_t i = 0; i < n; i++)
            {
                x[i * incx] *= alpha;
            }
        }

    // Rotates strided vectors using provided cosine and sine value
    template <class DataType>
        void rot(const uint64_t &n,
                 DataType *x, const uint64_t &incx,
                 DataType *y, const uint64_t &incy,
                 const decltype(std::abs(DataType{})) &c, const DataType &s)
        {
#pragma omp parallel for
            for (uint64_t i = 0; i < n; i++)
            {
                DataType t = c * x[i * incx] + s * y[i * incy];
                y[i * incy] = -std::conj(s) * x[i * incx] + c * y[i * incy];
                x[i * incx] = t;
            }
        }
    //Level 2: matrix-vector operations

    //General dense matrix rank-1 update with conjugation of second vector:
    //       A = \alpha x y^H + A
    template <class DataType>
        void gerc(const uint64_t &m, const uint64_t &n,
                const DataType &alpha,
                const DataType *x, const int64_t &incx,
                const DataType *y, const int64_t &incy,
                DataType *A, const uint64_t &ldA)
        {
#pragma omp parallel for
            for (uint64_t j = 0; j < n; j++)
            {
                const DataType alpha_y = alpha * std::conj(y[j * incy]);
                for (uint64_t i = 0; i < m; i++)
                {
                    A[i + j * ldA] += x[i * incx] * alpha_y;
                }
            }
        }

    //General dense matrix rank-1 update WITHOUT conjugation of second vector:
    //       A = \alpha x y^T + A
    template <class DataType>
    void geru(const uint64_t &m, const uint64_t &n,
              const DataType &alpha,
              const DataType *x, const int64_t &incx,
              const DataType *y, const int64_t &incy,
              DataType *A, const uint64_t &ldA)
    {
#pragma omp parallel for
        for (uint64_t j = 0; j < n; j++)
        {
            const DataType alpha_y = alpha * y[j * incy];
            for (uint64_t i = 0; i < m; i++)
            {
                A[i + j * ldA] += x[i * incx] * alpha_y;
            }
        }
    }
    
    //General dense matrix by vector multiplication:
    //       y = \alpha A x + \beta y
    template <class DataType>
        void gemv(const char &transpose,
                const uint64_t &m, const uint64_t &n,
                const DataType &alpha,
                const DataType *A, const uint64_t &ldA,
                const DataType *x, const int64_t &incx,
                const DataType beta,
                DataType *y, const int64_t &incy)
        {
#pragma omp parallel
            {
                if (transpose == 'n' || transpose == 'N')
                {
                    uint64_t P = omp_get_num_threads();
                    uint64_t p = omp_get_thread_num();
                    uint64_t sm = (m / P) * p + std::min(p, m % P);
                    uint64_t em = (m / P) * (p + 1) + std::min(p + 1, m % P);
                    if (beta == DataType{0.0})
                    {
                        for (uint64_t i = sm; i < em; i++)
                        {
                            y[i * incy] = DataType{0.0};
                        }
                    }
                    else
                    {
                        for (uint64_t i = sm; i < em; i++)
                        {
                            y[i * incy] *= beta;
                        }
                    }
                    for (uint64_t j = 0; j < n; j++)
                    {
                        for (uint64_t i = sm; i < em; i++)
                        {
                            y[i * incy] += alpha * A[i + ldA * j] * x[j * incx];
                        }
                    }
                }
                else if (transpose == 't' || transpose == 'T')
                {
                    if (beta == DataType{0.0})
                    {
#pragma omp for
                        for (uint64_t i = 0; i < n; i++)
                        {
                            y[i * incy] = DataType{0.0};
                        }
                    }
#pragma omp for
                    for (uint64_t j = 0; j < n; j++)
                    {
                        y[j * incy] *= beta;
                        for (uint64_t i = 0; i < m; i++)
                        {
                            y[j * incy] += alpha * A[i + j * m] * x[i * incx];
                        }
                    }
                }
                else if (transpose == 'c' || transpose == 'C')
                {
                    if (beta == DataType{0.0})
                    {
#pragma omp for
                        for (uint64_t i = 0; i < n; i++)
                        {
                            y[i * incy] = DataType{0.0};
                        }
                    }
#pragma omp for
                    for (uint64_t j = 0; j < n; j++)
                    {
                        y[j * incy] *= beta;
                        for (uint64_t i = 0; i < m; i++)
                        {
                            y[j * incy] += alpha * std::conj(A[i + j * m]) * x[i * incx];
                        }
                    }
                }
                else if (transpose == 'z' || transpose == 'Z')
                {
                    uint64_t P = omp_get_num_threads();
                    uint64_t p = omp_get_thread_num();
                    uint64_t sm = (m / P) * p + std::min(p, m % P);
                    uint64_t em = (m / P) * (p + 1) + std::min(p + 1, m % P);
                    if (beta == DataType{0.0})
                    {
                        for (uint64_t i = sm; i < em; i++)
                        {
                            y[i * incy] = DataType{0.0};
                        }
                    }
                    else
                    {
                        for (uint64_t i = sm; i < em; i++)
                        {
                            y[i * incy] *= beta;
                        }
                    }
                    for (uint64_t j = 0; j < n; j++)
                    {
                        for (uint64_t i = sm; i < em; i++)
                        {
                            y[i * incy] += alpha * std::conj(A[i + ldA * j]) * x[j * incx];
                        }
                    }
                }
            }
        }


    // Solution of the triangular system supplied in packed form
    template <class DataType>
        void tpsv(const char &UpperLower,
                  const char &transpose,
                  const char &diagonal,
                  const uint64_t &n,
                  const DataType *A,
                  DataType *x, const int64_t &incx)
        {
            bool unit_diagonal = (diagonal == 'U' || diagonal == 'u');
            if (UpperLower == 'u' || UpperLower == 'U')
            {
                if (transpose == 'n' || transpose == 'N')
                {
                    A += n * (n + 1) / 2;
                    DataType *ix = x + n * incx;
                    for (uint64_t i = n; i > 0; i--)
                    {
                        ix -= incx;
                        A--;
                        if (!unit_diagonal)
                        {
                            *ix /= *A;
                        }
                        DataType *jx = ix;
                        for (uint64_t j = i - 1; j > 0; j--)
                        {
                            A--;
                            jx -= incx;
                            *jx -= *A * *ix;
                        }
                    }
                }
                else if(transpose == 't' || transpose == 'T')
                {
                    DataType *ix = x;
                    for (uint64_t i = 0; i < n; i++)
                    {
                        const DataType *jx = x;
                        for (uint64_t j = 0; j < i; j++)
                        {
                            *ix -= *A * *jx;
                            A++;
                            jx += incx;
                        }
                        if (!unit_diagonal)
                        {
                            *ix /= *A;
                        }
                        A++;
                        ix += incx;
                    }
                }
                else
                {
                    DataType *ix = x;
                    for (uint64_t i = 0; i < n; i++)
                    {
                        const DataType *jx = x;
                        for (uint64_t j = 0; j < i; j++)
                        {
                            *ix -= std::conj(*A) * *jx;
                            A++;
                            jx += incx;
                        }
                        if (!unit_diagonal)
                        {
                            *ix /= std::conj(*A);
                        }
                        A++;
                        ix += incx;
                    }
                }
            }
            else
            {
                if (transpose == 'n' || transpose == 'N')
                {
                    DataType *ix = x;
                    for (uint64_t i = 0; i < n; i++)
                    {
                        if (!unit_diagonal)
                        {
                            *ix /= *A;
                        }
                        A++;
                        DataType *jx = ix + incx;
                        for (uint64_t j = i + 1; j < n; j++)
                        {
                            *jx -= *A * *ix;
                            A++;
                            jx += incx;
                        }
                    }
                }
                else if(transpose == 't' || transpose == 'T')
                {
                    A += n * (n + 1) / 2;
                    DataType *ix = x + n * incx;
                    for (uint64_t i = n; i > 0; i--)
                    {
                        ix -= incx;
                        const DataType *jx = x + n * incx;
                        for (uint64_t j = n; j > i; j--)
                        {
                            A--;
                            jx -= incx;
                            *ix -= *A * *jx;
                        }
                        A--;
                        if (!unit_diagonal)
                        {
                            *ix /= *A;
                        }
                    }
                }
                else if(transpose == 'c' || transpose == 'C')
                {
                    A += n * (n + 1) / 2;
                    DataType *ix = x + n * incx;
                    for (uint64_t i = n; i > 0; i--)
                    {
                        ix -= incx;
                        const DataType *jx = x + n * incx;
                        for (uint64_t j = n; j > i; j--)
                        {
                            A--;
                            jx -= incx;
                            *ix -= std::conj(*A) * *jx;
                        }
                        A--;
                        if (!unit_diagonal)
                        {
                            *ix /= std::conj(*A);
                        }
                    }
                }
            }
        }


    //Level 3: matrix-matrix operations

    // Matrix by matrix multiplication for general matrices
    template <class DataType>
        void gemm(const char &TransA,
                  const char &TransB,
                  const uint64_t &m, const uint64_t &n,
                  const uint64_t &k,
                  const DataType &alpha,
                  const DataType *A, const uint64_t &ldA,
                  const DataType *B, const uint64_t &ldB,
                  const DataType &beta,
                  DataType *C, const uint64_t &ldC)
        {
#pragma omp parallel for collapse(2)
            for (uint64_t j = 0; j < n; j++)
            {
                for (uint64_t i = 0; i < m; i++)
                {
                    C[i + j * ldC] *= beta;
                }
            }
            if ((TransA == 'N' || TransA == 'n') && (TransB == 'N' || TransB == 'n'))
            {
#pragma omp parallel for
                for (uint64_t j = 0; j < n; j++)
                {
                    for (uint64_t l = 0; l < k; l++)
                    {
                        for (uint64_t i = 0; i < m; i++)
                        {
                            C[i + j * ldC] += A[i + l * ldA] * B[l + j * ldB];
                        }
                    }
                }
            }
            else if ((TransA == 'T' || TransA == 't') && (TransB == 'N' || TransB == 'n'))
            {
#pragma omp parallel for collapse(2)
                for (uint64_t j = 0; j < n; j++)
                {
                    for (uint64_t i = 0; i < m; i++)
                    {
                        for (uint64_t l = 0; l < k; l++)
                        {
                            C[i + j * ldC] += A[i * ldA + l] * B[l + j * ldB];
                        }
                    }
                }
            }
            else if ((TransA == 'C' || TransA == 't') && (TransB == 'N' || TransB == 'n'))
            {
#pragma omp parallel for collapse(2)
                for (uint64_t j = 0; j < n; j++)
                {
                    for (uint64_t i = 0; i < m; i++)
                    {
                        for (uint64_t l = 0; l < k; l++)
                        {
                            C[i + j * ldC] += std::conj(A[i * ldA + l]) * B[l + j * ldB];
                        }
                    }
                }
            }
            else if ((TransA == 'N' || TransA == 'N') && (TransB == 'T' || TransB == 't'))
            {
#pragma omp parallel for collapse(2)
                for (uint64_t j = 0; j < n; j++)
                {
                    for (uint64_t l = 0; l < k; l++)
                    {
                        for (uint64_t i = 0; i < m; i++)
                        {
                            C[i + j * ldC] += A[i + l * ldA] * B[l * ldB + j];
                        }
                    }
                }
            }
            else if ((TransA == 'T' || TransA == 't') && (TransB == 'T' || TransB == 't'))
            {
#pragma omp parallel for collapse(2)
                for (uint64_t j = 0; j < n; j++)
                {
                    for (uint64_t i = 0; i < m; i++)
                    {
                        for (uint64_t l = 0; l < k; l++)
                        {
                            C[i + j * ldC] += A[i * ldA + l] * B[l * ldB + j];
                        }
                    }
                }
            }
            else if ((TransA == 'C' || TransA == 'c') && (TransB == 'T' || TransB == 't'))
            {
#pragma omp parallel for collapse(2)
                for (uint64_t j = 0; j < n; j++)
                {
                    for (uint64_t i = 0; i < m; i++)
                    {
                        for (uint64_t l = 0; l < k; l++)
                        {
                            C[i + j * ldC] += std::conj(A[i * ldA + l]) * B[l * ldB + j];
                        }
                    }
                }
            }
            else if ((TransA == 'N' || TransA == 'N') && (TransB == 'C' || TransB == 'c'))
            {
#pragma omp parallel for
                for (uint64_t j = 0; j < n; j++)
                {
                    for (uint64_t l = 0; l < k; l++)
                    {
                        for (uint64_t i = 0; i < m; i++)
                        {
                            C[i + j * ldC] += A[i + l * ldA] * std::conj(B[l * ldB + j]);
                        }
                    }
                }
            }
            else if ((TransA == 'T' || TransA == 't') && (TransB == 'C' || TransB == 'c'))
            {
#pragma omp parallel for collapse(2)
                for (uint64_t j = 0; j < n; j++)
                {
                    for (uint64_t i = 0; i < m; i++)
                    {
                        for (uint64_t l = 0; l < k; l++)
                        {
                            C[i + j * ldC] += A[i * ldA + l] * std::conj(B[l * ldB + j]);
                        }
                    }
                }
            }
            else if ((TransA == 'C' || TransA == 'c') && (TransB == 'C' || TransB == 'c'))
            {
#pragma omp parallel for collapse(2)
                for (uint64_t j = 0; j < n; j++)
                {
                    for (uint64_t i = 0; i < m; i++)
                    {
                        for (uint64_t l = 0; l < k; l++)
                        {
                            C[i + j * ldC] += std::conj(A[i * ldA + l]) * std::conj(B[l * ldB + j]);
                        }
                    }
                }
            }
        }

    template <class DataType>
    void trmm(const char &side, const char &uplo, const char &transa, const char &diag,
              const uint64_t &m, const uint64_t &n, const DataType &alpha, const DataType *A, const uint64_t &ldA,
              DataType *B, const uint64_t &dlB);

    template <class DataType>
    void trsm(const char &side, const char &uplo, const char &transa, const char &diag,
              const uint64_t &m, const uint64_t &n, const DataType &alpha, const DataType *A, const uint64_t &ldA,
              DataType *B, const uint64_t &dlB);


    template <class DataType>
    void trsm(const char &side, const char &uplo, const char &transa, const char &diag,
              const uint64_t &m, const uint64_t &n, const DataType &alpha, const DataType *A, const uint64_t &ldA,
              DataType *B, const uint64_t &dlB);



    namespace auxiliary {
    //General dense matrix by vector multiplication:
    //       y = \alpha A x + \beta y
    template <class DataType>
        void conj_no_transp_gemv(const char &transpose,
                const uint64_t &m, const uint64_t &n,
                const DataType &alpha,
                const DataType *A, const uint64_t &ldA,
                const DataType *x, const int64_t &incx,
                const DataType beta,
                DataType *y, const int64_t &incy)
        {
#pragma omp parallel
            {
                if (transpose == 'z' || transpose == 'Z')
                {
                    uint64_t P = omp_get_num_threads();
                    uint64_t p = omp_get_thread_num();
                    uint64_t sm = (m / P) * p + std::min(p, m % P);
                    uint64_t em = (m / P) * (p + 1) + std::min(p + 1, m % P);
                    if (beta == DataType{0.0})
                    {
                        for (uint64_t i = sm; i < em; i++)
                        {
                            y[i * incy] = DataType{0.0};
                        }
                    }
                    else
                    {
                        for (uint64_t i = sm; i < em; i++)
                        {
                            y[i * incy] *= beta;
                        }
                    }
                    for (uint64_t j = 0; j < n; j++)
                    {
                        for (uint64_t i = sm; i < em; i++)
                        {
                            y[i * incy] += alpha * std::conj(A[i + ldA * j]) * x[j * incx];
                        }
                    }
                }
                else 
                {
                    gemv(transpose, m, n, alpha, A, ldA, x, incx, beta, y, incy);
                }
            }
        }

        template<class DataType>
        void conj_no_transp_gemm(const char &TransA,
                  const char &TransB,
                  const uint64_t &m, const uint64_t &n,
                  const uint64_t &k,
                  const DataType &alpha,
                  const DataType *A, const uint64_t &ldA,
                  const DataType *B, const uint64_t &ldB,
                  const DataType &beta,
                  DataType *C, const uint64_t &ldC)
        {
#pragma omp parallel for collapse(2)
            for (uint64_t j = 0; j < n; j++)
            {
                for (uint64_t i = 0; i < m; i++)
                {
                    C[i + j * ldC] *= beta;
                }
            }

            // Non-standard options: 'Trans = Z' means conjugated but not transposed
            if ((TransA == 'Z' || TransA == 'z') && (TransB == 'N' || TransB == 'n'))
            {
#pragma omp parallel for
                for (uint64_t j = 0; j < n; j++)
                {
                    for (uint64_t l = 0; l < k; l++)
                    {
                        for (uint64_t i = 0; i < m; i++)
                        {
                            C[i + j * ldC] += std::conj(A[i + l * ldA]) * B[l + j * ldB];
                        }
                    }
                }
            }
            else if ((TransA == 'Z' || TransA == 'z') && (TransB == 'T' || TransB == 't'))
            {
#pragma omp parallel for collapse(2)
                for (uint64_t j = 0; j < n; j++)
                {
                    for (uint64_t l = 0; l < k; l++)
                    {
                        for (uint64_t i = 0; i < m; i++)
                        {
                            C[i + j * ldC] += std::conj(A[i + l * ldA]) * B[l * ldB + j];
                        }
                    }
                }
            }
            else if ((TransA == 'Z' || TransA == 'z') && (TransB == 'C' || TransB == 'c'))
            {
#pragma omp parallel for
                for (uint64_t j = 0; j < n; j++)
                {
                    for (uint64_t l = 0; l < k; l++)
                    {
                        for (uint64_t i = 0; i < m; i++)
                        {
                            C[i + j * ldC] += std::conj(A[i + l * ldA] * B[l * ldB + j]);
                        }
                    }
                }
            }
            else if ((TransA == 'Z' || TransA == 'z') && (TransB == 'Z' || TransB == 'z'))
            {
#pragma omp parallel for
                for (uint64_t j = 0; j < n; j++)
                {
                    for (uint64_t l = 0; l < k; l++)
                    {
                        for (uint64_t i = 0; i < m; i++)
                        {
                            C[i + j * ldC] += std::conj(A[i + l * ldA] * B[l + j * ldB]);
                        }
                    }
                }
            }
            

            if ((TransA == 'N' || TransA == 'n') && (TransB == 'Z' || TransB == 'z'))
            {
#pragma omp parallel for
                for (uint64_t j = 0; j < n; j++)
                {
                    for (uint64_t l = 0; l < k; l++)
                    {
                        for (uint64_t i = 0; i < m; i++)
                        {
                            C[i + j * ldC] += A[i + l * ldA] * std::conj(B[l + j * ldB]);
                        }
                    }
                }
            }
            else if ((TransA == 'T' || TransA == 't') && (TransB == 'Z' || TransB == 'z'))
            {
#pragma omp parallel for collapse(2)
                for (uint64_t j = 0; j < n; j++)
                {
                    for (uint64_t i = 0; i < m; i++)
                    {
                        for (uint64_t l = 0; l < k; l++)
                        {
                            C[i + j * ldC] += A[i * ldA + l] * std::conj(B[l + j * ldB]);
                        }
                    }
                }
            }
            else if ((TransA == 'C' || TransA == 't') && (TransB == 'Z' || TransB == 'z'))
            {
#pragma omp parallel for collapse(2)
                for (uint64_t j = 0; j < n; j++)
                {
                    for (uint64_t i = 0; i < m; i++)
                    {
                        for (uint64_t l = 0; l < k; l++)
                        {
                            C[i + j * ldC] += std::conj(A[i * ldA + l] * B[l + j * ldB]);
                        }
                    }
                }
            }
            else
            {
                gemm(TransA, TransB, m, n, k, alpha, A, ldA, B, ldB, beta, C, ldC);
            }
        }
    }

    //Implementations for standard types use the BLAS library

    // swap
    template <>
        void swap(const uint64_t &n,
                  const float *x, const int64_t &incx,
                  float *y, const int64_t &incy);
    template <>
        void swap(const uint64_t &n,
                  const double *x, const int64_t &incx,
                  double *y, const int64_t &incy);
    template <>
        void swap(const uint64_t &n,
                  const std::complex<float> *x, const int64_t &incx,
                  std::complex<float> *y, const int64_t &incy);
    template <>
        void swap(const uint64_t &n,
                  const std::complex<double> *x, const int64_t &incx,
                  std::complex<double> *y, const int64_t &incy);

    // axpy
    template <>
        void axpy(const uint64_t &n,
                  const float &alpha,
                  const float *x, const int64_t &incx,
                  float *y, const int64_t &incy);
    template <>
        void axpy(const uint64_t &n,
                  const double &alpha,
                  const double *x, const int64_t &incx,
                  double *y, const int64_t &incy);
    template <>
        void axpy(const uint64_t &n,
                  const std::complex<float> &alpha,
                  const std::complex<float> *x, const int64_t &incx,
                  std::complex<float> *y, const int64_t &incy);
    template <>
        void axpy(const uint64_t &n,
                  const std::complex<double> &alpha,
                  const std::complex<double> *x, const int64_t &incx,
                  std::complex<double> *y, const int64_t &incy);

    // copy
    template <>
        void copy(const uint64_t &n,
                  const float *x, const int64_t &incx,
                  float *y, const int64_t &incy);
    template <>
        void copy(const uint64_t &n,
                  const double *x, const int64_t &incx,
                  double *y, const int64_t &incy);
    template <>
        void copy(const uint64_t &n,
                  const std::complex<float> *x, const int64_t &incx,
                  std::complex<float> *y, const int64_t &incy);
    template <>
        void copy(const uint64_t &n,
                  const std::complex<double> *x, const int64_t &incx,
                  std::complex<double> *y, const int64_t &incy);

    // dotc
    template <>
        float dotc(const uint64_t &n,
                   const float *x, const int64_t &incx,
                   const float *y, const int64_t &incy);
    template <>
        double dotc(const uint64_t &n,
                    const double *x, const int64_t &incx,
                    const double *y, const int64_t &incy);
    template <>
        std::complex<float> dotc(const uint64_t &n,
                                 const std::complex<float> *x, const int64_t &incx,
                                 const std::complex<float> *y, const int64_t &incy);
    template <>
        std::complex<double> dotc(const uint64_t &n,
                                  const std::complex<double> *x, const int64_t &incx,
                                  const std::complex<double> *y, const int64_t &incy);

    // dotc
    template <>
    float dotu(const uint64_t &n,
               const float *x, const int64_t &incx,
               const float *y, const int64_t &incy);
    template <>
    double dotu(const uint64_t &n,
                const double *x, const int64_t &incx,
                const double *y, const int64_t &incy);
    template <>
    std::complex<float> dotu(const uint64_t &n,
                             const std::complex<float> *x, const int64_t &incx,
                             const std::complex<float> *y, const int64_t &incy);
    template <>
    std::complex<double> dotu(const uint64_t &n,
                              const std::complex<double> *x, const int64_t &incx,
                              const std::complex<double> *y, const int64_t &incy);

    // nrm2
    template <>
        float nrm2(const uint64_t &n,
                   const float *x, const int64_t &incx);

    template <>
        double nrm2(const uint64_t &n,
                    const double *x, const int64_t &incx);

    template <>
        float nrm2(const uint64_t &n,
                   const std::complex<float> *x, const int64_t &incx);

    template <>
        double nrm2(const uint64_t &n,
                   const std::complex<double> *x, const int64_t &incx);

    // rscal
    template <>
        void rscal(const uint64_t &n,
                  const float &alpha,
                  float *x, const int64_t &incx);
    template <>
        void rscal(const uint64_t &n,
                  const double &alpha,
                  double *x, const int64_t &incx);
    template <>
        void rscal(const uint64_t &n,
                  const float &alpha,
                  std::complex<float> *x, const int64_t &incx);
    template <>
        void rscal(const uint64_t &n,
                  const double &alpha,
                  std::complex<double> *x, const int64_t &incx);

    // scal
    template <>
        void scal(const uint64_t &n,
                  const float &alpha,
                  float *x, const int64_t &incx);
    template <>
        void scal(const uint64_t &n,
                  const double &alpha,
                  double *x, const int64_t &incx);
    template <>
        void scal(const uint64_t &n,
                  const std::complex<float> &alpha,
                  std::complex<float> *x, const int64_t &incx);
    template <>
        void scal(const uint64_t &n,
                  const std::complex<double> &alpha,
                  std::complex<double> *x, const int64_t &incx);

    //rot
    template<>
        void rot(const uint64_t &n,
                 float *x, const uint64_t &incx,
                 float *y, const uint64_t &incy,
                 const float &c, const float &s);

    template<>
        void rot(const uint64_t &n,
                 double *x, const uint64_t &incx,
                 double *y, const uint64_t &incy,
                 const double &c, const double &s);

    template<>
        void rot(const uint64_t &n,
                 std::complex<float> *x, const uint64_t &incx,
                 std::complex<float> *y, const uint64_t &incy,
                 const float &c, const std::complex<float> &s);

    template<>
        void rot(const uint64_t &n,
                 std::complex<double> *x, const uint64_t &incx,
                 std::complex<double> *y, const uint64_t &incy,
                 const double &c, const std::complex<double> &s);

    // gerc
    template <>
        void gerc(const uint64_t &m, const uint64_t &n,
                const float &alpha,
                const float *x, const int64_t &incx,
                const float *y, const int64_t &incy,
                float *A, const uint64_t &ldA);
    template <>
        void gerc(const uint64_t &m, const uint64_t &n,
                const double &alpha,
                const double *x, const int64_t &incx,
                const double *y, const int64_t &incy,
                double *A, const uint64_t &ldA);
    template <>
        void gerc(const uint64_t &m, const uint64_t &n,
                const std::complex<float> &alpha,
                const std::complex<float> *x, const int64_t &incx,
                const std::complex<float> *y, const int64_t &incy,
                std::complex<float> *A, const uint64_t &ldA);
    template <>
        void gerc(const uint64_t &m, const uint64_t &n,
                const std::complex<double> &alpha,
                const std::complex<double> *x, const int64_t &incx,
                const std::complex<double> *y, const int64_t &incy,
                std::complex<double> *A, const uint64_t &ldA);

    // geru
    template <>
    void geru(const uint64_t &m, const uint64_t &n,
              const float &alpha,
              const float *x, const int64_t &incx,
              const float *y, const int64_t &incy,
              float *A, const uint64_t &ldA);
    template <>
    void geru(const uint64_t &m, const uint64_t &n,
              const double &alpha,
              const double *x, const int64_t &incx,
              const double *y, const int64_t &incy,
              double *A, const uint64_t &ldA);
    template <>
    void geru(const uint64_t &m, const uint64_t &n,
              const std::complex<float> &alpha,
              const std::complex<float> *x, const int64_t &incx,
              const std::complex<float> *y, const int64_t &incy,
              std::complex<float> *A, const uint64_t &ldA);
    template <>
    void geru(const uint64_t &m, const uint64_t &n,
              const std::complex<double> &alpha,
              const std::complex<double> *x, const int64_t &incx,
              const std::complex<double> *y, const int64_t &incy,
              std::complex<double> *A, const uint64_t &ldA);

    // gemv
    template <>
        void gemv(const char &transpose,
                const uint64_t &m, const uint64_t &n,
                const float &alpha,
                const float *A, const uint64_t &ldA,
                const float *x, const int64_t &incx,
                const float beta,
                float *y, const int64_t &incy);
    template <>
        void gemv(const char &transpose,
                const uint64_t &m, const uint64_t &n,
                const double &alpha,
                const double *A, const uint64_t &ldA,
                const double *x, const int64_t &incx,
                const double beta,
                double *y, const int64_t &incy);
    template <>
        void gemv(const char &transpose,
                const uint64_t &m, const uint64_t &n,
                const std::complex<float> &alpha,
                const std::complex<float> *A, const uint64_t &ldA,
                const std::complex<float> *x, const int64_t &incx,
                const std::complex<float> beta,
                std::complex<float> *y, const int64_t &incy);
    template <>
        void gemv(const char &transpose,
                const uint64_t &m, const uint64_t &n,
                const std::complex<double> &alpha,
                const std::complex<double> *A, const uint64_t &ldA,
                const std::complex<double> *x, const int64_t &incx,
                const std::complex<double> beta,
                std::complex<double> *y, const int64_t &incy);

    // tpsv
    template <>
        void tpsv(const char &UpperLower,
                  const char &transpose,
                  const char &diagonal,
                  const uint64_t &n,
                  const float *A,
                  float *x, const int64_t &incx);

    template <>
        void tpsv(const char &UpperLower,
                  const char &transpose,
                  const char &diagonal,
                  const uint64_t &n,
                  const double *A,
                  double *x, const int64_t &incx);

    template <>
        void tpsv(const char &UpperLower,
                  const char &transpose,
                  const char &diagonal,
                  const uint64_t &n,
                  const std::complex<float> *A,
                  std::complex<float> *x, const int64_t &incx);

    template <>
        void tpsv(const char &UpperLower,
                  const char &transpose,
                  const char &diagonal,
                  const uint64_t &n,
                  const std::complex<double> *A,
                  std::complex<double> *x, const int64_t &incx);

    template <>
        void gemm(const char &TransA,
                  const char &TransB,
                  const uint64_t &m, const uint64_t &n,
                  const uint64_t &k,
                  const float &alpha,
                  const float *A, const uint64_t &ldA,
                  const float *B, const uint64_t &ldB,
                  const float &beta,
                  float *C, const uint64_t &ldC);
    template <>
        void gemm(const char &TransA,
                  const char &TransB,
                  const uint64_t &m, const uint64_t &n,
                  const uint64_t &k,
                  const double &alpha,
                  const double *A, const uint64_t &ldA,
                  const double *B, const uint64_t &ldB,
                  const double &beta,
                  double *C, const uint64_t &ldC);
    template <>
        void gemm(const char &TransA,
                  const char &TransB,
                  const uint64_t &m, const uint64_t &n,
                  const uint64_t &k,
                  const std::complex<float> &alpha,
                  const std::complex<float> *A, const uint64_t &ldA,
                  const std::complex<float> *B, const uint64_t &ldB,
                  const std::complex<float> &beta,
                  std::complex<float> *C, const uint64_t &ldC);
    template <>
        void gemm(const char &TransA,
                  const char &TransB,
                  const uint64_t &m, const uint64_t &n,
                  const uint64_t &k,
                  const std::complex<double> &alpha,
                  const std::complex<double> *A, const uint64_t &ldA,
                  const std::complex<double> *B, const uint64_t &ldB,
                  const std::complex<double> &beta,
                  std::complex<double> *C, const uint64_t &ldC);


    template <>
    void trmm(const char &side, const char &uplo, const char &transa, const char &diag,
              const uint64_t &m, const uint64_t &n, const float &alpha, const float *A, const uint64_t &ldA,
              float *B, const uint64_t &LDB);

    template <>
    void trmm(const char &side, const char &uplo, const char &transa, const char &diag,
              const uint64_t &m, const uint64_t &n, const double &alpha, const double *A, const uint64_t &ldA,
              double *B, const uint64_t &dlB);

    template <>
    void trmm(const char &side, const char &uplo, const char &transa, const char &diag,
              const uint64_t &m, const uint64_t &n, const std::complex<float> &alpha, const std::complex<float> *A,
              const uint64_t &ldA, std::complex<float> *B, const uint64_t &ldB);

    template <>
    void trmm(const char &side, const char &uplo, const char &transa, const char &diag,
              const uint64_t &m, const uint64_t &n, const std::complex<double> &alpha, const std::complex<double> *A,
              const uint64_t &ldA, std::complex<double> *B, const uint64_t &ldB);

    template <>
    void trsm(const char &side, const char &uplo, const char &transa, const char &diag,
              const uint64_t &m, const uint64_t &n, const float &alpha, const float *A, const uint64_t &ldA,
              float *B, const uint64_t &LDB);

    template <>
    void trsm(const char &side, const char &uplo, const char &transa, const char &diag,
              const uint64_t &m, const uint64_t &n, const double &alpha, const double *A, const uint64_t &ldA,
              double *B, const uint64_t &dlB);

    template <>
    void trsm(const char &side, const char &uplo, const char &transa, const char &diag,
              const uint64_t &m, const uint64_t &n, const std::complex<float> &alpha, const std::complex<float> *A,
              const uint64_t &ldA, std::complex<float> *B, const uint64_t &ldB);

    template <>
    void trsm(const char &side, const char &uplo, const char &transa, const char &diag,
              const uint64_t &m, const uint64_t &n, const std::complex<double> &alpha, const std::complex<double> *A,
              const uint64_t &ldA, std::complex<double> *B, const uint64_t &ldB);



    template <>
    void trsm(const char &side, const char &uplo, const char &transa, const char &diag,
              const uint64_t &m, const uint64_t &n, const float &alpha, const float *A, const uint64_t &ldA,
              float *B, const uint64_t &LDB);

    template <>
    void trsm(const char &side, const char &uplo, const char &transa, const char &diag,
              const uint64_t &m, const uint64_t &n, const double &alpha, const double *A, const uint64_t &ldA,
              double *B, const uint64_t &dlB);

    template <>
    void trsm(const char &side, const char &uplo, const char &transa, const char &diag,
              const uint64_t &m, const uint64_t &n, const std::complex<float> &alpha, const std::complex<float> *A,
              const uint64_t &ldA, std::complex<float> *B, const uint64_t &ldB);

    template <>
    void trsm(const char &side, const char &uplo, const char &transa, const char &diag,
              const uint64_t &m, const uint64_t &n, const std::complex<double> &alpha, const std::complex<double> *A,
              const uint64_t &ldA, std::complex<double> *B, const uint64_t &ldB);


}

#endif
