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


#include "BLAS.h"

#include "cpp_blas/BLAS.hpp"


namespace BLAS
{
    // swap
    template <>
        void swap(const uint64_t &n,
                  const float *x, const int64_t &incx,
                  float *y, const int64_t &incy)
        {
            sswap(n, x, incx, y, incy);
        }

    template <>
        void swap(const uint64_t &n,
                  const double *x, const int64_t &incx,
                  double *y, const int64_t &incy)
        {
            dswap(n, x, incx, y, incy);
        }

    template <>
        void swap(const uint64_t &n,
                  const std::complex<float> *x, const int64_t &incx,
                  std::complex<float> *y, const int64_t &incy)
        {
            cswap(n, x, incx, y, incy);
        }

    template <>
        void swap(const uint64_t &n,
                  const std::complex<double> *x, const int64_t &incx,
                  std::complex<double> *y, const int64_t &incy)
        {
            zswap(n, x, incx, y, incy);
        }

    // axpy: y = alpha * x + y
    template <>
        void axpy(const uint64_t &n,
                  const float &alpha,
                  const float *x, const int64_t &incx,
                  float *y, const int64_t &incy)
        {
            saxpy(n, alpha, x, incx, y, incy);
        }

    template <>
        void axpy(const uint64_t &n,
                  const double &alpha,
                  const double *x, const int64_t &incx,
                  double *y, const int64_t &incy)
        {
            daxpy(n, alpha, x, incx, y, incy);
        }

    template <>
        void axpy(const uint64_t &n,
                  const std::complex<float> &alpha,
                  const std::complex<float> *x, const int64_t &incx,
                  std::complex<float> *y, const int64_t &incy)
        {
            caxpy(n, alpha, x, incx, y, incy);
        }

    template <>
        void axpy(const uint64_t &n,
                  const std::complex<double> &alpha,
                  const std::complex<double> *x, const int64_t &incx,
                  std::complex<double> *y, const int64_t &incy)
        {
            zaxpy(n, alpha, x, incx, y, incy);
        }

    // copy
    template <>
        void copy(const uint64_t &n,
                  const float *x, const int64_t &incx,
                  float *y, const int64_t &incy)
        {
            scopy(n, x, incx, y, incy);
        }

    template <>
        void copy(const uint64_t &n,
                  const double *x, const int64_t &incx,
                  double *y, const int64_t &incy)
        {
            dcopy(n, x, incx, y, incy);
        }

    template <>
        void copy(const uint64_t &n,
                  const std::complex<float> *x, const int64_t &incx,
                  std::complex<float> *y, const int64_t &incy)
        {
            ccopy(n, x, incx, y, incy);
        }

    template <>
        void copy(const uint64_t &n,
                  const std::complex<double> *x, const int64_t &incx,
                  std::complex<double> *y, const int64_t &incy)
        {
            zcopy(n, x, incx, y, incy);
        }

    // dotc
    template <>
        float dotc(const uint64_t &n,
                   const float *x, const int64_t &incx,
                   const float *y, const int64_t &incy)
        {
            return sdot(n, x, incx, y, incy);
        }

    template <>
        double dotc(const uint64_t &n,
                    const double *x, const int64_t &incx,
                    const double *y, const int64_t &incy)
        {
            return ddot(n, x, incx, y, incy);
        }

    template <>
        std::complex<float> dotc(const uint64_t &n,
                                 const std::complex<float> *x, const int64_t &incx,
                                 const std::complex<float> *y, const int64_t &incy)
        {
            std::complex<float> r;
            cdotc_stackreturn(r, n, x, incx, y, incy);
            return r;
        }

    template <>
        std::complex<double> dotc(const uint64_t &n,
                                  const std::complex<double> *x, const int64_t &incx,
                                  const std::complex<double> *y, const int64_t &incy)
        {
            std::complex<double> r;
            zdotc_stackreturn(r, n, x, incx, y, incy);
            return r;
        }

// dotu
    template <>
    float dotu(const uint64_t &n,
               const float *x, const int64_t &incx,
               const float *y, const int64_t &incy)
    {
        return sdot(n, x, incx, y, incy);
    }

    template <>
    double dotu(const uint64_t &n,
                const double *x, const int64_t &incx,
                const double *y, const int64_t &incy)
    {
        return ddot(n, x, incx, y, incy);
    }

    template <>
    std::complex<float> dotu(const uint64_t &n,
                             const std::complex<float> *x, const int64_t &incx,
                             const std::complex<float> *y, const int64_t &incy)
    {
        std::complex<float> r;
        cdotu_stackreturn(r, n, x, incx, y, incy);
        return r;
    }

    template <>
    std::complex<double> dotu(const uint64_t &n,
                              const std::complex<double> *x, const int64_t &incx,
                              const std::complex<double> *y, const int64_t &incy)
    {
        std::complex<double> r;
        zdotu_stackreturn(r, n, x, incx, y, incy);
        return r;
    }


    // nrm2
    template <>
        float nrm2(const uint64_t &n,
                   const float *x, const int64_t &incx)
        {
            return snrm2(n, x, incx);
        }

    template <>
        double nrm2(const uint64_t &n,
                    const double *x, const int64_t &incx)
        {
            return dnrm2(n, x, incx);
        }

    template <>
        float nrm2(const uint64_t &n,
                   const std::complex<float> *x, const int64_t &incx)
        {
            return scnrm2(n, x, incx);
        }

    template <>
        double nrm2(const uint64_t &n,
                   const std::complex<double> *x, const int64_t &incx)
        {
            return dznrm2(n, x, incx);
        }

    // rscal
    template <>
        void rscal(const uint64_t &n,
                  const float &alpha,
                  float *x, const int64_t &incx)
        {
            sscal(n, alpha, x, incx);
        }
    template <>
        void rscal(const uint64_t &n,
                  const double &alpha,
                  double *x, const int64_t &incx)
        {
            dscal(n, alpha, x, incx);
        }

    template <>
        void rscal(const uint64_t &n,
                  const float &alpha,
                  std::complex<float> *x, const int64_t &incx)
        {
            csscal(n, alpha, x, incx);
        }
        
    template <>
        void rscal(const uint64_t &n,
                  const double &alpha,
                  std::complex<double> *x, const int64_t &incx)
        {
            zdscal(n, alpha, x, incx);
        }

    // scal
    template <>
        void scal(const uint64_t &n,
                  const float &alpha,
                  float *x, const int64_t &incx)
        {
            sscal(n, alpha, x, incx);
        }
    template <>
        void scal(const uint64_t &n,
                  const double &alpha,
                  double *x, const int64_t &incx)
        {
            dscal(n, alpha, x, incx);
        }

    template <>
        void scal(const uint64_t &n,
                  const std::complex<float> &alpha,
                  std::complex<float> *x, const int64_t &incx)
        {
            cscal(n, alpha, x, incx);
        }
        
    template <>
        void scal(const uint64_t &n,
                  const std::complex<double> &alpha,
                  std::complex<double> *x, const int64_t &incx)
        {
            zscal(n, alpha, x, incx);
        }

    // rot
    template<>
        void rot(const uint64_t &n,
                  float *x, const uint64_t &incx,
                  float *y, const uint64_t &incy,
                  const float &c, const float &s)
        {
            srot(n, x, incx, y, incy, c, s);
        }

    template<>
        void rot(const uint64_t &n,
                  double *x, const uint64_t &incx,
                  double *y, const uint64_t &incy,
                  const double &c, const double &s)
        {
            drot(n, x, incx, y, incy, c, s);
        }

    template<>
        void rot(const uint64_t &n,
                  std::complex<float> *x, const uint64_t &incx,
                  std::complex<float> *y, const uint64_t &incy,
                  const float &c, const std::complex<float> &s)
        {
            crot(n, x, incx, y, incy, c, s);
        }

    template<>
        void rot(const uint64_t &n,
                  std::complex<double> *x, const uint64_t &incx,
                  std::complex<double> *y, const uint64_t &incy,
                  const double &c, const std::complex<double> &s)
        {
            zrot(n, x, incx, y, incy, c, s);
        }

    // gerc
    template <>
        void gerc(const uint64_t &m, const uint64_t &n,
                const float &alpha,
                const float *x, const int64_t &incx,
                const float *y, const int64_t &incy,
                float *A, const uint64_t &ldA)
        {
            sger(m, n, alpha, x, incx, y, incy, A, ldA);
        }
    template <>
        void gerc(const uint64_t &m, const uint64_t &n,
                const double &alpha,
                const double *x, const int64_t &incx,
                const double *y, const int64_t &incy,
                double *A, const uint64_t &ldA)
        {
            dger(m, n, alpha, x, incx, y, incy, A, ldA);
        }
    template <>
        void gerc(const uint64_t &m, const uint64_t &n,
                const std::complex<float> &alpha,
                const std::complex<float> *x, const int64_t &incx,
                const std::complex<float> *y, const int64_t &incy,
                std::complex<float> *A, const uint64_t &ldA)
        {
            cgerc(m, n, alpha, x, incx, y, incy, A, ldA);
        }
    template <>
        void gerc(const uint64_t &m, const uint64_t &n,
                const std::complex<double> &alpha,
                const std::complex<double> *x, const int64_t &incx,
                const std::complex<double> *y, const int64_t &incy,
                std::complex<double> *A, const uint64_t &ldA)
        {
            zgerc(m, n, alpha, x, incx, y, incy, A, ldA);
        }

    // geru
    template <>
    void geru(const uint64_t &m, const uint64_t &n,
              const float &alpha,
              const float *x, const int64_t &incx,
              const float *y, const int64_t &incy,
              float *A, const uint64_t &ldA)
    {
        sger(m, n, alpha, x, incx, y, incy, A, ldA);
    }
    template <>
    void geru(const uint64_t &m, const uint64_t &n,
              const double &alpha,
              const double *x, const int64_t &incx,
              const double *y, const int64_t &incy,
              double *A, const uint64_t &ldA)
    {
        dger(m, n, alpha, x, incx, y, incy, A, ldA);
    }
    template <>
    void geru(const uint64_t &m, const uint64_t &n,
              const std::complex<float> &alpha,
              const std::complex<float> *x, const int64_t &incx,
              const std::complex<float> *y, const int64_t &incy,
              std::complex<float> *A, const uint64_t &ldA)
    {
        cgeru(m, n, alpha, x, incx, y, incy, A, ldA);
    }
    template <>
    void geru(const uint64_t &m, const uint64_t &n,
              const std::complex<double> &alpha,
              const std::complex<double> *x, const int64_t &incx,
              const std::complex<double> *y, const int64_t &incy,
              std::complex<double> *A, const uint64_t &ldA)
    {
        zgeru(m, n, alpha, x, incx, y, incy, A, ldA);
    }

    // gemv
    template<>
    void gemv(const char &transpose,
            const uint64_t &m, const uint64_t &n,
            const float &alpha,
            const float *A, const uint64_t &ldA,
            const float *x, const int64_t &incx,
            const float beta,
            float *y, const int64_t &incy)
    {
        char trans{transpose};
        if (transpose == 'Z' || transpose == 'z') {
            trans = 'N';
        }
        sgemv(trans, m, n, alpha, A, ldA, x, incx, beta, y, incy);
    }

    template<>
    void gemv(const char &transpose,
            const uint64_t &m, const uint64_t &n,
            const double &alpha,
            const double *A, const uint64_t &ldA,
            const double *x, const int64_t &incx,
            const double beta,
            double *y, const int64_t &incy)
    {
        char trans{transpose};
        if (transpose == 'Z' || transpose == 'z') {
            trans = 'N';
        }
        dgemv(trans, m, n, alpha, A, ldA, x, incx, beta, y, incy);
    }

    template<>
        void gemv(const char &transpose,
                const uint64_t &m, const uint64_t &n,
                const std::complex<float> &alpha,
                const std::complex<float> *A, const uint64_t &ldA,
                const std::complex<float> *x, const int64_t &incx,
                const std::complex<float> beta,
                std::complex<float> *y, const int64_t &incy)
        {
            if (transpose == 'Z' || transpose == 'z') {
                auxiliary::conj_no_transp_gemv(transpose, m, n, alpha, A, ldA, x, incx, beta, y, incy);
            } else {
                cgemv(transpose, m, n, alpha, A, ldA, x, incx, beta, y, incy);
            }
        }

    template<>
    void gemv(const char &transpose,
            const uint64_t &m, const uint64_t &n,
            const std::complex<double> &alpha,
            const std::complex<double> *A, const uint64_t &ldA,
            const std::complex<double> *x, const int64_t &incx,
            const std::complex<double> beta,
            std::complex<double> *y, const int64_t &incy)
    {
        if (transpose == 'Z' || transpose == 'z') {
            auxiliary::conj_no_transp_gemv(transpose, m, n, alpha, A, ldA, x, incx, beta, y, incy);
        } else {
            zgemv(transpose, m, n, alpha, A, ldA, x, incx, beta, y, incy);
        }
    }

    // tpsv
    template <>
        void tpsv(const char &UpperLower,
                  const char &transpose,
                  const char &diagonal,
                  const uint64_t &n,
                  const float *A,
                  float *x, const int64_t &incx)
        {
            stpsv(UpperLower, transpose, diagonal, n, A, x, incx);
        }

    template <>
        void tpsv(const char &UpperLower,
                  const char &transpose,
                  const char &diagonal,
                  const uint64_t &n,
                  const double *A,
                  double *x, const int64_t &incx)
        {
            dtpsv(UpperLower, transpose, diagonal, n, A, x, incx);
        }

    template <>
        void tpsv(const char &UpperLower,
                  const char &transpose,
                  const char &diagonal,
                  const uint64_t &n,
                  const std::complex<float> *A,
                  std::complex<float> *x, const int64_t &incx)
        {
            ctpsv(UpperLower, transpose, diagonal, n, A, x, incx);
        }

    template <>
        void tpsv(const char &UpperLower,
                  const char &transpose,
                  const char &diagonal,
                  const uint64_t &n,
                  const std::complex<double> *A,
                  std::complex<double> *x, const int64_t &incx)
        {
            ztpsv(UpperLower, transpose, diagonal, n, A, x, incx);
        }

    template <>
        void gemm(const char &TransA,
                  const char &TransB,
                  const uint64_t &m, const uint64_t &n,
                  const uint64_t &k,
                  const float &alpha,
                  const float *A, const uint64_t &ldA,
                  const float *B, const uint64_t &ldB,
                  const float &beta,
                  float *C, const uint64_t &ldC)
        {
            char ZTransA{TransA};
            if (TransA == 'Z' || TransA == 'z') {
                ZTransA = 'N';
            }
            char ZTransB{TransB};
            if (TransB == 'Z' || TransB == 'z') {
                ZTransB = 'N';
            }
            sgemm(ZTransA, ZTransB, m, n, k, alpha, A, ldA, B, ldB, beta, C, ldC);
        }
    template <>
        void gemm(const char &TransA,
                  const char &TransB,
                  const uint64_t &m, const uint64_t &n,
                  const uint64_t &k,
                  const double &alpha,
                  const double *A, const uint64_t &ldA,
                  const double *B, const uint64_t &ldB,
                  const double &beta,
                  double *C, const uint64_t &ldC)
        {
            char ZTransA{TransA};
            if (TransA == 'Z' || TransA == 'z') {
                ZTransA = 'N';
            }
            char ZTransB{TransB};
            if (TransB == 'Z' || TransB == 'z') {
                ZTransB = 'N';
            }
            dgemm(ZTransA, ZTransB, m, n, k, alpha, A, ldA, B, ldB, beta, C, ldC);
        }
    template <>
        void gemm(const char &TransA,
                  const char &TransB,
                  const uint64_t &m, const uint64_t &n,
                  const uint64_t &k,
                  const std::complex<float> &alpha,
                  const std::complex<float> *A, const uint64_t &ldA,
                  const std::complex<float> *B, const uint64_t &ldB,
                  const std::complex<float> &beta,
                  std::complex<float> *C, const uint64_t &ldC)
        {        
            if (TransA == 'Z' || TransA == 'z' || TransB == 'Z' || TransB == 'z') {
                auxiliary::conj_no_transp_gemm(TransA, TransB, m, n, k, alpha, A, ldA, B, ldB, beta, C, ldC);
            } else {
                cgemm(TransA, TransB, m, n, k, alpha, A, ldA, B, ldB, beta, C, ldC);
            }
        }
    template <>
        void gemm(const char &TransA,
                  const char &TransB,
                  const uint64_t &m, const uint64_t &n,
                  const uint64_t &k,
                  const std::complex<double> &alpha,
                  const std::complex<double> *A, const uint64_t &ldA,
                  const std::complex<double> *B, const uint64_t &ldB,
                  const std::complex<double> &beta,
                  std::complex<double> *C, const uint64_t &ldC)
        {
            if (TransA == 'Z' || TransA == 'z' || TransB == 'Z' || TransB == 'z') {
                auxiliary::conj_no_transp_gemm(TransA, TransB, m, n, k, alpha, A, ldA, B, ldB, beta, C, ldC);
            } else {
                zgemm(TransA, TransB, m, n, k, alpha, A, ldA, B, ldB, beta, C, ldC);
            }
        }

    template<>
    void trmm(const char &side, const char &uplo, const char &transa, const char &diag,
              const uint64_t &m, const uint64_t &n, const float &alpha,
              const float *A, const uint64_t &ldA, float *B, const uint64_t &ldB)
    {
        strmm(side, uplo, transa, diag, m, n, alpha, A, ldA, B, ldB);
    }

    template<>
    void trmm(const char &side, const char &uplo, const char &transa, const char &diag,
              const uint64_t &m, const uint64_t &n, const double &alpha,
              const double *A, const uint64_t &ldA, double *B, const uint64_t &ldB)
    {
        dtrmm(side, uplo, transa, diag, m, n, alpha, A, ldA, B, ldB);
    }

    template<>
    void trmm(const char &side, const char &uplo, const char &transa, const char &diag,
              const uint64_t &m, const uint64_t &n, const std::complex<float> &alpha,
              const std::complex<float> *A, const uint64_t &ldA, std::complex<float> *B, const uint64_t &ldB)
    {
        ctrmm(side, uplo, transa, diag, m, n, alpha, A, ldA, B, ldB);
    }

    template<>
    void trmm(const char &side, const char &uplo, const char &transa, const char &diag,
              const uint64_t &m, const uint64_t &n, const std::complex<double> &alpha,
              const std::complex<double> *A, const uint64_t &ldA, std::complex<double> *B, const uint64_t &ldB)
    {
        ztrmm(side, uplo, transa, diag, m, n, alpha, A, ldA, B, ldB);
    }


    template<>
    void trsm(const char &side, const char &uplo, const char &transa, const char &diag,
              const uint64_t &m, const uint64_t &n, const float &alpha,
              const float *A, const uint64_t &ldA, float *B, const uint64_t &ldB)
    {
        strsm(side, uplo, transa, diag, m, n, alpha, A, ldA, B, ldB);
    }

    template<>
    void trsm(const char &side, const char &uplo, const char &transa, const char &diag,
              const uint64_t &m, const uint64_t &n, const double &alpha,
              const double *A, const uint64_t &ldA, double *B, const uint64_t &ldB)
    {
        dtrsm(side, uplo, transa, diag, m, n, alpha, A, ldA, B, ldB);
    }

    template<>
    void trsm(const char &side, const char &uplo, const char &transa, const char &diag,
              const uint64_t &m, const uint64_t &n, const std::complex<float> &alpha,
              const std::complex<float> *A, const uint64_t &ldA, std::complex<float> *B, const uint64_t &ldB)
    {
        ctrsm(side, uplo, transa, diag, m, n, alpha, A, ldA, B, ldB);
    }

    template<>
    void trsm(const char &side, const char &uplo, const char &transa, const char &diag,
              const uint64_t &m, const uint64_t &n, const std::complex<double> &alpha,
              const std::complex<double> *A, const uint64_t &ldA, std::complex<double> *B, const uint64_t &ldB)
    {
        ztrsm(side, uplo, transa, diag, m, n, alpha, A, ldA, B, ldB);
    }

}
