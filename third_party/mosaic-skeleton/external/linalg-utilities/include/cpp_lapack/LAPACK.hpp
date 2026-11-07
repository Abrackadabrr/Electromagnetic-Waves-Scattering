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


#ifndef LINALG_UTILITIES_CPP_LAPACK_LAPACK_HPP
#define LINALG_UTILITIES_CPP_LAPACK_LAPACK_HPP


#include <cstdint>
#include <cmath>
#include <complex>
#include <type_traits>
#include <algorithm>
#include <limits>


namespace LAPACK
{
    //Auxiliary routines
    
    template <class DataType>
    void lartg(const DataType &a, const DataType &b, decltype(std::abs(DataType{})) &c, DataType &s, DataType &r)
    {
        typedef decltype(std::abs(DataType{})) RealType;
        RealType mx = std::max(std::abs(std::real(a)), std::abs(std::real(b)));
        if (!std::is_same<RealType, DataType>::value)
        {
            mx = std::max(mx, std::abs(std::imag(a)));
            mx = std::max(mx, std::abs(std::imag(b)));
        }
        if (mx == 0)
        {
            c = RealType{1.0};
            s = DataType{0.0};
            r = DataType{0.0};
            return;
        }
        int exp;
        std::frexp(mx, &exp);
        DataType A, B;
        if (std::is_same<RealType, DataType>::value)
        {
            A = std::ldexp(std::real(a), -exp);
            B = std::ldexp(std::real(b), -exp);
        }
        else
        {
            A = DataType{std::ldexp(std::real(a), -exp), std::ldexp(std::imag(a), -exp)};
            B = DataType{std::ldexp(std::real(b), -exp), std::ldexp(std::imag(b), -exp)};
        }
        if (std::abs(A) < std::numeric_limits<RealType>::epsilon())
        {
            c = RealType{0.0};
            s = std::conj(b) / std::abs(b);
            r = abs(b);
            return;
        }
        if (std::abs(A) < std::numeric_limits<RealType>::epsilon())
        {
            c = RealType{1.0};
            s = DataType{0.0};
            r = DataType{0.0};
            return;                
        }
        RealType R = std::sqrt(std::norm(A) + std::norm(B));
        DataType alpha = A / std::abs(A);
        RealType S = RealType{1.0} / R;
        c = std::abs(A) * S;
        s = std::conj(B) * alpha * S;
        r = alpha * std::ldexp(R, exp);
    }

    template<class DataType>
    void lacgv(const uint64_t &n, DataType *x, const uint64_t &incx) {
        for (uint64_t j = 0; j < n; j += incx) {
            x[j] = std::conj(x[j]);
        }
    }

    template<class DataType>
    void tptrs(const char &uplo, const char &transa, const char &diag, const uint64_t &n, const uint64_t &nrhs,
                const DataType *AP, const DataType *B, const uint64_t &ldB, const uint64_t &info);


    template<class DataType>
    int gesdd(const char &jobz, const uint64_t &m, const uint64_t &n, DataType *A, const uint64_t &ldA,
              decltype(std::abs(DataType{})) *S, DataType *U, const uint64_t &ldU, DataType *VT, const uint64_t &ldVT);

    template<class DataType>
    void geqrf(const uint64_t &m, const uint64_t &n, DataType *A, const uint64_t &ldA, DataType *tau,
               DataType *work, const uint64_t &lwork, int &info);

    template<class DataType>
    void larf(const char &side, const uint64_t &m, const uint64_t &n, const DataType *V, const uint64_t &incV,
               const DataType &tau, DataType *C, const uint64_t &ldC, DataType *work);

    template<class DataType>
    void larfb(const char &side, const char &trans, const char &direct, const char &storev,
               const uint64_t &m, const uint64_t &n, const uint64_t &k,
               const DataType *V, const uint64_t &ldV, const DataType *T, const uint64_t &ldT,
               DataType *C, const uint64_t &ldC, DataType *work, const uint64_t &ldwork);

    template <class DataType>
    uint64_t unmqr(const char &side, const char &trans, const uint64_t &m, const uint64_t &n, const uint64_t &k,
                   const DataType *A, const uint64_t &ldA, const DataType *tau, DataType *C, const uint64_t &ldC);

    // Reflection (I - \tau u u*)[alpha x]^T = [beta 0]^T, size(x) = n - 1
    // so that u = [1 ***], beta = sgn(Re(alpha)) * sqrt(|alpha|^2 + ||x||^2)
    template <class DataType>
    void larfg(const uint64_t &n, DataType &alpha, DataType *x, const uint64_t &incx, DataType &tau);


    template <>
    void lartg(const float &a, const float &b, float &c, float &s, float &r);

    template <>
    void lartg(const double &a, const double &b, double &c, double &s, double &r);

    template <>
    void lartg(const std::complex<float> &a, const std::complex<float> &b, float &c, std::complex<float> &s, std::complex<float> &r);

    template <>
    void lartg(const std::complex<double> &a, const std::complex<double> &b, double &c, std::complex<double> &s, std::complex<double> &r);

    template<>
    void lacgv(const uint64_t &n, float *x, const uint64_t &incx);

    template<>
    void lacgv(const uint64_t &n, double *x, const uint64_t &incx);

    template<>
    void lacgv(const uint64_t &n, std::complex<float> *x, const uint64_t &incx);

    template<>
    void lacgv(const uint64_t &n, std::complex<double> *x, const uint64_t &incx);

    template<>
    void tptrs(const char &uplo, const char &transa, const char &diag, const uint64_t &n, const uint64_t &nrhs,
                const float *AP, const float *B, const uint64_t &ldB, const uint64_t &info);

    template<>
    void tptrs(const char &uplo, const char &transa, const char &diag, const uint64_t &n, const uint64_t &nrhs,
                const double *AP, const double *B, const uint64_t &ldB, const uint64_t &info);

    template<>
    void tptrs(const char &uplo, const char &transa, const char &diag, const uint64_t &n, const uint64_t &nrhs,
                const std::complex<float> *AP, const std::complex<float> *B, const uint64_t &ldB, const uint64_t &info);

    template<>
    void tptrs(const char &uplo, const char &transa, const char &diag, const uint64_t &n, const uint64_t &nrhs,
                const std::complex<double> *AP, const std::complex<double> *B, const uint64_t &ldB, const uint64_t &info);

    template<>
    int gesdd(const char &jobz, const uint64_t &m, const uint64_t &n, float *A, const uint64_t &ldA,
              float *S, float *U, const uint64_t &ldU, float *VT, const uint64_t &ldVT);

    template<>
    int gesdd(const char &jobz, const uint64_t &m, const uint64_t &n, double *A, const uint64_t &ldA,
              double *S, double *U, const uint64_t &ldU, double *VT, const uint64_t &ldVT);

    template<>
    int gesdd(const char &jobz, const uint64_t &m, const uint64_t &n, std::complex<float> *A, const uint64_t &ldA,
              float *S, std::complex<float> *U, const uint64_t &ldU, std::complex<float> *VT, const uint64_t &ldVT);

    template<>
    int gesdd(const char &jobz, const uint64_t &m, const uint64_t &n, std::complex<double> *A, const uint64_t &ldA,
              double *S, std::complex<double> *U, const uint64_t &ldU, std::complex<double> *VT, const uint64_t &ldVT);

    template<>
    void geqrf(const uint64_t &m, const uint64_t &n, float *A, const uint64_t &ldA,
               float *tau, float *work, const uint64_t &lwork, int &info);

    template<>
    void geqrf(const uint64_t &m, const uint64_t &n, double *A, const uint64_t &ldA,
               double *tau, double *work, const uint64_t &lwork, int &info);

    template<>
    void geqrf(const uint64_t &m, const uint64_t &n, std::complex<float> *A, const uint64_t &ldA,
               std::complex<float> *tau, std::complex<float> *work, const uint64_t &lwork, int &info);

    template<>
    void geqrf(const uint64_t &m, const uint64_t &n, std::complex<double> *A, const uint64_t &ldA,
               std::complex<double> *tau, std::complex<double> *work, const uint64_t &lwork, int &info);

    template<>
    void larf(const char &side, const uint64_t &m, const uint64_t &n, const float *V, const uint64_t &incV,
              const float &tau, float *C, const uint64_t &ldC, float *work);

    template<>
    void larf(const char &side, const uint64_t &m, const uint64_t &n, const double *V, const uint64_t &incV,
              const double &tau, double *C, const uint64_t &ldC, double *work);

    template<>
    void larf(const char &side, const uint64_t &m, const uint64_t &n, const std::complex<float> *V, const uint64_t &incV,
              const std::complex<float> &tau, std::complex<float> *C, const uint64_t &ldC, std::complex<float> *work);

    template<>
    void larf(const char &side, const uint64_t &m, const uint64_t &n, const std::complex<double> *V, const uint64_t &incV,
              const std::complex<double> &tau, std::complex<double> *C, const uint64_t &ldC, std::complex<double> *work);

    template<>
    void larfb(const char &side, const char &trans, const char &direct, const char &storev,
               const uint64_t &m, const uint64_t &n, const uint64_t &k,
               const float *V, const uint64_t &ldV, const float *T, const uint64_t &ldT,
               float *C, const uint64_t &ldC, float *work, const uint64_t &ldwork);

    template<>
    void larfb(const char &side, const char &trans, const char &direct, const char &storev,
               const uint64_t &m, const uint64_t &n, const uint64_t &k,
               const double *V, const uint64_t &ldV, const double *T, const uint64_t &ldT,
               double *C, const uint64_t &ldC, double *work, const uint64_t &ldwork);

    template<>
    void larfb(const char &side, const char &trans, const char &direct, const char &storev,
               const uint64_t &m, const uint64_t &n, const uint64_t &k,
               const std::complex<float> *V, const uint64_t &ldV, const std::complex<float> *T, const uint64_t &ldT,
               std::complex<float> *C, const uint64_t &ldC, std::complex<float> *work, const uint64_t &ldwork);

    template<>
    void larfb(const char &side, const char &trans, const char &direct, const char &storev,
               const uint64_t &m, const uint64_t &n, const uint64_t &k,
               const std::complex<double> *V, const uint64_t &ldV, const std::complex<double> *T, const uint64_t &ldT,
               std::complex<double> *C, const uint64_t &ldC, std::complex<double> *work, const uint64_t &ldwork);

    template <>
    void larfg(const uint64_t &n, float &alpha, float *x, const uint64_t &incx, float &tau);

    template <>
    void larfg(const uint64_t &n, double &alpha, double *x, const uint64_t &incx, double &tau);

    template <>
    void larfg(const uint64_t &n, std::complex<float> &alpha, std::complex<float> *x, const uint64_t &incx, std::complex<float> &tau);

    template <>
    void larfg(const uint64_t &n, std::complex<double> &alpha, std::complex<double> *x, const uint64_t &incx, std::complex<double> &tau);


    template<>
    uint64_t unmqr(const char &side, const char &trans, const uint64_t &m, const uint64_t &n, const uint64_t &k,
                   const float *A, const uint64_t &ldA, const float *tau, float *C, const uint64_t &ldC);

    template<>
    uint64_t unmqr(const char &side, const char &trans, const uint64_t &m, const uint64_t &n, const uint64_t &k,
                   const double *A, const uint64_t &ldA, const double *tau, double *C, const uint64_t &ldC);

    template <>
    uint64_t unmqr(const char &side, const char &trans, const uint64_t &m, const uint64_t &n, const uint64_t &k,
                   const std::complex<float> *A, const uint64_t &ldA, const std::complex<float> *tau,
                   std::complex<float> *C, const uint64_t &ldC);

    template<>
    uint64_t unmqr(const char &side, const char &trans, const uint64_t &m, const uint64_t &n, const uint64_t &k,
                   const std::complex<double> *A, const uint64_t &ldA, const std::complex<double> *tau,
                   std::complex<double> *C, const uint64_t &ldC);

}
#endif
