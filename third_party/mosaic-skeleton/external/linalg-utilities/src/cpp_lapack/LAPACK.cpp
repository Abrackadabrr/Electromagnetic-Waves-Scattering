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


#include "LAPACK.h"

#include "cpp_lapack/LAPACK.hpp"
#include <complex>
#include <cstdint>


namespace LAPACK {

template <>
void lartg(const float &a, const float &b, float &c, float &s, float &r)
{
    slartg(a, b, c, s, r);
}

template <>
void lartg(const double &a, const double &b, double &c, double &s, double &r)
{
    dlartg(a, b, c, s, r);            
}

template <>
void lartg(const std::complex<float> &a, const std::complex<float> &b, float &c, std::complex<float> &s, std::complex<float> &r)
{
    clartg(a, b, c, s, r);
}

template <>
void lartg(const std::complex<double> &a, const std::complex<double> &b, double &c, std::complex<double> &s, std::complex<double> &r)
{
    zlartg(a, b, c, s, r);
}


template<>
void lacgv(const uint64_t &n, float *x, const uint64_t &incx) {
    (void) n;
    (void) x;
    (void) incx;
}

template<>
void lacgv(const uint64_t &n, double *x, const uint64_t &incx) {
    (void) n;
    (void) x;
    (void) incx;
}

template<>
void lacgv(const uint64_t &n, std::complex<float> *x, const uint64_t &incx) {
    clacgv(n, x, incx);
}

template<>
void lacgv(const uint64_t &n, std::complex<double> *x, const uint64_t &incx) {
    zlacgv(n, x, incx);
}

template<>
void tptrs(const char &uplo, const char &transa, const char &diag, const uint64_t &n, const uint64_t &nrhs,
           const float *AP, const float *B, const uint64_t &ldB, const uint64_t &info) {
    stptrs(uplo, transa, diag, n, nrhs, AP, B, ldB, info);
}

template<>
void tptrs(const char &uplo, const char &transa, const char &diag, const uint64_t &n, const uint64_t &nrhs,
           const double *AP, const double *B, const uint64_t &ldB, const uint64_t &info) {
    dtptrs(uplo, transa, diag, n, nrhs, AP, B, ldB, info);
}

template<>
void tptrs(const char &uplo, const char &transa, const char &diag, const uint64_t &n, const uint64_t &nrhs,
           const std::complex<float> *AP, const std::complex<float> *B, const uint64_t &ldB, const uint64_t &info){
    ctptrs(uplo, transa, diag, n, nrhs, AP, B, ldB, info);
}

template<>
void tptrs(const char &uplo, const char &transa, const char &diag, const uint64_t &n, const uint64_t &nrhs,
           const std::complex<double> *AP, const std::complex<double> *B, const uint64_t &ldB, const uint64_t &info) {
    ztptrs(uplo, transa, diag, n, nrhs, AP, B, ldB, info);
}

template<>
int gesdd(const char &jobz, const uint64_t &m, const uint64_t &n, float *A, const uint64_t &ldA,
          float *S, float *U, const uint64_t &ldU, float *VT, const uint64_t &ldVT) {
    const std::size_t min_mn = std::min(m, n);
    
    // Prepare working arrays

    // Integer work array
    const std::size_t iwork_size = 8 * min_mn;
    auto iwork = new int[iwork_size];

    // Ask for lwork for real work array
    int info{0};
    float q_work{0};
    sgesdd(jobz, m, n, nullptr, ldA, nullptr, nullptr, ldU, nullptr, ldVT, &q_work, -1, iwork, info);
    if (info != 0) {
        return info;
    }
    const std::size_t work_size = q_work + 1;
    auto work = new float[work_size];

    sgesdd(jobz, m, n, A, ldA, S, U, ldU, VT, ldVT, work, work_size, iwork, info);

    delete[] iwork;
    delete[] work;

    return info;
}

template<>
int gesdd(const char &jobz, const uint64_t &m, const uint64_t &n, double *A, const uint64_t &ldA,
          double *S, double *U, const uint64_t &ldU, double *VT, const uint64_t &ldVT) {
    const std::size_t min_mn = std::min(m, n);
    
    // Prepare working arrays

    // Integer work array
    const std::size_t iwork_size = 8 * min_mn;
    auto iwork = new int[iwork_size];

    // Ask for lwork for real work array
    int info{0};
    double q_work{0};
    dgesdd(jobz, m, n, nullptr, ldA, nullptr, nullptr, ldU, nullptr, ldVT, &q_work, -1, iwork, info);
    if (info != 0) {
        return info;
    }
    const std::size_t work_size = q_work + 1;
    auto work = new double[work_size];

    dgesdd(jobz, m, n, A, ldA, S, U, ldU, VT, ldVT, work, work_size, iwork, info);
    
    delete[] iwork;
    delete[] work;

    return info;
}

template<>
int gesdd(const char &jobz, const uint64_t &m, const uint64_t &n, std::complex<float> *A, const uint64_t &ldA,
          float *S, std::complex<float> *U, const uint64_t &ldU, std::complex<float> *VT, const uint64_t &ldVT) {
   const std::size_t min_mn = std::min(m, n);
    
    // Prepare working arrays

    // Integer work array
    const std::size_t iwork_size = 8 * min_mn;
    auto iwork = new int[iwork_size];

    // Ask for lwork for complex work array
    int info{0};
    std::complex<float> q_work{0};
    cgesdd(jobz, m, n, nullptr, ldA, nullptr, nullptr, ldU, nullptr, ldVT, &q_work, -1, nullptr, iwork, info);
    if (info != 0) {
        return info;
    }
    const std::size_t work_size = std::real(q_work) + 1;
    auto work = new std::complex<float>[work_size];

    // Real work array
    const std::size_t N_lrwork = 7 * min_mn;
    const std::size_t rect_lrwork = 5 * min_mn * (min_mn + 1);
    const std::size_t sq_lrwork = 2 * m * n + min_mn * (2 * min_mn + 1);
    std::size_t rwork_size{0};
    if (jobz == 'N') {
        rwork_size = N_lrwork;
    } else if (std::max(m, n) > 8 * min_mn) {
        rwork_size = rect_lrwork;
    } else {
        rwork_size = std::max(rect_lrwork, sq_lrwork);
    }
    auto rwork = new float[rwork_size];

    cgesdd(jobz, m, n, A, ldA, S, U, ldU, VT, ldVT, work, work_size, rwork, iwork, info);
    
    delete[] iwork;
    delete[] work;
    delete[] rwork;

    return info;
}

template<>
int gesdd(const char &jobz, const uint64_t &m, const uint64_t &n, std::complex<double> *A, const uint64_t &ldA,
          double *S, std::complex<double> *U, const uint64_t &ldU, std::complex<double> *VT, const uint64_t &ldVT) {
    const std::size_t min_mn = std::min(m, n);
    
    // Prepare working arrays

    // Integer work array
    const std::size_t iwork_size = 8 * min_mn;
    auto iwork = new int[iwork_size];

    // Ask for lwork for complex work array
    int info{0};
    std::complex<double> q_work{0};
    zgesdd(jobz, m, n, nullptr, ldA, nullptr, nullptr, ldU, nullptr, ldVT, &q_work, -1, nullptr, iwork, info);
    if (info != 0) {
        return info;
    }
    const std::size_t work_size = std::real(q_work) + 1;
    auto work = new std::complex<double>[work_size];

    // Real work array
    const std::size_t N_lrwork = 7 * min_mn;
    const std::size_t rect_lrwork = 5 * min_mn * (min_mn + 1);
    const std::size_t sq_lrwork = 2 * m * n + min_mn * (2 * min_mn + 1);
    std::size_t rwork_size{0};
    if (jobz == 'N') {
        rwork_size = N_lrwork;
    } else if (std::max(m, n) > 8 * min_mn) {
        rwork_size = rect_lrwork;
    } else {
        rwork_size = std::max(rect_lrwork, sq_lrwork);
    }
    auto rwork = new double[rwork_size];

    zgesdd(jobz, m, n, A, ldA, S, U, ldU, VT, ldVT, work, work_size, rwork, iwork, info);
    
    delete[] iwork;
    delete[] work;
    delete[] rwork;

    return info;
}

template<>
void geqrf(const uint64_t &m, const uint64_t &n, float *A, const uint64_t &ldA,
            float *tau, float *work, const uint64_t &lwork, int &info) {
    sgeqrf(m, n, A, ldA, tau, work, lwork, info);
}

template<>
void geqrf(const uint64_t &m, const uint64_t &n, double *A, const uint64_t &ldA,
            double *tau, double *work, const uint64_t &lwork, int &info) {
    dgeqrf(m, n, A, ldA, tau, work, lwork, info);
}

template<>
void geqrf(const uint64_t &m, const uint64_t &n, std::complex<float> *A, const uint64_t &ldA,
            std::complex<float> *tau, std::complex<float> *work, const uint64_t &lwork, int &info) {
    cgeqrf(m, n, A, ldA, tau, work, lwork, info);
}

template<>
void geqrf(const uint64_t &m, const uint64_t &n, std::complex<double> *A, const uint64_t &ldA,
            std::complex<double> *tau, std::complex<double> *work, const uint64_t &lwork, int &info) {
    zgeqrf(m, n, A, ldA, tau, work, lwork, info);
}

template<>
void larf(const char &side, const uint64_t &m, const uint64_t &n, const float *V, const uint64_t &incV,
            const float &tau, float *C, const uint64_t &ldC, float *work) {
    slarf(side, m, n, V, incV, tau, C, ldC, work);
}

template<>
void larf(const char &side, const uint64_t &m, const uint64_t &n, const double *V, const uint64_t &incV,
            const double &tau, double *C, const uint64_t &ldC, double *work) {
    dlarf(side, m, n, V, incV, tau, C, ldC, work);
}

template<>
void larf(const char &side, const uint64_t &m, const uint64_t &n, const std::complex<float> *V, const uint64_t &incV,
            const std::complex<float> &tau, std::complex<float> *C, const uint64_t &ldC, std::complex<float> *work) {
    clarf(side, m, n, V, incV, tau, C, ldC, work);
}

template<>
void larf(const char &side, const uint64_t &m, const uint64_t &n, const std::complex<double> *V, const uint64_t &incV,
            const std::complex<double> &tau, std::complex<double> *C, const uint64_t &ldC, std::complex<double> *work) {
    zlarf(side, m, n, V, incV, tau, C, ldC, work);
}

template<>
void larfb(const char &side, const char &trans, const char &direct, const char &storev,
            const uint64_t &m, const uint64_t &n, const uint64_t &k,
            const float *V, const uint64_t &ldV, const float *T, const uint64_t &ldT,
            float *C, const uint64_t &ldC, float *work, const uint64_t &ldwork) {
    slarfb(side, trans, direct, storev, m, n, k, V, ldV, T, ldT, C, ldC, work, ldwork);
}

template<>
void larfb(const char &side, const char &trans, const char &direct, const char &storev,
            const uint64_t &m, const uint64_t &n, const uint64_t &k,
            const double *V, const uint64_t &ldV, const double *T, const uint64_t &ldT,
            double *C, const uint64_t &ldC, double *work, const uint64_t &ldwork) {
    dlarfb(side, trans, direct, storev, m, n, k, V, ldV, T, ldT, C, ldC, work, ldwork);
}

template<>
void larfb(const char &side, const char &trans, const char &direct, const char &storev,
            const uint64_t &m, const uint64_t &n, const uint64_t &k,
            const std::complex<float> *V, const uint64_t &ldV, const std::complex<float> *T, const uint64_t &ldT,
            std::complex<float> *C, const uint64_t &ldC, std::complex<float> *work, const uint64_t &ldwork) {
    clarfb(side, trans, direct, storev, m, n, k, V, ldV, T, ldT, C, ldC, work, ldwork);
}

template<>
void larfb(const char &side, const char &trans, const char &direct, const char &storev,
            const uint64_t &m, const uint64_t &n, const uint64_t &k,
            const std::complex<double> *V, const uint64_t &ldV, const std::complex<double> *T, const uint64_t &ldT,
            std::complex<double> *C, const uint64_t &ldC, std::complex<double> *work, const uint64_t &ldwork) {
    zlarfb(side, trans, direct, storev, m, n, k, V, ldV, T, ldT, C, ldC, work, ldwork);
}

template<>
void larfg(const uint64_t &n, float &alpha, float *x, const uint64_t &incx, float &tau)
{
    slarfg(n, alpha, x, incx, tau);
}

template<>
void larfg(const uint64_t &n, double &alpha, double *x, const uint64_t &incx, double &tau)
{
    dlarfg(n, alpha, x, incx, tau);
}

template<>
void larfg(const uint64_t &n, std::complex<float> &alpha, std::complex<float> *x, const uint64_t &incx, std::complex<float> &tau)
{
    clarfg(n, alpha, x, incx, tau);
}

template<>
void larfg(const uint64_t &n, std::complex<double> &alpha, std::complex<double> *x, const uint64_t &incx, std::complex<double> &tau)
{
    zlarfg(n, alpha, x, incx, tau);
}


template<>
uint64_t unmqr(const char &side, const char &trans, const uint64_t &m, const uint64_t &n, const uint64_t &k,
                const float *A, const uint64_t &ldA, const float *tau, float *C, const uint64_t &ldC)
{
    LAPACK_INT info{};
    char t = trans;
    if (t == 'C') t = 'T';
    uint64_t lwork = (side == 'L') * n + (side == 'R') * m + 1;
    auto work = new float[lwork];
    sormqr(side, t, m, n, k, A, ldA, tau, C, ldC, work, lwork, info);
    delete []work;
    return info;
}

template<>
uint64_t unmqr(const char &side, const char &trans, const uint64_t &m, const uint64_t &n, const uint64_t &k,
                const double *A, const uint64_t &ldA, const double *tau, double *C, const uint64_t &ldC)
{
    LAPACK_INT info{};
    char t = trans;
    if (t == 'C') t = 'T';
    uint64_t lwork = (side == 'L') * n + (side == 'R') * m + 1;
    auto work = new double[lwork];
    dormqr(side, t, m, n, k, A, ldA, tau, C, ldC, work, lwork, info);
    delete []work;
    return info;
}

template<>
uint64_t unmqr(const char &side, const char &trans, const uint64_t &m, const uint64_t &n, const uint64_t &k,
                const std::complex<float> *A, const uint64_t &ldA, const std::complex<float> *tau,
                std::complex<float> *C, const uint64_t &ldC)
{
    uint64_t info{};
    uint64_t lwork = (side == 'L') * n + (side == 'R') * m + 1;
    auto work = new std::complex<float>[lwork];
    cunmqr(side, trans, m, n, k, A, ldA, tau, C, ldC, work, lwork, info);
    delete []work;
    return info;
}

template<>
uint64_t unmqr(const char &side, const char &trans, const uint64_t &m, const uint64_t &n, const uint64_t &k,
                const std::complex<double> *A, const uint64_t &ldA, const std::complex<double> *tau,
                std::complex<double> *C, const uint64_t &ldC)
{
    uint64_t info{};
    uint64_t lwork = (side == 'L') * n + (side == 'R') * m + 1;
    auto work = new std::complex<double>[lwork];
    zunmqr(side, trans, m, n, k, A, ldA, tau, C, ldC, work, lwork, info);
    delete []work;
    return info;
}


}
