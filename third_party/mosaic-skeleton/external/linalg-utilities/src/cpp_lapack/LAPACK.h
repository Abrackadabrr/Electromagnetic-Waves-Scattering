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


#ifndef LINALG_UTILITIES_CPP_LAPACK_LAPACK_H
#define LINALG_UTILITIES_CPP_LAPACK_LAPACK_H


#include "LAPACKFortranMangle.h"


#define slartg LAPACK_GLOBAL(slartg,SLARTG)
#define dlartg LAPACK_GLOBAL(dlartg,DLARTG)
#define clartg LAPACK_GLOBAL(clartg,CLARTG)
#define zlartg LAPACK_GLOBAL(zlartg,ZLARTG)
#define clacgv LAPACK_GLOBAL(clacgv,CLACGV)
#define zlacgv LAPACK_GLOBAL(zlacgv,ZLACGV)
#define stptrs LAPACK_GLOBAL(stptrs,STPTRS)
#define dtptrs LAPACK_GLOBAL(dtptrs,DTPTRS)
#define ctptrs LAPACK_GLOBAL(ctptrs,CTPTRS)
#define ztptrs LAPACK_GLOBAL(ztptrs,ZTPTRS)
#define sgesdd LAPACK_GLOBAL(sgesdd,SGESDD)
#define dgesdd LAPACK_GLOBAL(dgesdd,DGESDD)
#define cgesdd LAPACK_GLOBAL(cgesdd,CGESDD)
#define zgesdd LAPACK_GLOBAL(zgesdd,ZGESDD)
#define sgeqrf LAPACK_GLOBAL(sgeqrf,SGEQRF)
#define dgeqrf LAPACK_GLOBAL(dgeqrf,DGEQRF)
#define cgeqrf LAPACK_GLOBAL(cgeqrf,CGEQRF)
#define zgeqrf LAPACK_GLOBAL(zgeqrf,ZGEQRF)
#define slarf LAPACK_GLOBAL(slarf,SLARF)
#define dlarf LAPACK_GLOBAL(dlarf,DLARF)
#define clarf LAPACK_GLOBAL(clarf,CLARF)
#define zlarf LAPACK_GLOBAL(zlarf,ZLARF)
#define slarfb LAPACK_GLOBAL(slarfb,SLARFB)
#define dlarfb LAPACK_GLOBAL(dlarfb,DLARFB)
#define clarfb LAPACK_GLOBAL(clarfb,CLARFB)
#define zlarfb LAPACK_GLOBAL(zlarfb,ZLARFB)
#define slarfg LAPACK_GLOBAL(slarfg,SLARFG)
#define dlarfg LAPACK_GLOBAL(dlarfg,DLARFG)
#define clarfg LAPACK_GLOBAL(clarfg,CLARFG)
#define zlarfg LAPACK_GLOBAL(zlarfg,ZLARFG)
#define sormqr LAPACK_GLOBAL(sormqr,SORMQR)
#define dormqr LAPACK_GLOBAL(dormqr,DORMQR)
#define cunmqr LAPACK_GLOBAL(cunmqr,CUNMQR)
#define zunmqr LAPACK_GLOBAL(zunmqr,ZUNMQR)

#define LAPACK_INT int


#include <complex>


extern "C"
{
    void slartg(const float &, const float &, float &, float &, float &);

    void dlartg(const double &, const double &, double &, double &, double &);

    void clartg(const std::complex<float> &, const std::complex<float> &, float &, std::complex<float> &, std::complex<float> &);

    void zlartg(const std::complex<double> &, const std::complex<double> &, double &, std::complex<double> &, std::complex<double> &);

    void clacgv(const LAPACK_INT &, std::complex<float> *, const LAPACK_INT &);

    void zlacgv(const LAPACK_INT &, std::complex<double> *, const LAPACK_INT &);

    void stptrs(const char &, const char &, const char &, const LAPACK_INT &, const LAPACK_INT &,
                const float *, const float *, const LAPACK_INT &, const LAPACK_INT &);

    void dtptrs(const char &, const char &, const char &, const LAPACK_INT &, const LAPACK_INT &,
                const double *, const double *, const LAPACK_INT &, const LAPACK_INT &);

    void ctptrs(const char &, const char &, const char &, const LAPACK_INT &, const LAPACK_INT &,
                const std::complex<float> *, const std::complex<float> *, const LAPACK_INT &, const LAPACK_INT &);

    void ztptrs(const char &, const char &, const char &, const LAPACK_INT &, const LAPACK_INT &,
                const std::complex<double> *, const std::complex<double> *, const LAPACK_INT &, const LAPACK_INT &);

    void sgesdd(const char &, const LAPACK_INT &, const LAPACK_INT &, float *, const LAPACK_INT &,
                float *, float *, const LAPACK_INT &, float *, const LAPACK_INT &,
                float *, const LAPACK_INT &, LAPACK_INT *, LAPACK_INT &);

    void dgesdd(const char &, const LAPACK_INT &, const LAPACK_INT &, double *, const LAPACK_INT &,
                double *, double *, const LAPACK_INT &, double *, const LAPACK_INT &,
                double *, const LAPACK_INT &, LAPACK_INT *, LAPACK_INT &);

    void cgesdd(const char &, const LAPACK_INT &, const LAPACK_INT &, std::complex<float> *, const LAPACK_INT &,
                float *, std::complex<float> *, const LAPACK_INT &, std::complex<float> *, const LAPACK_INT &,
                std::complex<float> *, const LAPACK_INT &, float *, LAPACK_INT *, LAPACK_INT &);

    void zgesdd(const char &, const LAPACK_INT &, const LAPACK_INT &, std::complex<double> *, const LAPACK_INT &,
                double *, std::complex<double> *, const LAPACK_INT &, std::complex<double> *, const LAPACK_INT &,
                std::complex<double> *, const LAPACK_INT &, double *, LAPACK_INT *, LAPACK_INT &);

    void sgeqrf(const LAPACK_INT &m, const LAPACK_INT &n, float *A, const LAPACK_INT &ldA,
                float *tau, float *work, const LAPACK_INT &lwork, LAPACK_INT &info);

    void dgeqrf(const LAPACK_INT &m, const LAPACK_INT &n, double *A, const LAPACK_INT &ldA,
                double *tau, double *work, const LAPACK_INT &lwork, LAPACK_INT &info);

    void cgeqrf(const LAPACK_INT &m, const LAPACK_INT &n, std::complex<float> *A, const LAPACK_INT &ldA,
                std::complex<float> *tau, std::complex<float> *work, const LAPACK_INT &lwork, LAPACK_INT &info);

    void zgeqrf(const LAPACK_INT &m, const LAPACK_INT &n, std::complex<double> *A, const LAPACK_INT &ldA,
                std::complex<double> *tau, std::complex<double> *work, const LAPACK_INT &lwork, LAPACK_INT &info);

    void slarf(const char &side, const LAPACK_INT &m, const LAPACK_INT &n, const float *V, const LAPACK_INT &incV,
               const float &tau, float *C, const LAPACK_INT &ldC, float *work);

    void dlarf(const char &side, const LAPACK_INT &m, const LAPACK_INT &n, const double *V, const LAPACK_INT &incV,
               const double &tau, double *C, const LAPACK_INT &ldC, double *work);

    void clarf(const char &side, const LAPACK_INT &m, const LAPACK_INT &n, const std::complex<float> *V, const LAPACK_INT &incV,
               const std::complex<float> &tau, std::complex<float> *C, const LAPACK_INT &ldC, std::complex<float> *work);

    void zlarf(const char &side, const LAPACK_INT &m, const LAPACK_INT &n, const std::complex<double> *V, const LAPACK_INT &incV,
               const std::complex<double> &tau, std::complex<double> *C, const LAPACK_INT &ldC, std::complex<double> *work);

    void slarfb(const char &side, const char &trans, const char &direct, const char &storev,
                const LAPACK_INT &m, const LAPACK_INT &n, const LAPACK_INT &k,
                const float *V, const LAPACK_INT &ldV, const float *T, const LAPACK_INT &ldT,
                float *C, const LAPACK_INT &ldC, float *work, const LAPACK_INT &ldwork);

    void dlarfb(const char &side, const char &trans, const char &direct, const char &storev,
                const LAPACK_INT &m, const LAPACK_INT &n, const LAPACK_INT &k,
                const double *V, const LAPACK_INT &ldV, const double *T, const LAPACK_INT &ldT,
                double *C, const LAPACK_INT &ldC, double *work, const LAPACK_INT &ldwork);

    void clarfb(const char &side, const char &trans, const char &direct, const char &storev,
                const LAPACK_INT &m, const LAPACK_INT &n, const LAPACK_INT &k,
                const std::complex<float> *V, const LAPACK_INT &ldV, const std::complex<float> *T, const LAPACK_INT &ldT,
                std::complex<float> *C, const LAPACK_INT &ldC, std::complex<float> *work, const LAPACK_INT &ldwork);

    void zlarfb(const char &side, const char &trans, const char &direct, const char &storev,
                const LAPACK_INT &m, const LAPACK_INT &n, const LAPACK_INT &k,
                const std::complex<double> *V, const LAPACK_INT &ldV, const std::complex<double> *T, const LAPACK_INT &ldT,
                std::complex<double> *C, const LAPACK_INT &ldC, std::complex<double> *work, const LAPACK_INT &ldwork);

    void slarfg(const LAPACK_INT &, float &, float *, const LAPACK_INT &, float &);

    void dlarfg(const LAPACK_INT &, double &, double *, const LAPACK_INT &, double &);

    void clarfg(const LAPACK_INT &, std::complex<float> &, std::complex<float> *, const LAPACK_INT &, std::complex<float> &);

    void zlarfg(const LAPACK_INT &, std::complex<double> &, std::complex<double> *, const LAPACK_INT &, std::complex<double> &);


    void sormqr(const char &, const char &, const LAPACK_INT &, const LAPACK_INT &, const LAPACK_INT &,
                const float *, const LAPACK_INT &, const float *, float *, const LAPACK_INT &,
                float *, const LAPACK_INT &, const LAPACK_INT &);

    void dormqr(const char &, const char &, const LAPACK_INT &, const LAPACK_INT &, const LAPACK_INT &,
                const double *, const LAPACK_INT &, const double *, double *, const LAPACK_INT &,
                double *, const LAPACK_INT &, const LAPACK_INT &);

    void cunmqr(const char &, const char &, const LAPACK_INT &, const LAPACK_INT &, const LAPACK_INT &,
                const std::complex<float> *, const LAPACK_INT &, const std::complex<float> *,
                std::complex<float> *, const LAPACK_INT &, std::complex<float> *, const LAPACK_INT &, const LAPACK_INT &);

    void zunmqr(const char &, const char &, const LAPACK_INT &, const LAPACK_INT &, const LAPACK_INT &,
                const std::complex<double> *, const LAPACK_INT &, const std::complex<double> *,
                std::complex<double> *, const LAPACK_INT &, std::complex<double> *, const LAPACK_INT &, const LAPACK_INT &);

}

#endif
