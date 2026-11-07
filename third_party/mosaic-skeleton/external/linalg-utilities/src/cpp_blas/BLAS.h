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


#ifndef LINALG_UTILITIES_CPP_BLAS_BLAS_H
#define LINALG_UTILITIES_CPP_BLAS_BLAS_H


#include "BLASFortranMangle.h"


#define sswap BLAS_GLOBAL(sswap,SSWAP)
#define dswap BLAS_GLOBAL(dswap,DSWAP)
#define cswap BLAS_GLOBAL(cswap,CSWAP)
#define zswap BLAS_GLOBAL(zswap,ZSWAP)
#define saxpy BLAS_GLOBAL(saxpy,SAXPY)
#define daxpy BLAS_GLOBAL(daxpy,DAXPY)
#define caxpy BLAS_GLOBAL(caxpy,CAXPY)
#define zaxpy BLAS_GLOBAL(zaxpy,ZAXPY)
#define scopy BLAS_GLOBAL(scopy,SCOPY)
#define dcopy BLAS_GLOBAL(dcopy,DCOPY)
#define ccopy BLAS_GLOBAL(ccopy,CCOPY)
#define zcopy BLAS_GLOBAL(zcopy,ZCOPY)
#define snrm2 BLAS_GLOBAL(snrm2,SNRM2)
#define dnrm2 BLAS_GLOBAL(dnrm2,DNRM2)
#define scnrm2 BLAS_GLOBAL(scnrm2,sCNRM2)
#define dznrm2 BLAS_GLOBAL(dznrm2,DZNRM2)
#define sdot BLAS_GLOBAL(sdot,SDOT)
#define ddot BLAS_GLOBAL(ddot,DDOT)
#define cdotc_stackreturn BLAS_GLOBAL(cdotc_stackreturn,CDOTC_STACKRETURN)
#define cdotu_stackreturn BLAS_GLOBAL(cdotu_stackreturn,CDOTU_STACKRETURN)
#define zdotc_stackreturn BLAS_GLOBAL(zdotc_stackreturn,ZDOTC_STACKRETURN)
#define zdotu_stackreturn BLAS_GLOBAL(zdotu_stackreturn,ZDOTU_STACKRETURN)
#define sscal BLAS_GLOBAL(sscal,SSCAL)
#define dscal BLAS_GLOBAL(dscal,DSCAL)
#define cscal BLAS_GLOBAL(cscal,CSCAL)
#define csscal BLAS_GLOBAL(csscal,CSSCAL)
#define zscal BLAS_GLOBAL(zscal,ZSCAL)
#define zdscal BLAS_GLOBAL(zdscal,ZDSCAL)
#define srot BLAS_GLOBAL(srot,SROT)
#define drot BLAS_GLOBAL(drot,DROT)
#define crot BLAS_GLOBAL(crot,CROT)
#define zrot BLAS_GLOBAL(zrot,ZROT)
#define sger BLAS_GLOBAL(sger,SGER)
#define dger BLAS_GLOBAL(dger,DGER)
#define cgerc BLAS_GLOBAL(cgerc,CGERC)
#define cgeru BLAS_GLOBAL(cgeru,CGERU)
#define zgerc BLAS_GLOBAL(zgerc,ZGERC)
#define zgeru BLAS_GLOBAL(zgeru,ZGERU)
#define sgemv BLAS_GLOBAL(sgemv,SGEMV)
#define dgemv BLAS_GLOBAL(dgemv,DGEMV)
#define cgemv BLAS_GLOBAL(cgemv,CGEMV)
#define zgemv BLAS_GLOBAL(zgemv,ZGEMV)
#define stpsv BLAS_GLOBAL(stpsv,STPSV)
#define dtpsv BLAS_GLOBAL(dtpsv,DTPSV)
#define ctpsv BLAS_GLOBAL(ctpsv,CTPSV)
#define ztpsv BLAS_GLOBAL(ztpsv,ZTPSV)
#define sgemm BLAS_GLOBAL(sgemm,SGEMM)
#define dgemm BLAS_GLOBAL(dgemm,DGEMM)
#define cgemm BLAS_GLOBAL(cgemm,CGEMM)
#define zgemm BLAS_GLOBAL(zgemm,ZGEMM)
#define strmm BLAS_GLOBAL(strmm,STRMM)
#define dtrmm BLAS_GLOBAL(dtrmm,DTRMM)
#define ctrmm BLAS_GLOBAL(ctrmm,CTRMM)
#define ztrmm BLAS_GLOBAL(ztrmm,ZTRMM)
#define strsm BLAS_GLOBAL(strsm,STRSM)
#define dtrsm BLAS_GLOBAL(dtrsm,DTRSM)
#define ctrsm BLAS_GLOBAL(ctrsm,CTRSM)
#define ztrsm BLAS_GLOBAL(ztrsm,ZTRSM)

#define BLAS_INT int


#include <complex>


extern "C"
{
    void sswap(const BLAS_INT &,
               const float *, const BLAS_INT &,
               float *, const BLAS_INT &);

    void dswap(const BLAS_INT &,
               const double *, const BLAS_INT &,
               double *, const BLAS_INT &);

    void cswap(const BLAS_INT &,
               const std::complex<float> *, const BLAS_INT &,
               std::complex<float> *, const BLAS_INT &);

    void zswap(const BLAS_INT &,
               const std::complex<double> *, const BLAS_INT &,
               std::complex<double> *, const BLAS_INT &);

    void saxpy(const BLAS_INT &,
               const float &,
               const float *, const BLAS_INT &,
               float *, const BLAS_INT &);

    void daxpy(const BLAS_INT &,
               const double &,
               const double *, const BLAS_INT &,
               double *, const BLAS_INT &);

    void caxpy(const BLAS_INT &,
               const std::complex<float> &,
               const std::complex<float> *, const BLAS_INT &,
               std::complex<float> *, const BLAS_INT &);

    void zaxpy(const BLAS_INT &,
               const std::complex<double> &,
               const std::complex<double> *, const BLAS_INT &,
               std::complex<double> *, const BLAS_INT &);

    void scopy(const BLAS_INT &,
               const float *, const BLAS_INT &,
               float *, const BLAS_INT &);

    void dcopy(const BLAS_INT &,
               const double *, const BLAS_INT &,
               double *, const BLAS_INT &);

    void ccopy(const BLAS_INT &,
               const std::complex<float> *, const BLAS_INT &,
               std::complex<float> *, const BLAS_INT &);

    void zcopy(const BLAS_INT &,
               const std::complex<double> *, const BLAS_INT &,
               std::complex<double> *, const BLAS_INT &);

    float sdot(const BLAS_INT &,
               const float *, const BLAS_INT &,
               const float *, const BLAS_INT &);

    double ddot(const BLAS_INT &,
               const double *, const BLAS_INT &,
               const double *, const BLAS_INT &);

    void cdotu_stackreturn(std::complex<float> &, const BLAS_INT &,
               const std::complex<float> *, const BLAS_INT &,
               const std::complex<float> *, const BLAS_INT &);

    void cdotc_stackreturn(std::complex<float> &, const BLAS_INT &,
               const std::complex<float> *, const BLAS_INT &,
               const std::complex<float> *, const BLAS_INT &);

    void zdotu_stackreturn(std::complex<double> &, const BLAS_INT &,
               const std::complex<double> *, const BLAS_INT &,
               const std::complex<double> *, const BLAS_INT &);

    void zdotc_stackreturn(std::complex<double> &, const BLAS_INT &,
               const std::complex<double> *, const BLAS_INT &,
               const std::complex<double> *, const BLAS_INT &);

    float snrm2(const BLAS_INT &,
                const float *, const BLAS_INT &);

    double dnrm2(const BLAS_INT &,
                  const double *, const BLAS_INT &);

    float scnrm2(const BLAS_INT &,
                 const std::complex<float> *, const BLAS_INT &);

    double dznrm2(const BLAS_INT &,
                  const std::complex<double> *, const BLAS_INT &);

    void sscal(const BLAS_INT &,
               const float &,
               const float *, const BLAS_INT &);

    void dscal(const BLAS_INT &,
               const double &,
               const double *, const BLAS_INT &);

    void cscal(const BLAS_INT &,
               const std::complex<float> &,
               const std::complex<float> *, const BLAS_INT &);

    void csscal(const BLAS_INT &,
                const float &,
                const std::complex<float> *, const BLAS_INT &);

    void zscal(const BLAS_INT &,
               const std::complex<double> &,
               const std::complex<double> *, const BLAS_INT &);

    void zdscal(const BLAS_INT &,
                const double &,
                const std::complex<double> *, const BLAS_INT &);

    void srot(const BLAS_INT &,
              float *, const BLAS_INT &,
              float *, const BLAS_INT &,
              const float &, const float &);

    void drot(const BLAS_INT &,
              double *, const BLAS_INT &,
              double *, const BLAS_INT &,
              const double &, const double &);

    void crot(const BLAS_INT &,
              std::complex<float> *, const BLAS_INT &,
              std::complex<float> *, const BLAS_INT &,
              const float &, const std::complex<float> &);

    void zrot(const BLAS_INT &,
              std::complex<double> *, const BLAS_INT &,
              std::complex<double> *, const BLAS_INT &,
              const double &, const std::complex<double> &);

    void sger(const BLAS_INT &, const BLAS_INT &,
              const float &,
              const float *, const BLAS_INT &,
              const float *, const BLAS_INT &,
              float *, const BLAS_INT &);

    void dger(const BLAS_INT &, const BLAS_INT &,
              const double &,
              const double *, const BLAS_INT &,
              const double *, const BLAS_INT &,
              double *, const BLAS_INT &);

    void cgerc(const BLAS_INT &, const BLAS_INT &,
               const std::complex<float> &,
               const std::complex<float> *, const BLAS_INT &,
               const std::complex<float> *, const BLAS_INT &,
               std::complex<float> *, const BLAS_INT &);

    void zgerc(const BLAS_INT &, const BLAS_INT &,
               const std::complex<double> &,
               const std::complex<double> *, const BLAS_INT &,
               const std::complex<double> *, const BLAS_INT &,
               std::complex<double> *, const BLAS_INT &);

    void cgeru(const BLAS_INT &, const BLAS_INT &,
               const std::complex<float> &,
               const std::complex<float> *, const BLAS_INT &,
               const std::complex<float> *, const BLAS_INT &,
               std::complex<float> *, const BLAS_INT &);

    void zgeru(const BLAS_INT &, const BLAS_INT &,
               const std::complex<double> &,
               const std::complex<double> *, const BLAS_INT &,
               const std::complex<double> *, const BLAS_INT &,
               std::complex<double> *, const BLAS_INT &);

    void sgemv(const char &,
               const BLAS_INT &, const BLAS_INT &,
               const float &,
               const float *, const BLAS_INT &,
               const float *, const BLAS_INT &,
               const float &, 
               float *, const BLAS_INT &);

    void dgemv(const char &,
               const BLAS_INT &, const BLAS_INT &,
               const double &,
               const double *, const BLAS_INT &,
               const double *, const BLAS_INT &,
               const double &, 
               double *, const BLAS_INT &);

    void cgemv(const char &,
               const BLAS_INT &, const BLAS_INT &,
               const std::complex<float> &,
               const std::complex<float> *, const BLAS_INT &,
               const std::complex<float> *, const BLAS_INT &,
               const std::complex<float> &, 
               std::complex<float> *, const BLAS_INT &);

    void zgemv(const char &,
               const BLAS_INT &, const BLAS_INT &,
               const std::complex<double> &,
               const std::complex<double> *, const BLAS_INT &,
               const std::complex<double> *, const BLAS_INT &,
               const std::complex<double> &, 
               std::complex<double> *, const BLAS_INT &);

    void stpsv(const char &,
               const char &,
               const char &,
               const BLAS_INT &,
               const float *,
               float *, const BLAS_INT &);

    void dtpsv(const char &,
               const char &,
               const char &,
               const BLAS_INT &,
               const double *,
               double *, const BLAS_INT &);

    void ctpsv(const char &,
               const char &,
               const char &,
               const BLAS_INT &,
               const std::complex<float> *,
               std::complex<float> *, const BLAS_INT &);

    void ztpsv(const char &,
               const char &,
               const char &,
               const BLAS_INT &,
               const std::complex<double> *,
               std::complex<double> *, const BLAS_INT &);

    void sgemm(const char &,
               const char &,
               const BLAS_INT &, const BLAS_INT &,
               const BLAS_INT &,
               const float &,
               const float *, const BLAS_INT &,
               const float *, const BLAS_INT &,
               const float &,
               float *, const BLAS_INT &);

    void dgemm(const char &,
               const char &,
               const BLAS_INT &, const BLAS_INT &,
               const BLAS_INT &,
               const double &,
               const double *, const BLAS_INT &,
               const double *, const BLAS_INT &,
               const double &,
               double *, const BLAS_INT &);

    void cgemm(const char &,
               const char &,
               const BLAS_INT &, const BLAS_INT &,
               const BLAS_INT &,
               const std::complex<float> &,
               const std::complex<float> *, const BLAS_INT &,
               const std::complex<float> *, const BLAS_INT &,
               const std::complex<float> &,
               std::complex<float> *, const BLAS_INT &);

    void zgemm(const char &,
               const char &,
               const BLAS_INT &, const BLAS_INT &,
               const BLAS_INT &,
               const std::complex<double> &,
               const std::complex<double> *, const BLAS_INT &,
               const std::complex<double> *, const BLAS_INT &,
               const std::complex<double> &,
               std::complex<double> *, const BLAS_INT &);

    void strmm(const char &, const char &, const char &, const char &, const BLAS_INT &, const BLAS_INT &,
               const float &, const float *, const BLAS_INT &,
               float *, const BLAS_INT &);

    void dtrmm(const char &, const char &, const char &, const char &, const BLAS_INT &, const BLAS_INT &,
               const double &, const double *, const BLAS_INT &,
               double *, const BLAS_INT &);

    void ctrmm(const char &, const char &, const char &, const char &, const BLAS_INT &, const BLAS_INT &,
               const std::complex<float> &, const std::complex<float> *, const BLAS_INT &,
               std::complex<float> *, const BLAS_INT &);

    void ztrmm(const char &, const char &, const char &, const char &, const BLAS_INT &, const BLAS_INT &,
               const std::complex<double> &, const std::complex<double> *, const BLAS_INT &,
               std::complex<double> *, const BLAS_INT &);


    void strsm(const char &, const char &, const char &, const char &, const BLAS_INT &, const BLAS_INT &,
                const float &, const float *, const BLAS_INT &,
                float *, const BLAS_INT &);
 
     void dtrsm(const char &, const char &, const char &, const char &, const BLAS_INT &, const BLAS_INT &,
                const double &, const double *, const BLAS_INT &,
                double *, const BLAS_INT &);
 
     void ctrsm(const char &, const char &, const char &, const char &, const BLAS_INT &, const BLAS_INT &,
                const std::complex<float> &, const std::complex<float> *, const BLAS_INT &,
                std::complex<float> *, const BLAS_INT &);
 
     void ztrsm(const char &, const char &, const char &, const char &, const BLAS_INT &, const BLAS_INT &,
                const std::complex<double> &, const std::complex<double> *, const BLAS_INT &,
                std::complex<double> *, const BLAS_INT &);

}

#endif
