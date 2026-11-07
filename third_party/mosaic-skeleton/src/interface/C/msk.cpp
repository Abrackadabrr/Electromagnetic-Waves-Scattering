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
//
// Created by bulat_v on 26.07.23.
//


#include "msk/interface/C/msk.h"

#include "msk/interface/C/msk.hpp"


// Additional helping fucntions
MSk::oracle::Parameters convert_parameters(mskParameters *c_parameters) {
    return MSk::oracle::Parameters{
        c_parameters->min_rows,
        c_parameters->min_cols,
        c_parameters->rho,
        c_parameters->sphere_norm,
        MSk::oracle::DividePolicy::dual,
        c_parameters->force_lowrank,
        false,
        MSk::log::build_logger("MSkOracle", c_parameters->log_file, MSk::log::string_to_level(c_parameters->log_level)),
        {c_parameters->abs_tol, c_parameters->rel_tol, c_parameters->max_rank},
        false,
        MSk::log::build_logger("MSkApprox", c_parameters->log_file, MSk::log::string_to_level(c_parameters->log_level))
    };
}

template<class DataType>
void set_matrix_data(MSk::MosaicSkeleton<DataType> *msk_ptr, mskMosaicMatrix *msk_matrix) {

    msk_matrix->_matrix = msk_ptr;

    std::tie(msk_matrix->row_size, msk_matrix->col_size) = msk_ptr->size();

    msk_matrix->frob_norm = msk_ptr->frob_norm();

    msk_matrix->memory = msk_ptr->memory();

#ifdef MSK_USE_MPI
    msk_matrix->local_row_size = msk_ptr->local_row_size();
    msk_matrix->local_col_size = msk_ptr->local_col_size();
    msk_matrix->local_row_off = msk_ptr->local_row_off();
    msk_matrix->local_col_off = msk_ptr->local_col_off();

    msk_matrix->local_frob_norm = msk_ptr->frob_norm(false);

    msk_matrix->local_memory = msk_ptr->memory(false);
#endif
}


int msk_build_basic_s3(mskMosaicMatrix *msk_matrix,
                       mskParameters *parameters,
                       uint64_t m, uint64_t n,
                       float (*f_elem)(uint64_t, uint64_t, void *), void *data,
                       const float *row_pts, const float *col_pts,
                       uint64_t n_jobs
#ifdef MSK_USE_MPI
                        , MPI_Comm comm
#endif
)
{
    using namespace std::placeholders;
    std::function<float(uint64_t, uint64_t)> f = std::bind(f_elem, _1, _2, data);
    MSk::oracle::Parameters msk_params = convert_parameters(parameters);
    MSk::MosaicSkeleton<float> *msk_ptr{nullptr};
    MSk::C::Result result = MSk::C::build_basic<float, float, 3>(
            &msk_ptr, m, n, f, row_pts, col_pts, msk_params, n_jobs
#ifdef MSK_USE_MPI
            , comm
#endif
    );
    if (result != MSk::C::Result::OK) {
        return result;
    }
    set_matrix_data(msk_ptr, msk_matrix);
    return MSk::C::Result::OK;
}

int msk_build_basic_d3(mskMosaicMatrix *msk_matrix,
                        mskParameters *parameters,
                        uint64_t m, uint64_t n,
                        double (*f_elem)(uint64_t, uint64_t, void *), void *data,
                        const double *row_pts, const double *col_pts,
                        uint64_t n_jobs
#ifdef MSK_USE_MPI
        , MPI_Comm comm
#endif
)
{
    using namespace std::placeholders;
    std::function<double(uint64_t, uint64_t)> f = std::bind(f_elem, _1, _2, data);
    MSk::oracle::Parameters msk_params = convert_parameters(parameters);
    MSk::MosaicSkeleton<double> *msk_ptr{nullptr};
    MSk::C::Result result = MSk::C::build_basic<double, double, 3>(
            &msk_ptr, m, n, f, row_pts, col_pts, msk_params, n_jobs
#ifdef MSK_USE_MPI
            , comm
#endif
    );
    if (result != MSk::C::Result::OK) {
        return result;
    }
    set_matrix_data(msk_ptr, msk_matrix);
    return MSk::C::Result::OK;
}

int msk_build_basic_c3(mskMosaicMatrix *msk_matrix,
                        mskParameters *parameters,
                        uint64_t m, uint64_t n,
                        float _Complex (*f_elem)(uint64_t, uint64_t, void *), void *data,
                        const float *row_pts, const float *col_pts,
                        uint64_t n_jobs
#ifdef MSK_USE_MPI
        , MPI_Comm comm
#endif
)
{
    using namespace std::placeholders;
    std::function<std::complex<float>(uint64_t, uint64_t)> f = std::bind(f_elem, _1, _2, data);
    MSk::oracle::Parameters msk_params = convert_parameters(parameters);
    MSk::MosaicSkeleton<std::complex<float>> *msk_ptr{nullptr};
    MSk::C::Result result = MSk::C::build_basic<std::complex<float>, float, 3>(
            &msk_ptr, m, n, f, row_pts, col_pts, msk_params, n_jobs
#ifdef MSK_USE_MPI
                                                                   , comm
#endif
                                                                   );
    if (result != MSk::C::Result::OK) {
        return result;
    }
    set_matrix_data(msk_ptr, msk_matrix);
    return MSk::C::Result::OK;
}

int msk_build_basic_z3(mskMosaicMatrix *msk_matrix,
                        mskParameters *parameters,
                        uint64_t m, uint64_t n,
                        double _Complex (*f_elem)(uint64_t, uint64_t, void *), void *data,
                        const double *row_pts, const double *col_pts,
                        uint64_t n_jobs
#ifdef MSK_USE_MPI
        , MPI_Comm comm
#endif
)
{
    using namespace std::placeholders;
    std::function<std::complex<double>(uint64_t, uint64_t)> f = std::bind(f_elem, _1, _2, data);
    MSk::oracle::Parameters msk_params = convert_parameters(parameters);
    MSk::MosaicSkeleton<std::complex<double>> *msk_ptr{nullptr};
    MSk::C::Result result = MSk::C::build_basic<std::complex<double>, double, 3>(
            &msk_ptr, m, n, f, row_pts, col_pts, msk_params, n_jobs
#ifdef MSK_USE_MPI
            , comm
#endif
    );
    if (result != MSk::C::Result::OK) {
        return result;
    }
    set_matrix_data(msk_ptr, msk_matrix);
    return MSk::C::Result::OK;
}


int msk_build_block_s3(mskMosaicMatrix *msk_matrix,
                    mskParameters *parameters,
                    uint64_t m, uint64_t m_blocks, const uint64_t *m_block_sizes,
                    uint64_t n, uint64_t n_blocks, const uint64_t *n_block_sizes,
                    float (*f_elem)(uint64_t, uint64_t, void *), void *data,
                    const float *row_pts, const float *col_pts,
                    uint64_t n_jobs
#ifdef MSK_USE_MPI
                    , MPI_Comm comm
#endif
)
{
    using namespace std::placeholders;
    std::function<float(uint64_t, uint64_t)> f = std::bind(f_elem, _1, _2, data);
    MSk::oracle::Parameters msk_params = convert_parameters(parameters);
    MSk::MosaicSkeleton<float> *msk_ptr{nullptr};
    MSk::C::Result result = MSk::C::build_block<float, float, 3>(
            &msk_ptr,
            m, m_blocks, m_block_sizes,
            n, n_blocks, n_block_sizes,
            f, row_pts, col_pts, msk_params, n_jobs
#ifdef MSK_USE_MPI
            , comm
#endif
    );
    if (result != MSk::C::Result::OK) {
        return result;
    }
    set_matrix_data(msk_ptr, msk_matrix);
    return MSk::C::Result::OK;
}

int msk_build_block_d3(mskMosaicMatrix *msk_matrix,
                        mskParameters *parameters,
                        uint64_t m, uint64_t m_blocks, const uint64_t *m_block_sizes,
                        uint64_t n, uint64_t n_blocks, const uint64_t *n_block_sizes,
                        double (*f_elem)(uint64_t, uint64_t, void *), void *data,
                        const double *row_pts, const double *col_pts,
                        uint64_t n_jobs
#ifdef MSK_USE_MPI
        , MPI_Comm comm
#endif
)
{
    using namespace std::placeholders;
    std::function<double(uint64_t, uint64_t)> f = std::bind(f_elem, _1, _2, data);
    MSk::oracle::Parameters msk_params = convert_parameters(parameters);
    MSk::MosaicSkeleton<double> *msk_ptr{nullptr};
    MSk::C::Result result = MSk::C::build_block<double, double, 3>(
            &msk_ptr,
            m, m_blocks, m_block_sizes,
            n, n_blocks, n_block_sizes,
            f, row_pts, col_pts, msk_params, n_jobs
#ifdef MSK_USE_MPI
            , comm
#endif
    );
    if (result != MSk::C::Result::OK) {
        return result;
    }
    set_matrix_data(msk_ptr, msk_matrix);
    return MSk::C::Result::OK;
}

int msk_build_block_c3(mskMosaicMatrix *msk_matrix,
                        mskParameters *parameters,
                        uint64_t m, uint64_t m_blocks, const uint64_t *m_block_sizes,
                        uint64_t n, uint64_t n_blocks, const uint64_t *n_block_sizes,
                        float _Complex (*f_elem)(uint64_t, uint64_t, void *), void *data,
                        const float *row_pts, const float *col_pts,
                        uint64_t n_jobs
#ifdef MSK_USE_MPI
        , MPI_Comm comm
#endif
)
{
    using namespace std::placeholders;
    std::function<std::complex<float>(uint64_t, uint64_t)> f = std::bind(f_elem, _1, _2, data);
    MSk::oracle::Parameters msk_params = convert_parameters(parameters);
    MSk::MosaicSkeleton<std::complex<float>> *msk_ptr{nullptr};
    MSk::C::Result result = MSk::C::build_block<std::complex<float>, float, 3>(
            &msk_ptr,
            m, m_blocks, m_block_sizes,
            n, n_blocks, n_block_sizes, f, row_pts, col_pts, msk_params, n_jobs
#ifdef MSK_USE_MPI
            , comm
#endif
    );
    if (result != MSk::C::Result::OK) {
        return result;
    }
    set_matrix_data(msk_ptr, msk_matrix);
    return MSk::C::Result::OK;
}

int msk_build_block_z3(mskMosaicMatrix *msk_matrix,
                        mskParameters *parameters,
                        uint64_t m, uint64_t m_blocks, const uint64_t *m_block_sizes,
                        uint64_t n, uint64_t n_blocks, const uint64_t *n_block_sizes,
                        double _Complex (*f_elem)(uint64_t, uint64_t, void *), void *data,
                        const double *row_pts, const double *col_pts,
                        uint64_t n_jobs
#ifdef MSK_USE_MPI
        , MPI_Comm comm
#endif
)
{
    using namespace std::placeholders;
    std::function<std::complex<double>(uint64_t, uint64_t)> f = std::bind(f_elem, _1, _2, data);
    MSk::oracle::Parameters msk_params = convert_parameters(parameters);
    MSk::MosaicSkeleton<std::complex<double>> *msk_ptr{nullptr};
    MSk::C::Result result = MSk::C::build_block<std::complex<double>, double, 3>(
            &msk_ptr, 
            m, m_blocks, m_block_sizes,
            n, n_blocks, n_block_sizes,
            f, row_pts, col_pts, msk_params, n_jobs
#ifdef MSK_USE_MPI
            , comm
#endif
    );
    if (result != MSk::C::Result::OK) {
        return result;
    }
    set_matrix_data(msk_ptr, msk_matrix);
    return MSk::C::Result::OK;
}


int msk_recompress_s3(mskMosaicMatrix *msk_matrix, double rel_tol) {
    return MSk::C::recompress<float>(msk_matrix->_matrix, rel_tol);
}

int msk_recompress_d3(mskMosaicMatrix *msk_matrix, double rel_tol) {
    return MSk::C::recompress<double>(msk_matrix->_matrix, rel_tol);
}

int msk_recompress_c3(mskMosaicMatrix *msk_matrix, double rel_tol) {
    return MSk::C::recompress<std::complex<float>>(msk_matrix->_matrix, rel_tol);
}

int msk_recompress_z3(mskMosaicMatrix *msk_matrix, double rel_tol) {
    return MSk::C::recompress<std::complex<double>>(msk_matrix->_matrix, rel_tol);
}


int msk_destroy_s3(mskMosaicMatrix *msk_matrix) {
    return MSk::C::destroy<float>(msk_matrix->_matrix);
}

int msk_destroy_d3(mskMosaicMatrix *msk_matrix) {
    return MSk::C::destroy<double>(msk_matrix->_matrix);
}

int msk_destroy_c3(mskMosaicMatrix *msk_matrix) {
    return MSk::C::destroy<std::complex<float>>(msk_matrix->_matrix);
}

int msk_destroy_z3(mskMosaicMatrix *msk_matrix) {
    return MSk::C::destroy<std::complex<double>>(msk_matrix->_matrix);
}


int msk_matvec_s3(const mskMosaicMatrix *msk_matrix, const char trans, const float *x, float *y) {
    return MSk::C::matvec<float>(msk_matrix->_matrix, trans, x, y);
}

int msk_matvec_d3(const mskMosaicMatrix *msk_matrix, const char trans, const double *x, double *y) {
    return MSk::C::matvec<double>(msk_matrix->_matrix, trans, x, y);
}

int msk_matvec_c3(const mskMosaicMatrix *msk_matrix, const char trans, const float _Complex *x, float _Complex *y) {
    return MSk::C::matvec<std::complex<float>>(msk_matrix->_matrix, trans,
            reinterpret_cast<const std::complex<float> *>(x),
            reinterpret_cast<std::complex<float> *>(y)
    );
}

int msk_matvec_z3(const mskMosaicMatrix *msk_matrix, const char trans, const double _Complex *x, double _Complex *y) {
    return MSk::C::matvec<std::complex<double>>(msk_matrix->_matrix, trans,
            reinterpret_cast<const std::complex<double> *>(x),
            reinterpret_cast<std::complex<double> *>(y)
    );
}


#ifdef MSK_USE_MPI
int msk_to_dense_s3(const mskMosaicMatrix *msk_matrix, float *ptr, uint64_t ldim, int root) {
    return MSk::C::to_dense<float>(msk_matrix->_matrix, ptr, ldim, root);
}

int msk_to_dense_d3(const mskMosaicMatrix *msk_matrix, double *ptr, uint64_t ldim, int root) {
    return MSk::C::to_dense<double>(msk_matrix->_matrix, ptr, ldim, root);
}

int msk_to_dense_c3(const mskMosaicMatrix *msk_matrix, float _Complex *ptr, uint64_t ldim, int root) {
    return MSk::C::to_dense<std::complex<float>>(msk_matrix->_matrix, reinterpret_cast<std::complex<float> *>(ptr), ldim, root);
}

int msk_to_dense_z3(const mskMosaicMatrix *msk_matrix, double _Complex *ptr, uint64_t ldim, int root) {
    return MSk::C::to_dense<std::complex<double>>(msk_matrix->_matrix, reinterpret_cast<std::complex<double> *>(ptr), ldim, root);
}
#else
int msk_to_dense_s3(const mskMosaicMatrix *msk_matrix, float *ptr, uint64_t ldim) {
    return MSk::C::to_dense<float>(msk_matrix->_matrix, ptr, ldim);
}

int msk_to_dense_d3(const mskMosaicMatrix *msk_matrix, double *ptr, uint64_t ldim) {
    return MSk::C::to_dense<double>(msk_matrix->_matrix, ptr, ldim);
}

int msk_to_dense_c3(const mskMosaicMatrix *msk_matrix, float _Complex *ptr, uint64_t ldim) {
    return MSk::C::to_dense<std::complex<float>>(msk_matrix->_matrix, reinterpret_cast<std::complex<float> *>(ptr), ldim);
}

int msk_to_dense_z3(const mskMosaicMatrix *msk_matrix, double _Complex *ptr, uint64_t ldim) {
    return MSk::C::to_dense<std::complex<double>>(msk_matrix->_matrix, reinterpret_cast<std::complex<double> *>(ptr), ldim);
}
#endif


#ifdef MSK_USE_MPI
int msk_permute_col_s3(const mskMosaicMatrix *msk_matrix, const float * x_in, float *x_out, bool global) {
    return MSk::C::permute_col<float>(msk_matrix->_matrix, x_in, x_out, global);
}

int msk_permute_col_d3(const mskMosaicMatrix *msk_matrix, const double *x_in, double *x_out, bool global) {
    return MSk::C::permute_col<double>(msk_matrix->_matrix, x_in, x_out, global);
}

int msk_permute_col_c3(const mskMosaicMatrix *msk_matrix, const float _Complex *x_in, float _Complex *x_out, bool global) {
    return MSk::C::permute_col<std::complex<float>>(msk_matrix->_matrix,
                                          reinterpret_cast<const std::complex<float> *>(x_in),
                                          reinterpret_cast<std::complex<float> *>(x_out),
                                          global);
}

int msk_permute_col_z3(const mskMosaicMatrix *msk_matrix, const double _Complex *x_in, double _Complex *x_out, bool global) {
    return MSk::C::permute_col<std::complex<double>>(msk_matrix->_matrix,
                                           reinterpret_cast<const std::complex<double> *>(x_in),
                                           reinterpret_cast<std::complex<double> *>(x_out),
                                           global);
}


int msk_permute_col_inplace_s3(const mskMosaicMatrix *msk_matrix, float *x, bool global) {
    return MSk::C::permute_col<float>(msk_matrix->_matrix, x, global);
}

int msk_permute_col_inplace_d3(const mskMosaicMatrix *msk_matrix, double *x, bool global) {
    return MSk::C::permute_col<double>(msk_matrix->_matrix, x, global);
}

int msk_permute_col_inplace_c3(const mskMosaicMatrix *msk_matrix, float _Complex *x, bool global) {
    return MSk::C::permute_col<std::complex<float>>(msk_matrix->_matrix, reinterpret_cast<std::complex<float> *>(x), global);
}

int msk_permute_col_inplace_z3(const mskMosaicMatrix *msk_matrix, double _Complex *x, bool global) {
    return MSk::C::permute_col<std::complex<double>>(msk_matrix->_matrix, reinterpret_cast<std::complex<double> *>(x), global);
}


int msk_unpermute_col_s3(const mskMosaicMatrix *msk_matrix, const float * x_in, float *x_out, bool global) {
    return MSk::C::unpermute_col<float>(msk_matrix->_matrix, x_in, x_out, global);
}

int msk_unpermute_col_d3(const mskMosaicMatrix *msk_matrix, const double *x_in, double *x_out, bool global) {
    return MSk::C::unpermute_col<double>(msk_matrix->_matrix, x_in, x_out, global);
}

int msk_unpermute_col_c3(const mskMosaicMatrix *msk_matrix, const float _Complex *x_in, float _Complex *x_out, bool global) {
    return MSk::C::unpermute_col<std::complex<float>>(msk_matrix->_matrix,
                                            reinterpret_cast<const std::complex<float> *>(x_in),
                                            reinterpret_cast<std::complex<float> *>(x_out),
                                            global);
}

int msk_unpermute_col_z3(const mskMosaicMatrix *msk_matrix, const double _Complex *x_in, double _Complex *x_out, bool global) {
    return MSk::C::unpermute_col<std::complex<double>>(msk_matrix->_matrix,
                                             reinterpret_cast<const std::complex<double> *>(x_in),
                                             reinterpret_cast<std::complex<double> *>(x_out),
                                             global);
}


int msk_unpermute_col_inplace_s3(const mskMosaicMatrix *msk_matrix, float *x, bool global) {
    return MSk::C::unpermute_col<float>(msk_matrix->_matrix, x, global);
}

int msk_unpermute_col_inplace_d3(const mskMosaicMatrix *msk_matrix, double *x, bool global) {
    return MSk::C::unpermute_col<double>(msk_matrix->_matrix, x, global);
}

int msk_unpermute_col_inplace_c3(const mskMosaicMatrix *msk_matrix, float _Complex *x, bool global) {
    return MSk::C::unpermute_col<std::complex<float>>(msk_matrix->_matrix, reinterpret_cast<std::complex<float> *>(x), global);
}

int msk_unpermute_col_inplace_z3(const mskMosaicMatrix *msk_matrix, double _Complex *x, bool global) {
    return MSk::C::unpermute_col<std::complex<double>>(msk_matrix->_matrix, reinterpret_cast<std::complex<double> *>(x), global);
}


int msk_permute_row_s3(const mskMosaicMatrix *msk_matrix, const float * x_in, float *x_out, bool global) {
    return MSk::C::permute_row<float>(msk_matrix->_matrix, x_in, x_out, global);
}

int msk_permute_row_d3(const mskMosaicMatrix *msk_matrix, const double *x_in, double *x_out, bool global) {
    return MSk::C::permute_row<double>(msk_matrix->_matrix, x_in, x_out, global);
}

int msk_permute_row_c3(const mskMosaicMatrix *msk_matrix, const float _Complex *x_in, float _Complex *x_out, bool global) {
    return MSk::C::permute_row<std::complex<float>>(msk_matrix->_matrix,
                                          reinterpret_cast<const std::complex<float> *>(x_in),
                                          reinterpret_cast<std::complex<float> *>(x_out),
                                          global);
}

int msk_permute_row_z3(const mskMosaicMatrix *msk_matrix, const double _Complex *x_in, double _Complex *x_out, bool global) {
    return MSk::C::permute_row<std::complex<double>>(msk_matrix->_matrix,
                                           reinterpret_cast<const std::complex<double> *>(x_in),
                                           reinterpret_cast<std::complex<double> *>(x_out),
                                           global);
}


int msk_permute_row_inplace_s3(const mskMosaicMatrix *msk_matrix, float *x, bool global) {
    return MSk::C::permute_row<float>(msk_matrix->_matrix, x, global);
}

int msk_permute_row_inplace_d3(const mskMosaicMatrix *msk_matrix, double *x, bool global) {
    return MSk::C::permute_row<double>(msk_matrix->_matrix, x, global);
}

int msk_permute_row_inplace_c3(const mskMosaicMatrix *msk_matrix, float _Complex *x, bool global) {
    return MSk::C::permute_row<std::complex<float>>(msk_matrix->_matrix, reinterpret_cast<std::complex<float> *>(x), global);
}

int msk_permute_row_inplace_z3(const mskMosaicMatrix *msk_matrix, double _Complex *x, bool global) {
    return MSk::C::permute_row<std::complex<double>>(msk_matrix->_matrix, reinterpret_cast<std::complex<double> *>(x), global);
}


int msk_unpermute_row_s3(const mskMosaicMatrix *msk_matrix, const float * x_in, float *x_out, bool global) {
    return MSk::C::unpermute_row<float>(msk_matrix->_matrix, x_in, x_out, global);
}

int msk_unpermute_row_d3(const mskMosaicMatrix *msk_matrix, const double *x_in, double *x_out, bool global) {
    return MSk::C::unpermute_row<double>(msk_matrix->_matrix, x_in, x_out, global);
}

int msk_unpermute_row_c3(const mskMosaicMatrix *msk_matrix, const float _Complex *x_in, float _Complex *x_out, bool global) {
    return MSk::C::unpermute_row<std::complex<float>>(msk_matrix->_matrix,
                                            reinterpret_cast<const std::complex<float> *>(x_in),
                                            reinterpret_cast<std::complex<float> *>(x_out),
                                            global);
}

int msk_unpermute_row_z3(const mskMosaicMatrix *msk_matrix, const double _Complex *x_in, double _Complex *x_out, bool global) {
    return MSk::C::unpermute_row<std::complex<double>>(msk_matrix->_matrix,
                                             reinterpret_cast<const std::complex<double> *>(x_in),
                                             reinterpret_cast<std::complex<double> *>(x_out),
                                             global);
}


int msk_unpermute_row_inplace_s3(const mskMosaicMatrix *msk_matrix, float *x, bool global) {
    return MSk::C::unpermute_row<float>(msk_matrix->_matrix, x, global);
}

int msk_unpermute_row_inplace_d3(const mskMosaicMatrix *msk_matrix, double *x, bool global) {
    return MSk::C::unpermute_row<double>(msk_matrix->_matrix, x, global);
}

int msk_unpermute_row_inplace_c3(const mskMosaicMatrix *msk_matrix, float _Complex *x, bool global) {
    return MSk::C::unpermute_row<std::complex<float>>(msk_matrix->_matrix, reinterpret_cast<std::complex<float> *>(x), global);
}

int msk_unpermute_row_inplace_z3(const mskMosaicMatrix *msk_matrix, double _Complex *x, bool global) {
    return MSk::C::unpermute_row<std::complex<double>>(msk_matrix->_matrix, reinterpret_cast<std::complex<double> *>(x), global);
}
#else
int msk_permute_col_s3(const mskMosaicMatrix *msk_matrix, const float * x_in, float *x_out) {
    return MSk::C::permute_col<float>(msk_matrix->_matrix, x_in, x_out);
}

int msk_permute_col_d3(const mskMosaicMatrix *msk_matrix, const double *x_in, double *x_out) {
    return MSk::C::permute_col<double>(msk_matrix->_matrix, x_in, x_out);
}

int msk_permute_col_c3(const mskMosaicMatrix *msk_matrix, const float _Complex *x_in, float _Complex *x_out) {
    return MSk::C::permute_col<std::complex<float>>(msk_matrix->_matrix,
            reinterpret_cast<const std::complex<float> *>(x_in),
            reinterpret_cast<std::complex<float> *>(x_out)
    );
}

int msk_permute_col_z3(const mskMosaicMatrix *msk_matrix, const double _Complex *x_in, double _Complex *x_out) {
    return MSk::C::permute_col<std::complex<double>>(msk_matrix->_matrix,
            reinterpret_cast<const std::complex<double> *>(x_in),
            reinterpret_cast<std::complex<double> *>(x_out)
    );
}


int msk_permute_col_inplace_s3(const mskMosaicMatrix *msk_matrix, float *x) {
    return MSk::C::permute_col<float>(msk_matrix->_matrix, x);
}

int msk_permute_col_inplace_d3(const mskMosaicMatrix *msk_matrix, double *x) {
    return MSk::C::permute_col<double>(msk_matrix->_matrix, x);
}

int msk_permute_col_inplace_c3(const mskMosaicMatrix *msk_matrix, float _Complex *x) {
    return MSk::C::permute_col<std::complex<float>>(msk_matrix->_matrix, reinterpret_cast<std::complex<float> *>(x));
}

int msk_permute_col_inplace_z3(const mskMosaicMatrix *msk_matrix, double _Complex *x) {
    return MSk::C::permute_col<std::complex<double>>(msk_matrix->_matrix, reinterpret_cast<std::complex<double> *>(x));
}


int msk_unpermute_col_s3(const mskMosaicMatrix *msk_matrix, const float * x_in, float *x_out) {
    return MSk::C::unpermute_col<float>(msk_matrix->_matrix, x_in, x_out);
}

int msk_unpermute_col_d3(const mskMosaicMatrix *msk_matrix, const double *x_in, double *x_out) {
    return MSk::C::unpermute_col<double>(msk_matrix->_matrix, x_in, x_out);
}

int msk_unpermute_col_c3(const mskMosaicMatrix *msk_matrix, const float _Complex *x_in, float _Complex *x_out) {
    return MSk::C::unpermute_col<std::complex<float>>(msk_matrix->_matrix,
                                            reinterpret_cast<const std::complex<float> *>(x_in),
                                            reinterpret_cast<std::complex<float> *>(x_out));
}

int msk_unpermute_col_z3(const mskMosaicMatrix *msk_matrix, const double _Complex *x_in, double _Complex *x_out) {
    return MSk::C::unpermute_col<std::complex<double>>(msk_matrix->_matrix,
                                             reinterpret_cast<const std::complex<double> *>(x_in),
                                             reinterpret_cast<std::complex<double> *>(x_out));
}


int msk_unpermute_col_inplace_s3(const mskMosaicMatrix *msk_matrix, float *x) {
    return MSk::C::unpermute_col<float>(msk_matrix->_matrix, x);
}

int msk_unpermute_col_inplace_d3(const mskMosaicMatrix *msk_matrix, double *x) {
    return MSk::C::unpermute_col<double>(msk_matrix->_matrix, x);
}

int msk_unpermute_col_inplace_c3(const mskMosaicMatrix *msk_matrix, float _Complex *x) {
    return MSk::C::unpermute_col<std::complex<float>>(msk_matrix->_matrix, reinterpret_cast<std::complex<float> *>(x));
}

int msk_unpermute_col_inplace_z3(const mskMosaicMatrix *msk_matrix, double _Complex *x) {
    return MSk::C::unpermute_col<std::complex<double>>(msk_matrix->_matrix, reinterpret_cast<std::complex<double> *>(x));
}


int msk_permute_row_s3(const mskMosaicMatrix *msk_matrix, const float * x_in, float *x_out) {
    return MSk::C::permute_row<float>(msk_matrix->_matrix, x_in, x_out);
}

int msk_permute_row_d3(const mskMosaicMatrix *msk_matrix, const double *x_in, double *x_out) {
    return MSk::C::permute_row<double>(msk_matrix->_matrix, x_in, x_out);
}

int msk_permute_row_c3(const mskMosaicMatrix *msk_matrix, const float _Complex *x_in, float _Complex *x_out) {
    return MSk::C::permute_row<std::complex<float>>(msk_matrix->_matrix,
                                          reinterpret_cast<const std::complex<float> *>(x_in),
                                          reinterpret_cast<std::complex<float> *>(x_out));
}

int msk_permute_row_z3(const mskMosaicMatrix *msk_matrix, const double _Complex *x_in, double _Complex *x_out) {
    return MSk::C::permute_row<std::complex<double>>(msk_matrix->_matrix,
                                           reinterpret_cast<const std::complex<double> *>(x_in),
                                           reinterpret_cast<std::complex<double> *>(x_out));
}


int msk_permute_row_inplace_s3(const mskMosaicMatrix *msk_matrix, float *x) {
    return MSk::C::permute_row<float>(msk_matrix->_matrix, x);
}

int msk_permute_row_inplace_d3(const mskMosaicMatrix *msk_matrix, double *x) {
    return MSk::C::permute_row<double>(msk_matrix->_matrix, x);
}

int msk_permute_row_inplace_c3(const mskMosaicMatrix *msk_matrix, float _Complex *x) {
    return MSk::C::permute_row<std::complex<float>>(msk_matrix->_matrix, reinterpret_cast<std::complex<float> *>(x));
}

int msk_permute_row_inplace_z3(const mskMosaicMatrix *msk_matrix, double _Complex *x) {
    return MSk::C::permute_row<std::complex<double>>(msk_matrix->_matrix, reinterpret_cast<std::complex<double> *>(x));
}


int msk_unpermute_row_s3(const mskMosaicMatrix *msk_matrix, const float * x_in, float *x_out) {
    return MSk::C::unpermute_row<float>(msk_matrix->_matrix, x_in, x_out);
}

int msk_unpermute_row_d3(const mskMosaicMatrix *msk_matrix, const double *x_in, double *x_out) {
    return MSk::C::unpermute_row<double>(msk_matrix->_matrix, x_in, x_out);
}

int msk_unpermute_row_c3(const mskMosaicMatrix *msk_matrix, const float _Complex *x_in, float _Complex *x_out) {
    return MSk::C::unpermute_row<std::complex<float>>(msk_matrix->_matrix,
                                            reinterpret_cast<const std::complex<float> *>(x_in),
                                            reinterpret_cast<std::complex<float> *>(x_out));
}

int msk_unpermute_row_z3(const mskMosaicMatrix *msk_matrix, const double _Complex *x_in, double _Complex *x_out) {
    return MSk::C::unpermute_row<std::complex<double>>(msk_matrix->_matrix,
                                             reinterpret_cast<const std::complex<double> *>(x_in),
                                             reinterpret_cast<std::complex<double> *>(x_out));
}


int msk_unpermute_row_inplace_s3(const mskMosaicMatrix *msk_matrix, float *x) {
    return MSk::C::unpermute_row<float>(msk_matrix->_matrix, x);
}

int msk_unpermute_row_inplace_d3(const mskMosaicMatrix *msk_matrix, double *x) {
    return MSk::C::unpermute_row<double>(msk_matrix->_matrix, x);
}

int msk_unpermute_row_inplace_c3(const mskMosaicMatrix *msk_matrix, float _Complex *x) {
    return MSk::C::unpermute_row<std::complex<float>>(msk_matrix->_matrix, reinterpret_cast<std::complex<float> *>(x));
}

int msk_unpermute_row_inplace_z3(const mskMosaicMatrix *msk_matrix, double _Complex *x) {
    return MSk::C::unpermute_row<std::complex<double>>(msk_matrix->_matrix, reinterpret_cast<std::complex<double> *>(x));
}
#endif

#ifdef MSK_USE_MPI
int msk_check_error_s3(mskMosaicMatrix *msk_matrix,
                        float (*f_elem)(uint64_t, uint64_t, void *), void *data,
                        bool relative, bool global, double *error) {
    using namespace std::placeholders;
    std::function<float(uint64_t, uint64_t)> f = std::bind(f_elem, _1, _2, data);
    return MSk::C::check_error<float>(msk_matrix->_matrix, f, relative, global, error);
}

int msk_check_error_d3(mskMosaicMatrix *msk_matrix,
                        double (*f_elem)(uint64_t, uint64_t, void *), void *data,
                        bool relative, bool global, double *error) {
    using namespace std::placeholders;
    std::function<double(uint64_t, uint64_t)> f = std::bind(f_elem, _1, _2, data);
    return MSk::C::check_error<double>(msk_matrix->_matrix, f, relative, global, error);
}

int msk_check_error_c3(mskMosaicMatrix *msk_matrix,
                        float _Complex (*f_elem)(uint64_t, uint64_t, void *), void *data,
                        bool relative, bool global, double *error) {
    using namespace std::placeholders;
    std::function<std::complex<float>(uint64_t, uint64_t)> f = std::bind(f_elem, _1, _2, data);
    return MSk::C::check_error<std::complex<float>>(msk_matrix->_matrix, f, relative, global, error);
}

int msk_check_error_z3(mskMosaicMatrix *msk_matrix,
                        double _Complex (*f_elem)(uint64_t, uint64_t, void *), void *data,
                        bool relative, bool global, double *error) {
    using namespace std::placeholders;
    std::function<std::complex<double>(uint64_t, uint64_t)> f = std::bind(f_elem, _1, _2, data);
    return MSk::C::check_error<std::complex<double>>(msk_matrix->_matrix, f, relative, global, error);
}

#else

int msk_check_error_s3(mskMosaicMatrix *msk_matrix,
                        float (*f_elem)(uint64_t, uint64_t, void *), void *data,
                        bool relative, double *error) {
    using namespace std::placeholders;
    std::function<float(uint64_t, uint64_t)> f = std::bind(f_elem, _1, _2, data);
    return MSk::C::check_error<float>(msk_matrix->_matrix, f, relative, error);
}

int msk_check_error_d3(mskMosaicMatrix *msk_matrix,
                        double (*f_elem)(uint64_t, uint64_t, void *), void *data,
                        bool relative, double *error) {
    using namespace std::placeholders;
    std::function<double(uint64_t, uint64_t)> f = std::bind(f_elem, _1, _2, data);
    return MSk::C::check_error<double>(msk_matrix->_matrix, f, relative, error);
}

int msk_check_error_c3(mskMosaicMatrix *msk_matrix,
                        float _Complex (*f_elem)(uint64_t, uint64_t, void *), void *data,
                        bool relative, double *error) {
    using namespace std::placeholders;
    std::function<std::complex<float>(uint64_t, uint64_t)> f = std::bind(f_elem, _1, _2, data);
    return MSk::C::check_error<std::complex<float>>(msk_matrix->_matrix, f, relative, error);
}

int msk_check_error_z3(mskMosaicMatrix *msk_matrix,
                        double _Complex (*f_elem)(uint64_t, uint64_t, void *), void *data,
                        bool relative, double *error) {
    using namespace std::placeholders;
    std::function<std::complex<double>(uint64_t, uint64_t)> f = std::bind(f_elem, _1, _2, data);
    return MSk::C::check_error<std::complex<double>>(msk_matrix->_matrix, f, relative, error);
}

#endif


int msk_approx_block_solve_s3(uint64_t m, uint64_t m_blocks, const uint64_t *m_block_sizes,
                               uint64_t n, uint64_t n_blocks, const uint64_t *n_block_sizes,
                               float (*f_elem)(uint64_t, uint64_t, void *), void *data,
                               const float *row_pts, const float *col_pts,
                               mskParameters *parameters,
                               uint64_t n_rhs, float *B, float *X,
                               mskSolveParameters *solve_params,
                               uint64_t n_jobs
#ifdef MSK_USE_MPI
                               , MPI_Comm comm = MPI_COMM_WORLD
#endif
) {
    using namespace std::placeholders;
    std::function<float(uint64_t, uint64_t)> f = std::bind(f_elem, _1, _2, data);
    MSk::oracle::Parameters msk_params = convert_parameters(parameters);
    MSk::solvers::linsys::Parameters<float> gmres_params(solve_params->rel_tol, solve_params->abs_tol, 0, solve_params->iter_limit,
                                          MSk::solvers::linsys::gmres::Orthogonalization::ClassicalGramSchmidt,
                                          MSk::solvers::linsys::gmres::Reorthogonalization::AsNeeded);
    return MSk::C::approx_block_solve<float, float, 3>(
            m, m_blocks, m_block_sizes,
            n, n_blocks, n_block_sizes,
            f, row_pts, col_pts, msk_params,
            n_rhs, B, X,
            gmres_params,
            n_jobs
#ifdef MSK_USE_MPI
            , comm
#endif
    );
}

int msk_approx_block_solve_d3(uint64_t m, uint64_t m_blocks, const uint64_t *m_block_sizes,
                               uint64_t n, uint64_t n_blocks, const uint64_t *n_block_sizes,
                               double (*f_elem)(uint64_t, uint64_t, void *), void *data,
                               const double *row_pts, const double *col_pts,
                               mskParameters *parameters,
                               uint64_t n_rhs, double *B, double *X,
                               mskSolveParameters *solve_params,
                               uint64_t n_jobs
#ifdef MSK_USE_MPI
                               , MPI_Comm comm = MPI_COMM_WORLD
#endif
) {
    using namespace std::placeholders;
    std::function<double(uint64_t, uint64_t)> f = std::bind(f_elem, _1, _2, data);
    MSk::oracle::Parameters msk_params = convert_parameters(parameters);
    MSk::solvers::linsys::Parameters<double> gmres_params(solve_params->rel_tol, solve_params->abs_tol, 0, solve_params->iter_limit,
                                           MSk::solvers::linsys::gmres::Orthogonalization::ClassicalGramSchmidt,
                                           MSk::solvers::linsys::gmres::Reorthogonalization::AsNeeded);
    return MSk::C::approx_block_solve<double, double, 3>(
            m, m_blocks, m_block_sizes,
            n, n_blocks, n_block_sizes,
            f, row_pts, col_pts, msk_params,
            n_rhs, B, X,
            gmres_params,
            n_jobs
#ifdef MSK_USE_MPI
            , comm
#endif
    );
}

int msk_approx_block_solve_c3(uint64_t m, uint64_t m_blocks, const uint64_t *m_block_sizes,
                               uint64_t n, uint64_t n_blocks, const uint64_t *n_block_sizes,
                               float _Complex (*f_elem)(uint64_t, uint64_t, void *), void *data,
                               const float *row_pts, const float *col_pts,
                               mskParameters *parameters,
                               uint64_t n_rhs, float _Complex *B, float _Complex *X,
                               mskSolveParameters *solve_params,
                               uint64_t n_jobs
#ifdef MSK_USE_MPI
                               , MPI_Comm comm = MPI_COMM_WORLD
#endif
) {
    using namespace std::placeholders;
    std::function<std::complex<float>(uint64_t, uint64_t)> f = std::bind(f_elem, _1, _2, data);
    MSk::oracle::Parameters msk_params = convert_parameters(parameters);
    MSk::solvers::linsys::Parameters<float> gmres_params(solve_params->rel_tol, solve_params->abs_tol, 0, solve_params->iter_limit,
                                          MSk::solvers::linsys::gmres::Orthogonalization::ClassicalGramSchmidt,
                                          MSk::solvers::linsys::gmres::Reorthogonalization::AsNeeded);
    return MSk::C::approx_block_solve<std::complex<float>, float, 3>(
            m, m_blocks, m_block_sizes,
            n, n_blocks, n_block_sizes,
            f, row_pts, col_pts, msk_params,
            n_rhs, reinterpret_cast<std::complex<float> *>(B), reinterpret_cast<std::complex<float> *>(X),
            gmres_params,
            n_jobs
#ifdef MSK_USE_MPI
            , comm
#endif
    );
}

int msk_approx_block_solve_z3(uint64_t m, uint64_t m_blocks, const uint64_t *m_block_sizes,
                               uint64_t n, uint64_t n_blocks, const uint64_t *n_block_sizes,
                               double _Complex (*f_elem)(uint64_t, uint64_t, void *), void *data,
                               const double *row_pts, const double *col_pts,
                               mskParameters *parameters,
                               uint64_t n_rhs, double _Complex *B, double _Complex *X,
                               mskSolveParameters *solve_params,
                               uint64_t n_jobs
#ifdef MSK_USE_MPI
                               , MPI_Comm comm = MPI_COMM_WORLD
#endif
) {
    using namespace std::placeholders;
    std::function<std::complex<double>(uint64_t, uint64_t)> f = std::bind(f_elem, _1, _2, data);
    MSk::oracle::Parameters msk_params = convert_parameters(parameters);
    MSk::solvers::linsys::Parameters<double> gmres_params(solve_params->rel_tol, solve_params->abs_tol, 0, solve_params->iter_limit,
                                           MSk::solvers::linsys::gmres::Orthogonalization::ClassicalGramSchmidt,
                                           MSk::solvers::linsys::gmres::Reorthogonalization::AsNeeded);
    return MSk::C::approx_block_solve<std::complex<double>, double, 3>(
            m, m_blocks, m_block_sizes,
            n, n_blocks, n_block_sizes,
            f, row_pts, col_pts, msk_params,
            n_rhs, reinterpret_cast<std::complex<double> *>(B), reinterpret_cast<std::complex<double> *>(X),
            gmres_params,
            n_jobs
#ifdef MSK_USE_MPI
            , comm
#endif
    );
}
