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
// Created by bulat_v on 10.11.2021.
//


/**
 * @file
 *
 * @brief Contains the orthogonalization routine for multiple-rhs GMRES algorithm defined in GMRES.hpp.
 */


#ifndef MOSAIC_SKELETON_SOLVERS_LINSYS_ORTHOGONALIZATION_HPP
#define MOSAIC_SKELETON_SOLVERS_LINSYS_ORTHOGONALIZATION_HPP


#ifdef MSK_USE_MPI
#include "cpp_mpi/MPI.hpp"
#endif

#include "parameters.hpp"


namespace MSk {
namespace solvers {
namespace linsys {
namespace gmres {

template<typename DataType>
void orthogonalize(
#ifdef MSK_USE_MPI
            MPI_Comm comm,
#endif
            uint64_t basis_size, DataType *Basis, uint64_t basis_ld, DataType *q_new,
            DataType *coefs, int64_t inc_coefs = 1,
            gmres::Orthogonalization orth_method = gmres::Orthogonalization::ModifiedGramSchmidt,
            gmres::Reorthogonalization reorth = gmres::Reorthogonalization::Never,
            uint64_t reorth_number = 0) {

    using RealType = decltype(std::abs(DataType()));
    uint64_t orth_number = 1 + (reorth == gmres::Reorthogonalization::Always ? reorth_number : 0);
    const RealType dot_tolerance = 5 * std::numeric_limits<RealType>::epsilon();

    for (uint64_t i = 0; i < orth_number; i++) {
        if (orth_method == gmres::Orthogonalization::ClassicalGramSchmidt) {
            SmartArray<DataType> Bc_q(basis_size);
            BLAS::gemv('C', basis_ld, basis_size, DataType{1.0}, Basis, basis_ld,
                       q_new, 1, DataType{0.0}, Bc_q.data(), 1);
#ifdef MSK_USE_MPI
            MPI::Allreduce_inplace_sum(Bc_q.data(), basis_size, comm);
#endif
            BLAS::gemv('N', basis_ld, basis_size, DataType{-1.0}, Basis, basis_ld,
                       Bc_q.data(), 1, DataType{1.0}, q_new, 1);
            BLAS::axpy(basis_size, DataType{1.0}, Bc_q.data(), 1, coefs, inc_coefs);
        } else { //Modified Gram-Schmidt
            const DataType *bj = Basis;
            DataType *coefj = coefs;
            for (uint64_t j = 0; j < basis_size; j++) {
                DataType alpha = BLAS::dotc(basis_ld, bj, 1, q_new, 1);
#ifdef MSK_USE_MPI
                MPI::Allreduce_inplace_sum(&alpha, 1, comm);
#endif
                *coefj += alpha;
                BLAS::axpy(basis_ld, -alpha, bj, 1, q_new, 1);
                coefj += inc_coefs;
                bj += basis_ld;
            }
        }
    }
    if (reorth == gmres::Reorthogonalization::AsNeeded) {
        // Check if reorthgonalization is necessary

        //Compute scalar product
        SmartArray<DataType> Bc_q(basis_size);
        BLAS::gemv('C', basis_ld, basis_size, DataType{1.0}, Basis, basis_ld,
                   q_new, 1, DataType{0.0}, Bc_q.data(), 1);
#ifdef MSK_USE_MPI
        MPI::Allreduce_inplace_sum(Bc_q.data(), basis_size, comm);
#endif
        // Find max dot product
        RealType max_dot = RealType{0.0};
        for (uint64_t i = 0; i < basis_size; ++i) {
            max_dot = std::max(std::abs(Bc_q[i]), max_dot);
        }
        RealType q_norm = BLAS::nrm2(basis_ld, q_new, 1);
#ifdef MSK_USE_MPI
        q_norm *= q_norm;
        MPI::Allreduce_inplace_sum(&q_norm, 1, comm);
        q_norm = std::sqrt(q_norm);
#endif
        // If maximal scalar product is too large perform reorthogonalization
        if (max_dot > dot_tolerance * q_norm) {
            if (orth_method == gmres::Orthogonalization::ClassicalGramSchmidt) {
                BLAS::gemv('N', basis_ld, basis_size, DataType{-1.0}, Basis, basis_ld,
                           Bc_q.data(), 1, DataType{1.0}, q_new, 1);
                BLAS::axpy(basis_size, DataType{1.0}, Bc_q.data(), 1, coefs, inc_coefs);
            } else { // Modified Gram-Schmidt
                const DataType *bj = Basis;
                DataType *coefj = coefs;
                for (uint64_t j = 0; j < basis_size; ++j) {
                    DataType alpha = BLAS::dotc(basis_ld, bj, 1, q_new, 1);
#ifdef MSK_USE_MPI
                    MPI::Allreduce_inplace_sum(&alpha, 1, comm);
#endif
                    *coefj += alpha;
                    BLAS::axpy(basis_ld, -alpha, bj, 1, q_new, 1);
                    bj += basis_ld;
                    coefj += inc_coefs;
                }
            }
        }
    }
}


}
}
}
}

#endif //MOSAIC_SKELETON_SOLVERS_LINSYS_ORTHOGONALIZATION_HPP
