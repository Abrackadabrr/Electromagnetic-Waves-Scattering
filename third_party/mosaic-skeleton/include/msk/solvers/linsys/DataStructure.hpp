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
 * @brief Contains internal data structures for multiple-rhs GMRES algorithm defined in GMRES.hpp.
 */


#ifndef MOSAIC_SKELETON_SOLVERS_LINSYS_DATASTRUCTURE_HPP
#define MOSAIC_SKELETON_SOLVERS_LINSYS_DATASTRUCTURE_HPP


#ifdef MSK_USE_MPI
#include "cpp_mpi/MPI.hpp"
#endif

#include "parameters.hpp"
#include "Orthogonalization.hpp"


namespace MSk {
namespace solvers {
namespace linsys {
namespace gmres {


struct BasisStat {
    uint64_t basis_size;
    uint64_t sol_basis_size;
    uint64_t rhs_num;
};

template <class DataType>
struct Rotation {
    using RealType = decltype(std::abs(DataType()));
    RealType c{1.0};
    DataType s{0.0};
    uint64_t row1{0};
    uint64_t row2{0};
};

template<typename DataType>
class RhsData;

template<typename DataType>
class BasisData {
public:
    friend class RhsData<DataType>;
    using RealType = decltype(std::abs(DataType()));

private:
    DataType zero{0.0};

    // Q is matrix containing Basis of Krylov subspace
    SmartArray<DataType> Q;
    // R is in the packed form -- for order n is vector of size n (n + 1) / 2
    SmartArray<DataType> R;
    // G^H R is Hessenberg matrix, G is set of Givens rotation
    SmartArray<Rotation<DataType>> G;
    // Matrix C contains coefs (before q_last's) of P vectors - to be multiplied by A (in packed form)
    // C: A Q[1..n+p-1] C = Q[1..n+p] H
    SmartArray<DataType> c_last;
    SmartArray<DataType> C_special;
    // Position (offset) of last vector of rhs' Basis in total Basis matrix
    SmartArray<uint64_t> q_lasts_pos;

    const uint64_t row_number;
    uint64_t full_row_number;
    uint64_t basis_size, sol_basis_size, rhs_num;
    uint64_t old_G_size, G_size, C_special_size;
    gmres::Orthogonalization orth_method;
    gmres::Reorthogonalization reorth_sched;
    uint64_t reorth_number;

#ifdef MSK_USE_MPI
    MPI_Comm matrix_comm;
    int proc_id{0}, proc_size{0};
#endif

public:

    BasisData(uint64_t row_number, uint64_t full_row_number,
#ifdef MSK_USE_MPI
              MPI_Comm comm,
#endif
              gmres::Orthogonalization orth_m = gmres::Orthogonalization::ModifiedGramSchmidt,
              gmres::Reorthogonalization reorth_s = gmres::Reorthogonalization::Never, uint64_t reorth_num = 0) :
            row_number{row_number}, full_row_number{full_row_number},
            basis_size{0}, sol_basis_size{0}, rhs_num{static_cast<uint64_t>(-1)},
            old_G_size{0}, G_size{0}, C_special_size{0},
            orth_method{orth_m}, reorth_sched{reorth_s}, reorth_number{reorth_num}
#ifdef MSK_USE_MPI
    , matrix_comm{comm}
#endif
    {
#ifdef MSK_USE_MPI
        MPI_Comm_size(matrix_comm, &proc_size);
        MPI_Comm_rank(matrix_comm, &proc_id);
        //full_row_number = row_number * proc_size;
#endif
    }
    ~BasisData() {}

    uint64_t get_size() const {
        return Q.byte_size() + R.byte_size() + G.byte_size() + c_last.byte_size() +
               C_special.byte_size() + q_lasts_pos.byte_size();
    }

    uint64_t get_rhs_num() const {
        return rhs_num;
    }

    DataType *get_basis() {
        return Q.data();
    }

    uint64_t get_reorth_number() const {
        return reorth_number;
    }

    gmres::Orthogonalization get_orth_method() const {
        return orth_method;
    }

    gmres::Reorthogonalization get_reorth() const {
        return reorth_sched;
    }

    void get_last_dir(DataType *p) {
        BLAS::gemv('N', row_number, basis_size, DataType{1.0}, Q.data(), row_number, c_last.data(), 1,
                   DataType{0.0}, p, 1);
    }

    DataType *get_space() {
        Q.reserve((basis_size + 1) * row_number);
        return Q.data() + row_number * basis_size;
    }

    DataType *update() {
        DataType *q_last = Q.data() + basis_size * row_number;
        RealType R_size = ((sol_basis_size + 1) * (sol_basis_size + 2)) / 2 + (basis_size - sol_basis_size);
        R.reserve(R_size);
        // Current column of R. Its size == basis_size + 1 a.k.a. new basis_size
        DataType *r_last = R.data() + (sol_basis_size * (sol_basis_size + 1)) / 2;

        // Orthogonalize current q to previous ones, put coefficients in r_last
        BLAS::copy(basis_size, &zero, 0, r_last, 1);
        orthogonalize(
#ifdef MSK_USE_MPI
                matrix_comm,
#endif
                      basis_size, Q.data(), row_number, q_last, r_last, 1,
                      orth_method, reorth_sched, reorth_number);

        // Norm of new q after orthogonalization and normalizing it
        RealType q_norm = BLAS::nrm2(row_number, q_last, 1);

#ifdef MSK_USE_MPI
        q_norm *= q_norm;
        MPI::Allreduce_inplace_sum(&q_norm, 1, matrix_comm);
        q_norm = std::sqrt(q_norm);
#endif
        r_last[basis_size] = q_norm;

        BLAS::rscal(row_number, RealType{1.0} / q_norm, q_last, 1);
        // New basis vector is finally computed
        ++basis_size, ++sol_basis_size, ++q_lasts_pos[rhs_num];

        if (basis_size > full_row_number) {
            basis_size = full_row_number;
            --q_lasts_pos[rhs_num];
        }

        // Apply previously computed Givens rotations to r_last
        for (uint64_t i = 0; i < G_size; ++i) {
            BLAS::rot(1, r_last + G[i].row1, 1, r_last + G[i].row2, 1, G[i].c, G[i].s);
        }
        // Compute new rotations
        old_G_size = G_size;
        G_size += basis_size - sol_basis_size;
        G.reserve(G_size);

        uint64_t base_row = basis_size - 2; // sol_basis_size - 1;
        uint64_t next_row = basis_size - 1; // base_row + 1;
        for (uint64_t i = old_G_size; i < G_size; ++i, --next_row, --base_row) {
            Rotation<DataType> &rot = G[i];
            LAPACK::lartg(r_last[base_row], r_last[next_row], rot.c, rot.s, r_last[base_row]);
            r_last[next_row] = 0;
            rot.row1 = base_row;
            rot.row2 = next_row;
        }

        c_last.reserve(basis_size);
        BLAS::copy(basis_size, &zero, 0, c_last.data(), 1);
        c_last[sol_basis_size - 1] = DataType{1.0};
        for (uint64_t i = G_size; i > 0; --i) {
            BLAS::rot(1, c_last.data() + G[i - 1].row1, 1, c_last.data() + G[i - 1].row2, 1, G[i - 1].c, -G[i - 1].s);
        }
        return q_last;
    }

    DataType *expand(DataType *res_orth, RealType orth_norm, DataType *Qc_r0, DataType *GQc_r0, uint64_t ldim) {
        ++rhs_num;
        q_lasts_pos.reserve(rhs_num + 1);
        Q.reserve((basis_size + 1) * row_number);
        DataType *q_last = Q.data() + basis_size * row_number;

        BLAS::copy(row_number, res_orth, 1, q_last, 1);
        q_lasts_pos[rhs_num] = basis_size;

        gmres::orthogonalize(
#ifdef MSK_USE_MPI
                matrix_comm,
#endif
                             basis_size, Q.data(), row_number, q_last, Qc_r0, ldim,
                             orth_method, reorth_sched, reorth_number);
        (void) orth_norm;
        RealType q_norm = BLAS::nrm2(row_number, q_last, 1); // orth_norm;
#ifdef MSK_USE_MPI
        q_norm *= q_norm;
        MPI::Allreduce_inplace_sum(&q_norm, 1, matrix_comm);
        q_norm = std::sqrt(q_norm);
#endif
        Qc_r0[basis_size * ldim] = q_norm;
        BLAS::rscal(row_number, 1 / q_norm, q_last, 1);
        ++basis_size;
        // renew GQc_r0
        BLAS::copy(basis_size, Qc_r0, ldim, GQc_r0, ldim);
        for (uint64_t i = 0; i < G_size; ++i) {
            Rotation<DataType> &rot = G[i];
            BLAS::rot(1, GQc_r0 + rot.row1 * ldim, 1, GQc_r0 + rot.row2 * ldim, 1, rot.c, rot.s);
        }

        if (basis_size > full_row_number) {
            basis_size = full_row_number;
            --q_lasts_pos[rhs_num];
        }

        //RealType residual_norm = BLAS::nrm2(basis_size - sol_basis_size, GQc_r0 + sol_basis_size * ldim, ldim);

        c_last.reserve(basis_size);

        // Save coefs to C
        BLAS::copy(sol_basis_size, &zero, 0, c_last.data(), 1);
        BLAS::copy(basis_size - sol_basis_size, GQc_r0 + sol_basis_size * ldim, ldim, c_last.data() + sol_basis_size, 1);
        // Normalize column of C
        RealType c_norm = BLAS::nrm2(basis_size - sol_basis_size, c_last.data() + sol_basis_size, 1);
        BLAS::rscal(basis_size - sol_basis_size, RealType{1.0} / c_norm, c_last.data() + sol_basis_size, 1);
        // Apply previous rotations
        for (uint64_t i = G_size; i > 0; --i) {
            BLAS::rot(1, c_last.data() + G[i - 1].row1, 1, c_last.data() + G[i - 1].row2, 1, G[i - 1].c, -G[i - 1].s);
        }
        C_special.reserve(C_special_size + basis_size);
        BLAS::copy(basis_size, c_last.data(), 1, C_special.data() + C_special_size, 1);
        C_special_size += basis_size;

        return q_last;
    }

    void restart() { // cold restart
        basis_size = sol_basis_size = 0;
        rhs_num = -1;
        old_G_size = G_size = C_special_size = 0;
        q_lasts_pos[0] = 0;
    }

    std::pair<uint64_t, Rotation<DataType> *> get_last_rotations() {
        return std::make_pair<uint64_t, Rotation<DataType> *>(basis_size - sol_basis_size, G.data() + old_G_size);
    }

    std::pair<uint64_t, Rotation<DataType> *> get_rotations() {
        return std::make_pair<uint64_t, Rotation<DataType> *>(static_cast<uint64_t>(G_size), G.data());
    }

    uint64_t get_basis_size() const {
        return basis_size;
    }

    uint64_t get_sol_basis_size() const {
        return sol_basis_size;
    }
};


template<typename DataType>
class RhsData {
public:
    using RealType = decltype(std::abs(DataType()));

private:
    const DataType zero{0.0};

    // not swapped, stored as given
    DataType *x;
    uint64_t ldx;

    // always stored as [ lead, [uncoverged], [converged] ]
    DataType *res_orth;
    uint64_t ldres;

    const uint64_t row_number;
    //uint64_t full_row_number;
    uint64_t unconv_cnt, /*basis_size, sol_basis_size,*/ Nrhs;

    SmartArray<uint64_t> id_at_pos;
    SmartArray<uint64_t> pos_of_id;

    SmartArray<RealType> init_res_norm;
    SmartArray<RealType> aim_res_norm;
    SmartArray<RealType> res_norm;
    SmartArray<RealType> orth_norm;

    SmartArray<DataType> GQc_R0_T;
    SmartArray<DataType> Qc_R0_T;

    BasisData<DataType> basis;
#ifdef MSK_USE_MPI
    MPI_Comm matrix_comm;
    int proc_id{0}, proc_size{0};
#endif
private:
    void swap_rhs(uint64_t i, uint64_t j) {
        if (i == j) {
            return;
        }
        uint64_t basis_size = basis.basis_size;
        BLAS::swap(row_number, res_orth + i * ldres, 1, res_orth + j * ldres, 1);
        BLAS::swap(basis_size, GQc_R0_T.data() + i, Nrhs, GQc_R0_T.data() + j, Nrhs);
        BLAS::swap(basis_size, Qc_R0_T.data() + i, Nrhs, Qc_R0_T.data() + j, Nrhs);
        std::swap(init_res_norm[i], init_res_norm[j]);
        std::swap(aim_res_norm[i], aim_res_norm[j]);
        std::swap(res_norm[i], res_norm[j]);
        std::swap(orth_norm[i], orth_norm[j]);
        std::swap(id_at_pos[i], id_at_pos[j]);
        pos_of_id[id_at_pos[i]] = i;
        pos_of_id[id_at_pos[j]] = j;
    }

    // apply new rotations to unconv rhs'
    void new_rotations(uint64_t rot_cnt, Rotation<DataType> *rots) {
        for (uint64_t i = 0; i < rot_cnt; ++i) {
            Rotation<DataType> &rot = rots[i];
            BLAS::rot(unconv_cnt, GQc_R0_T.data() + rot.row1 * Nrhs, 1,
                      GQc_R0_T.data() + rot.row2 * Nrhs, 1, rot.c, rot.s);
        }
    }

    bool build(uint64_t pos_start, uint64_t to_solve = 1u, bool force_sol = false) {
        uint64_t basis_size = basis.basis_size;
        uint64_t sol_basis_size = basis.sol_basis_size;
        if (basis_size == 0) {
            return false;
        }

        uint64_t cur_rhs = basis.get_rhs_num();
        //            logger->debug("Build: order = {}, basis_size = {}, sol_basis_size = {}", cur_rhs, basis_size, sol_basis_size);
        SmartArray<DataType> sys_coefs(sol_basis_size * to_solve); // Solutions of triangular systems
        for (uint64_t i = 0; i < to_solve; ++i) {
            BLAS::copy(sol_basis_size, GQc_R0_T.data() + pos_start + i, Nrhs, sys_coefs.data() + i * sol_basis_size, 1);
        }

        uint64_t info{0};
        LAPACK::tptrs('U', 'N', 'N', sol_basis_size, to_solve, basis.R.data(), sys_coefs.data(), sol_basis_size, info);

        // Rebuild coefs for Q: z = C y
        uint64_t image_basis_size = std::max(basis_size - 1, sol_basis_size);
        SmartArray<DataType> solution_coefs(image_basis_size * to_solve);
        // Faster solution computation via C_special
        BLAS::copy(image_basis_size * to_solve, &zero, 0, solution_coefs.data(), 1);
        uint64_t c_pos{0}, c_len{1}, c_off{0};

        for (uint64_t i = 0; i <= cur_rhs; ++i) {
            BLAS::geru(c_len, to_solve, DataType{1.0},
                       basis.C_special.data() + c_off, 1,
                       sys_coefs.data() + c_pos, sol_basis_size,
                       solution_coefs.data(), image_basis_size);
            c_pos = basis.q_lasts_pos[i] - i;
            c_off += c_len;
            c_len = std::min(basis.q_lasts_pos[i] + 2, basis_size);
        }

        for (uint64_t i = 0; i < basis.G_size - (basis_size - sol_basis_size); ++i) {
            auto &rot = basis.G[i];
            BLAS::rot(to_solve, solution_coefs.data() + rot.row1, image_basis_size,
                      solution_coefs.data() + rot.row2, image_basis_size, rot.c, rot.s);
        }
        for (uint64_t i = 0; i < to_solve; ++i) {
            BLAS::axpy(sol_basis_size - 1, DataType{1.0}, sys_coefs.data() + 1 + i * sol_basis_size, 1,
                       solution_coefs.data() + i * image_basis_size, 1);
        }
        for (uint64_t i = 0, y_off = basis.q_lasts_pos[0]; i < cur_rhs; ++i) {
            BLAS::axpy(to_solve, DataType{-1.0}, sys_coefs.data() + y_off, sol_basis_size,
                       solution_coefs.data() + y_off - 1, image_basis_size);
            y_off += basis.q_lasts_pos[i + 1] - basis.q_lasts_pos[i] - 1;
        }

        for (uint64_t i = basis.G_size - (basis_size - sol_basis_size); i > 0; --i) {
            auto &rot = basis.G[i - 1];
            BLAS::rot(to_solve, solution_coefs.data()+ rot.row1, image_basis_size,
                      solution_coefs.data() + rot.row2, image_basis_size, rot.c, -rot.s);
        }

        // Gain x = x_0 + Qz = x_0 + QCy
        SmartArray<DataType> cur_x(row_number * to_solve);
        for (uint64_t i = 0; i < to_solve; ++i) {
            BLAS::copy(row_number, x + id_at_pos[pos_start + i] * ldx, 1, cur_x.data() + i * row_number, 1);
        }

        BLAS::gemm('N', 'N', row_number, to_solve, image_basis_size, DataType{1.0}, basis.Q.data(), row_number,
                   solution_coefs.data(), image_basis_size, DataType{1.0}, cur_x.data(), row_number);

        SmartArray<DataType> Hy(to_solve * basis_size);
        //assert(basis_size == sol_basis_size + (cur_rhs + 1) && "Basis sizes mismatch");
        if (pos_start == 0 && to_solve == Nrhs) {
            BLAS::copy(to_solve * sol_basis_size, GQc_R0_T.data(), 1, Hy.data(), 1);
            BLAS::copy(to_solve * (basis_size - sol_basis_size), &zero, 0, Hy.data() + to_solve * sol_basis_size, 1);
        } else {
            for (uint64_t i = 0; i < to_solve; ++i) {
                BLAS::copy(sol_basis_size, GQc_R0_T.data() + pos_start + i, Nrhs, Hy.data() + i, to_solve);
            }
            BLAS::copy(to_solve * (basis_size - sol_basis_size), &zero, 0, Hy.data() + sol_basis_size * to_solve, 1);
        }
        for (uint64_t i = basis.G_size; i > 0; --i) {
            auto rot = basis.G[i - 1];
            BLAS::rot(to_solve, Hy.data()+ rot.row1 * to_solve, 1, Hy.data() + rot.row2 * to_solve, 1, rot.c, -rot.s);
        }

        // Calculate r0
        if (pos_start == 0 && to_solve == Nrhs) {
            BLAS::axpy(to_solve * basis_size, DataType{-1.0}, Qc_R0_T.data(), 1, Hy.data(), 1);
        } else {
            for (uint64_t i = 0; i < to_solve; ++i) {
                BLAS::axpy(basis_size, DataType{-1.0}, Qc_R0_T.data() + pos_start + i, Nrhs, Hy.data(), to_solve);
            }
        }

        DataType *res = res_orth + pos_start * ldres;
        BLAS::gemm('N', 'T', row_number, to_solve, basis_size, DataType{-1.0},
                   basis.Q.data(), row_number, Hy.data(), to_solve, DataType{1.0}, res, ldres);

        for (uint64_t i = 0; i < to_solve; ++i) {
            res_norm[pos_start + i] = BLAS::nrm2(row_number, res + i * ldres, 1);
        }

        bool converged{false};
        for (uint64_t i = 0; i < to_solve; ++i) {
            converged = res_norm[pos_start + i] <= aim_res_norm[pos_start + i];
            if (converged || force_sol) {
                BLAS::copy(row_number, cur_x.data() + i * row_number, 1, x + id_at_pos[pos_start + i] * ldx, 1);
            }
        }
        return true;
    }

public:
    RhsData(uint64_t rhs_cnt, uint64_t row_num, uint64_t full_row_num, DataType *X0, uint64_t ldX,
            DataType *res, uint64_t ldr, Parameters<RealType> params
#ifdef MSK_USE_MPI
            , MPI_Comm comm
#endif
            ) :
            x{X0}, ldx{ldX}, res_orth{res}, ldres{ldr},
            row_number{row_num}, unconv_cnt{rhs_cnt}, Nrhs{rhs_cnt},
            id_at_pos{rhs_cnt}, pos_of_id{rhs_cnt},
            init_res_norm{rhs_cnt}, aim_res_norm{rhs_cnt}, res_norm{rhs_cnt}, orth_norm{rhs_cnt},
            GQc_R0_T{rhs_cnt}, Qc_R0_T{rhs_cnt},
            basis{row_num, full_row_num,
#ifdef MSK_USE_MPI
                  comm,
#endif
                  params.get_orthogonalization_method(), params.get_reorthogonalization_schedule(),
                  (params.get_reorthogonalization_schedule() == gmres::Reorthogonalization::Always) ? 1u : 0}
#ifdef MSK_USE_MPI
    , matrix_comm{comm}
#endif
    {

        RealType abs_tol = params.get_absolute_tolerance(), rel_tol = params.get_relative_tolerance();
#ifdef MSK_USE_MPI
        MPI_Comm_size(matrix_comm, &proc_size);
        MPI_Comm_rank(matrix_comm, &proc_id);
        //full_row_number = row_number * proc_size;

        for (uint64_t i = 0; i < Nrhs; ++i) {
            id_at_pos[i] = pos_of_id[i] = i;
            orth_norm[i] = BLAS::nrm2(row_number, res_orth + i * ldres, 1);
            orth_norm[i] *= orth_norm[i];
        }
        MPI::Allreduce_inplace_sum(orth_norm.data(), Nrhs, matrix_comm);
        for (uint64_t i = 0; i < Nrhs; ++i) {
            orth_norm[i] = std::sqrt(orth_norm[i]);
            init_res_norm[i] = res_norm[i] = orth_norm[i];
            aim_res_norm[i] = init_res_norm[i] * rel_tol + abs_tol;
        }
#else
        for (uint64_t i = 0; i < Nrhs; ++i) {
            id_at_pos[i] = pos_of_id[i] = i;
            init_res_norm[i] = res_norm[i] = orth_norm[i] = BLAS::nrm2(row_number, res_orth + i * ldres, 1);
            aim_res_norm[i] = init_res_norm[i] * rel_tol + abs_tol;
        }
#endif
    }

    ~RhsData() {
        for (uint64_t i = 0; i < Nrhs; ++i) {
            swap_rhs(i, pos_of_id[i]);
        }
    }

    uint64_t get_size() const {
        return  basis.get_size() + GQc_R0_T.byte_size() + Qc_R0_T.byte_size() +
                Nrhs * (2 * sizeof(uint64_t) + 4 * sizeof(RealType));
    }

    BasisStat get_basis_stat () const {
        return BasisStat{basis.get_basis_size(), basis.get_sol_basis_size(), basis.get_rhs_num()};
    }

    void get_last_dir(DataType *p) {
        basis.get_last_dir(p);
    }

    DataType *get_vector_space() {
        return basis.get_space();
    }

    // Orth-ze all unconv rhs' (except leading) to new basis vector
    void new_basis_vector(uint64_t orthogonalizations_number, bool is_r0 = false) {
        uint64_t basis_size = basis.basis_size;

        Qc_R0_T.reserve(Nrhs * (basis_size + 1));
        GQc_R0_T.reserve(Nrhs * (basis_size + 1));

        DataType *Qc_r0_last = Qc_R0_T.data() + Nrhs * basis_size;
        BLAS::copy(Nrhs, &zero, 0, Qc_r0_last, 1);
        DataType *q{nullptr};
        if (is_r0) {
            q = basis.expand(res_orth, orth_norm[0], Qc_R0_T.data(), GQc_R0_T.data(), Nrhs);
            Qc_r0_last[0] = orth_norm[0];
            orth_norm[0] = 0;
            BLAS::copy(row_number, &zero, 0, res_orth, 1);
        } else {
            q = basis.update();
        }
        basis_size = basis.basis_size;

        SmartArray<DataType> coefs(unconv_cnt - 1);
        for (uint64_t i = 0; i < orthogonalizations_number; i++) {
            BLAS::gemv('C', row_number, unconv_cnt - 1, DataType{1.0}, res_orth + 1 * ldres, ldres,
                       q, 1, DataType{0.0}, coefs.data(), 1);
            LAPACK::lacgv(unconv_cnt - 1, coefs.data(), 1);
#ifdef MSK_USE_MPI
            MPI::Allreduce_inplace_sum(coefs.data(), unconv_cnt - 1, matrix_comm);
#endif
            BLAS::geru(row_number, unconv_cnt - 1, DataType{-1.0}, q, 1, coefs.data(), 1,
                       res_orth + 1 * ldres, ldres);
            BLAS::axpy(unconv_cnt - 1, DataType{1.0}, coefs.data(), 1, Qc_r0_last + 1, 1);
        }
#ifdef MSK_USE_MPI
        for (uint64_t j = 1; j < unconv_cnt; ++j) {
            orth_norm[j] = BLAS::nrm2(row_number, res_orth + j * ldres, 1);
            orth_norm[j] *= orth_norm[j];
        }
        MPI::Allreduce_inplace_sum(orth_norm.data() + 1, unconv_cnt - 1, matrix_comm);
        for (uint64_t j = 1; j < unconv_cnt; ++j) {
            orth_norm[j] = std::sqrt(orth_norm[j]);
        }
#else
        for (uint64_t j = 1; j < unconv_cnt; ++j) {
            orth_norm[j] = BLAS::nrm2(row_number, res_orth + j * ldres, 1);
        }
#endif
        BLAS::copy(Nrhs, Qc_r0_last, 1, GQc_R0_T.data() + Nrhs * (basis_size - 1), 1);

        if (!is_r0) {
            auto rots = basis.get_last_rotations();
            new_rotations(rots.first, rots.second);
        }
    }

    // Swap: (lead, ..., new, ...) -> (new, ..., lead, ...)
    // Returns true if leading changed
    bool check_update(bool force_sol = false) {
        uint64_t basis_size = basis.basis_size;
        uint64_t sol_basis_size = basis.sol_basis_size;
        RealType max_ratio = 0;
        uint64_t max_ratio_pos = Nrhs;

        // check all preconverged
        bool lead_upd{false};
        for (uint64_t i = unconv_cnt; i > 0; --i) {
            uint64_t pos = i - 1;

            DataType *GQc_r0 = GQc_R0_T.data() + pos;
            DataType *Qc_r0 = Qc_R0_T.data() + pos;
            RealType lst_sq_err = BLAS::nrm2(basis_size - sol_basis_size, GQc_r0 + sol_basis_size * Nrhs, Nrhs);
            res_norm[pos] = std::sqrt(lst_sq_err * lst_sq_err + orth_norm[pos] * orth_norm[pos]);
            bool conv{false};
            if ((res_norm[pos] <= aim_res_norm[pos]) || force_sol) {
                // If converged by approx residual norm,
                // then recalculate it and double check
                // No orthogonalization, just reorth if needed
                gmres::orthogonalize(
#ifdef MSK_USE_MPI
                        matrix_comm,
#endif
                                     basis_size, basis.Q.data(), row_number, res_orth + pos * ldres,
                                     Qc_r0, Nrhs, basis.orth_method, gmres::Reorthogonalization::AsNeeded, -1);
                BLAS::copy(basis_size, Qc_r0, Nrhs, GQc_r0, Nrhs);
                for (uint64_t j = 0; j < basis.G_size; ++j) {
                    Rotation<DataType> &rot = basis.G[j];
                    BLAS::rot(1, GQc_r0 + rot.row1 * Nrhs, 1, GQc_r0 + rot.row2 * Nrhs, 1, rot.c, rot.s);
                }
                lst_sq_err = BLAS::nrm2(basis_size - sol_basis_size, GQc_r0 + sol_basis_size * Nrhs, Nrhs);
                orth_norm[pos] = BLAS::nrm2(row_number, res_orth + pos * ldres, 1);
#ifdef MSK_USE_MPI
                orth_norm[pos] *= orth_norm[pos];
                MPI::Allreduce_inplace_sum(&orth_norm[pos], 1, matrix_comm);
                res_norm[pos] = std::sqrt(lst_sq_err * lst_sq_err + orth_norm[pos]);
#else
                res_norm[pos] = std::sqrt(lst_sq_err * lst_sq_err + orth_norm[pos] * orth_norm[pos]);
#endif
                if ((res_norm[pos] <= aim_res_norm[pos]) || force_sol) {
                    BLAS::copy(basis_size - sol_basis_size, &zero, 0, GQc_r0 + sol_basis_size * Nrhs, Nrhs);
                    if (res_norm[pos] <= aim_res_norm[pos]) {
                        conv = true;
                    }
                }
            }
            if (conv) {
                if (pos == unconv_cnt - 1) {
                    // if convergence continues from tail, just drop the last one
                    --unconv_cnt;
                }
                if (pos == 0) { // Current working residual has converged
                    lead_upd = true;
                }
            } else {
                // This residual could become the next working
                RealType ratio = res_norm[pos] / init_res_norm[pos];
                if (ratio > max_ratio) {
                    max_ratio_pos = pos;
                    max_ratio = ratio;
                }
            }
        }

        // All rhs converged
        if (unconv_cnt == 0 || force_sol) {
            build(0, Nrhs, force_sol);
            return true;
        }
        // Now position (unconv_cnt - 1) marks the last unconverged rhs
        // Put the new working rhs
        if (lead_upd) {
            swap_rhs(0, max_ratio_pos);
            GQc_R0_T.reserve((basis_size + 1) * Nrhs);
            Qc_R0_T.reserve((basis_size + 1) * Nrhs);
        }
        // Check convergence of rhs in the whole list, 
        // not only from its tail
        // There is at least 1 unconverged rhs,
        // so unconv_cnt cannot become zero anyway
        while (unconv_cnt > 0 && res_norm[unconv_cnt - 1] <= aim_res_norm[unconv_cnt - 1]) {
            --unconv_cnt;
        }
        for (uint64_t i = 0; i < unconv_cnt - 1; ++i) {
            if (res_norm[i] > aim_res_norm[i]) {
                continue;
            }
            swap_rhs(i, unconv_cnt - 1);
            do {
                --unconv_cnt;
            } while (res_norm[unconv_cnt - 1] <= aim_res_norm[unconv_cnt - 1]);
        }
        return lead_upd;
    }

    void restart() { // cold restart
        check_update(true);
        orth_norm[0] = BLAS::nrm2(row_number, res_orth, 1);
#ifdef MSK_USE_MPI
        orth_norm[0] *= orth_norm[0];
        MPI::Allreduce_inplace_sum(orth_norm.data(), 1, matrix_comm);
        orth_norm[0] = std::sqrt(orth_norm[0]);
#endif
        basis.restart();
    }

    DataType *get_res_by_id(uint64_t id) {
        return res_orth + pos_of_id[id] * ldres;
    }

    DataType *get_res_by_pos(uint64_t pos)  {
        return res_orth + pos * ldres;
    }

    DataType *get_init_res_by_id(uint64_t id, DataType *result) const {
        return get_init_res_by_pos(pos_of_id[id], basis, result);
    }

    DataType *get_init_res_by_pos(uint64_t pos, DataType *result) const {
        uint64_t basis_size = basis.basis_size;
        if (pos == 0) {
            BLAS::gemv('N', row_number, basis_size, DataType{1.0}, basis.Q.data(), row_number,
                       Qc_R0_T.data() + pos, Nrhs, DataType{0.0}, result, 1);
        } else {
            BLAS::copy(row_number, res_orth + pos * ldres, 1, result, 1);
            BLAS::gemv('N', row_number, basis_size, DataType{1.0}, basis.Q.data(), row_number,
                       Qc_R0_T.data() + pos, Nrhs, DataType{1.0}, result, 1);
        }
        return result;
    }

    DataType *get_GQc_r0_start_by_id(uint64_t id) {
        return GQc_R0_T.data() + pos_of_id[id];
    }

    DataType *get_GQc_r0_start_by_pos(uint64_t pos) {
        return GQc_R0_T.data() + pos;
    }

    uint64_t get_pos_by_id(uint64_t id) const {
        return pos_of_id[id];
    }

    uint64_t get_id_by_pos(uint64_t pos) const {
        return id_at_pos[pos];
    }

    RealType get_init_res_norm_by_id(uint64_t id) const {
        return init_res_norm[pos_of_id[id]];
    }

    RealType get_init_res_norm_by_pos(uint64_t pos) const {
        return init_res_norm[pos];
    }

    RealType get_aim_res_norm_by_id(uint64_t id) const {
        return aim_res_norm[pos_of_id[id]];
    }

    RealType get_aim_res_norm_by_pos(uint64_t pos) const {
        return aim_res_norm[pos];
    }

    RealType get_res_norm_by_id(uint64_t id) const {
        return res_norm[pos_of_id[id]];
    }

    RealType get_res_norm_by_pos(uint64_t pos) const {
        return res_norm[pos];
    }

    RealType get_orth_res_norm_by_id(uint64_t id) const {
        return orth_norm[pos_of_id[id]];
    }

    RealType get_orth_res_norm_by_pos(uint64_t pos) const {
        return orth_norm[pos];
    }

    uint64_t get_unconv_cnt() const {
        return unconv_cnt;
    }

    uint64_t get_conv_cnt() const {
        return Nrhs - unconv_cnt;
    }

    uint64_t get_basis_size() const {
        return basis.basis_size;
    }
};

}
}
}
}

#endif //MOSAIC_SKELETON_SOLVERS_LINSYS_DATASTRUCTURE_HPP
