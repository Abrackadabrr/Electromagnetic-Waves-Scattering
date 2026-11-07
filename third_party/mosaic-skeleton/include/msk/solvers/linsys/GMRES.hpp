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


/**
 * @file
 *
 * @brief Defines basic GMRES method and its parameters.
 */

#ifndef MOSAIC_SKELETON_SOLVERS_LINSYS_GMRES_HPP
#define MOSAIC_SKELETON_SOLVERS_LINSYS_GMRES_HPP


#include <complex>
#include <iostream>
#include <fstream>
#include <limits>

#include "cpp_blas/BLAS.hpp"
#include "cpp_lapack/LAPACK.hpp"
#include "smart_array/smart_array.hpp"

#ifdef MSK_USE_MPI
#include "cpp_mpi/MPI.hpp"
#endif

#include "parameters.hpp"
#include "DataStructure.hpp"


namespace MSk {
namespace solvers {
namespace linsys {

/**
 * @brief Solve a system of linear equations using GMRES method.
 * @tparam Matrix Abstract matrix type, required to provide matvec method and size getters.
 * @tparam DataType Matrix element type.
 * @param[in] A Matrix object.
 * @param[in] b Right-hand side.
 * @param[in, out] x Initial guess on input, solution on output.
 * @param parameters Solving parameters.
 */
template <class Matrix, class DataType>
void GMRES(const Matrix &A, const DataType *b, DataType *x, Parameters<decltype(std::real(DataType()))> parameters = Parameters<decltype(std::abs(DataType()))>())
{
    if (A.get_row_number() != A.get_column_number())
    {
        throw std::runtime_error("Number of matrix rows and columns is not equal in GMRES!");
    }

    typedef decltype(std::abs(DataType())) RealType; //RealType is needed if DataType is complex

    const uint64_t row_number = A.get_row_number();

    // If reorthogonalization is always performed just do orthogonalization twice
    const uint64_t orthogonalizations_number = (parameters.get_reorthogonalization_schedule() == gmres::Reorthogonalization::Always) ? 2 : 1;

    // Tolerance for reorthogonalization
    const RealType dot_tolerance = 5 * std::numeric_limits<RealType>::epsilon();


    std::cout << "Starting GMRES for matrix of order " << row_number << std::endl;
    std::cout << "Relative tolerance = " << parameters.get_relative_tolerance() << std::endl;
    std::cout << "Absolute tolerance = " << parameters.get_absolute_tolerance() << std::endl;
    if (parameters.get_orthogonalization_method() == gmres::Orthogonalization::ClassicalGramSchmidt)
    {
        std::cout << "Using classical Gram-Schmidt method for orthogonalization" << std::endl;
    }
    else
    {
        std::cout << "Using modified Gram-Schmidt method for orthogonalization" << std::endl;
    }

    if (parameters.get_reorthogonalization_schedule() == gmres::Reorthogonalization::Never)
    {
        std::cout << "No reorthogonalizations are performed" << std::endl;
    }
    else if (parameters.get_reorthogonalization_schedule() == gmres::Reorthogonalization::AsNeeded)
    {
        std::cout << "Reorthogonalizations are performed as needed" << std::endl;
    }
    else
    {
        std::cout << "Reorthogonalizations are always performed" << std::endl;
    }

    // Q is matrix containing basis of Krylov subspace
    SmartArray<DataType> Q(row_number);

    // G^H R is Hessenberg matrix, G is set of Givens rotation, k-th rotation is for rows k and k + 1
    // R is in the packed form -- for order n is vector of size n (n + 1) / 2
    SmartArray<DataType> R;
    SmartArray<std::pair<RealType, DataType> > G;

    // G Q^H ^ r_0 = \|r_0\|_2 G e_0, where r_0 is initial residual
    SmartArray<DataType> GQCr(1);


    //Initial residual computation
    BLAS::copy(row_number, b, 1, Q.data(), 1);//Copy of rhs
    if (BLAS::nrm2(row_number, x, 1) > 0) //If x nonzero, subtract Ax
    {
        SmartArray<DataType> Ax(row_number);
        A.matvec(x, Ax.data());
        BLAS::axpy(row_number, DataType{-1.0}, Ax.data(), 1, Q.data(), 1);
    }

    // Initial residual norm
    RealType initial_residual_norm = BLAS::nrm2(row_number, Q.data(), 1);
    if (initial_residual_norm == 0) //Fast exit
    {
        return;
    }

    RealType residual_norm = initial_residual_norm;
    uint64_t iteration = 0;
    do
    {
        // First iteration after restart
        uint64_t restart_iteration = 0;
        // No rotations so far, projection vector is just \|r_0\|_2 e_0
        GQCr[0] = residual_norm;
        // Computation of q_0
        BLAS::rscal(row_number, RealType{1.0} / residual_norm, Q.data(), 1);

        while (residual_norm > initial_residual_norm * parameters.get_relative_tolerance() + parameters.get_absolute_tolerance()
                && (parameters.get_restart() == 0 || restart_iteration < parameters.get_restart()))
        {
            iteration++;
            restart_iteration++;

            // Size checks and resizes
            uint64_t Q_size = (restart_iteration + 1) * row_number;
            uint64_t R_size = (restart_iteration + 1) * restart_iteration / 2;
            uint64_t G_size = restart_iteration;
            uint64_t GQCr_size = restart_iteration + 1;
            if (Q_size > Q.size())
            {
                Q.resize(Q_size);
            }
            if (R_size > R.size())
            {
                R.resize(R_size);
            }
            if (G_size > G.size())
            {
                G.resize(G_size);
            }
            if (GQCr_size > GQCr.size())
            {
                GQCr.resize(GQCr_size);
            }

            // New element of projection vector haven't participated in rotations yet, so it is 0
            GQCr[restart_iteration] = 0;

            // Current column of R
            DataType *r = R.data() + (restart_iteration - 1) * restart_iteration / 2;
            for (uint64_t i = 0; i < restart_iteration; i++)
            {
                r[i] = 0; //put zero
            }

            // Current q is lat column of Q
            DataType *q = Q.data() + row_number * restart_iteration;

            //Multiply previous q by matrix and put result in current q
            A.matvec(q - row_number, q);

            // Orthogonalize current q to previous ones, put coefficients in r
            for (uint64_t i = 0; i < orthogonalizations_number; i++)
            {
                if (parameters.get_orthogonalization_method() == gmres::Orthogonalization::ClassicalGramSchmidt) // Classical Gram-Schmidt method
                {
                    SmartArray<DataType> QCq(restart_iteration);
                    BLAS::gemv('C', row_number, restart_iteration, DataType{1.0}, Q.data(), row_number, q, 1, DataType{0.0}, QCq.data(), 1);
                    BLAS::axpy(restart_iteration, DataType{1.0}, QCq.data(), 1, r, 1);
                    BLAS::gemv('N', row_number, restart_iteration, DataType{-1.0}, Q.data(), row_number, QCq.data(), 1, DataType{1.0}, q, 1);
                }
                else //Modified Gram-Schmidt
                {
                    const DataType *qj = Q.data();
                    for (uint64_t j = 0; j < restart_iteration; j++)
                    {
                        DataType alpha = BLAS::dotc(row_number, qj, 1, q, 1);
                        r[j] += alpha;
                        //                                std::cout << "alpha[" << j << "] = " << alpha << std::endl;
                        BLAS::axpy(row_number, -alpha, qj, 1, q, 1);
                        qj += row_number;
                    }
                }
            }
            if (parameters.get_reorthogonalization_schedule() == gmres::Reorthogonalization::AsNeeded)
            {
                // Check if reortagonalization nececcary
                SmartArray<DataType> QCq(restart_iteration);
                //Compute scalar product
                BLAS::gemv('C', row_number, restart_iteration, DataType{1.0}, Q.data(), row_number, q, 1, DataType{0.0}, QCq.data(), 1);

                // Find max
                RealType max_dot = RealType{0.0};
                for (uint64_t i = 0; i < restart_iteration; i++)
                {
                    max_dot = std::max(std::abs(QCq[i]), max_dot);
                }
                RealType q_norm = BLAS::nrm2(row_number, q, 1);

                // If maximal scalar product is too large perform reothogonalization
                if (max_dot > dot_tolerance * q_norm)
                {
                    if (parameters.get_orthogonalization_method() == gmres::Orthogonalization::ClassicalGramSchmidt)
                    {
                        BLAS::axpy(restart_iteration, DataType{1.0}, QCq.data(), 1, r, 1);
                        BLAS::gemv('N', row_number, restart_iteration, DataType{-1.0}, Q.data(), row_number, QCq.data(), 1, DataType{1.0}, q, 1);
                    }
                    else
                    {
                        const DataType *qj = Q.data();
                        for (uint64_t j = 0; j < restart_iteration; j++)
                        {
                            DataType alpha = BLAS::dotc(row_number, qj, 1, q, 1);
                            r[j] += alpha;
                            BLAS::axpy(row_number, -alpha, qj, 1, q, 1);
                            qj += row_number;
                        }
                    }
                }
            }

            // Norm of current q after orthogonalization
            RealType q_norm = BLAS::nrm2(row_number, q, 1);

            // Normalize q
            BLAS::rscal(row_number, RealType{1.0} / q_norm, q, 1);

            // Apply previously computed Givens rotations to r
            for (uint64_t i = 0; i < restart_iteration - 1; i++)
            {
                BLAS::rot(1, r + i, 1, r + i + 1, 1, G[i].first, G[i].second);
            }

            // Compute new rotation -- with last element of r and norm of q after orthogonalization
            LAPACK::lartg(r[restart_iteration - 1], DataType(q_norm), G[restart_iteration - 1].first, G[restart_iteration - 1].second, r[restart_iteration - 1]);

            // Apply new rotation to projection vector
            BLAS::rot(1, GQCr.data() + restart_iteration - 1, 1, GQCr.data() + restart_iteration, 1, G[restart_iteration - 1].first, G[restart_iteration - 1].second);

            // Residual norm is equal to absolute value last component of projection vector
            residual_norm = std::abs(GQCr[restart_iteration]);
        }

        // R^{-1} by projection vector
        BLAS::tpsv('U', 'N', 'N', restart_iteration, R.data(), GQCr.data(), 1);

        // x += Q R^{-1} G Q^H r_0
        BLAS::gemv('N', row_number, restart_iteration, DataType{1.0}, Q.data(), row_number, GQCr.data(), 1, DataType{1.0}, x, 1);

        // Compute real residual
        BLAS::copy(row_number, b, 1, Q.data(), 1); // copy vector to Q
        if (BLAS::nrm2(row_number, x, 1) > 0) // Subtract Ax
        {
            SmartArray<DataType> Ax(row_number);
            A.matvec(x, Ax.data());
            BLAS::axpy(row_number, DataType{-1.0}, Ax.data(), 1, Q.data(), 1);
        }

        //Compute real residual norm
        residual_norm = BLAS::nrm2(row_number, Q.data(), 1);

    }
    while (residual_norm > initial_residual_norm * parameters.get_relative_tolerance() + parameters.get_absolute_tolerance());
    std::ofstream stat;
    stat.open("stats.txt", std::ios_base::app);
    stat << iteration + 1 << " basis,\t";
    stat.close();
}


/**
     * @brief Solve a system of linear equations with multiple right-hand sides using GMRES method.
     *
     * @tparam Matrix Abstract matrix type, required to provide matvec method and size getters.
     * @tparam DataType Matrix element type.
     * @param[in] A Matrix object.
     * @param[in] Nrhs Number of right-hand sides.
     * @param[in, out] B Array of right-hand sides on input, residuals on output.
     * @param[in] ldB Leading dimension of RHS array.
     * @param[in, out] X Array of initial guess vectors on input, solutions on output.
     * @param[in] ldX Leading dimension of solutions array.
     * @param[in] parameters Solving parametes.
     * @return Error/success code.
     */
template <class Matrix, class DataType>
uint64_t GMRES(const Matrix &A, const uint64_t Nrhs, DataType *B, uint64_t ldB, DataType *X, uint64_t ldX,
               Parameters<decltype(std::real(DataType()))> parameters =
               Parameters<decltype(std::abs(DataType()))>()) {

    int id{0}, size{1};
#ifdef MSK_USE_MPI
    MPI_Comm comm = A.get_comm();
    MPI_Comm_rank(comm, &id);
    MPI_Comm_size(comm, &size);
#endif

    if (Nrhs == 0) {
        throw std::runtime_error("Nonpositive number of rhs in GMRES!");
    }

    const uint64_t row_number = A.get_row_number();
    //const uint64_t full_row_number = size * row_number;
    const uint64_t full_row_number = A.get_global_row_number();
    //        if (full_row_number != A.get_column_number()) {
    //            throw std::runtime_error("Number of matrix rows and columns is not equal in GMRES!");
    //        }
    using RealType = decltype(std::abs(DataType())); //RealType is needed if DataType is complex

    // If reorthogonalization is always performed just do orthogonalization twice
    uint64_t orth_num = 1;
    if (parameters.get_reorthogonalization_schedule() == gmres::Reorthogonalization::Always) {
        orth_num = 2;
    }

    double matvec_time{0}, new_basis_vec_time{0}, check_upd_time{0};

    if (id == 0) {
        std::cout << "Starting GMRES for matrix of order " << full_row_number << " with " << Nrhs << " rhs\n";
        std::cout << "Relative tolerance = " << parameters.get_relative_tolerance() << std::endl;
        std::cout << "Absolute tolerance = " << parameters.get_absolute_tolerance() << std::endl;
        if (parameters.get_orthogonalization_method() == gmres::Orthogonalization::ClassicalGramSchmidt) {
            std::cout << "Using classical Gram-Schmidt method for orthogonalization\n";
        } else {
            std::cout << "Using modified Gram-Schmidt method for orthogonalization\n";
        }

        if (parameters.get_reorthogonalization_schedule() == gmres::Reorthogonalization::Never) {
            std::cout << "No reorthogonalizations are performed\n";
        } else if (parameters.get_reorthogonalization_schedule() == gmres::Reorthogonalization::AsNeeded) {
            std::cout << "Reorthogonalizations are performed as needed\n";
        } else {
            std::cout << "Reorthogonalizations are always performed\n";
        }
    }


    // Precomputations of residuals
    SmartArray<RealType> x0_norm_sq(Nrhs);
    for (uint64_t i = 0; i < Nrhs; ++i) {
        x0_norm_sq[i] = BLAS::nrm2(row_number, X + i * ldX, 1);
        x0_norm_sq[i] *= x0_norm_sq[i];
    }
#ifdef MSK_USE_MPI
    MPI::Allreduce_inplace_sum(x0_norm_sq.data(), Nrhs, A.get_comm());
#endif
    SmartArray<DataType> Ax(row_number);


    for (uint64_t i = 0; i < Nrhs; ++i) {
        DataType *b = B + i * ldB;
        DataType *x = X + i * ldX;
        if (x0_norm_sq[i] > 0) { //If x nonzero, subtract Ax
            auto mv_t_start = std::chrono::steady_clock::now();
            A.matvec(x, Ax.data());
            auto mv_t_end = std::chrono::steady_clock::now();
            std::chrono::duration<double> t_delta = mv_t_end - mv_t_start;
            matvec_time += t_delta.count();
            BLAS::axpy(row_number, DataType{-1.0}, Ax.data(), 1, b, 1);
        }
    }

    gmres::RhsData<DataType> Rhs(Nrhs, row_number, full_row_number, X, ldX, B, ldB, parameters
#ifdef MSK_USE_MPI
                          , comm
#endif
                          );

    const DataType zero{0.0};
    auto t_start = std::chrono::steady_clock::now();

    Rhs.new_basis_vector(orth_num, true);

    bool choose_rhs{false};
    uint64_t iterations{0}, max_iter{parameters.get_iterations_limit()};
    if (max_iter == 0) {
        max_iter = -1;
    }
    gmres::BasisStat stat = Rhs.get_basis_stat();
    auto out_stat = [&]() {
        std::cout << "It: " << iterations;
        std::cout << "; rhs: ord = " << stat.rhs_num << ", id = " << Rhs.get_id_by_pos(0);
        std::cout << "; res = " << Rhs.get_res_norm_by_pos(0) << ", aim = " << Rhs.get_aim_res_norm_by_pos(0);
        std::cout << "; unconv = " << Rhs.get_unconv_cnt() << std::endl;
    };

    while (stat.sol_basis_size < full_row_number) {
        if (Rhs.get_unconv_cnt() == 0) {
            break;
        }
        if (iterations >= max_iter) {
            Rhs.check_update(true);
            out_stat();
            break;
        }

        if (choose_rhs) {
            auto nbv_start = std::chrono::steady_clock::now();
            Rhs.new_basis_vector(orth_num, true);
            auto nbv_end = std::chrono::steady_clock::now();
            std::chrono::duration<double> t_nbv = nbv_end - nbv_start;
            new_basis_vec_time += t_nbv.count();
        }

        SmartArray<DataType> p(row_number);
        Rhs.get_last_dir(p.data());

        auto mv_t_start = std::chrono::steady_clock::now();
        DataType *y = Rhs.get_vector_space();
        BLAS::copy(row_number, &zero, 0, y, 1);
        A.matvec(p.data(), Rhs.get_vector_space());
        auto mv_t_end = std::chrono::steady_clock::now();
        std::chrono::duration<double> t_matvec = mv_t_end - mv_t_start;
        matvec_time += t_matvec.count();

        auto nbv_start = std::chrono::steady_clock::now();
        Rhs.new_basis_vector(orth_num, false);
        auto nbv_end = std::chrono::steady_clock::now();
        std::chrono::duration<double> t_nbv = nbv_end - nbv_start;
        new_basis_vec_time += t_nbv.count();

        stat = Rhs.get_basis_stat();

        auto chupd_start = std::chrono::steady_clock::now();
        // If need restart
        if (parameters.get_restart() != 0 && stat.sol_basis_size >= parameters.get_restart()) {
            Rhs.restart();
            choose_rhs = true;
            if (id == 0) {
                out_stat();
                std::cout << "RESTART by Basis overflow" << std::endl;
            }
        } else {
            choose_rhs = Rhs.check_update(false);
            if (id == 0 && (iterations % 50 == 0 || Rhs.get_unconv_cnt() == 0 || choose_rhs)) {
                out_stat();
            }
        }
        auto chupd_end = std::chrono::steady_clock::now();
        std::chrono::duration<double> t_chupd = chupd_end - chupd_start;
        check_upd_time += t_chupd.count();
        ++iterations;
    }


    auto t_end = std::chrono::steady_clock::now();
    std::chrono::duration<double> t_delta = t_end - t_start;

    auto KMTb_mem = [](long double mem) {
        std::string suffix;
        if (mem < (1u << 30)) {
            mem /= (1u << 20);
            suffix = " MB";
        } else if (mem < (1lu << 40)) {
            mem /= (1u << 30);
            suffix = " GB";
        } else {
            mem /= (1lu << 40);
            suffix = " TB";
        }
        return std::to_string(mem) + suffix;
    };

    long double mem = Rhs.get_size() + Ax.size();
    long double mem_full = mem * size;
    if (id == 0) {
        stat = Rhs.get_basis_stat();
        std::cout << "aGMRES full time = " << t_delta.count() << "s, clean matvec time = " << matvec_time << "s\n";
        std::cout << "Converged: " << Rhs.get_conv_cnt() << " / " << Nrhs << std::endl;
        std::cout << "Full basis_size = " << stat.basis_size << ", solution basis_size = ";
        std::cout << stat.sol_basis_size << ", iterations_cnt = " << iterations << std::endl;
        std::cout << "Memory usage: per_proc = " << KMTb_mem(mem) << ", full = " << KMTb_mem(mem_full) << std::endl;
        std::cout << "New_bv_time = " << new_basis_vec_time << ", check_upd_time = " << check_upd_time << ", matvec_time = " << matvec_time << std::endl;
    }
    return Rhs.get_unconv_cnt();
}


}

}

}


#endif
