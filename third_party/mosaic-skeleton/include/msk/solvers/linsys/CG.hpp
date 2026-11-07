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
 * @brief Defines Conjugate Gradient method and its utilities.
 */


#ifndef MOSAIC_SKELETON_SOLVERS_LINSYS_CG_HPP
#define MOSAIC_SKELETON_SOLVERS_LINSYS_CG_HPP


#include <cstdint>
#include <limits>
#include <string>
#include <type_traits>
#include <complex>
#include <iostream>
#include <fstream>

#include "cpp_blas/BLAS.hpp"
#include "smart_array/smart_array.hpp"

#ifdef MSK_USE_MPI
#include "cpp_mpi/MPI.hpp"
#endif

#include "parameters.hpp"


namespace MSk {

namespace solvers {

namespace linsys {



//--------------------CG for symmetric matrix, A.matvec is required--------------------------------
/**
 * @brief Preconditioned Conjugate Gradient method.
 *
 * @details Used for iterative solving a system of linear equations with symmetric positive definite matrix.
 * Works also for multiprocessing systems using MPI.
 *
 * @tparam DataType Elements type.
 * @tparam Matrix Abstract matrix object, required toe have method matvec(x, y), performing matrix-vector product.
 * @tparam Precond Abstract preconditioner, required to have method apply(x, y), performing preconditioning.
 * @param[in] A Matrix.
 * @param[in] b Right-hand side.
 * @param[in, out] x Initial guess (input) and solution (output).
 * @param[in] parameters Solving parameters.
 * @param[in] precond Preconditioner.
 * @return Success/error code.
 */
template <class DataType, class Matrix, class Precond>
uint64_t CG(const Matrix &A, const DataType *b, DataType *x,
            Parameters<decltype(std::real(DataType()))> parameters,
            const Precond &precond
#ifdef MSK_USE_MPI
        , MPI_Comm comm = MPI_COMM_WORLD
#endif
) {
    using RealType = decltype(std::abs(DataType())); // RealType is needed if DataType is complex
    int proc_id{0};
#ifdef MSK_USE_MPI
    int n_proc{1};
    MPI_Comm_rank(comm, &proc_id);
    MPI_Comm_size(comm, &n_proc);
#endif
    const uint64_t row_number = A.get_row_number();
    uint64_t global_row_number = row_number;
#ifdef MSK_USE_MPI
    MPI_Allreduce(MPI_IN_PLACE, &global_row_number, 1, MPI_UINT64_T, MPI_SUM, comm);
#endif
    int verb = parameters.get_verbosity();
    if (proc_id == 0 && verb >= 0) {
        std::cout << "Starting CG for matrix of order " << global_row_number << std::endl;
        parameters.print_parameters();
    }

    SmartArray<DataType> r(row_number);	// r_i, current residual, orthogonal basis of Krylov space
    SmartArray<DataType> z(row_number); // z_i, preconditioned residual
    SmartArray<DataType> p(row_number);	// p_{i+1}, current direction, A-orthogonal basis of Krylov space
    SmartArray<DataType> Ap(row_number);// Ap_{i+1}, current direction multiplied by the matrix

    // Cumpute initial residual
    BLAS::copy(row_number, b, 1, r.data(), 1); // Copy of rhs
    if (BLAS::nrm2(row_number, x, 1) > 0) { // If x nonzero, subtract Ax
        A.matvec(x, Ap.data()); // Ap = Ax
        BLAS::axpy(row_number, DataType{-1.0}, Ap.data(), 1, r.data(), 1);
    }

    // Get initial residual norm
    RealType initial_residual_norm = BLAS::nrm2(row_number, r.data(), 1);
#ifdef MSK_USE_MPI
    initial_residual_norm = initial_residual_norm * initial_residual_norm;
    MPI::Allreduce_inplace_sum(&initial_residual_norm, 1, comm);
    initial_residual_norm = std::sqrt(initial_residual_norm);
#endif
    // Get initial residual norm according to parameters
    RealType target_residual = initial_residual_norm * parameters.get_relative_tolerance() + parameters.get_absolute_tolerance();
    if (initial_residual_norm == 0) { // Fast exit
        if (verb >= 0 && proc_id == 0) {
            std::cout << "CG converged at start. Residual norm = " << initial_residual_norm
                      << ", target residual norm = " << target_residual << std::endl;
        }
        return 0;
    }

    // Current residual norm, squared for precision of calculations
    RealType residual_norm = initial_residual_norm * initial_residual_norm;

    precond.apply(r.data(), z.data());
    BLAS::copy(row_number, z.data(), 1, p.data(), 1); // p_1 = z_0

    // Cycle until the iteration limit is reached
    uint64_t iteration_limit = parameters.get_iterations_limit();
    if (iteration_limit == 0) {
        iteration_limit = global_row_number;
    }

    for (uint64_t iteration = 1; iteration <= iteration_limit; iteration++) {
        // Matvec
        A.matvec(p.data(), Ap.data());	// calc Ap_i
        DataType pAp = std::real(BLAS::dotc(row_number, Ap.data(), 1, p.data(), 1)); // (p_i, Ap_i)
        DataType rCz = BLAS::dotc(row_number, r.data(), 1, z.data(), 1); // (z_i, r_i)
#ifdef MSK_USE_MPI
        DataType tmp_dot[2] = {pAp, rCz};
        MPI::Allreduce_inplace_sum(tmp_dot, 2, comm);
        pAp = tmp_dot[0];
        rCz = tmp_dot[1];
#endif

        // Update residual and solution
        // alpha = (z_i-1, r_i-1) / (p_i, Ap_i)
        DataType alpha = rCz / pAp;
        // x_i = x_i-1 + alpha * p_i
        BLAS::axpy(row_number, alpha, p.data(), 1, x, 1);
        // r_i = r_i-1 - alpha * Ap_i
        BLAS::axpy(row_number, -alpha, Ap.data(), 1, r.data(), 1);

        residual_norm = std::real(BLAS::dotc(row_number, r.data(), 1, r.data(), 1));
#ifdef MSK_USE_MPI
        MPI::Allreduce_inplace_sum(&residual_norm, 1, comm);
#endif
        if (proc_id == 0 && verb > 0 && iteration % verb == 0) {
            std::cout << "Iteration = " << iteration << ", residual norm = " << std::sqrt(residual_norm)
                      << ", target residual norm = " << target_residual << std::endl;
        }
        // Check convergence criteria
        if (sqrt(residual_norm) <= target_residual) {
            if (proc_id == 0 && verb >= 0) {
                std::cout << "CG converged. Iteration = " << iteration << ", residual norm = " << std::sqrt(residual_norm)
                          << ", target residual norm = " << target_residual << std::endl;
            }
            return 0;
        }

        // Get new direction
        // beta = (z_i, r_i) / (z_i-1, r_i-1)
        RealType beta = rCz;
        precond.apply(r.data(), z.data());
        rCz = BLAS::dotc(row_number, r.data(), 1, z.data(), 1); // (z_i, r_i)
#ifdef MSK_USE_MPI
        MPI::Allreduce_inplace_sum(&rCz, 1, comm);
#endif
        beta = rCz / beta;
        // p_i+1 = r_i + beta p_i
        BLAS::rscal(row_number, beta, p.data(), 1);
        BLAS::axpy(row_number, DataType{1.0}, z.data(), 1, p.data(), 1);
    }
    if (verb >= 0 && proc_id == 0) {
        std::cout << "CG did not converge at iterations limit = " << iteration_limit << ", residual norm = " << std::sqrt(residual_norm)
                  << ", target residual norm = " << target_residual << std::endl;
    }
    return -1;
}


//--------------------Regularized preconditioned CGNR--------------------------------
/**
 * @brief Regularized Preconditioned Conjugate Gradient Normal Residual method.
 *
 * @details Used for iterative solving a system of linear equations in normal form:
 * (A^* A + lambda^2 I) C^-1 y = A^* b, where M^-1 is a preconditioner,
 * lambda is the regularization parameter.
 * Works also for multiprocessing systems using MPI.
 *
 * @tparam DataType Elements type.
 * @tparam Matrix Abstract matrix object, required toe have method matvec(x, y), performing matrix-vector product.
 * @tparam Precond Abstract preconditioner, required to have method apply(x, y), performing preconditioning.
 * @param[in] A Matrix.
 * @param[in] C Preconditioner.
 * @param[in] b Right-hand side.
 * @param[in, out] x Initial guess (input) and solution (output).
 * @param[in] parameters Solving parameters.
 * @param[in] lambda2 Regularization parameter.
 * @return Success/error code.
 */
template <class DataType, class Matrix, class Precond>
int64_t rCGNR(const Matrix &A, const Precond &C, const DataType *b, DataType *x,
              Parameters<decltype(std::real(DataType()))> parameters,
              decltype(std::real(DataType())) lambda2
#ifdef MSK_USE_MPI
        , MPI_Comm comm = MPI_COMM_WORLD
#endif
) {
    using RealType = decltype(std::abs(DataType())); //RealType is needed if DataType is complex
    int proc_id{0};
#ifdef MSK_USE_MPI
    int n_proc{1};
    MPI_Comm_rank(comm, &proc_id);
    MPI_Comm_size(comm, &n_proc);
#endif
    const uint64_t row_number = A.get_row_number();
    const uint64_t column_number = A.get_column_number();
    uint64_t global_row_number = row_number, global_column_number{column_number};
#ifdef MSK_USE_MPI
    MPI_Allreduce(MPI_IN_PLACE, &global_row_number, 1, MPI_UINT64_T, MPI_SUM, comm);
    MPI_Allreduce(MPI_IN_PLACE, &global_column_number, 1, MPI_UINT64_T, MPI_SUM, comm);
#endif
    int verb = parameters.get_verbosity();
    if (proc_id == 0 && verb >= 0) {
        std::cout << "Starting rCGNR for matrix of size " << global_row_number << " x ";
        std::cout << global_column_number << std::endl;
        parameters.print_parameters();
    }

    SmartArray<DataType> r(row_number);	// r_i, current residual, orthogonal basis of Krylov space
    SmartArray<DataType> ner(column_number); // ner_i, current residual of normal system
    SmartArray<DataType> Cner(column_number); // C ner_i, current preconditioned residual of normal system
    SmartArray<DataType> p(column_number);	// p_{i+1}, current direction, A-orthogonal basis of Krylov space
    SmartArray<DataType> Ap(row_number); // Ap_{i+1}
    SmartArray<DataType> AHAp(column_number); // (A^H)Ap_{i+1}

    // Compute initial residual
    BLAS::copy(row_number, b, 1, r.data(), 1); // Copy of rhs
    if (BLAS::nrm2(column_number, x, 1) > 0) { // If x nonzero, subtract Ax
        A.matvec(x, Ap.data()); // Ap = Ax
        BLAS::axpy(row_number, DataType{-1.0}, Ap.data(), 1, r.data(), 1);
    }
    // Compute initial residual norm
#ifdef MSK_USE_MPI
    RealType init_res_norm = std::real(BLAS::dotc(row_number, r.data(), 1, r.data(), 1));
    MPI::Allreduce_inplace_sum(&init_res_norm, 1, comm);
    RealType res_norm_sq = init_res_norm;
    init_res_norm = std::sqrt(init_res_norm);
#else
    RealType init_res_norm = BLAS::nrm2(row_number, r.data(), 1);
    RealType res_norm_sq = init_res_norm * init_res_norm;
#endif
    // Get target residual norm
    RealType target_res_norm = init_res_norm * parameters.get_relative_tolerance() + parameters.get_absolute_tolerance();
    if (init_res_norm <= target_res_norm) { // Fast exit
        if (verb >= 0 && proc_id == 0) {
            std::cout << "rCGNR converged at start by ";
            std::cout << "residual = " << init_res_norm << ", target residual = " << target_res_norm;
            std::cout << std::endl;
        }
        return 0;
    }

    // Get initial residual of normal system
    A.conj_matvec(r.data(), ner.data());
    if (lambda2 > 0) {
        BLAS::axpy(column_number, -DataType{lambda2}, x, 1, ner.data(), 1);
    }
    RealType init_ne_res_norm = BLAS::nrm2(column_number, ner.data(), 1);
#ifdef MSK_USE_MPI
    init_ne_res_norm = init_ne_res_norm * init_ne_res_norm;
    MPI::Allreduce_inplace_sum(&init_ne_res_norm, 1, comm);
    init_ne_res_norm = std::sqrt(init_ne_res_norm);
#endif
    // Initial residual tolerance is chosen close to machine precision
    RealType target_ne_res_norm = std::numeric_limits<RealType>::epsilon() * 5 * init_ne_res_norm;
    if (init_ne_res_norm <= target_ne_res_norm) { // Fast exit
        if (verb >= 0 && proc_id == 0) {
            std::cout << "rCGNR converged at start by ";
            std::cout << "NE residual = " << init_ne_res_norm << ", target NE residual = " << target_ne_res_norm;
            std::cout << std::endl;
        }
        return 0;
    }
    RealType ne_res_norm = init_ne_res_norm;

    // Compute preconditioned residual
    C.apply(ner.data(), Cner.data());
    BLAS::copy(column_number, Cner.data(), 1, p.data(), 1); // p_1 = Cner_0

    DataType ner_Cner = BLAS::dotc(column_number, ner.data(), 1, Cner.data(), 1);
#ifdef MSK_USE_MPI
    MPI::Allreduce_inplace_sum(&ner_Cner, 1, comm);
#endif

    // Cycle until iteration limit reached
    uint64_t iteration_limit = parameters.get_iterations_limit();
    if (iteration_limit == 0) {
        iteration_limit = global_column_number;
    }

    RealType sol_norm{BLAS::nrm2(column_number, x, 1)};
    RealType reg_res_norm = std::sqrt(res_norm_sq + lambda2 * sol_norm * sol_norm);
    const RealType init_reg_res_norm = reg_res_norm;

    RealType target_reg_res_norm = init_reg_res_norm * parameters.get_relative_tolerance() + parameters.get_absolute_tolerance();
    if (init_reg_res_norm <= target_reg_res_norm) { // Fast exit
        if (verb >= 0 && proc_id == 0) {
            std::cout << "rCGNR converged at start by ";
            std::cout << "REG residual = " << init_reg_res_norm << ", target REG residual = " << target_reg_res_norm;
            std::cout << std::endl;
        }
        return 0;
    }
    
    constexpr int64_t n_recalc_res_every{20};
    int64_t n_last_restart{0};
    int64_t n_converge_restart = 4; // how many times to recompute true residual and direction vector

    auto message = [&](uint64_t iteration, const std::string &pre_msg = "") {
        if (proc_id == 0 && verb >= 0) {
            std::cout << pre_msg;
            std::cout << "It = " << iteration;
            std::cout << ", residual = " << std::sqrt(res_norm_sq) << ", target residual = " << target_res_norm;
            std::cout << ", NE residual = " << ne_res_norm << ", target NE residual = " << target_ne_res_norm;
            std::cout << ", REG residual = " << reg_res_norm << ", target residual = " << target_reg_res_norm;
            std::cout << std::endl;
        }
    };

    for (uint64_t iteration = 1; iteration <= iteration_limit; iteration++) {
        // Perform matvec
        A.matvec(p.data(), Ap.data()); // calc Ap_i
        DataType p_norm_sq = std::real(BLAS::dotc(column_number, p.data(), 1, p.data(), 1)); // ||p_i||^2
        DataType Ap_norm_sq = std::real(BLAS::dotc(row_number, Ap.data(), 1, Ap.data(), 1)); // ||Ap_i||^2
#ifdef MSK_USE_MPI
        DataType tmp_dot[2] = {p_norm_sq, Ap_norm_sq};
        MPI::Allreduce_inplace_sum(tmp_dot, 2, comm);
        p_norm_sq = tmp_dot[0];
        Ap_norm_sq = tmp_dot[1];
#endif

        // Update residual and solution
        // alpha = (ner_{i-1}, Cner_{i-1}) / (||Ap_i||^2 + lambda2 * ||p_i||^2)
        DataType alpha = ner_Cner / (Ap_norm_sq + lambda2 * p_norm_sq);
        // x_i = x_i-1 + alpha * p_i
        BLAS::axpy(column_number, alpha, p.data(), 1, x, 1);
        // r_i = r_i-1 - alpha * Ap_i
        BLAS::axpy(row_number, -alpha, Ap.data(), 1, r.data(), 1);

        RealType old_res_norm_sq{res_norm_sq};
        res_norm_sq = std::real(BLAS::dotc(row_number, r.data(), 1, r.data(), 1));
#ifdef MSK_USE_MPI
        MPI::Allreduce_inplace_sum(&res_norm_sq, 1, comm);
#endif
        // Update normal residual
        // n_{i+1} = n_i - alpha * (A^* Ap_i + lambda2 p_i)
        A.conj_matvec(Ap.data(), AHAp.data());
        BLAS::axpy(column_number, -alpha, AHAp.data(), 1, ner.data(), 1);
        if (lambda2 > 0) {
            BLAS::axpy(column_number, -lambda2 * alpha, p.data(), 1, ner.data(), 1);
        }
        // RealType old_ne_res_norm = ne_res_norm;
#ifdef MSK_USE_MPI
        ne_res_norm = BLAS::dotc(column_number, ner.data(), 1, ner.data(), 1);
        MPI::Allreduce_inplace_sum(&ne_res_norm, 1, comm);
        ne_res_norm = std::sqrt(std::real(ne_res_norm));
#else
        ne_res_norm = BLAS::nrm2(column_number, ner.data(), 1);
#endif
        RealType old_reg_res_norm{reg_res_norm};
        sol_norm = BLAS::nrm2(column_number, x, 1);
        reg_res_norm = std::sqrt(res_norm_sq + lambda2 * sol_norm * sol_norm);

        bool cold_restart{false}; // Flag if to fully restart the method
        RealType reg_res_drop = old_reg_res_norm / reg_res_norm - 1;
        const bool low_res_def{reg_res_drop < 2 * std::sqrt(std::numeric_limits<RealType>::epsilon())};
        // If residual has increased then try to restart iterations
        if (iteration - n_last_restart >= n_recalc_res_every || low_res_def) {
            if (n_converge_restart <= 0) { // Can't restart, finish
                message(iteration, "rCGNR stopped by REG residual. ");
                return iteration;
            }
            // Restart
            if (proc_id == 0 && verb >= 0) {
                std::cout << "Restart at it " << iteration;
                std::cout << ", REG residual changed from " << old_reg_res_norm / init_reg_res_norm << " to " << reg_res_norm / init_reg_res_norm;
                std::cout << std::endl;
            }
            // Perform residuals recalculation
            A.matvec(x, Ap.data());
            BLAS::copy(row_number, b, 1, r.data(), 1);
            BLAS::axpy(row_number, DataType{-1}, Ap.data(), 1, r.data(), 1); // correct residual
            A.conj_matvec(r.data(), ner.data()); // correct normal residual
            if (lambda2 > 0) {
                BLAS::axpy(column_number, -DataType{lambda2}, x, 1, ner.data(), 1);
            }

            res_norm_sq = std::real(BLAS::dotc(row_number, r.data(), 1, r.data(), 1));
#ifdef MSK_USE_MPI
            MPI::Allreduce_inplace_sum(&res_norm_sq, 1, comm);
#endif
            reg_res_norm = std::sqrt(res_norm_sq + lambda2 * sol_norm * sol_norm);
#ifdef MSK_USE_MPI
            ne_res_norm = BLAS::dotc(column_number, ner.data(), 1, ner.data(), 1);
            MPI::Allreduce_inplace_sum(&ne_res_norm, 1, comm);
            ne_res_norm = std::sqrt(std::real(ne_res_norm));
#else
            ne_res_norm = BLAS::nrm2(column_number, ner.data(), 1);
#endif
            // Update possible restarts count
            if (low_res_def) {
                // Sequential low residual change -> convergence
                if (n_last_restart == iteration - 1) {
                    --n_converge_restart;
                    cold_restart = n_converge_restart <= 0;
                } else {
                    n_converge_restart = 2;
                }
            }
            n_last_restart = iteration;
        }

        // Check convergence criteria
        if (reg_res_norm <= target_reg_res_norm) {
            message(iteration, "rCGNR converged by REG residual. ");
            return iteration;
        }
        if (ne_res_norm <= target_ne_res_norm) { // if converged by normal residual
            message(iteration, "rCGNR converged by NE residual. ");
            return iteration;
        }

        if (proc_id == 0 && verb > 0 && iteration % verb == 0) {
            std::cout << "It " << iteration;
            std::cout << ": \tres = " << std::sqrt(res_norm_sq) / init_res_norm;
            std::cout << ", \tREG res = " << reg_res_norm;
            std::cout << ", \tNE res = " << ne_res_norm / init_ne_res_norm << std::endl;
        }

        // Get new direction
        // beta = (ner_i, Cner_i) / (ner_{i-1}, Cner_{i-1})
        // beta is if RealType if Prec is symmetric, but not in general case
        DataType beta = ner_Cner;
        // Precondition normal residual
        C.apply(ner.data(), Cner.data());
        ner_Cner = BLAS::dotc(column_number, ner.data(), 1, Cner.data(), 1); // (ner_i, Cner_i)
#ifdef MSK_USE_MPI
        MPI::Allreduce_inplace_sum(&ner_Cner, 1, comm);
#endif
        beta = ner_Cner / beta;
        if (cold_restart) { // Fully restart the method from new residual: p_{i+1} = Cner_i
            BLAS::copy(column_number, Cner.data(), 1, p.data(), 1);
        } else {
            // p_{i+1} = Cner_i + beta * p_i
            BLAS::scal(column_number, beta, p.data(), 1);
            BLAS::axpy(column_number, DataType{1.0}, Cner.data(), 1, p.data(), 1);
        }
    }
    message(iteration_limit, "rCGNR did not converge at iterations limit. ");
    return -1;
}

}

}

}

#endif
