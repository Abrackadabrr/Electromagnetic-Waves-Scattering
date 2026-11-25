//
// Created by evgen on 25.11.2026.
//

#include "math/matrix/ToeplitzFullyTemplated.hpp"
#include "math/matrix/iterative_solvers_coverage/DiagonalPreconditioner.hpp"
#include "math/matrix/iterative_solvers_coverage/MatrixTraits.hpp"
#include "math/matrix/iterative_solvers_coverage/ToeplitzTwoLevelGaussSeidel.hpp"

#include <gtest/gtest.h>
#include <unsupported/Eigen/IterativeSolvers>

#include <algorithm>
#include <complex>
#include <iomanip>
#include <iostream>
#include <random>

using namespace EMW;
using namespace EMW::Math::LinAgl::Matrix;

using scalar_t = Types::complex_d;
using dense_t = Types::MatrixX<scalar_t>;
using vector_t = Types::VectorX<scalar_t>;
using inner_t = ToeplitzStructure<scalar_t, dense_t>;
using outer_t = ToeplitzStructure<scalar_t, inner_t>;

// Adapter: the current MatrixReplacement API requires solve(rhs) returning
// a vector, whereas TwoLevelBlockLowerPreconditioner uses solve(rhs, x, ws).
template <typename Matrix> using ToeplitzGaussSeidelEigenAdapter = Preconditioning::ToeplitzGaussSeidelEigen<Matrix>;

using wrapper_no_t = Wrappers::MatrixReplacement<outer_t>;
using wrapper_gs_t = Wrappers::MatrixReplacement<outer_t, ToeplitzGaussSeidelEigenAdapter<outer_t>>;
using wrapper_jac_t =
    Wrappers::MatrixReplacement<outer_t, Preconditioning::BlockDiagonalPreconditioner<scalar_t, outer_t>>;

// Shared runner for both matrix wrappers. MatrixReplacement applies RIGHT
// preconditioning: wrapped * y = A * P^{-1} * y.
struct GMRESRunResult {
    Eigen::Index iterations = 0;
    double estimated_error = 0.0;
    double relative_residual = 0.0;
    double relative_solution_error = 0.0;
    Eigen::ComputationInfo info = Eigen::InvalidInput;
};

std::ostream &operator<<(std::ostream &os, const GMRESRunResult &result) {
    os << std::left << std::setw(16) << std::right << std::setw(12) << result.iterations << std::scientific
       << std::setprecision(3) << std::setw(20) << result.estimated_error << std::setw(20) << result.relative_residual
       << std::setw(20) << result.relative_solution_error << '\n';
    return os;
}

template <typename matrix_t> using GMRES = Eigen::GMRES<matrix_t, Eigen::IdentityPreconditioner>;

template <typename WrappedMatrix>
GMRESRunResult run_gmres(const WrappedMatrix &wrapped, const vector_t &rhs, const dense_t &original_matrix,
                         const vector_t &x_exact, int max_iterations, double tolerance) {

    // Do not apply a second preconditioner inside Eigen::GMRES.
    GMRES<WrappedMatrix> gmres;
    gmres.setTolerance(tolerance);
    gmres.setMaxIterations(max_iterations);
    gmres.set_restart(max_iterations); // no restart within the iteration budget
    gmres.compute(wrapped);

    // The Krylov method solves A * P^{-1} * y = rhs.
    const vector_t y = gmres.solve(rhs);
    // Recover x = P^{-1} * y (identity for the unpreconditioned wrapper).
    const vector_t x = wrapped.get_preconditioner().solve(y);

    return {
        .iterations = gmres.iterations(),
        .estimated_error = gmres.error(),
        .relative_residual = (original_matrix * x - rhs).norm() / rhs.norm(),
        .relative_solution_error = (x - x_exact).norm() / x_exact.norm(),
        .info = gmres.info(),
    };
}

// Gaussian-free deterministic data generation, independent from Eigen::Random.
dense_t random_block(Eigen::Index b, std::mt19937_64 &gen, double scale) {
    std::uniform_real_distribution<double> uniform(-1.0, 1.0);
    dense_t block(b, b);
    for (Eigen::Index col = 0; col < b; ++col) {
        for (Eigen::Index row = 0; row < b; ++row) {
            block(row, col) = scale * scalar_t(uniform(gen), uniform(gen));
        }
    }
    return block;
}

TEST(TwoLevelBlockLowerPreconditioner, GMRESIterationCount) {
    constexpr Types::index n_outer = 5;
    constexpr Types::index n_inner = 5;
    constexpr Eigen::Index b = 8;
    constexpr int max_iterations = 150;
    constexpr double relative_tolerance = 1e-9;

    std::mt19937_64 gen(20261125);

    Containers::vector<inner_t> outer_blocks;
    outer_blocks.reserve(outer_t::get_size_of_container(n_outer, n_outer));

    // Unique Toeplitz storage order: first column bottom-to-top,
    // then first row left-to-right. Diagonal index on each level is n-1.
    for (Types::index oi = 0; oi < outer_t::get_size_of_container(n_outer, n_outer); ++oi) {
        Containers::vector<dense_t> inner_blocks;
        inner_blocks.reserve(inner_t::get_size_of_container(n_inner, n_inner));

        for (Types::index ii = 0; ii < inner_t::get_size_of_container(n_inner, n_inner); ++ii) {
            const auto d_outer = std::abs(static_cast<double>(oi) - static_cast<double>(n_outer - 1));
            const auto d_inner = std::abs(static_cast<double>(ii) - static_cast<double>(n_inner - 1));

            const double distance = std::hypot(d_outer, d_inner);
            const bool diagonal = (distance == 0.0);

            dense_t block = random_block(b, gen, diagonal ? 0.2 : 0.5);

            if (diagonal) {
                block.diagonal().array() += scalar_t(8.0);
            } else {
                block /= distance;
            }

            inner_blocks.emplace_back(std::move(block));
        }

        outer_blocks.emplace_back(n_inner, n_inner, std::move(inner_blocks));
    }
    const outer_t A(n_outer, n_outer, std::move(outer_blocks));
    const dense_t dense_A = A.to_dense(); // Test reference only (200 x 200).
    ASSERT_EQ(dense_A.rows(), dense_A.cols());

    const Eigen::Index n = dense_A.rows();
    const vector_t x_exact = vector_t::Ones(n);
    const vector_t rhs = dense_A * x_exact;

    // One and the same A, rhs, iteration budget, tolerance, restart policy.
    const wrapper_no_t A_without_prec(A);
    const wrapper_gs_t A_with_gs(A);
    const wrapper_jac_t A_with_jac(A);

    // Check the matrix-free operators before starting GMRES.
    const vector_t probe = vector_t::LinSpaced(n, scalar_t(0.1), scalar_t(0.9));
    const vector_t expected_no = dense_A * probe;
    const vector_t actual_no = A_without_prec * probe;
    ASSERT_LT((actual_no - expected_no).norm() / expected_no.norm(), 1e-12);

    const vector_t z_probe = A_with_gs.get_preconditioner().solve(probe);
    const vector_t expected_gs = dense_A * z_probe;
    const vector_t actual_gs = A_with_gs * probe;
    ASSERT_LT((actual_gs - expected_gs).norm() / expected_gs.norm(), 1e-12);

    // Same GMRES configuration for both wrappers; only the right
    // preconditioner embedded in MatrixReplacement differs.
    const auto no = run_gmres(A_without_prec, rhs, dense_A, x_exact, max_iterations, relative_tolerance);
    const auto gs = run_gmres(A_with_gs, rhs, dense_A, x_exact, max_iterations, relative_tolerance);
    const auto jac = run_gmres(A_with_jac, rhs, dense_A, x_exact, max_iterations, relative_tolerance);

    std::cout << std::endl;
    std::cout << "no prec" << no << std::endl;
    std::cout << "gs prec" << gs << std::endl;
    std::cout << "jac prec" << jac << std::endl;

    EXPECT_EQ(no.info, Eigen::Success);
    EXPECT_EQ(gs.info, Eigen::Success);
    EXPECT_EQ(jac.info, Eigen::Success);
    EXPECT_LT(no.relative_residual, 1e-8);
    EXPECT_LT(gs.relative_residual, 1e-8);
    EXPECT_LT(jac.relative_residual, 1e-8);
    EXPECT_LT(no.relative_solution_error, 1e-7);
    EXPECT_LT(gs.relative_solution_error, 1e-7);
    EXPECT_LT(jac.relative_solution_error, 1e-7);
}
