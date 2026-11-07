/*
 * Copyright 2024-2026 Bulat Valiakhmetov
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
// Created by b.valiakhmetov on 2/22/24.
//

#ifndef MOSAIC_SKELETON_EXAMPLES_COMMON_COMMON_HPP
#define MOSAIC_SKELETON_EXAMPLES_COMMON_COMMON_HPP


#include "msk_examples_config.h"

#include "toml.hpp"

#include "msk/oracle/oracle.hpp"
#include "msk/utility/logging.hpp"
#include "msk/solvers/linsys/GMRES.hpp"


namespace MSk {


/**
 * @brief Create Logger from toml table.
 *
 * @details Toml table may contain fields "name", "file" and "name".
 * If filename equals "stdout" then logging will be directed to standard output.
 *
 * @param[in] log_params Logging parameters.
 * @param[in] default_name Default logger name.
 * @param[in] default_file Default logger filename.
 * @param[in] default_level Default logging level.
 * @return Logger.
 */
std::shared_ptr<log::Logger> build_logger(const toml::value &log_params,
                                          const std::string &default_name = "Log",
                                          const std::string &default_file = "stdout",
                                          const std::string &default_level = "info") {
    const std::string name = toml::find_or<std::string>(log_params, "name", default_name);
    const std::string file = toml::find_or<std::string>(log_params, "file", default_file);
    const std::string level = toml::find_or<std::string>(log_params, "level", default_level);
    return log::build_logger(name, file, MSk::log::string_to_level(level));
}


namespace oracle {


/**
 * @brief Convert string to policy.
 *
 * @details In case of incorrect name DividePolicy::dual is reterned.
 * @param[in] policy Policy name.
 * @return Policy value.
 */
DividePolicy string_to_policy(const std::string &policy) {
    if (policy == "single") {
        return DividePolicy::single;
    }
    return DividePolicy::dual;
}


/**
 * @brief Construct oracle::Parameters from toml table.
 *
 * @param[in] table Toml table.
 */
Parameters read_parameters(toml::value &oracle_params, toml::value &approx_params) {
    Parameters params;
    // auto &oracle_params = toml::find(table, "oracle");

    params.min_rows = toml::find_or<toml::integer>(oracle_params, "min_rows", 128);
    params.min_cols = toml::find_or<toml::integer>(oracle_params, "min_cols", 128);

    auto v_rho = toml::find(oracle_params, "rho");
    params.rho = v_rho.is_floating() ? v_rho.as_floating(std::nothrow) : static_cast<double>(v_rho.as_integer());
    params.sphere_norm = toml::find_or<toml::integer>(oracle_params, "sphere_norm", 0);

    params.policy = string_to_policy(toml::find_or<toml::string>(oracle_params, "policy", "dual"));

    params.force_lowrank = toml::find_or<toml::boolean>(oracle_params, "force_lowrank", false);

    auto oracle_log_params = toml::find_or<toml::table>(oracle_params, "log", {{"name", "log"},
                                                                              {"file", "stdout"},
                                                                              {"level", "off"}});
    params.log = build_logger(oracle_log_params, "Oracle", "./log/oracle.txt", "info");

    // auto &approx_params = toml::find(table, "approximation");
    double abs_tol = toml::find_or<toml::floating>(approx_params, "abs_tol", 0.0);
    double rel_tol = toml::find_or<toml::floating>(approx_params, "rel_tol", 1e-3);
    uint64_t max_rank = toml::find_or<toml::integer>(approx_params, "max_rank", 0);
    params.block_comp = {abs_tol, rel_tol, max_rank};

    auto block_log_params = toml::find_or<toml::table>(approx_params, "log", {{"name", "log"},
                                                                             {"file", "stdout"},
                                                                             {"level", "off"}});
    params.block_log = build_logger(block_log_params, "Approx", "./log/approximation.txt", "info");

    return params;
}

}


/**
 * @brief Fill GMRES::Parameters structure from toml table.
 *
 * @details Default values for parameters that table does not provide:
 * absolute tolerance 0.0, relative tolerance 1e-6, no restarts, no iterations limis,
 * classical Gram-Schmidt orthogonalization and reorthogonalization as needed.
 *
 * @tparam RealType Real parameters type.
 * @param[in] table Toml table.
 * @param[out] params Parameters structure.
 */
template<typename RealType>
void set_gmres_parameters(const toml::value &table, MSk::solvers::linsys::Parameters<RealType> &params) {
    RealType sol_abs_tol = toml::find_or<toml::floating>(table, "abs_tol", 0.0);
    RealType sol_rel_tol = toml::find_or<toml::floating>(table, "rel_tol", 1e-6);
    uint64_t restart = toml::find_or<toml::integer>(table, "restart", 0);
    uint64_t iter_limit = toml::find_or<toml::integer>(table, "iterations_limit", 0);
    std::string orth = toml::find_or<toml::string>(table, "orthogonalization", "GS");
    std::string reorth = toml::find_or<toml::string>(table, "reorthogonalization", "as_needed");

    params.set_absolute_tolerance(sol_abs_tol);
    params.set_relative_tolerance(sol_rel_tol);
    params.set_iterations_limit(iter_limit);
    params.set_restart(restart);
    if (orth == "MGS") {
        params.set_orthogonalization_method(MSk::solvers::linsys::gmres::Orthogonalization::ModifiedGramSchmidt);
    } else {
        params.set_orthogonalization_method(MSk::solvers::linsys::gmres::Orthogonalization::ClassicalGramSchmidt);
    }
    if (reorth == "always" || reorth == "Always") {
        params.set_reorthogonalization_schedule(MSk::solvers::linsys::gmres::Reorthogonalization::Always);
    } else if (reorth == "never" || reorth == "Never") {
        params.set_reorthogonalization_schedule(MSk::solvers::linsys::gmres::Reorthogonalization::Never);
    } else {
        params.set_reorthogonalization_schedule(MSk::solvers::linsys::gmres::Reorthogonalization::AsNeeded);
    }
}


/**
 * @brief Calculate the 3-dimensional Helmholtz equation kernel for two points.
 *
 * @details Returns exp(i k r) / r, r = |x - y|.
 *
 * @tparam RealType Points coordinate type.
 * @tparam dim Points dimension.
 * @param[in] point_x First point.
 * @param[in] point_y Second point.
 * @param[in] k Wave number.
 * @return Kernel value.
 */
template<typename RealType, uint64_t dim>
std::complex<RealType> helmgoltz_fundamental(const geom::Point<RealType, dim> &point_x,
                                             const geom::Point<RealType, dim> &point_y,
                                             const std::complex<RealType> k = {1.0}) {
    std::complex<RealType> res{};
    RealType dist = (point_x - point_y).norm();
    if (dist < 10 * std::numeric_limits<double>::epsilon()) {
        res = 0;
    } else {
        res = std::exp(std::complex<RealType>(0.0, 1.0) * k * dist) / dist;
    }
    return res;
}


/**
 * @brief Get the euclidean distance between two points or its inverse value.
 *
 * @details For equal points zero in returned instea of infinite inverse distance.
 *
 * @tparam RealType Points coordinate type.
 * @tparam dim Points dimension.
 * @param[in] point_x First point.
 * @param[in] point_y Second point.
 * @param[in] inverse If to inverse the distance.
 * @return Distance or its inverse.
 */
template<typename RealType, uint64_t dim>
RealType distance(const geom::Point<RealType, dim> &point_x,
                  const geom::Point<RealType, dim> &point_y,
                  bool inverse = false) {
    RealType dist{(point_x - point_y).norm()};
    if (!inverse) {
        return dist;
    }
    if (dist < 10 * std::numeric_limits<RealType>::epsilon()) {
        return 0;
    } else {
        return RealType{1.0} / dist;
    }
}


/**
 * @brief Matrix Hilbert function.
 *
 * @details Indices are counted from zero, returned value is 1 / (i + j + 1).
 *
 * @tparam DataType Matrix elemet type.
 * @param i Row index.
 * @param j Column index.
 * @return Hilbert function value.
 */
template<class DataType>
DataType hilbert(uint64_t i, uint64_t j) {
    return DataType{1.0} / (1.0 + i + j);
}

}


#endif //MOSAIC_SKELETON_EXAMPLES_COMMON_COMMON_HPP
