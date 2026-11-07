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
// Created by b.valiakhmetov on 10/28/24.
//

#ifndef MOSAIC_SKELETON_SOLVERS_LINSYS_PARAMETERS_HPP
#define MOSAIC_SKELETON_SOLVERS_LINSYS_PARAMETERS_HPP


#include <cstdint>


namespace MSk {

namespace solvers {

namespace linsys {

namespace gmres {

enum class Orthogonalization {
    ClassicalGramSchmidt,
    ModifiedGramSchmidt
};

enum class Reorthogonalization {
    Never,
    AsNeeded,
    Always
};

}


/**
 * @brief Parameters of linear system solving using an iterative method.
 *
 * @tparam RealType Real data type.
 */
template <class RealType>
class Parameters {
private:
    RealType absolute_tolerance{0.0}; ///< Error absolute tolerance
    RealType relative_tolerance{1e-4}; ///< Error relative tolerance
    uint64_t iterations_limit{0}; ///< Limit on number of iterations
    uint64_t restart{0}; ///< Number of iterations prior restart
    int verbosity{0}; ///< Verbosity level
    gmres::Orthogonalization orthogonalization_method{gmres::Orthogonalization::ClassicalGramSchmidt}; ///< ClassicalGramSchmidt or ModifiedGramSchmidt
    gmres::Reorthogonalization reorthogonalization_schedule{gmres::Reorthogonalization::AsNeeded}; ///< Never, or AsNeeded, or Always

public:
    Parameters() = default;
    /**
     * @brief Initialize parameters.
     *
     * @param rel_tol Error relative tolerance.
     * @param abs_tol Error absolute tolerance.
     * @param iter_lim Iterations limit.
     * @param restart_lim Restart iterations at the number.
     * @param verb_every Verbose at iterations: if -1 - never, if 0 - start and finish, if k>0 - every k.
     */
    explicit Parameters(RealType rel_tol, RealType abs_tol = RealType{0.0}, uint64_t iter_lim = 0, uint64_t restart_lim = 0,
               gmres::Orthogonalization orth_method = gmres::Orthogonalization::ClassicalGramSchmidt,
               gmres::Reorthogonalization reorth_sched = gmres::Reorthogonalization::AsNeeded,
               int verb_every = 0) :
            absolute_tolerance{abs_tol},
            relative_tolerance{rel_tol},
            iterations_limit{iter_lim},
            restart{restart_lim},
            orthogonalization_method{orth_method},
            reorthogonalization_schedule{reorth_sched},
            verbosity{verb_every} {}


    const RealType &get_absolute_tolerance() const {return absolute_tolerance;}
    const RealType &get_relative_tolerance() const {return relative_tolerance;}
    const uint64_t &get_iterations_limit() const {return iterations_limit;}
    const uint64_t &get_restart() const {return restart;}
    const int &get_verbosity() const {return verbosity;}
    const gmres::Orthogonalization &get_orthogonalization_method() const {return orthogonalization_method;}
    const gmres::Reorthogonalization &get_reorthogonalization_schedule() const {return reorthogonalization_schedule;}

    void set_absolute_tolerance(const RealType &atol) {absolute_tolerance = atol;}
    void set_relative_tolerance(const RealType &rtol) {relative_tolerance = rtol;}
    void set_iterations_limit(uint64_t r) {iterations_limit = r;}
    void set_restart(uint64_t r) {restart = r;}
    void set_verbosity(int verb_every) {verbosity = verb_every;}
    void set_orthogonalization_method(const gmres::Orthogonalization &m) {orthogonalization_method = m;}
    void set_reorthogonalization_schedule(const gmres::Reorthogonalization &s) {reorthogonalization_schedule = s;}

    void print_parameters() {
        std::cout << "Relative tolerance = " << relative_tolerance << std::endl;
        std::cout << "Absolute tolerance = " << absolute_tolerance << std::endl;
        std::cout << "Iterations limit = " << iterations_limit << std::endl;
        std::cout << "Verbosity = " << verbosity << std::endl;
    }
};

}

}

}



#endif //MOSAIC_SKELETON_SOLVERS_LINSYS_PARAMETERS_HPP
