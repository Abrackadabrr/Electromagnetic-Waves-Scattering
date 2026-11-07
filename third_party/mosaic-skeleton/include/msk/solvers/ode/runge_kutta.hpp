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

#ifndef MOSAIC_SKELETON_SOLVERS_ODE_RUNGE_KUTTA_HPP
#define MOSAIC_SKELETON_SOLVERS_ODE_RUNGE_KUTTA_HPP


#include <iostream>

#include "smart_array/smart_array.hpp"
#include "cpp_blas/BLAS.hpp"


namespace MSk {
namespace solvers {
namespace ode {

template<typename RealType, class Operator>
void runge_kutta2(uint64_t size, Operator &A, RealType *y, RealType T, uint64_t n_steps) {
    SmartArray<RealType> y_mid(size), A_y_mid(size);
    const RealType step{T / n_steps};
    for (uint64_t i = 0; i < n_steps; ++i) {
        // mid-step
        A.apply(y, y_mid.data());
        BLAS::rscal(size, step / 2.0, y_mid.data(), 1);
        BLAS::axpy(size, RealType(1.0), y, 1, y_mid.data(), 1);

        // step
        A.apply(y_mid.data(), A_y_mid.data());
        BLAS::axpy(size, step, A_y_mid.data(), 1, y, 1);
    }
}

template<typename RealType, class Operator>
void runge_kutta4(uint64_t size, Operator &A, RealType *y, RealType T, uint64_t n_steps) {
    SmartArray<RealType> y_step(size), A_y_step(size), dy(size);
    const RealType step{T / n_steps};
    for (uint64_t i = 0; i < n_steps; ++i) {
        if (i % (n_steps / 10) == 0) {
            RealType norm{BLAS::nrm2(size, y, 1)};
            std::cout << i << " :: t = " << step * i << ", norm = " << norm << std::endl;
            if (std::isnan(norm)) {
                break;
            }
        }
        // step 1
        A.apply(y, A_y_step.data()); // k1
        BLAS::copy(size, A_y_step.data(), 1, dy.data(), 1); // dy = k1

        // step 2
        BLAS::copy(size, y, 1, y_step.data(), 1);
        BLAS::axpy(size, step / 2, A_y_step.data(), 1, y_step.data(), 1); // y + h/2 k1
        A.apply(y_step.data(), A_y_step.data()); // k2
        BLAS::axpy(size, RealType{2}, A_y_step.data(), 1, dy.data(), 1); // dy = k1 + 2 k2

        // step 3
        BLAS::copy(size, y, 1, y_step.data(), 1);
        BLAS::axpy(size, step / 2, A_y_step.data(), 1, y_step.data(), 1); // y + h/2 k2
        A.apply(y_step.data(), A_y_step.data()); // k3
        BLAS::axpy(size, RealType{2}, A_y_step.data(), 1, dy.data(), 1); // dy = k1 + 2 k2 + 2 k3

        // step 4
        BLAS::copy(size, y, 1, y_step.data(), 1);
        BLAS::axpy(size, step, A_y_step.data(), 1, y_step.data(), 1); // y + h k3
        A.apply(y_step.data(), A_y_step.data()); // k4
        BLAS::axpy(size, RealType{1}, A_y_step.data(), 1, dy.data(), 1); // dy = k1 + 2 k2 + 2 k3 + k4

        BLAS::axpy(size, step / 6, dy.data(), 1, y, 1);
    }
}


template<typename RealType, class Operator>
void runge_kutta4_adaptive(uint64_t size, Operator &A, RealType *y,
                           RealType T, RealType init_step, RealType tol = 1e-4) {
    SmartArray<RealType> y_step(size), A_y_step(size), dy(size), dy_check(size);
    RealType t{0}, step{init_step}, next_t_verbose{T / 100};
    uint64_t steps_total{0}, steps_hit{0};
    while (t < T) {
        if (t >= next_t_verbose) {
            RealType norm{BLAS::nrm2(size, y, 1)};
            std::cout << "t = " << t << " :: step = " << step << ", norm = " << norm;
            std::cout << "; steps: total = " << steps_total << ", successful = " << steps_hit << std::endl;
            if (std::isnan(norm)) {
                break;
            }
            next_t_verbose += T / 100;
        }
        // step 1
        A.apply(y, A_y_step.data()); // k1
        BLAS::copy(size, A_y_step.data(), 1, dy.data(), 1); // dy = k1
        BLAS::copy(size, A_y_step.data(), 1, dy_check.data(), 1); // dy = k1

        // step 2
        BLAS::copy(size, y, 1, y_step.data(), 1);
        BLAS::axpy(size, step / 2, A_y_step.data(), 1, y_step.data(), 1); // y + h/2 k1
        A.apply(y_step.data(), A_y_step.data()); // k2
        BLAS::axpy(size, RealType{2}, A_y_step.data(), 1, dy.data(), 1); // dy = k1 + 2 k2

        // step 3
        BLAS::copy(size, y, 1, y_step.data(), 1);
        BLAS::axpy(size, step / 2, A_y_step.data(), 1, y_step.data(), 1); // y + h/2 k2
        A.apply(y_step.data(), A_y_step.data()); // k3
        BLAS::axpy(size, RealType{2}, A_y_step.data(), 1, dy.data(), 1); // dy = k1 + 2 k2 + 2 k3

        // step 4
        BLAS::copy(size, y, 1, y_step.data(), 1);
        BLAS::axpy(size, step, A_y_step.data(), 1, y_step.data(), 1); // y + h k3
        A.apply(y_step.data(), A_y_step.data()); // k4
        BLAS::axpy(size, RealType{1}, A_y_step.data(), 1, dy.data(), 1); // dy = k1 + 2 k2 + 2 k3 + k4

        // check step tolerance
        BLAS::axpy(size, RealType{-1}, A_y_step.data(), 1, dy_check.data(), 1); // k1 - k4
        RealType dy_err = step / 6 * BLAS::nrm2(size, dy_check.data(), 1);
        if (dy_err < tol) { // make step
            t += step;
            BLAS::axpy(size, step / 6, dy.data(), 1, y, 1);
            ++steps_hit;
        }
        // update step size
        step *= std::min(RealType{4}, std::max(RealType{0.25}, std::pow(tol / (2 * dy_err), 0.25)));
        if (t + step > T) {
            step = T - t;
        }
        ++steps_total;
    }
    std::cout << "Total steps: " << steps_total << ", successful steps: " << steps_hit << std::endl;
}


}
}
}
#endif //MOSAIC_SKELETON_SOLVERS_ODE_RUNGE_KUTTA_HPP
