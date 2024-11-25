//
// Created by evgen on 17.11.2024.
//

#ifndef SIMPLEITERATIONMETHOD_HPP
#define SIMPLEITERATIONMETHOD_HPP
#include "types/Types.hpp"

#include <iostream>
#include <math/MathConstants.hpp>

namespace EMW::IterativeSolvers {

enum INFO_TYPE {
    SUCCESS = 0,
    FAILURE = 1
};

template<typename Operator_t, typename Vector_t>
struct SimpleIterationMethod {
    Types::scalar tolerance = 1e-2;
    Types::index max_iterations = 5000;
    Types::index iterations;
    Types::scalar error;
    Types::scalar tau;

    Vector_t solve(const Operator_t& A, const Operator_t& K, const Vector_t& b, const Vector_t& initial) {
        Vector_t resudial = b - K * initial;
        const Types::scalar b_norm = b.norm();
        error = resudial.norm() / b_norm;
        Vector_t result = initial;
        iterations = 0;
        while (error > tolerance && iterations < max_iterations) {
            result += A * resudial;
            resudial = b - K * result;
            error = resudial.norm() / b_norm;
            iterations++;
            std::cout << iterations << ' ' << error << std::endl;
        }
        return result;
    }

    INFO_TYPE info() const {
        if (iterations < max_iterations) return INFO_TYPE::SUCCESS;
        return INFO_TYPE::FAILURE;
    }
};

}

#endif //SIMPLEITERATIONMETHOD_HPP
