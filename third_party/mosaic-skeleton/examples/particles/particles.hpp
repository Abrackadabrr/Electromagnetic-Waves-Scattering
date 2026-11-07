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
// Created by b.valiakhmetov on 1/29/24.
//


/**
 * @file
 *
 * @brief Contains some examples of matrix functions for particles simulation theory.
 */


#ifndef MOSAIC_SKELETON_MATRIX_PARTICLES_HPP
#define MOSAIC_SKELETON_MATRIX_PARTICLES_HPP


#include <cstdint>
#include <cmath>


namespace MSk {


/**
 * @brief Kernel for atmosphere related Smoluchowsky process.
 *
 * @tparam RealType Kernel element type.
 * @param i Row number, starting from 0.
 * @param j Column number, starting from 0.
 * @return Kernel element value.
 */
template<class RealType>
RealType ker_atmos(uint64_t i, uint64_t j) {
    ++i, ++j;
    if (i == j) {
        return 4;
    }
    RealType num = (i + j) * std::pow(std::pow(i, 1.0 / 3.0) + std::pow(j, 1.0 / 3.0), 2.0 / 3.0);
    RealType denom = std::pow(i * j, 5.0 / 9.0) * std::abs(std::pow(i, 2.0 / 3.0) - std::pow(j, 2.0 / 3.0));
    return num / denom;
}


/**
 * @brief Kernel for liquid stream related Smoluchowsky process.
 *
 * @tparam RealType Kernel element type.
 * @param i Row number, starting from 0.
 * @param j Column number, starting from 0.
 * @return Kernel element value.
 */
template<class RealType>
RealType ker_stream(uint64_t i, uint64_t j) {
    ++i, ++j;
    if (i == j) {
        return 4;
    }
    RealType lr_part = std::pow(i, 1.0 / 3.0) + std::pow(j, 1.0 / 3.0);
    lr_part = lr_part * lr_part;
    RealType diff = std::abs(std::pow(i, 2.0 / 3.0) - std::pow(j, 2.0 / 3.0));
    return lr_part * diff;
}


/**
 * @brief Real fusion kernel.
 *
 * @tparam RealType Kernel element type.
 * @param i Row number, starting from 0.
 * @param j Column number, starting from 0.

 * @return Kernel element value.
 */
template<class RealType>
RealType ker_real_fusion(uint64_t i, uint64_t j) {
    ++i, ++j;
    if (i == j) {
      return 0;
    }
    RealType diff = std::abs(std::pow(i, 2.0 / 3.0) - std::pow(j, 2.0 / 3.0));
    RealType sum = std::sqrt(std::pow(i, 1.0 / 3.0) + std::pow(j, 1.0 / 3.0));
    return diff / sum * std::erf(diff / sum);
}


/**
 * @brief Fusion kernel.
 *
 * @tparam RealType Kernel element type.
 * @param i Row number, starting from 0.
 * @param j Column number, starting from 0.
 * @param alpha Temperature parameter.
 * @return Kernel element value.
 */
template<class RealType>
RealType ker_fusion(uint64_t i, uint64_t j, RealType alpha = 1.0) {
    ++i, ++j;
    if (i == j) {
        return 0;
    }
    RealType diff = i > j ? i - j : j - i;
    RealType sum = i + j;
    return diff * std::erf(alpha * diff / sum);
}


/**
 * @brief Difference module kernel.
 *
 * @tparam RealType Kernel element type.
 * @param i Row number, starting from 0.
 * @param j Column number, starting from 0.
 * @return Kernel element value.
 */
template<class RealType>
RealType ker_module(uint64_t i, uint64_t j) {
    if (i == j) {
        return 4;
    }
    return i > j ? i - j : j - i;
}

template<class RealType>
RealType ker_modified(uint64_t i, uint64_t j, RealType c = 1.0) {
    ++i, ++j;
    return ((i + j) * (i + j)) / (i * j + c);
}

template<class RealType>
RealType ker_emulsion(uint64_t i, uint64_t j, RealType c1 = {1.0}, RealType c2 = {1.0}) {
    ++i, ++j;
    RealType m1 = std::pow(i, (2.0 / 3.0)) + std::pow(j, (2.0 / 3.0));
    RealType m2 = std::sqrt(std::pow(i, (2.0 / 9.0)) + std::pow(j, (2.0 / 9.0)));
    RealType d = std::pow(i * j, 1.0 / 3.0) / (std::pow(i, 1.0 / 3.0) + std::pow(j, 1.0 / 3.0));
    return c1 * m1 * m2 * std::exp(-c2 * d * d * d * d);
}

/**
 * @brief Modified difference module kernel.
 *
 * @tparam RealType Kernel element type.
 * @param i Row number, starting from 0.
 * @param j Column number, starting from 0.
 * @param p Power of i and j under the module.
 * @return Kernel element value.
 */
template<class RealType>
RealType ker_module_powered(uint64_t i, uint64_t j, RealType p) {
    ++i, ++j;
    if (i == j) {
        return 4;
    }
    return std::abs(std::pow(i, p) - std::pow(j, p));
}


};


#endif //MOSAIC_SKELETON_MATRIX_PARTICLES_HPP
