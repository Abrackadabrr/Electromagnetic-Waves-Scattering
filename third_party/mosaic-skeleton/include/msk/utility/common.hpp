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
// Created by b.valiakhmetov on 14/4/25.
//



/**
 * @file
 *
 * @brief Containg logging interface utilities to spdlog library.
 */

#ifndef MOSAIC_SKELETON_UTILITY_COMMON_HPP
#define MOSAIC_SKELETON_UTILITY_COMMON_HPP


#include "msk_config.h"

#include <complex>
#include <cstdint>
#include <type_traits>
 
 
namespace MSk {


template<typename T>
typename std::enable_if<!std::is_same<T, decltype(std::abs(T{}))>::value, T>::type conjugate(const T &x) {
    return std::conj(x);
}

template<typename T>
typename std::enable_if<std::is_same<T, decltype(std::abs(T{}))>::value, T>::type conjugate(const T &x) {
    return x;
}

}

#endif // MOSAIC_SKELETON_UTILITY_COMMON_HPP
