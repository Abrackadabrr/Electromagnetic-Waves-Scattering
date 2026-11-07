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
// Created by b.valiakhmetov on 4/9/24.
//


#include "msk/matrix/basic.hpp"


template class MSk::matrix::Elementary<float>;

template class MSk::matrix::Elementary<double>;

template class MSk::matrix::Elementary<std::complex<float>>;

template class MSk::matrix::Elementary<std::complex<double>>;


template class MSk::matrix::BiPoint<float, float, 2>;

template class MSk::matrix::BiPoint<double, double, 2>;

template class MSk::matrix::BiPoint<std::complex<float>, float, 2>;

template class MSk::matrix::BiPoint<std::complex<double>, double, 2>;


template class MSk::matrix::BiPoint<float, float, 3>;

template class MSk::matrix::BiPoint<double, double, 3>;

template class MSk::matrix::BiPoint<std::complex<float>, float, 3>;

template class MSk::matrix::BiPoint<std::complex<double>, double, 3>;


template class MSk::matrix::Dense<float, 2>;

template class MSk::matrix::Dense<double, 2>;

template class MSk::matrix::Dense<std::complex<float>, 2>;

template class MSk::matrix::Dense<std::complex<double>, 2>;


template class MSk::matrix::Dense<float, 3>;

template class MSk::matrix::Dense<double, 3>;

template class MSk::matrix::Dense<std::complex<float>, 3>;

template class MSk::matrix::Dense<std::complex<double>, 3>;


template class MSk::matrix::BiLinear<float>;

template class MSk::matrix::BiLinear<double>;

template class MSk::matrix::BiLinear<std::complex<float>>;

template class MSk::matrix::BiLinear<std::complex<double>>;
