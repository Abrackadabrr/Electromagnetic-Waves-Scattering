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


#include "msk/block/manager.hpp"
#include "msk/evaluation/queue.hpp"


template class MSk::eval::Queue<float>;

template class MSk::eval::Queue<double>;

template class MSk::eval::Queue<std::complex<float>>;

template class MSk::eval::Queue<std::complex<double>>;



template void MSk::eval::Queue<float>::start(matrix::Elementary<float> &matrix, exec::ThreadPool *pool);
template void MSk::eval::Queue<float>::start(matrix::Dense<float, 2> &matrix, exec::ThreadPool *pool);
template void MSk::eval::Queue<float>::start(matrix::Dense<float, 3> &matrix, exec::ThreadPool *pool);
template void MSk::eval::Queue<float>::start(matrix::BiPoint<float, float, 2> &matrix, exec::ThreadPool *pool);
template void MSk::eval::Queue<float>::start(matrix::BiPoint<float, float, 3> &matrix, exec::ThreadPool *pool);

template void MSk::eval::Queue<double>::start(matrix::Elementary<double> &matrix, exec::ThreadPool *pool);
template void MSk::eval::Queue<double>::start(matrix::Dense<double, 2> &matrix, exec::ThreadPool *pool);
template void MSk::eval::Queue<double>::start(matrix::Dense<double, 3> &matrix, exec::ThreadPool *pool);
template void MSk::eval::Queue<double>::start(matrix::BiPoint<double, double, 2> &matrix, exec::ThreadPool *pool);
template void MSk::eval::Queue<double>::start(matrix::BiPoint<double, double, 3> &matrix, exec::ThreadPool *pool);

template void MSk::eval::Queue<std::complex<float>>::start(matrix::Elementary<std::complex<float>> &matrix, exec::ThreadPool *pool);
template void MSk::eval::Queue<std::complex<float>>::start(matrix::Dense<std::complex<float>, 2> &matrix, exec::ThreadPool *pool);
template void MSk::eval::Queue<std::complex<float>>::start(matrix::Dense<std::complex<float>, 3> &matrix, exec::ThreadPool *pool);
template void MSk::eval::Queue<std::complex<float>>::start(matrix::BiPoint<std::complex<float>, float, 2> &matrix, exec::ThreadPool *pool);
template void MSk::eval::Queue<std::complex<float>>::start(matrix::BiPoint<std::complex<float>, float, 3> &matrix, exec::ThreadPool *pool);

template void MSk::eval::Queue<std::complex<double>>::start(matrix::Elementary<std::complex<double>> &matrix, exec::ThreadPool *pool);
template void MSk::eval::Queue<std::complex<double>>::start(matrix::Dense<std::complex<double>, 2> &matrix, exec::ThreadPool *pool);
template void MSk::eval::Queue<std::complex<double>>::start(matrix::Dense<std::complex<double>, 3> &matrix, exec::ThreadPool *pool);
template void MSk::eval::Queue<std::complex<double>>::start(matrix::BiPoint<std::complex<double>, double, 2> &matrix, exec::ThreadPool *pool);
template void MSk::eval::Queue<std::complex<double>>::start(matrix::BiPoint<std::complex<double>, double, 3> &matrix, exec::ThreadPool *pool);


template void MSk::eval::Queue<float>::run_extract(matrix::Elementary<float> &matrix, exec::ThreadPool *pool);
template void MSk::eval::Queue<float>::run_extract(matrix::Dense<float, 2> &matrix, exec::ThreadPool *pool);
template void MSk::eval::Queue<float>::run_extract(matrix::Dense<float, 3> &matrix, exec::ThreadPool *pool);
template void MSk::eval::Queue<float>::run_extract(matrix::BiPoint<float, float, 2> &matrix, exec::ThreadPool *pool);
template void MSk::eval::Queue<float>::run_extract(matrix::BiPoint<float, float, 3> &matrix, exec::ThreadPool *pool);

template void MSk::eval::Queue<double>::run_extract(matrix::Elementary<double> &matrix, exec::ThreadPool *pool);
template void MSk::eval::Queue<double>::run_extract(matrix::Dense<double, 2> &matrix, exec::ThreadPool *pool);
template void MSk::eval::Queue<double>::run_extract(matrix::Dense<double, 3> &matrix, exec::ThreadPool *pool);
template void MSk::eval::Queue<double>::run_extract(matrix::BiPoint<double, double, 2> &matrix, exec::ThreadPool *pool);
template void MSk::eval::Queue<double>::run_extract(matrix::BiPoint<double, double, 3> &matrix, exec::ThreadPool *pool);

template void MSk::eval::Queue<std::complex<float>>::run_extract(matrix::Elementary<std::complex<float>> &matrix, exec::ThreadPool *pool);
template void MSk::eval::Queue<std::complex<float>>::run_extract(matrix::Dense<std::complex<float>, 2> &matrix, exec::ThreadPool *pool);
template void MSk::eval::Queue<std::complex<float>>::run_extract(matrix::Dense<std::complex<float>, 3> &matrix, exec::ThreadPool *pool);
template void MSk::eval::Queue<std::complex<float>>::run_extract(matrix::BiPoint<std::complex<float>, float, 2> &matrix, exec::ThreadPool *pool);
template void MSk::eval::Queue<std::complex<float>>::run_extract(matrix::BiPoint<std::complex<float>, float, 3> &matrix, exec::ThreadPool *pool);

template void MSk::eval::Queue<std::complex<double>>::run_extract(matrix::Elementary<std::complex<double>> &matrix, exec::ThreadPool *pool);
template void MSk::eval::Queue<std::complex<double>>::run_extract(matrix::Dense<std::complex<double>, 2> &matrix, exec::ThreadPool *pool);
template void MSk::eval::Queue<std::complex<double>>::run_extract(matrix::Dense<std::complex<double>, 3> &matrix, exec::ThreadPool *pool);
template void MSk::eval::Queue<std::complex<double>>::run_extract(matrix::BiPoint<std::complex<double>, double, 2> &matrix, exec::ThreadPool *pool);
template void MSk::eval::Queue<std::complex<double>>::run_extract(matrix::BiPoint<std::complex<double>, double, 3> &matrix, exec::ThreadPool *pool);


template void MSk::eval::Queue<float>::evaluate_batch(matrix::Elementary<float> &matrix, std::vector<int> eval_req_ids);
template void MSk::eval::Queue<float>::evaluate_batch(matrix::Dense<float, 2> &matrix, std::vector<int> eval_req_ids);
template void MSk::eval::Queue<float>::evaluate_batch(matrix::Dense<float, 3> &matrix, std::vector<int> eval_req_ids);
template void MSk::eval::Queue<float>::evaluate_batch(matrix::BiPoint<float, float, 2> &matrix, std::vector<int> eval_req_ids);
template void MSk::eval::Queue<float>::evaluate_batch(matrix::BiPoint<float, float, 3> &matrix, std::vector<int> eval_req_ids);

template void MSk::eval::Queue<double>::evaluate_batch(matrix::Elementary<double> &matrix, std::vector<int> eval_req_ids);
template void MSk::eval::Queue<double>::evaluate_batch(matrix::Dense<double, 2> &matrix, std::vector<int> eval_req_ids);
template void MSk::eval::Queue<double>::evaluate_batch(matrix::Dense<double, 3> &matrix, std::vector<int> eval_req_ids);
template void MSk::eval::Queue<double>::evaluate_batch(matrix::BiPoint<double, double, 2> &matrix, std::vector<int> eval_req_ids);
template void MSk::eval::Queue<double>::evaluate_batch(matrix::BiPoint<double, double, 3> &matrix, std::vector<int> eval_req_ids);

template void MSk::eval::Queue<std::complex<float>>::evaluate_batch(matrix::Elementary<std::complex<float>> &matrix, std::vector<int> eval_req_ids);
template void MSk::eval::Queue<std::complex<float>>::evaluate_batch(matrix::Dense<std::complex<float>, 2> &matrix, std::vector<int> eval_req_ids);
template void MSk::eval::Queue<std::complex<float>>::evaluate_batch(matrix::Dense<std::complex<float>, 3> &matrix, std::vector<int> eval_req_ids);
template void MSk::eval::Queue<std::complex<float>>::evaluate_batch(matrix::BiPoint<std::complex<float>, float, 2> &matrix, std::vector<int> eval_req_ids);
template void MSk::eval::Queue<std::complex<float>>::evaluate_batch(matrix::BiPoint<std::complex<float>, float, 3> &matrix, std::vector<int> eval_req_ids);

template void MSk::eval::Queue<std::complex<double>>::evaluate_batch(matrix::Elementary<std::complex<double>> &matrix, std::vector<int> eval_req_ids);
template void MSk::eval::Queue<std::complex<double>>::evaluate_batch(matrix::Dense<std::complex<double>, 2> &matrix, std::vector<int> eval_req_ids);
template void MSk::eval::Queue<std::complex<double>>::evaluate_batch(matrix::Dense<std::complex<double>, 3> &matrix, std::vector<int> eval_req_ids);
template void MSk::eval::Queue<std::complex<double>>::evaluate_batch(matrix::BiPoint<std::complex<double>, double, 2> &matrix, std::vector<int> eval_req_ids);
template void MSk::eval::Queue<std::complex<double>>::evaluate_batch(matrix::BiPoint<std::complex<double>, double, 3> &matrix, std::vector<int> eval_req_ids);
