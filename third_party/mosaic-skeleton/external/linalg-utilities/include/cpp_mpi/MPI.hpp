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
// Created by b.valiakhmetov on 31.03.23.
//


#ifndef LINALG_UTILITIES_CPP_MPI_MPI_HPP
#define LINALG_UTILITIES_CPP_MPI_MPI_HPP


#include <mpi.h>

#include <cstdint>
#include <complex>


namespace MPI {
    template<class DataType>
    MPI_Datatype data_type() {
        if (std::is_same<DataType, float>())
        {
            return MPI_FLOAT;
        }
        if (std::is_same<DataType, double>())
        {
            return MPI_DOUBLE;
        }
        if (std::is_same<DataType, std::complex<float>>())
        {
            return MPI_C_FLOAT_COMPLEX;
        }
        if (std::is_same<DataType, std::complex<double>>())
        {
            return MPI_C_DOUBLE_COMPLEX;
        }
        if (std::is_same<DataType, uint64_t>())
        {
            return MPI_UINT64_T;
        }
        if (std::is_same<DataType, uint32_t>())
        {
            return MPI_UINT32_T;
        }
        if (std::is_same<DataType, int>())
        {
            return MPI_INT;
        }
        throw std::runtime_error("Unexpected datatype");
    }

    template<typename DataType>
    int Ireduce_sum(const DataType *out_buf, DataType *in_buf, uint64_t count, uint64_t root,
                    MPI_Comm comm, MPI_Request *request);

    template<class DataType>
    int Allreduce_inplace_sum(DataType *sendrecvbuf, uint64_t count, const MPI_Comm &comm);

    template<typename DataType>
    int Allgatherv_inplace(DataType *send_recv_buf, uint64_t count, const int *sizes, const int *offs,
                           const MPI_Comm &comm);


    // template resolutions
    template<>
    int Ireduce_sum(const float *out_buf, float *in_buf, uint64_t count, uint64_t root,
                    MPI_Comm comm, MPI_Request *request);

    template<>
    int Ireduce_sum(const double *out_buf, double *in_buf, uint64_t count, uint64_t root,
                    MPI_Comm comm, MPI_Request *request);

    template< >
    int Ireduce_sum(const std::complex<float> *out_buf, std::complex<float> *in_buf, uint64_t count, uint64_t root,
                    MPI_Comm comm, MPI_Request *request);

    template< >
    int Ireduce_sum(const std::complex<double> *out_buf, std::complex<double> *in_buf, uint64_t count, uint64_t root,
                    MPI_Comm comm, MPI_Request *request);


    template<>
    int Allreduce_inplace_sum(float *sendrecvbuf, uint64_t count, const MPI_Comm &comm);

    template<>
    int Allreduce_inplace_sum(double *sendrecvbuf, uint64_t count, const MPI_Comm &comm);

    template<>
    int Allreduce_inplace_sum(std::complex<float> *sendrecvbuf, uint64_t count, const MPI_Comm &comm);

    template<>
    int Allreduce_inplace_sum(std::complex<double> *sendrecvbuf, uint64_t count, const MPI_Comm &comm);


    template<>
    int Allgatherv_inplace(float *send_recv_buf, uint64_t count, const int *sizes, const int *offs,
                           const MPI_Comm &comm);

    template<>
    int Allgatherv_inplace(double *send_recv_buf, uint64_t count, const int *sizes, const int *offs,
                           const MPI_Comm &comm);

    template<>
    int Allgatherv_inplace(std::complex<float> *send_recv_buf, uint64_t count, const int *sizes, const int *offs,
                           const MPI_Comm &comm);

    template<>
    int Allgatherv_inplace(std::complex<double> *send_recv_buf, uint64_t count, const int *sizes, const int *offs,
                           const MPI_Comm &comm);
}

#endif //LINALG_UTILITIES_CPP_MPI_MPI_HPP
