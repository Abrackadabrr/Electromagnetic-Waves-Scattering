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


#include "cpp_mpi/MPI.hpp"


namespace MPI {
    template<>
    int Ireduce_sum(const float *out_buf, float *in_buf, uint64_t count, uint64_t root,
                    MPI_Comm comm, MPI_Request *request) {
        return MPI_Ireduce(out_buf, in_buf, count, MPI_FLOAT, MPI_SUM, root, comm, request);
    }

    template<>
    int Ireduce_sum(const double *out_buf, double *in_buf, uint64_t count, uint64_t root,
                    MPI_Comm comm, MPI_Request *request) {
        return MPI_Ireduce(out_buf, in_buf, count, MPI_DOUBLE, MPI_SUM, root, comm, request);
    }

    template< >
    int Ireduce_sum(const std::complex<float> *out_buf, std::complex<float> *in_buf, uint64_t count, uint64_t root,
                    MPI_Comm comm, MPI_Request *request) {
        return MPI_Ireduce(out_buf, in_buf, count, MPI_C_FLOAT_COMPLEX, MPI_SUM, root, comm, request);
    }

    template< >
    int Ireduce_sum(const std::complex<double> *out_buf, std::complex<double> *in_buf, uint64_t count, uint64_t root,
                    MPI_Comm comm, MPI_Request *request) {
        return MPI_Ireduce(out_buf, in_buf, count, MPI_C_DOUBLE_COMPLEX, MPI_SUM, root, comm, request);
    }


    template<>
    int Allreduce_inplace_sum(float *sendrecvbuf, uint64_t count, const MPI_Comm &comm)
    {
        int res = 0;
        for(uint64_t i = 0; i < count; i+=(1<<30))
        {
            int num = std::min(count - i, uint64_t(1) << 30);
            res += MPI_Allreduce(MPI_IN_PLACE, sendrecvbuf+i, num, MPI_FLOAT, MPI_SUM, comm);
        }
        return res;
    }


    template<>
    int Allreduce_inplace_sum(double *sendrecvbuf, uint64_t count, const MPI_Comm &comm) {
        int res = 0;
        for(uint64_t i = 0; i < count; i+=(1<<30))
        {
            int num = std::min(count - i, uint64_t(1) << 30);
            res += MPI_Allreduce(MPI_IN_PLACE, sendrecvbuf+i, num, MPI_DOUBLE, MPI_SUM, comm);
        }
        return res;
    }

    template<>
    int Allreduce_inplace_sum(std::complex<float> *sendrecvbuf, uint64_t count, const MPI_Comm &comm) {
        int res = 0;
        for(uint64_t i = 0; i < count; i+=(1<<30))
        {
            int num = std::min(count - i, uint64_t(1) << 30);
            res += MPI_Allreduce(MPI_IN_PLACE, sendrecvbuf+i, num, MPI_C_FLOAT_COMPLEX, MPI_SUM, comm);
        }
        return res;
    }

    template<>
    int Allreduce_inplace_sum(std::complex<double> *sendrecvbuf, uint64_t count, const MPI_Comm &comm) {
        int res = 0;
        for(uint64_t i = 0; i < count; i+=(1<<30))
        {
            int num = std::min(count - i, uint64_t(1) << 30);
            res += MPI_Allreduce(MPI_IN_PLACE, sendrecvbuf+i, num, MPI_C_DOUBLE_COMPLEX, MPI_SUM, comm);
        }
        return res;
    }


    template<>
    int Allgatherv_inplace(float *send_recv_buf, uint64_t count, const int *sizes, const int *offs,
                           const MPI_Comm &comm) {
        return MPI_Allgatherv(MPI_IN_PLACE, count, MPI_FLOAT,
                              send_recv_buf, sizes, offs, MPI_FLOAT, comm);
    }

    template<>
    int Allgatherv_inplace(double *send_recv_buf, uint64_t count, const int *sizes, const int *offs,
                           const MPI_Comm &comm) {
        return MPI_Allgatherv(MPI_IN_PLACE, count, MPI_DOUBLE,
                              send_recv_buf, sizes, offs, MPI_DOUBLE, comm);
    }

    template<>
    int Allgatherv_inplace(std::complex<float> *send_recv_buf, uint64_t count, const int *sizes, const int *offs,
                           const MPI_Comm &comm) {
        return MPI_Allgatherv(MPI_IN_PLACE, count, MPI_C_FLOAT_COMPLEX,
                              send_recv_buf, sizes, offs, MPI_C_FLOAT_COMPLEX, comm);
    }

    template<>
    int Allgatherv_inplace(std::complex<double> *send_recv_buf, uint64_t count, const int *sizes, const int *offs,
                           const MPI_Comm &comm) {
        return MPI_Allgatherv(MPI_IN_PLACE, count, MPI_C_DOUBLE_COMPLEX,
                              send_recv_buf, sizes, offs, MPI_C_DOUBLE_COMPLEX, comm);
    }
}
