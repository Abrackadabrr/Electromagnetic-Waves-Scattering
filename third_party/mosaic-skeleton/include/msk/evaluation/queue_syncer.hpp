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
// Created by b.valiakhmetov on 10/14/22.
//


#ifndef MOSAIC_SKELETON_EVALUATION_QUEUE_SYNCER_HPP
#define MOSAIC_SKELETON_EVALUATION_QUEUE_SYNCER_HPP


#include "msk_config.h"

#ifdef MSK_USE_MPI
#include "cpp_mpi/MPI.hpp"

#include <sstream>

#include "msk/utility/logging.hpp"
#include "msk/evaluation/queue.hpp"


constexpr int Q_TAG_INFO = 0;
constexpr int Q_TAG_MAKE = 1;
constexpr int Q_TAG_DATA = 2;

constexpr int Q_MSG_RECV = 0;
constexpr int Q_MSG_SEND = 1;
constexpr int Q_MSG_FNSH = 2;
constexpr int Q_MSG_NONE = 3;


namespace MSk {

namespace eval {


template<class DataType>
class QueueSyncer {
public:

private:
    Queue<DataType> *eq_ptr;
    /* tag:
     * Q_TAG_INFO - inform master (ask/want to give/ elements)
     * Q_TAG_MAKE - master forces to exchange elements
     * Q_TAG_DATA - data send/receive
     *
     * INFO: id -> 0: (ask to receive(0)/send(1)/finish(2), id, cur_n_elems)
     * MAKE: 0 -> id: (should receive(0)/send(1)/finish(2), partner_id, n_elems)
     * DATA: id1 -> id2: ...
     */

    enum class Fill {
        empty,
        low,
        normal,
        high
    };

    std::atomic<Fill> fill{Fill::empty};
    std::atomic<bool> local_finish{false}, to_inform{false};
    uint64_t low_wm{2 * 2048}, high_wm{10 * 2048}; // watermarks of queue's low or high load

    // communication of queues
    MPI_Comm comm{MPI_COMM_WORLD}; // communicator for all queues running the same matrix
    int id, n_procs;
    static constexpr int msg_len{3};
    int info_msg[msg_len], make_msg[msg_len];
    MPI_Request make_req;
    std::vector<MPI_Request> q_reqs; // master-queue gather info from other queues'
    SmartArray<int> master_msg_data; // master: messages from all procs, non-master: 1 message from master
    SmartArray<int> all_q_levels, all_q_queries; // info gathered on master

    std::shared_ptr<log::Logger> log{nullptr};

public:
    QueueSyncer(Queue<DataType> *ptr, uint64_t low, uint64_t high,
                std::shared_ptr<log::Logger> logger = nullptr) :
            eq_ptr{ptr}, low_wm{low}, high_wm{high}, log{logger} {
        MPI_Comm_size(comm, &n_procs);
        MPI_Comm_rank(comm, &id);
        info_msg[0] = info_msg[1] = info_msg[2] = -2;
        MPI_Irecv(&make_msg[0], msg_len, MPI_INT, 0, Q_TAG_MAKE, comm, &make_req);
        if (id == 0) {
            master_msg_data.resize(n_procs * msg_len);
            q_reqs.resize(n_procs);
            all_q_levels.resize(n_procs);
            all_q_queries.resize(n_procs);
            for (int i = 0; i < n_procs; ++i) {
                all_q_levels[i] = all_q_queries[i] = -1;
                MPI_Irecv(master_msg_data.data() + msg_len * i, msg_len, MPI_INT,
                          i, Q_TAG_INFO, comm, q_reqs.data() + i);
            }
        }
    }

    void run() {
        if (log) log->info("Start EQSyncer run");
        int iter{0};
        using State = typename Queue<DataType>::State;
        State eq_state{eq_ptr->get_state()};
        while (eq_state != State::stopping && eq_state != State::stopped) {
            if (log)
                log->debug("Iter {}, Q state:: running: {}, sleeping: {}, stopping: {}, stopped {}", iter++,
                           eq_state == State::running, eq_state == State::sleeping,
                           eq_state == State::stopping, eq_state == State::stopped);
            check_sync();
            std::this_thread::sleep_for(std::chrono::microseconds(100));
            eq_state = eq_ptr->get_state();
        }
        if (log) log->info("End EQSyncer run");
    }

    void inform_master() {
        to_inform = true;
    }

    void finish() {
        local_finish = true;
        to_inform = true;
    }

private:
    int _inform_master() {
        uint64_t e_left{eq_ptr->elements_left()};
        if (local_finish) {
            info_msg[0] = Q_MSG_FNSH;
            info_msg[1] = id, info_msg[2] = e_left;
            if (log) log->debug("{}: inform master: FINISH", id);
            MPI_Send(&info_msg[0], msg_len, MPI_INT, 0, Q_TAG_INFO, comm);
            return info_msg[0];
        }
        Fill new_fill;
        if (e_left == 0) {
            new_fill = Fill::empty;
        } else if (e_left < low_wm) {
            new_fill = Fill::low;
        } else if (e_left <= high_wm) {
            new_fill = Fill::normal;
        } else {
            new_fill = Fill::high;
        }
        if (new_fill != fill && new_fill != Fill::normal) {
            fill = new_fill;
            if (fill == Fill::high) {
                info_msg[0] = Q_MSG_SEND;
                if (log) log->debug("{}: inform master: SEND", id);
            } else {
                info_msg[0] = Q_MSG_RECV;
                if (log) log->debug("{}: inform master: RECEIVE", id);
            }
            info_msg[1] = id, info_msg[2] = e_left;
            MPI_Send(&info_msg[0], msg_len, MPI_INT, 0, Q_TAG_INFO, comm);
        }
        to_inform = false;
        return info_msg[0];
    }

    void check_master() {
        SmartArray<int> ready_idx(n_procs);
        int num_ready{};
        SmartArray<MPI_Status> stats(n_procs);
        MPI_Testsome(n_procs, q_reqs.data(), &num_ready, ready_idx.data(), stats.data());
        if (log) log->debug("Master {} / {} reqs ready", num_ready, n_procs);
        if (num_ready == MPI_UNDEFINED) {
            throw std::runtime_error("EQ requests are not set");
        } else if (num_ready == 0) {
            return;
        }
        for (int i = 0; i < num_ready; ++i) {
            int idx = ready_idx[i];
            int deal{master_msg_data[msg_len * idx]};
            int partner{master_msg_data[msg_len * idx + 1]};
            int n_elems{master_msg_data[msg_len * idx + 2]};
            all_q_levels[idx] = n_elems;
            if (deal == Q_MSG_RECV) { // ask
                all_q_queries[idx] = 0;
            } else if (deal == Q_MSG_SEND) { // wants to give
                all_q_queries[idx] = 1;
            } else if (deal == Q_MSG_FNSH) {
                all_q_queries[idx] = 2;
            } else {
                throw std::runtime_error("Invalid deal tag");
            }
            if (log) log->debug("Master received from {}: deal {}, partner {}, {} elems", idx, deal, partner, n_elems);
            MPI_Irecv(master_msg_data.data() + msg_len * idx, msg_len, MPI_INT,
                      idx, Q_TAG_INFO, comm, q_reqs.data() + idx);
        }
        int i0 = 0, i1 = 0;
        while (true) {
            for (; i0 < n_procs && all_q_queries[i0] == 1; ++i0);
            for (; i1 < n_procs && all_q_queries[i1] != 1; ++i1);
            if (i0 == n_procs || i1 == n_procs) {
                break;
            }
            SmartArray<int> answers(2 * msg_len);
            int n_elems = (all_q_levels[i1] - all_q_levels[i0]) / 2;
//                n_elems = std::min(100, all_q_levels[i1]);
            if (log) log->debug("master connect: {} batch, from {} to {}", 1, i1, i0);
            answers[0] = Q_MSG_RECV, answers[1] = i1, answers[2] = n_elems;
            answers[3] = Q_MSG_SEND, answers[4] = i0, answers[5] = n_elems;
            MPI_Send(&answers[0], msg_len, MPI_INT, i0, Q_TAG_MAKE, comm);
            MPI_Send(&answers[3], msg_len, MPI_INT, i1, Q_TAG_MAKE, comm);
            all_q_queries[i0] = all_q_queries[i1] = -1;
            ++i0, ++i1;
        }
        int n_finished{0};
        for (int i = 0; i < n_procs; ++i) {
            if (all_q_queries[i] == 2) {
                ++n_finished;
            }
        }
        if (log) log->debug("Finished: {} / {}", n_finished, n_procs);
        if (n_finished == n_procs) {
            for (int i = 0; i < n_procs; ++i) {
                if (log) log->debug("Master allows to finish: {}", i);
                SmartArray<int> answer(msg_len);
                answer[0] = Q_MSG_FNSH, answer[1] = i, answer[2] = 0;
                MPI_Send(&answer[0], msg_len, MPI_INT, i, Q_TAG_MAKE, comm);
            }
        }
    }

    void check_sync() {
        if (to_inform) {
            _inform_master();
        }
        if (id == 0) {
            check_master();
        }
        MPI_Status stat;
        int flag{};
        MPI_Test(&make_req, &flag, &stat);
        if (flag) { // master informs what to do
            if (log) log->debug("Message from master {}", make_msg[0]);
            if (make_msg[0] == Q_MSG_FNSH) {
                eq_ptr->stop(true);
                return;
            }
            exchange_elements(make_msg[0], make_msg[1], make_msg[2]);
            MPI_Irecv(&make_msg[0], msg_len, MPI_INT,
                      0, Q_TAG_MAKE, comm, &make_req);
        }
    }

    void exchange_elements(int deal, int partner_id, int n_elems) {
        if (n_elems <= 0) {
            return;
        }
        if (deal == Q_MSG_RECV) { // receive
            MPI_Status stat;
            if (log) log->debug("{}: wait {} elems from {}", id, n_elems, partner_id);
            MPI_Probe(partner_id, Q_TAG_DATA, comm, &stat);
            int data_len;
            MPI_Get_count(&stat, MPI_UINT8_T, &data_len);
            SmartArray<uint8_t> data(data_len);
            MPI_Recv(data.data(), data_len, MPI_UINT8_T, partner_id, Q_TAG_DATA, comm, &stat);
            uint8_t *ptr = data.data();
            uint64_t n_reqs = reinterpret_cast<uint64_t *>(ptr)[0];
            ptr += sizeof(n_reqs);
            if (log) log->debug("{} received: {} reqs, {} elems", id, n_reqs);
            if (n_reqs == 0) {
                return;
            }
            for (uint64_t i = 0; i < n_reqs; ++i) {
                Request<DataType> request;
                uint64_t p_size = request.unpack(ptr);
                eq_ptr->add_elements(std::move(request));
                ptr += p_size;
            }
            eq_ptr->awaken();
        } else if (deal == Q_MSG_SEND) { // send
            std::vector<std::pair<uint64_t, uint64_t>> eval_elems;
            std::vector<Request<DataType>>
            eval_reqs;
            eq_ptr->extract_batch(eval_reqs);
            inform_master();
            uint64_t n_reqs = eval_reqs.size();
            int data_len = sizeof(n_reqs);
            for (uint64_t i = 0; i < n_reqs; ++i) {
                data_len += eval_reqs[i].packed_size();
            }

            SmartArray<uint8_t> data(data_len);
            uint8_t *ptr = data.data();
            reinterpret_cast<uint64_t *>(ptr)[0] = n_reqs;
            ptr += sizeof(n_reqs);
            for (uint64_t i = 0; i < n_reqs; ++i) {
                auto &req = eval_reqs[i];
                uint64_t p_size = req.pack(ptr);
                ptr += p_size;
            }
            if (log) log->debug("{}: send {} elems to {}", id, n_elems, partner_id);
            MPI_Send(data.data(), data_len, MPI_UINT8_T, partner_id, Q_TAG_DATA, comm);
            if (log) log->debug("{}: send success", id);
        } else {
            throw std::runtime_error("Invalid deal tag");
        }
    }
};

}
}

#endif

#endif //MOSAIC_SKELETON_EVALUATION_QUEUE_SYNCER_HPP
