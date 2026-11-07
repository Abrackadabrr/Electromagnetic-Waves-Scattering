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
// Created by b.valiakhmetov on 4/4/24.
//


#ifndef MOSAIC_SKELETON_ORACLE_ELEMENTARY_HPP
#define MOSAIC_SKELETON_ORACLE_ELEMENTARY_HPP


#include "msk/oracle/oracle.hpp"
#include "msk/computation/dense.hpp"
#include "msk/computation/lowrank.hpp"
#include "msk/block/dense.hpp"
#include "msk/block/lowrank.hpp"


namespace MSk {

namespace oracle {


/**
 * @brief Mosaic-skeleton format generating Oracle for functional matrices with linear equidistannt geometry.
 *
 * @details A matrix is considered as an actual function of two integer indices,
 * that do not (may not) correspond to points in space.
 *
 * @tparam DataType Matrix element type.
 */
template<class DataType>
class Elementary : public Base<DataType> {
private:
    using Idx = typename Base<DataType>::Idx;
    using TBlock = typename Base<DataType>::TBlock;
    using IdxIt = typename Base<DataType>::IdxIt;

    /**
     * if rho == 1.0, then all not identical blocks are low-rank (like in HSS);
     * if rho == 2.0, then all not identical or not neighbouring blocks are low-rank.
     */
    double rho{1.0}; ///< Distance parameter

public:
    /**
     * @brief Set block division parameters.
     *
     * @param[in] params Oracle Parameters.
     */
    explicit Elementary(const Parameters &params) : Base<DataType>(params), rho{params.rho} {}

    ~Elementary() override = default;

    /**
     * @brief Get the action shat should be performed with a block of the matrix.
     *
     * @details It is the main Oracle purpose: to decide if a block should be computed in some format
     * (low-rank or dense) or should be divided into several subblocks.
     *
     * @param[in] block_it Pointer to a block to process.
     * @return ActionAnswer for the block.
     */
    ActionAnswer<DataType> get_action(TBlock *block_it) override {
        std::shared_ptr<log::Logger> log = this->log;
        if (block_it == nullptr) {
            if (log) log->error("Action on nullptr Block");
            return ActionAnswer<DataType>{};
        }
        if (log) log->debug("get action for block {}", block_it->id());
        uint64_t id{block_it->id()}, i0, i1, j0, j1, m, n;
        std::tie(i0, j0) = block_it->corner();
        std::tie(m, n) = block_it->size();
        i1 = i0 + m, j1 = j0 + n;
        if (log) log->debug("Action on Block {} of size {} x {}", id, m, n);

        // Check if the block is compressible
        bool lrank_cond{false};
        if (rho == 1.0) { // Check if blocks do not overlap (HSS mosaic)
            lrank_cond = (i1 <= j0) || (j1 <= i0);
        } else if (rho == 2.0) { // Check if blocks are not neighbouring
            lrank_cond = (i1 < j0) || (j1 < i0);
        } else {
            // Calculate centers (doubled) of row and column ranges
            uint64_t row_center2{i0 + i1}, col_center2{j0 + j1};
            // Calculate the distance between centers (doubled)
            uint64_t dist_centers2{row_center2 > col_center2 ? row_center2 - col_center2 : col_center2 - row_center2};
            lrank_cond = dist_centers2 >= rho * 2 * (m + n);
        }
        // Mark block as low rank
        const uint64_t &min_row_size{this->min_row_size}, &min_col_size{this->min_col_size};
        if (lrank_cond || (this->_force_lrank && m <= min_row_size && n <= min_col_size)) {
            if (log) log->debug("LOWRANK Block {}: ({}, {}) x ({}, {})", id, i0, i1, j0, j1);
            auto db = new block::LowRank<DataType>(m, n);
            comp::Computer<DataType> *comp = nullptr;
            if (std::min(m, n) > 128) {
                comp =  new comp::Cross<DataType>(db, this->block_params, this->block_recomp, this->block_log);
            } else {
                comp = new comp::SVD<DataType>(db, this->block_params, true);
            }
            block_it->set_data(db);
            block_it->set_computer(comp);
            return {db, comp};
        }

        // Mark block as dense
        if (m <= min_row_size && n <= min_col_size) {
            if (log) log->debug("DENSE Block {}: ({}, {}) x ({}, {})", id, i0, i1, j0, j1);
            auto db = new block::Dense<DataType>(m, n);
            auto comp = new comp::Dense<DataType>(db);
            block_it->set_data(db);
            block_it->set_computer(comp);
            return {db, comp};
        }

        // Divide block into subblocks
        Idx rows = block_it->rows();
        Idx cols = block_it->columns();

        ActionAnswer<DataType> result{};
        std::vector<TBlock> children(4);
        uint64_t up_blocks{0};
        if (m > min_row_size && n > min_col_size) { // Slice in half by both axes
            Idx rows1, rows2, cols1, cols2;
            IdxIt row_mid = rows.begin() + m / 2;
            IdxIt col_mid = cols.begin() + n / 2;
            std::tie(rows1, rows2) = rows.split(row_mid);
            std::tie(cols1, cols2) = cols.split(col_mid);

            up_blocks = 4;
            children[0] = std::move(TBlock(rows1, cols1));
            children[1] = std::move(TBlock(rows1, cols2));
            children[2] = std::move(TBlock(rows2, cols1));
            children[3] = std::move(TBlock(rows2, cols2));
        } else if (m > min_row_size) { // Slice by rows
            Idx rows1, rows2;
            IdxIt row_mid = rows.begin() + m / 2;
            std::tie(rows1, rows2) = rows.split(row_mid);

            up_blocks = 2;
            children[0] = std::move(TBlock(rows1, cols));
            children[1] = std::move(TBlock(rows2, cols));
        } else { // Slice by columns, n > min_block_size
            Idx cols1, cols2;
            IdxIt col_mid = cols.begin() + n / 2;
            std::tie(cols1, cols2) = cols.split(col_mid);

            up_blocks = 2;
            children[0] = std::move(TBlock(rows, cols1));
            children[1] = std::move(TBlock(rows, cols2));
        }
        if (this->_symmetry) {
            for (uint64_t i = 0; i < up_blocks; ++i) {
                uint64_t i0, j0;
                std::tie(i0, j0) = children[i].corner();
                if (i0 > j0) {
                    std::swap(children[i], children[up_blocks - 1]);
                    --up_blocks;
                }
            }
        }
        result = {up_blocks, std::move(children)};
        // Disable bermutation for output
        if (log) log->debug("DIVIDE Block {}: ({}, {}) x ({}, {}) into {}", id, i0, i1, j0, j1, result.n_child);
        return result;
    }

};


extern template class Elementary<float>;

extern template class Elementary<double>;

extern template class Elementary<std::complex<float>>;

extern template class Elementary<std::complex<double>>;


}
}


#endif //MOSAIC_SKELETON_ORACLE_ELEMENTARY_HPP
