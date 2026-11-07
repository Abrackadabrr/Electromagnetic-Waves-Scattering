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


#ifndef MOSAIC_SKELETON_ORACLE_DISTANT_HPP
#define MOSAIC_SKELETON_ORACLE_DISTANT_HPP


#include <cstdint>
#include <fstream>

#include "msk/oracle/oracle.hpp"
#include "msk/computation/dense.hpp"
#include "msk/computation/lowrank.hpp"
#include "msk/block/dense.hpp"
#include "msk/block/lowrank.hpp"
#include "msk/geometry/metrics.hpp"


namespace MSk {

namespace oracle {


/**
 * @brief Tree of coordinate line multilevel division.
 *
 * @tparam Iterator Segment dividers (cuts) type.
 */
template<class Iterator>
class DivisionTree {
private:
    std::vector<Iterator> _dividers; ///< Heap of dividers
    std::vector<bool> _is_set; ///< Marks of heap's items existence

public:
    using iterator_type = typename std::vector<Iterator>::iterator; ///< Tree iterator type

public:
    /**
     * @brief Initialize empty tree - no dividers on the line.
     */
    DivisionTree() : _dividers{}, _is_set{} {}

    /**
     * @brief Find if a range has already been divided.
     *
     * @param[in] range_b Range starting iterator.
     * @param[in] range_e Range ending iterator.
     * @return Segment dividing iterator or range_e if segment has no division yet.
     */
    Iterator find_division(Iterator range_b, Iterator range_e) {
        uint64_t idx{0};
        if (range_e - range_b == 0) {
            return range_e;
        }
        // Start binary search of the first divider DIV such that: range_b < DIV < range_e
        while (idx < _dividers.size()) {
            if (!_is_set[idx]) { // If no divider on the suitable place
                return range_e;
            }
            if (range_e <= _dividers[idx]) { // Current DIV: range_e <= DIV, go left
                idx = 2 * idx + 1;
            } else if (_dividers[idx] <= range_b) { // Current DIV: DIV <= range_b, go right
                idx = 2 * idx + 2;
            } else { // Suitable divider is found
                return _dividers[idx];
            }
        }
        return range_e;
    }

    /**
     * @brief Insert new divider to the tree
     *
     * @param[in] div Divider.
     */
    void divide(Iterator div) {
        if (_dividers.empty()) { // The first divider
            _is_set.resize(1);
            _is_set[0] = true;
            _dividers.resize(1);
            _dividers[0] = div;
            return;
        }
        // Find the smallest segment that the given divider is inside
        uint64_t idx{0};
        while (idx < _dividers.size()) { // Found suitable place, insert
            if (!_is_set[idx]) {
                _is_set[idx] = true;
                _dividers[idx] = div;
                return;
            }
            if (div < _dividers[idx]) { // Go left
                idx = 2 * idx + 1;
            } else if (_dividers[idx] < div) { // Go right
                idx = 2 * idx + 2;
            } else { // div == _dividers[idx] but it has been already set
                return;
            }
        }
        // Heap is full, resize it and set the divider
        _is_set.resize(idx + 1, false);
        _is_set[idx] = true;
        _dividers.resize(idx + 1);
        _dividers[idx] = div;
    }

    /**
     * @brief Dump dividers tree to a file.
     *
     * @details Data is appended to a file. The heap is written directly,
     * unset positions are replaced with '_' symbol.
     *
     * @param[in] filename File to write to.
     * @param[in] pre_msg Message to insert right before the tree.
     */
    void dump_dividers(const std::string &filename, const std::string &pre_msg = "") {
        std::ofstream out(filename, std::ios_base::app);
        out << pre_msg << _dividers.size() << "\n";
        for (uint64_t idx = 0; idx < _dividers.size(); ++idx) {
            if (_is_set[idx]) {
                out << *_dividers[idx] << " ";
            } else {
                out << "_" << " ";
            }
        }
        out << std::endl;
        out.close();
    }
};




/**
 * @brief Mosaic-skeleton format generating Oracle for functional matrices
 * with arbitrary spatial points geometry.
 *
 * @tparam DataType Matrix element type.
 * @tparam Matrix Matrix type.
 */
template<class DataType, class Matrix>
class Distant final : public Base<DataType> {
private:
    /**
     * @brief Policy for a matrix block division.
     */
    enum class BoundPolicy {
        rect, ///< Bound points with rectangles
        mean_sphere, ///< Bound points with a sphere centered in the centroid
        med_sphere, ///<  Bound points with a sphere centered in the bounding box center
    };
private:
    using obj_type = typename Matrix::obj_type;
    using Idx = typename Base<DataType>::Idx;
    using TBlock = typename Base<DataType>::TBlock;
    using IdxIt = typename Base<DataType>::IdxIt;

    const Matrix &_matrix;

    double rho{1.0}; ///< Distance parameter
    uint64_t p_dim; ///< Hoelder norm parameter for distance calculation
    BoundPolicy bnd_policy{BoundPolicy::rect}; ///< Distance and size measurement policy

    DivisionTree<IdxIt> row_tree, col_tree; ///< Indices division tree structure

private:
    /**
     * @brief Calculate centroid (mass center) and the farthest point for an index range of points.
     *
     * @param[in] it_begin Points index range begining iterator.
     * @param[in] it_end Points index range ending iterator.
     * @param[in] is_row Is it row points or column.
     * @return Center and radial point.
     */
    std::pair<obj_type, obj_type> center_radius(const IdxIt &it_begin, const IdxIt &it_end, const bool is_row) {
        obj_type center{};
        obj_type radius{};
        std::function<obj_type(IdxIt)> point_getter;
        if (is_row) {
            point_getter = [this](IdxIt it) {
                return _matrix.row_point(*it);
            };
        } else {
            point_getter = [this](IdxIt it) {
                return _matrix.col_point(*it);
            };
        }
        if (bnd_policy == BoundPolicy::rect) {
            obj_type low, high;
            std::tie(low, high) = geom::bound_box(point_getter, it_begin, it_end);
            center = (low + high) * 0.5;
            radius = (high - low) * 0.5;
        } else {
            if (bnd_policy == BoundPolicy::mean_sphere) {
                center = geom::centroid(point_getter, it_begin, it_end);
            } else { // if (bnd_policy == BoundPolicy::med_sphere)
                obj_type low, high;
                std::tie(low, high) = geom::bound_box(point_getter, it_begin, it_end);
                center = (low + high) * 0.5;
            }
            radius = geom::bound_sphere(point_getter, it_begin, it_end, center, p_dim);
        }
        return {center, radius};
    }

    /**
     * @brief Divide a range of points into two sets with the hyperplane.
     *
     * @details The hyperplane is built so it passes throught the center in the direction
     * orthogonal to farthest poit vector. In that way the set of points is sliced
     * orthiginal to the direction it is "elongated".
     *
     * @param[in] it_begin Points index range begining iterator.
     * @param[in] it_end Points index range ending iterator.
     * @param[in] is_row Is it row points or column.
     * @param[in] center Centroid of points range.
     * @param[in] radius Radially farthest point the center.
     * @return Dividing point index iterator.
     */
    IdxIt divide(const IdxIt &it_begin, const IdxIt &it_end,
                 const bool is_row, const obj_type &center, const obj_type &radius) {
        IdxIt mid;

        if (is_row) {
            mid = row_tree.find_division(it_begin, it_end);
        } else {
            mid = col_tree.find_division(it_begin, it_end);
        }
        if (mid == it_end) { // get new division
            obj_type normal{};
            std::function<obj_type(IdxIt)> point_getter;
            if (is_row) {
                point_getter = [this](IdxIt it) {
                    return _matrix.row_point(*it);
                };
            } else {
                point_getter = [this](IdxIt it) {
                    return _matrix.col_point(*it);
                };
            }
            if (bnd_policy == BoundPolicy::rect) {
                typename obj_type::value_type max_width{0};
                uint64_t max_dir{0};
                for (uint64_t i = 0; i < obj_type().size(); ++i) {
                    if (std::abs(radius[i]) > max_width) {
                        max_width = std::abs(radius[i]);
                        max_dir = i;
                    }
                }
                normal = obj_type(0);
                normal[max_dir] = 1;
            }/* else if (bnd_policy == BoundPolicy::mean_sphere) {
                normal = (1.0 / radius.norm()) * radius;
            }*/ else {
                normal = geom::iter_elong(point_getter, it_begin, it_end, center+ radius);
            }

            // Perform division
            mid = it_begin;
            for (IdxIt it = it_begin; it != it_end; ++it) {
                // Permute points respectively
                if (point_getter(it).hplane_side(normal, center)) {
                    transpose(mid, it);
                    ++mid;
                }
            }
            // remember division
            if (is_row) {
                row_tree.divide(mid);
            } else {
                col_tree.divide(mid);
            }
        }
        return mid;
    }

    bool lowrank_criteria(const obj_type &row_center, const obj_type &row_radius,
                          const obj_type &col_center, const obj_type &col_radius) const {

        if (bnd_policy == BoundPolicy::rect) {
            obj_type center = (row_center - col_center) * 0.5, radius = row_radius + col_radius;
            obj_type low = center - radius, high = center + radius;
            obj_type zero_dist{0};
            for (uint64_t i = 0; i < low.size(); ++i) {
                if (low[i] * high[i] > 0) {
                    zero_dist[i] = std::min(std::abs(high[i]), std::abs(low[i]));
                }
            }
            return zero_dist.norm(p_dim) > (rho - 1) * std::min(row_radius.norm(p_dim), col_radius.norm(p_dim));//(rho - 1) * radius.norm(p_dim);
        }
        // Spherical
        return (row_center - col_center).norm(p_dim) >= rho * (row_radius.norm(p_dim) + col_radius.norm(p_dim));
    }


public:
    /**
     * @brief Set block division parameters.
     *
     * @param[in] matrix Matrix proving points.
     * @param[in] params Oracle Parameters.
     */
    Distant(const Matrix &matrix, const Parameters &params) : Base<DataType>(params),
            _matrix{matrix}, rho{params.rho}, p_dim{params.sphere_norm}, row_tree{}, col_tree{} {}

    ~Distant() override = default;

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

        if (!this->_manual_actions.empty()) {
            for (auto &preset : this->_manual_actions) {
                if (preset.block == *block_it) {
                    auto ans = std::move(preset.answer);
                    if (ans.n_child > 0) {
                        for (auto &block : ans.child) {
                            block.rows().permute(block_it->rows().permutation());
                            block.columns().permute(block_it->columns().permutation());
                        }
                    }
                    return ans;
                }
            }
        }

        uint64_t id{block_it->id()}, i0, i1, j0, j1, m, n;
        std::tie(i0, j0) = block_it->corner();
        std::tie(m, n) = block_it->size();
        i1 = i0 + m, j1 = j0 + n;
        if (log) log->debug("Action on Block {} of size {} x {}", id, m, n);

        const Idx &rows = block_it->rows();
        const Idx &cols = block_it->columns();

        obj_type row_center{};
        obj_type row_radius{};
        std::tie(row_center, row_radius) = center_radius(rows.begin(true), rows.end(true), true);

        obj_type col_center{};
        obj_type col_radius{};
        std::tie(col_center, col_radius) = center_radius(cols.begin(true), cols.end(true), false);

        // Check if compressible
        const uint64_t &min_row_size{this->min_row_size}, &min_col_size{this->min_col_size};
        if (lowrank_criteria(row_center, row_radius, col_center, col_radius) ||
            (this->_force_lrank && m <= min_row_size && n <= min_col_size)) {
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

        // If Dense
        if (m <= min_row_size && n <= min_col_size) {
            if (log) log->debug("DENSE Block {}: ({}, {}) x ({}, {})", id, i0, i1, j0, j1);
            auto db = new block::Dense<DataType>(m, n);
            auto comp = new comp::Dense<DataType>(db);
            block_it->set_data(db);
            block_it->set_computer(comp);
            return {db, comp};
        }

        bool divide_rows{m > min_row_size}, divide_cols{n > min_col_size};
        // Split by relatively longest side if policy is single
        if (this->policy == DividePolicy::single && divide_rows && divide_cols) {
            if (m * min_col_size > n * min_row_size) {
                divide_cols = false;
            } else {
                divide_rows = false;
            }
        }
        ActionAnswer<DataType> result{};
        uint64_t up_blocks{0};
        std::vector<TBlock> children(4);
        if (divide_rows && divide_cols) {
            Idx rows1, rows2, cols1, cols2;
            IdxIt row_mid = divide(rows.begin(true), rows.end(true), true, row_center, row_radius);
            IdxIt col_mid = divide(cols.begin(true), cols.end(true), false, col_center, col_radius);
            std::tie(rows1, rows2) = rows.split(row_mid);
            std::tie(cols1, cols2) = cols.split(col_mid);

            up_blocks = 4;
            children[0] = std::move(TBlock(rows1, cols1));
            children[1] = std::move(TBlock(rows1, cols2));
            children[2] = std::move(TBlock(rows2, cols1));
            children[3] = std::move(TBlock(rows2, cols2));
        } else if (divide_rows) {
            Idx rows1, rows2;
            IdxIt row_mid = divide(rows.begin(true), rows.end(true), true, row_center, row_radius);
            std::tie(rows1, rows2) = rows.split(row_mid);

            up_blocks = 2;
            children[0] = std::move(TBlock(rows1, cols));
            children[1] = std::move(TBlock(rows2, cols));
        } else { // divide_cols
            Idx cols1, cols2;
            IdxIt col_mid = divide(cols.begin(true), cols.end(true), false, col_center, col_radius);
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
                    --i;
                }
            }
        }
        children.resize(up_blocks);
        result = {up_blocks, std::move(children)};
        if (log) log->debug("DIVIDE Block {}: ({}, {}) x ({}, {}) into {}", id, i0, i1, j0, j1, result.n_child);
        return result;
    }

    /**
     * @brief Get the Oracle distance factor.
     *
     * @details If for a block inequality:
     * distance(C_r, C_c) >= RHO * (radius(PTS_r) + radius(PTS_c))
     * holds, where PTS_r are rows related points, PTS_c are columns related points,
     * C_r is the centroid os PTS_r, C_c is the ceintroid of PTS_c and RHO is the distance factor,
     * then it is marked as low-rank.
     *
     * @return Distance factor.
     */
    double dist_factor() const {
        return rho;
    }

    /**
     * @brief Set the Oracle distance factor.
     *
     * @param[in] new_factor New factor value.
     */
    void dist_factor(double new_factor) {
        rho = new_factor;
    }

    /**
     * @brief Write rows and columns division trees to a file.
     *
     * @details For details see DivisionTree::dump_dividers().
     *
     * @param[in] filename File name.
     */
    void dump_trees(const std::string &filename) {
        row_tree.dump_dividers(filename, "row\n");
        col_tree.dump_dividers(filename, "column\n");
    }

};


extern template class DivisionTree<index::IndexIterator<>>;


}
}


#endif //MOSAIC_SKELETON_ORACLE_DISTANT_HPP
