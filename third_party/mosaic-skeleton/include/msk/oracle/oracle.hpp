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
// Created by bulat_v on 24.02.2022.
//


/**
 * @file
 *
 * @brief Contains basic classes implementing Oracle entity.
 */


#ifndef MOSAIC_SKELETON_ORACLE_ORACLE_HPP
#define MOSAIC_SKELETON_ORACLE_ORACLE_HPP


#include <memory>
#include <utility>
#include <vector>

#include "msk/block/block.hpp"
#include "msk/computation/computer.hpp"
#include "msk/block/data.hpp"
#include "msk/utility/logging.hpp"


namespace MSk {

namespace oracle {


/**
 * @brief Policy for a matrix block division.
 */
enum class DividePolicy {
    single, ///< Divide only rows or columns
    dual, ///< Divide both rows and columns
};


/**
 * @brief Strructure with the information about block::Block type and children.
 *
 * @details If no children then the block is treated as ready to be computed
 * (with data and its computer). Otherwise it is to be divided to a several subblocks
 * that should be brocessed independently.
 *
 * @tparam DataType Block element type.
 */
template<class DataType>
struct ActionAnswer {
public:
    using TBlock = block::Block<DataType>; ///< Block type
    uint64_t n_child{0}; ///< Number of block's children
    std::vector<TBlock> child{}; ///< Array of child blocks

    comp::Computer<DataType> *computer{nullptr}; ///< Computer of the block
    block::Data<DataType> *data{nullptr}; ///< Data representation of the block

public:
    explicit ActionAnswer() : n_child{0}, computer{nullptr}, data{nullptr} {}

    /**
     * @brief Initialize struct data.
     *
     * @param[in] n Children number.
     * @param[in,out] children Block's children array to be moved.
     */
    ActionAnswer(uint64_t n, std::vector<TBlock> &&children) : n_child{n}, child{std::move(children)} {}

    /**
     * @brief Initialize struct data.
     *
     * @param[in] db Data representation of the block.
     * @param[in] comp Computer for the data representation.
     */
    ActionAnswer(block::Data<DataType> *db, comp::Computer<DataType> *comp) : computer{comp}, data{db} {}

    ActionAnswer(const ActionAnswer &) = delete;
    ActionAnswer(ActionAnswer &&) = default;
    ActionAnswer &operator=(const ActionAnswer &) = delete;
    ActionAnswer &operator=(ActionAnswer &&) = default;
};


/**
 * @brief Oracle operational parameters.
 */
struct Parameters {
    uint64_t min_rows{128}, min_cols{128}; ///< Sizes of block that shouldn't be divided
    double rho{1.0}; ///< Distance parameter
    uint64_t sphere_norm{0}; ///< Hoelder norm parameter for distance calculation
    DividePolicy policy{DividePolicy::dual}; ///< Blocks division policy
    bool force_lowrank{false}; ///< If to mark undividable blocks as low-rank anyway
    bool symmetry{false}; ///< Matrix is assumed symmetric or hermitian
    std::shared_ptr<log::Logger> log{nullptr}; ///< Logger
    comp::Parameters block_comp{0.0, 1e-3, 0}; ///< Blocks computation parameters
    bool block_recomp{true}; ///< If to recompress blocks right after computation
    std::shared_ptr<log::Logger> block_log{nullptr}; ///< Logger for blocks

    Parameters() = default;

    Parameters(uint64_t min_r, uint64_t min_c, double rho_dist, uint64_t p_norm, DividePolicy div_policy,
               bool force_lrank, bool symmetric, std::shared_ptr<log::Logger> oracle_log,
               comp::Parameters comp_params, bool recompres_block, std::shared_ptr<log::Logger> approx_log) :
               min_rows{min_r}, min_cols{min_c}, rho{rho_dist}, sphere_norm{p_norm},
               policy{div_policy}, force_lowrank{force_lrank}, symmetry{symmetric},
               log{std::move(oracle_log)},
               block_comp{comp_params}, block_recomp{recompres_block}, block_log{std::move(approx_log)} {}
};


/**
 * @brief Mosaic-skeleton format generating Oracle for functional matrices with linear equidistannt geometry.
 *
 * @details A matrix is considered as an actual function of two integer indices,
 * that do not (may not) correspond to points in space.
 *
 * @tparam DataType Matrix element type.
 */
template<class DataType>
class Base {
protected:
    using Idx = index::Index<uint64_t>;
    using TBlock = block::Block<DataType, Idx>;
    using IdxIt = typename Idx::iterator_type;

    uint64_t min_row_size{0}, min_col_size{0}; ///< Sizes of block that shouldn't be divided
    DividePolicy policy{DividePolicy::dual}; ///< Blocks division policy
    bool _force_lrank{false}; ///< If to mark undividable blocks as low-rank anyway

    bool _symmetry{false};

    struct Preset {
        TBlock block;
        ActionAnswer<DataType> answer{};
    };
    std::vector<Preset> _manual_actions{};

    std::shared_ptr<log::Logger> log{nullptr}; ///< Logger

    comp::Parameters block_params; ///< Parameters for block
    bool block_recomp; ///< If to recompress blocks right after computation
    std::shared_ptr<log::Logger> block_log{nullptr}; ///< Logger for blocks

public:
    /**
     * @brief Set block division parameters.
     *
     * @param[in] params Oracle Parameters.
     */
    explicit Base(const Parameters &params) :
            min_row_size{params.min_rows}, min_col_size{params.min_cols}, policy{params.policy},
            _force_lrank{params.force_lowrank}, log{params.log},
            _symmetry{params.symmetry},
            block_params{params.block_comp}, block_recomp{params.block_recomp}, block_log{params.block_log} {}

    virtual ~Base() = default;

    /**
     * @brief Get the action shat should be performed with a block of the matrix.
     *
     * @details It is the main Oracle purpose: to decide if a block should be computed in some format
     * (low-rank or dense) or should be divided into several subblocks.
     *
     * @param[in] block_it Pointer to a block to process.
     * @return ActionAnswer for the block.
     */
    virtual ActionAnswer<DataType> get_action(TBlock *block_it) = 0;

    void add_manual(TBlock &&block, ActionAnswer<DataType> &&answer) {
        Preset preset;
        preset.block = std::move(block);
        preset.answer = std::move(answer);
        _manual_actions.push_back(std::move(preset));
    }

    /**
     * @brief Update low-rank blocks forcing flag.
     *
     * @details If the flag is true then blocks that should be divided because of
     * not fitting the low-rank criteria but are too small will also be marked as low-rank.
     * Normally they are marked as dense.
     *
     * @param[in] force New flag value.
     * @return Old flag value.
     */
    bool force_lrank(bool force = true) {
        bool old = _force_lrank;
        _force_lrank = force;
        return old;
    }

    /**
     * @brief Get the undividable number of columns.
     *
     * @details If the number of block's columns is less or equal to the number then
     * it will not be sliced along that axis while hierarchical division process.
     *
     * @return Undividable column number.
     */
    uint64_t min_cols() const {
        return min_col_size;
    }

    /**
     * @brief Get the undividable number of rows.
     *
     * @details If the number of block's rows is less or equal to the number then
     * it will not be sliced along that axis while hierarchical division process.
     *
     * @return Undividable row number.
     */
    uint64_t min_rows() const {
        return min_row_size;
    }

    /**
     * @brief Set the undividable number of rows.
     *
     * @copydetails min_rows()
     *
     * @param[in] new_min New number of rows.
     */
    void min_rows(uint64_t new_min) {
        min_row_size = new_min;
    }

    /**
     * @brief Set the undividable number of columns.
     *
     * @copydetails min_cols()
     *
     * @param[in] new_min New number of columns.
     */
    void min_cols(uint64_t new_min) {
        min_col_size = new_min;
    }

    /**
     * @brief Get the division policy of the Oracle.
     *
     * @return DividePolicy.
     */
    DividePolicy divide_policy() const {
        return policy;
    }

    /**
     * @brief Set the division policy of the Oracle.
     *
     * @param[in] new_policy New policy to apply.
     */
    void divide_policy(DividePolicy new_policy) {
        policy = new_policy;
    }

    /**
     * @brief Get the relative error tolerance that is set for blocks approximation.
     *
     * @return Relative error tolerance.
     */
    double block_tolerance() const {
        return block_params.rel_tol;
    }

    /**
     * @brief Set the relative error tolerance that is set for blocks approximation.
     *
     * @param[in] new_tol New relative error tolerance value.
     */
    void block_tolerance(double new_tol) {
        block_params.rel_tol = new_tol;
    }

    bool symmetry() const {
        return _symmetry;
    }

    void symmetry(bool sym) {
        _symmetry = sym;
    }
};


extern template class ActionAnswer<float>;

extern template class ActionAnswer<double>;

extern template class ActionAnswer<std::complex<float>>;

extern template class ActionAnswer<std::complex<double>>;


extern template class Base<float>;

extern template class Base<double>;

extern template class Base<std::complex<float>>;

extern template class Base<std::complex<double>>;


}
} //namespace MSk

#endif //MOSAIC_SKELETON_ORACLE_ORACLE_HPP
