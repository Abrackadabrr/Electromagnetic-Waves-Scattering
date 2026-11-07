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
// Created by b.valiakhmetov on 10/28/24.
//

#ifndef MOSAIC_SKELETON_SOLVERS_LINSYS_PRECONDITION_HPP
#define MOSAIC_SKELETON_SOLVERS_LINSYS_PRECONDITION_HPP


#include <cstdint>

#include "smart_array/smart_array.hpp"
#include "cpp_blas/BLAS.hpp"


namespace MSk {
namespace solvers {
namespace linsys {
namespace precond {

/**
 * @brief Empty (identical) preconditioner for system solving iterative methods.
 *
 * @details Can be used as a plug for solving withous preconditioning.
 *
 * @tparam DataType Vector element type.
 */
template<class DataType>
class None {
private:
    uint64_t _size; ///< Matrix size
public:
    /**
     * @brief Initialize.
     *
     * @param[in] size Matrix size.
     */
    None(uint64_t size) : _size{size} {}

    /**
     * @brief Perform preconditioning: simple copy.
     *
     * @param[in] x Input vector.
     * @param[out] y Output vector.
     */
    void apply(const DataType *x, DataType *y) const {
        BLAS::copy(_size, x, 1, y, 1);
    }
};


/**
 * @brief Diagonal preconditioner.
 *
 * @tparam DataType Element type.
 */
template<class DataType>
class Diag {
private:
    uint64_t _size; ///< Matrix size
    SmartArray<DataType> _inv_diag; ///< Inverse diagonal elements

public:
    /**
     * @brief Initialize diagonal elements by copy.
     *
     * @param[in] size Matrix size.
     * @param[in] diag Diagonal elements.
     * @param[in] inverse Flag: if elements are already inversed.
     */
    Diag(uint64_t size, const SmartArray<DataType> &diag, bool inverse = true) : _size{size}, _inv_diag{diag} {
        if (!inverse) {
            for (uint64_t i = 0; i < _inv_diag.size(); ++i) {
                _inv_diag[i] = DataType{1.0} / _inv_diag[i];
            }
        }
    }

    /**
     * @brief Initialize diagonal elements by value.
     *
     * @param[in] size Matrix size.
     * @param[in, out] diag Diagonal elements.
     * @param[in] inverse Flag: if elements are already inversed.
     */
    Diag(uint64_t size, SmartArray<DataType> &&diag, bool inverse = true) : _size{size}, _inv_diag{std::move(diag)} {
        if (!inverse) {
            for (uint64_t i = 0; i < _inv_diag.size(); ++i) {
                _inv_diag[i] = DataType{1.0} / _inv_diag[i];
            }
        }
    }

    /**
     * @brief Perform diagonal preconditioning.
     *
     * @param[in] x Input vector.
     * @param[out] y Output vector.
     */
    void apply(const DataType *x, DataType *y) const {
        for (uint64_t i = 0; i < _size; ++i) {
            y[i] = _inv_diag[i] * x[i];
        }
    }
};



template<typename DataType>
class RegNormal {
private:
    using RealType = decltype(std::abs(DataType()));
    const uint64_t _n, _r;
    const DataType *_basis;
    const RealType *_sing;
    const RealType _lam2{0.0};

    mutable SmartArray<DataType> _coef, _outer, _tmp_coef;

    uint64_t _n_reorth{2};

public:
    RegNormal(const uint64_t n, const uint64_t r, const DataType *V, const RealType *sing, const RealType lambda = 0) :
        _n{n}, _r{r}, _basis{V}, _sing{sing}, _lam2{lambda * lambda} {
        _coef.resize(_r);
        if (_lam2 > 0) {
            _tmp_coef.resize(_r);
            _outer.resize(_n);
        }
    }

    void apply(const DataType *x, DataType *y) const {
        BLAS::gemv('C', _n, _r, DataType{1}, _basis, _n, x, 1, DataType{0}, _coef.data(), 1);
        if (_lam2 > 0) {
            BLAS::copy(_r, _coef.data(), 1, _tmp_coef.data(), 1);
        }
        for (uint64_t i = 0; i < _r; ++i) {
            _coef[i] /= _sing[i] * _sing[i] + _lam2;
        }
        BLAS::gemv('N', _n, _r, DataType{1}, _basis, _n, _coef.data(), 1, DataType{0}, y, 1);
        if (_lam2 > 0) {
            BLAS::copy(_n, x, 1, _outer.data(), 1);
            BLAS::gemv('N', _n, _r, DataType{-1}, _basis, _n, _tmp_coef.data(), 1, DataType{1}, _outer.data(), 1);
            for (uint64_t i = 0; i < _n_reorth; ++i) {
                BLAS::gemv('C', _n, _r, DataType{1}, _basis, _n, _outer.data(), 1, DataType{0}, _tmp_coef.data(), 1);
                BLAS::gemv('N', _n, _r, DataType{-1}, _basis, _n, _tmp_coef.data(), 1, DataType{1}, _outer.data(), 1);
            }
            BLAS::axpy(_n, DataType{1} / _lam2, _outer.data(), 1, y, 1);
        }
    }
};

}
}
}
}















#endif //MOSAIC_SKELETON_SOLVERS_LINSYS_PRECONDITION_HPP
