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
// Created by bulat_v on 04.02.2022.
//


#include "msk_examples_config.h"

#include <cmath>
#include <fstream>
#include <random>
#include <complex>

#include "msk/solvers/linsys/CG.hpp"
#include "msk/block/block.hpp"
#include "msk/block/dense.hpp"
#include "msk/utility/svd.hpp"


template<typename DataType>
class RegNormalPrecond {
private:
    using RealType = decltype(std::abs(DataType()));
    const uint64_t _n, _r;
    const DataType *_basis;
    const RealType *_sing;
    const RealType _lam2{0.0};

    mutable SmartArray<DataType> _coef, _outer, _tmp_coef;

    uint64_t _n_reorth{2};

public:
    RegNormalPrecond(const uint64_t n, const uint64_t r, const DataType *V, const RealType *sing, const RealType lambda = 0) :
        _n{n}, _r{r}, _basis{V}, _sing{sing}, _lam2{lambda * lambda} {
        _coef.resize(_r);
        _tmp_coef.resize(_r);
        if (_lam2 > 0) {
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


template<class DataType>
void read_vector(const std::string &fname, uint64_t n, DataType *x) {
    std::ifstream fin(fname, std::ios_base::in | std::ios_base::binary);
    fin.read(reinterpret_cast<char *>(x), sizeof(DataType) * n);
}

template<class DataType>
void write_vector(const std::string &fname, uint64_t n, const DataType *x) {
    std::ofstream fin(fname, std::ios_base::out | std::ios_base::binary);
    fin.write(reinterpret_cast<const char *>(x), sizeof(DataType) * n);
}


template<typename DataType>
class MyMatrix {
private:
    const MSk::block::Block<DataType> &data;
    uint64_t m, n;

public:
    MyMatrix(const MSk::block::Block<DataType> &block) : data{block} {
        std::tie(m, n) = block.size();
    }

    void matvec(const DataType *x, DataType *y) const {
        constexpr DataType zero{0};
        BLAS::copy(m, &zero, 0, y, 1);
        data.matvec(x, y);
    }

    void conj_matvec(const DataType *x, DataType *y) const {
        constexpr DataType zero{0};
        BLAS::copy(n, &zero, 0, y, 1);
        data.conj_matvec(x, y);
    }

    uint64_t get_row_number() const {
        return m;
    }

    uint64_t get_column_number() const {
        return n;
    }
};



int main(int argc, char *argv[]) {
    using DataType = double;
    using RealType = decltype(std::abs(DataType{}));

    const uint64_t n_rows = 1000, n_cols = 500;
    RealType noise = 1e-4;

    std::normal_distribution<double> normal(0, 1);
    std::default_random_engine eng(42);
    SmartArray<DataType> data{n_rows * n_cols};

    /*
    for (uint64_t i = 0; i < n_rows * n_cols; ++i) {
        data[i] = normal(eng);
    }
    */
    const uint64_t pre_rank = 400;
    SmartArray<DataType> Q{n_rows * pre_rank}, R{n_cols * pre_rank};
    for (uint64_t i = 0; i < n_rows * pre_rank; ++i) {
        Q[i] = normal(eng);
    }
    for (uint64_t i = 0; i < n_cols * pre_rank; ++i) {
        R[i] = normal(eng);
    }
    BLAS::gemm('N', 'T', n_rows, n_cols, pre_rank, DataType{1}, Q.data(), n_rows,
        R.data(), n_cols, DataType{0}, data.data(), n_rows);

    if (noise > 0) {
        RealType matrix_norm = BLAS::nrm2(n_rows * n_cols, data.data(), 1);
        for (uint64_t i = 0; i < n_rows * n_cols; ++i) {
            data[i] += noise * normal(eng) * matrix_norm;
        }
    }
    SmartArray<DataType> dense_matrix(data);

    MSk::index::Index<uint64_t> I(0, n_rows), J(0, n_cols);
    auto data_block = new MSk::block::Dense<DataType>{n_rows, n_cols};
    data_block->set_data(std::move(data));
    MSk::block::Block<DataType> full_block{I, J, data_block, nullptr};
    MyMatrix<DataType> matrix{full_block};

    SmartArray<DataType> sol_true{n_cols}, rhs{n_rows};
    for (uint64_t i = 0; i < n_cols; ++i) {
        sol_true[i] = normal(eng);
    }
    matrix.matvec(sol_true.data(), rhs.data());



    RealType lambda = 0;
    /*
    CG::NonePrecond<DataType> prec{n_cols};
    */

    SmartArray<DataType> U, V;
    SmartArray<RealType> S;
    auto res = MSk::perform_svd(n_rows, n_cols, dense_matrix, U, S, V, {0, 1e-8});
    uint64_t n_sing = res.memory / (n_rows + n_cols);
    std::cout << "S from " << S[0] << " to " << S[n_sing - 1] << std::endl;

    const RealType alpha = 1e-1;

    RealType tolerance{alpha * std::abs(S[0])};
    uint64_t rank{0};
    for (; rank < n_sing && std::abs(S[rank]) > tolerance; ++rank) {}
    lambda = S[rank - 1] / 2;
    std::cout << "For alpha " << alpha << ": rank = " << rank  << ", lam = " << lambda << std::endl;
    RegNormalPrecond<DataType> prec{n_cols, rank, V.data(), S.data(), lambda};

    SmartArray<DataType> Vcx{rank}, Xin{n_cols};
    BLAS::gemv('C', n_cols, rank, DataType{1}, V.data(), n_cols, sol_true.data(), 1, DataType{0}, Vcx.data(), 1);
    RealType x_norm = BLAS::nrm2(n_cols, sol_true.data(), 1);
    RealType x_in_norm = BLAS::nrm2(rank, Vcx.data(), 1) / x_norm;
    std::cout << "X norm = " << x_norm << ", in = " << x_in_norm << ", out = " << std::sqrt(1.0 - x_in_norm * x_in_norm) << std::endl;


    SmartArray<DataType> X{n_cols};
    MSk::solvers::linsys::Parameters<RealType> params{1e-6, 0.0, 10};
    params.set_verbosity(1);
    MSk::solvers::linsys::rCGNR(matrix, prec, rhs.data(), X.data(), params, lambda * lambda);

    BLAS::axpy(n_cols, DataType{-1}, sol_true.data(), 1, X.data(), 1);
    RealType err = BLAS::nrm2(n_cols, X.data(), 1) / BLAS::nrm2(n_cols, sol_true.data(), 1);
    std::cout << "err = " << err << std::endl;
    return 0;
}


