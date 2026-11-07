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


/**
 * @file
 *
 * @brief Tools for indexing and storing permutation.
 */


#ifndef MOSAIC_SKELETON_INDEX_INDEX_HPP
#define MOSAIC_SKELETON_INDEX_INDEX_HPP


#include <cstdint>
#include <vector>
#include <mutex>

#include "smart_array/smart_array.hpp"


namespace MSk {

namespace index {

/**
 * @brief A class for storing permutation of indices in range [0, N).
 * 
 * @tparam T Integer index type
 */
template<typename T = uint64_t>
class Permutation {
private:
    const uint64_t _size{0};
    SmartArray<T> _perm;

    mutable bool _identical{true};

    // switchable synchronization
    bool _locked{true};
    mutable std::mutex _m;
    //std::vector<bool> _in_use;

private:
    template<typename E>
    bool _permute(const E *input, E *output) const {
        for (uint64_t i = 0; i < _size; ++i) {
            output[i] = input[_perm[i]];
        }
        return true;
    }

    template<typename E>
    bool _unpermute(const E *input, E *output) const {
        for (uint64_t i = 0; i < _size; ++i) {
            output[_perm[i]] = input[i];
        }
        return true;
    }

public:
    /**
     * @brief Construct a new Permutation object with identity permutation.
     * 
     * @param[in] size Sequence length.
     */
    explicit Permutation(const uint64_t size) : _size{size}, _perm(size)/*, _in_use(size, false)*/ {
        T idx{0};
        for (uint64_t i = 0; i < size; ++i, ++idx) {
            _perm[i] = idx;
        }
    }

    ~Permutation() = default;

    Permutation(const Permutation &) = delete;

    Permutation(Permutation &&other) noexcept = delete;

    Permutation &operator=(const Permutation &) = delete;

    Permutation &operator=(Permutation &&other) noexcept = delete;

    /**
     * @brief Size getter.
     * 
     * @return Size of the permutation array.
     */
    uint64_t size() const {
        return _size;
    }

    /**
     * @brief Enable/disable lock flag.
     *
     * @details Enable/disable locking the permutaion array while accesing.
     * 
     * @param[in] to New lock flag value.
     * @return Old lock flag value.
     */
    bool lock(bool to) {
        std::lock_guard<std::mutex> lock{_m};
        bool was = _locked;
        _locked = to;
        return was;
    }

    /**
     * @brief Get if permutation is identical a priori.
     *
     * @param[in] check If to scan the inner array.
     * @return Identity flag value.
     */
    bool identical(bool check = false) {
        std::unique_lock<std::mutex> lock(_m, std::defer_lock);
        if (_locked || check) {
            lock.lock();
        }
        if (check) {
            _identical = true;
            for (uint64_t i = 0; i < _size; ++i) {
                if (_perm[i] != i) {
                    _identical = false;
                    break;
                }
            }
        }
        return _identical;
    }

    /**
    * @brief Get number at the given index.
    * 
    * @param[in] idx Position of the number.
    * @return Number at position idx.
    */
    T operator[](uint64_t idx) const {
        std::unique_lock<std::mutex> lock(_m, std::defer_lock);
        if (_locked) {
            lock.lock();
        }
        if (_identical) {
            return idx;
        }
        return _perm[idx];
    }

    /**
     * @brief Transpose two numbers in the array.
     * 
     * @param[in] i First number position.
     * @param[in] j Second number position.
     * @return Trasposition success.
     */
    bool transpose(uint64_t i, uint64_t j) {
        if (i >= _size || j >= _size) {
            return false;
        }
        std::lock_guard<std::mutex> lock{_m};
        _identical = false;
        T tmp{_perm[i]};
        _perm[i] = _perm[j];
        _perm[j] = tmp;
        return true;
    }

    /**
     * @brief Update current permutation by applying another one.
     * 
     * @param[in] p Permutation to apply.
     * @param[in] i0 Index to start applying from.
     * @return Update success.
     */
    bool update(const Permutation &p, uint64_t i0 = 0) {
        uint64_t n_perm{p.size()};
        if (i0 + n_perm > _size) {
            return false;
        }
        std::lock_guard<std::mutex> lock{_m};
        _identical = false;
        return p.permute(&_perm[i0]);
    }

    /**
     * @brief Permute the given vector of arbitrary elements and write the result to another.
     * 
     * @tparam E Vector element type.
     * @param[in] elements Elements to permute.
     * @param[out] result Resulting vector.
     * @return Permutation success.
     *
     * @note Resulting vector is resisized if required.
     */
    template<typename E>
    bool permute(const std::vector<E> &elements, std::vector<E> &result) const {
        if (elements.size() != _size) {
            return false;
        }
        std::unique_lock<std::mutex> lock(_m, std::defer_lock);
        if (_locked) {
            lock.lock();
        }
        if (_identical) {
            result = elements;
            return true;
        }
        result.resize(_size);
        return _permute(elements.data(), result.data());
    }

    /**
     * @brief Unpermute (inverse permute) the given vector of arbitrary elements and write the result to another.
     * 
     * @tparam E Vector element type.
     * @param[in] elements Elements to unpermute.
     * @param[out] result Resulting vector.
     * @return Unpermutation success.
     *
     * @note Resulting vector is resisized if required.
     */
    template<typename E>
    bool unpermute(const std::vector<E> &elements, std::vector<E> &result) const {
        if (elements.size() != _size) {
            return false;
        }
        std::unique_lock<std::mutex> lock(_m, std::defer_lock);
        if (_locked) {
            lock.lock();
        }
        if (_identical) {
            result = elements;
            return true;
        }
        result.resize(_size);
        return _unpermute(elements.data(), result.data());
    }

    /**
     * @brief Permute the given vector of arbitrary elements inplace.
     * 
     * @tparam E Vector element type.
     * @param[in, out] elements Elements to permute.
     * @return Permutation success.
     */
    template<typename E>
    bool permute(std::vector<E> &elements) const {
        if (elements.size() != _size) {
            return false;
        }
        std::unique_lock<std::mutex> lock(_m, std::defer_lock);
        if (_locked) {
            lock.lock();
        }
        if (_identical) { // Fast exit with no permutation
            return true;
        }
        std::vector<E> old_elements(std::move(elements));
        elements.resize(_size);
        return _permute(old_elements.data(), elements.data());
    }

    /**
     * @brief Unpermute (inverse permute) the given vector of arbitrary elements inplace.
     * 
     * @tparam E Vector element type.
     * @param[in, out] elements Elements to unpermute.
     * @return Unpermutation success.
     */
    template<typename E>
    bool unpermute(std::vector<E> &elements) const {
        if (elements.size() != _size) {
            return false;
        }
        std::unique_lock<std::mutex> lock(_m, std::defer_lock);
        if (_locked) {
            lock.lock();
        }
        if (_identical) { // Fast exit with no unpermutation
            return true;
        }
        std::vector<E> old_elements(std::move(elements));
        elements.resize(_size);
        return _unpermute(old_elements.data(), elements.data());
    }


    /**
     * @brief Permute the given array of arbitrary elements inplace.
     * 
     * @tparam E Array element type.
     * @param[in, out] elements Elements to permute.
     * @return Permutation success.
     */
    template<typename E>
    bool permute(E *elements) const {
        std::unique_lock<std::mutex> lock(_m, std::defer_lock);
        if (_locked) {
            lock.lock();
        }
        if (_identical) { // Fast exit with no permutation
            return true;
        }
        SmartArray<E> old_elements(_size);
        std::copy(elements, elements + _size, old_elements.data());
        return _permute(old_elements.data(), elements);
    }

    /**
     * @brief Unpermute (inverse permute) the given array of arbitrary elements inplace.
     * 
     * @tparam E Array element type.
     * @param[in, out] elements Elements to unpermute.
     * @return Unpermutation success.
     */
    template<typename E>
    bool unpermute(E *elements) const {
        std::unique_lock<std::mutex> lock(_m, std::defer_lock);
        if (_locked) {
            lock.lock();
        }
        if (_identical) { // Fast exit with no unpermutation
            return true;
        }
        SmartArray<E> old_elements(_size);
        std::copy(elements, elements + _size, old_elements.data());
        return _unpermute(old_elements.data(), elements);
    }

    /**
     * @brief Permute the given array of arbitrary elements and write the result to another.
     * 
     * @tparam E Array element type.
     * @param[in] elements Elements to permute.
     * @param[out] result Resulting array.
     * @return Permutation success.
     */
    template<typename E>
    bool permute(const E *elements, E *result) const {
        if (result == nullptr) {
            return false;
        }
        if (elements == result) {
            return permute(result);
        }
        std::unique_lock<std::mutex> lock(_m, std::defer_lock);
        if (_locked) {
            lock.lock();
        }
        if (_identical) {
            for (uint64_t i = 0; i < _size; ++i) {
                result[i] = elements[i];
            }
            return true;
        }
        return _permute(elements, result);
    }

    /**
     * @brief Unpermute (inverse permute) the given array of arbitrary elements and write the result to another.
     * 
     * @tparam E Array element type.
     * @param[in] elements Elements to unpermute.
     * @param[out] result Resulting array.
     * @return Unpermutation success.
     */
    template<typename E>
    bool unpermute(const E *elements, E *result) const {
        if (result == nullptr) {
            return false;
        }
        if (elements == result) {
            return unpermute(result);
        }
        std::unique_lock<std::mutex> lock(_m, std::defer_lock);
        if (_locked) {
            lock.lock();
        }
        if (_identical) {
            for (uint64_t i = 0; i < _size; ++i) {
                result[i] = elements[i];
            }
            return true;
        }
        return _unpermute(elements, result);
    }

    /**
     * @brief Permute columns of a given matrix and write the result to another.
     * 
     * @tparam E Matrix element type.
     * @param[in] col_size Block column size.
     * @param[in] cols Pointer to the corner of a block
     * @param[in] cols_ldim Leading dimension of a block.
     * @param[out] result Resulting matrix.
     * @param[in] result_ldim Leading dimension of the result.
     * @return Permutation success,
     */
    template<typename E>
    bool permute_columns(const uint64_t col_size, const E *cols, const uint64_t cols_ldim, 
                         E *result, const uint64_t result_ldim) const {
        std::unique_lock<std::mutex> lock{_m, std::defer_lock};
        if (_locked) {
            lock.lock();
        }
        if (_identical) {
            for (uint64_t j = 0; j < _size; ++j) {
                std::copy(cols + cols_ldim * j, cols + cols_ldim * j + col_size,
                          result + result_ldim * j);
            }
        } else {
            for (uint64_t j = 0; j < _size; ++j) {
                std::copy(cols + cols_ldim * _perm[j], cols + cols_ldim * _perm[j] + col_size,
                          result + result_ldim * j);
            }
        }
        return true;
    }

    /**
     * @brief Permute columns of a given matrix inplace.
     * 
     * @tparam E Matrix element type.
     * @param[in] col_size Block column size.
     * @param[in, out] cols Pointer to the corner of a block
     * @param[in] cols_ldim Leading dimension of a block.
     * @return Permutation success
     */
    template<typename E>
    bool permute_columns(const uint64_t col_size, E *cols, const uint64_t cols_ldim) const {
        std::unique_lock<std::mutex> lock{_m, std::defer_lock};
        if (_locked) {
            lock.lock();
        }
        if (_identical) { // Fast exit with no permutation
            return true;
        }
        SmartArray<E> old_cols{col_size * _size};
        for (uint64_t j = 0; j < _size; ++j) {
            std::copy(cols + cols_ldim * j, cols + cols_ldim * j + col_size,
                        old_cols.data() + col_size * j);
        }
        for (uint64_t j = 0; j < _size; ++j) {
            std::copy(old_cols.data() + col_size * _perm[j], old_cols.data() + col_size * (_perm[j] + 1),
                      cols + cols_ldim * j);
        }
        return true;
    }

    /**
     * @brief Unpermute (inverse permute) columns of a given matrix and write the result to another.
     * 
     * @tparam E Matrix element type.
     * @param[in] col_size Block column size.
     * @param[in] cols Pointer to the corner of a block
     * @param[in] cols_ldim Leading dimension of a block.
     * @param[out] result Resulting matrix.
     * @param[in] result_ldim Leading dimension of the result.
     * @return Unpermutation success,
     */
    template<typename E>
    bool unpermute_columns(const uint64_t col_size, const E *cols, const uint64_t cols_ldim, 
                           E *result, const uint64_t result_ldim) const {
        std::unique_lock<std::mutex> lock{_m, std::defer_lock};
        if (_locked) {
            lock.lock();
        }
        if (_identical) { // Fast exit with no unpermutation
            for (uint64_t j = 0; j < _size; ++j) {
                std::copy(cols + cols_ldim * j, cols + cols_ldim * j + col_size,
                          result + result_ldim * j);
            }
        } else {
            for (uint64_t j = 0; j < _size; ++j) {
                std::copy(cols + cols_ldim * j, cols + cols_ldim * j + col_size,
                          result + result_ldim * _perm[j]);
            }
        }
        return true;
    }

    /**
     * @brief Unpermute (inverse permute) columns of a given matrix inplace.
     * 
     * @tparam E Matrix element type.
     * @param[in] col_size Block column size.
     * @param[in, out] cols Pointer to the corner of a block
     * @param[in] cols_ldim Leading dimension of a block.
     * @return Unpermutation success
     */
    template<typename E>
    bool unpermute_columns(uint64_t col_size, E *cols, uint64_t cols_ldim) const {
        std::unique_lock<std::mutex> lock{_m, std::defer_lock};
        if (_locked) {
            lock.lock();
        }
        if (_identical) { // Fast exit with no permutation
            return true;
        }
        SmartArray<E> old_cols{col_size * _size};
        for (uint64_t j = 0; j < _size; ++j) {
            std::copy(cols + cols_ldim * j, cols + cols_ldim * j + col_size, 
                      old_cols.data() + col_size * j);
        }
        for (uint64_t j = 0; j < _size; ++j) {
            std::copy(old_cols.data() + col_size * j, old_cols.data() + col_size * (j + 1), 
                      cols + cols_ldim * _perm[j]);
        }
        return true;
    }

};


template<typename T>
class IndexIterator;


/**
 * @brief Indexing over an integer segment with optional permutation.
 *
 * @tparam T Integer index type.
 */
template<typename T = uint64_t>
class Index {
public:
    using idx_type = T;
    using iterator_type = IndexIterator<T>;

private:
    T _begin{}, _end{};
    uint64_t _size{};
    Permutation<T> *_perm;

public:
    /**
     * @brief Construct Index for the integer range [0, size).
     *
     * @param[in] size Segment size.
     * @param[in] perm Optional Permutation.
     */
    Index(uint64_t size = 0, Permutation<T> *perm = nullptr) : _begin{0}, _end{size}, _size{size}, _perm{perm} {
        if (_perm != nullptr && _size > _perm->size()) {
            throw std::logic_error("Index out of permutation range");
        }
    }

    /**
     * @brief Construct Index for the integer range [begin, end).
     *
     * @param[in] idx_begin Segment beginning (included).
     * @param[in] idx_end Segment ending (excluded).
     * @param[in] perm Possible Permutation.
     */
    Index(T idx_begin, T idx_end, Permutation<T> *perm = nullptr) :
            _begin{idx_begin}, _end{idx_end}, _size{_end - _begin}, _perm{perm} {
        if (_begin > _end) {
            throw std::logic_error("Invalid Index bounds");
        }
        if (_perm != nullptr && _size > _perm->size()) {
            throw std::logic_error("Index out of permutation range");
        }
    }

    /**
     * @brief Construct Index from beginning and ending iterators.
     *
     * @param[in] it_begin Segment beginning (included).
     * @param[in] it_end Segment ending (excluded).
     */
    Index(const IndexIterator<T> &it_begin, const IndexIterator<T> &it_end) :
            _begin{it_begin._idx}, _end{it_end._idx}, _size{_end - _begin}, _perm{it_begin._perm} {
        if (_begin > _end) {
            throw std::logic_error("Invalid Index bounds");
        }
        if (it_begin._perm != it_end._perm) {
            throw std::logic_error("Different Permutations are passed to Index");
        }
    }

    Index(const Index &) = default;

    Index(Index &&) = default;

    Index &operator=(const Index &) = default;

    Index &operator=(Index &&) = default;

    ~Index() = default;

    uint64_t size() const {
        return _size;
    }

    bool empty() const {
        return _size == 0;
    }

    T operator[](uint64_t idx) const {
        if (idx >= _size) {
            throw std::runtime_error("Index index out of range");
        }
        if (_perm != nullptr) {
            return (*_perm)[_begin + idx];
        } else {
            return _begin + idx;
        }
    }

    void permute(Permutation<T> *perm) {
        _perm = perm;
    }

    Permutation<T> *permutation() {
        return _perm;
    }

    Index intersect(const Index &other) const {
        T inter_begin = _begin > other._begin ? _begin : other._begin;
        T inter_end = _end < other._end ? _end : other._end;
        if (inter_begin >= inter_end) {
            return {0, _perm};
        }
        return {inter_begin, inter_end, _perm};
    }

    std::pair<Index, Index> split(IndexIterator<T> mid) const {
        mid.permute(false);
        T _mid = *mid;
        _mid = std::max(_begin, std::min(_end, _mid));
        return {Index(_begin, _mid, _perm), Index(_mid, _end, _perm)};
    }

    friend bool operator==(const Index &r1, const Index &r2) {
        bool range{r1._begin == r2._begin && r1._end == r2._end};
        bool perm{r1._perm == nullptr || r2._perm == nullptr || r1._perm == r2._perm};
        return range && perm;
    }

    friend bool operator!=(const Index &r1, const Index &r2) {
        bool no_range{r1._begin != r2._begin || r1._end != r2._end};
        bool no_perm{r1._perm != nullptr && r2._perm != nullptr && r1._perm != r2._perm};
        return no_range || no_perm;
    }

    // Iterators
    /**
     * @brief General (begin) IndexIterator with optional permutation usage.
     *
     * @param[in] use_perm Flag to use permatation.
     * @return Beginning IndexIterator.
     */
    IndexIterator<T> begin(bool use_perm = true) const {
        return {_begin, _perm, use_perm};
    }

    /**
     * @brief General (end) IndexIterator with optional permutation usage.
     *
     * @param[in] use_perm Flag to use permatation.
     * @return Past-the-end IndexIterator.
     */
    IndexIterator<T> end(bool use_perm = true) const {
        return {_end, _perm, use_perm};
    }

    /**
     * @brief (begin) IndexIterator that uses permutation.
     *
     * @return Beginning IndexIterator.
     */
    IndexIterator<T> pbegin() const {
        return {_begin, _perm, true};
    }

    /**
     * @brief (end) IndexIterator that uses permutation.
     *
     * @return Past-the-end IndexIterator.
     */
    IndexIterator<T> pend() const {
        return {_end, _perm, true};
    }

    /**
     * @brief (begin) IndexIterator without permutation.
     *
     * @return Beginning IndexIterator.
     */
    IndexIterator<T> dbegin() const {
        return {_begin, _perm, false};
    }

    /**
     * @brief (end) IndexIterator without permutation.
     *
     * @return Past-the-end IndexIterator.
     */
    IndexIterator<T> dend() const {
        return {_end, _perm, false};
    }
};


/**
 * @brief Iterator to run over Index.
 *
 * @tparam T Integer index type.
 */
template<typename T = uint64_t>
class IndexIterator {
    friend class Index<T>;

private:
    uint64_t _idx;
    Permutation<T> *_perm;
    bool _use_perm;

private:
    /**
     * @brief Raw IndexIterator constructor.
     *
     * @param[in] idx Starting index.
     * @param[in] perm Permutation if available.
     * @param[in] use_perm Flag to use permutation.
     */
    IndexIterator(uint64_t idx, Permutation<T> *perm = nullptr, bool use_perm = false) :
            _idx{idx}, _perm{perm}, _use_perm{use_perm && perm != nullptr} {}

public:
    IndexIterator() : _idx{0}, _perm{nullptr}, _use_perm{false} {}

    IndexIterator(const IndexIterator &other) = default;

    /**
     * @brief See if permutation is used.
     *
     * @return bool flag value.
     */
    bool is_perm() const {
        return _use_perm;
    }

    /**
     * @brief Enable/disable permutation usage.
     *
     * @param[in] use_perm Flag to switch to.
     * @return Resulting flag value.
     */
    bool permute(bool use_perm = true) {
        _use_perm = use_perm && _perm != nullptr;
        return use_perm;
    }

    T operator*() const {
        if (_use_perm) {
            return (*_perm)[_idx];
        } else {
            return _idx;
        }
    }

    IndexIterator &operator++() {
        ++_idx;
        return *this;
    }

    IndexIterator &operator--() {
        --_idx;
        return *this;
    }

    IndexIterator &operator+=(int shift) {
        _idx += shift;
        return *this;
    }

    IndexIterator &operator-=(int shift) {
        _idx -= shift;
        return *this;
    }

    friend bool operator==(const IndexIterator &it1, const IndexIterator &it2) {
        return it1._idx == it2._idx;
    }

    friend bool operator!=(const IndexIterator &it1, const IndexIterator &it2) {
        return it1._idx != it2._idx;
    }

    friend bool operator<=(const IndexIterator &it1, const IndexIterator &it2) {
        return it1._idx <= it2._idx;
    }

    friend bool operator>=(const IndexIterator &it1, const IndexIterator &it2) {
        return it1._idx >= it2._idx;
    }

    friend bool operator<(const IndexIterator &it1, const IndexIterator &it2) {
        return it1._idx < it2._idx;
    }

    friend bool operator>(const IndexIterator &it1, const IndexIterator &it2) {
        return it1._idx > it2._idx;
    }

    friend IndexIterator operator+(const IndexIterator &it, int shift) {
        return {it._idx + shift, it._perm, it._use_perm};
    }

    friend IndexIterator operator-(const IndexIterator &it, int shift) {
        return {it._idx - shift, it._perm, it._use_perm};
    }

    friend int operator-(IndexIterator it1, IndexIterator it2) {
        return it1._idx - it2._idx;
    }

    /**
     * @brief Permute elements between two iterators.
     *
     * @param[in] it_begin Sequence beginning.
     * @param[in] it_end Sequence ending.
     * @param[in] perm Permutation to apply.
     * @return Permutation success.
     */
    friend bool permute(IndexIterator it_begin, IndexIterator it_end, const Permutation<T> &perm) {
        if (it_begin._perm != it_end._perm || it_begin._perm == nullptr) {
            return false;
        }
        if (perm.size() != it_end._idx - it_begin._idx) {
            return false;
        }
        return it_begin._perm->update(perm, it_begin._idx);
    }

    /**
     * @brief Transpose (swap) two elements.
     *
     * @param[in] it1 First element.
     * @param[in] it2 Second element.
     * @return Transposition success.
     */
    friend bool transpose(IndexIterator it1, IndexIterator it2) {
        if (it1._perm != it2._perm) {
            return false;
        }
        std::swap(it1._idx, it2._idx);
        if (it1._perm == nullptr) {
            return true;
        }
        return it1._perm->transpose(it2._idx, it1._idx);
    }

};


extern template class Permutation<>;

extern template class Index<>;

extern template class IndexIterator<>;


}
}


#endif //MOSAIC_SKELETON_INDEX_INDEX_HPP
