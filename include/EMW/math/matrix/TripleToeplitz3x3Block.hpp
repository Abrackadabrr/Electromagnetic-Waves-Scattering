//
// Created by evgen on 09.07.2026.
//

#ifndef TRIPLETOEPLITZ3X3BLOCK_HPP
#define TRIPLETOEPLITZ3X3BLOCK_HPP

#include "types/Types.hpp"

#include <cassert>
#include <vector>

namespace EMW::Math::LinAgl::Matrix {

template <typename T> class TripleToeplitz3x3Block {
  public:
    using scalar_type = T;
    using block_type = Eigen::Matrix<T, 3, 3>;

  private:
    struct block_storage {
        block_type value = block_type::Zero();
    };

    size_t first_layer_size_ = 0;
    size_t second_layer_size_ = 0;
    size_t third_layer_size_ = 0;

    size_t first_layer_toeplitz_size_ = 0;
    size_t second_layer_toeplitz_size_ = 0;
    size_t third_layer_toeplitz_size_ = 0;

    std::vector<block_storage> blocks_;

    [[nodiscard]] static size_t toeplitz_storage_size(size_t size) noexcept {
        assert(size > 0);
        return 2 * size - 1;
    }

    [[nodiscard]] static size_t toeplitz_index(size_t row, size_t col, size_t size) noexcept {
        assert(row < size);
        assert(col < size);
        return (col >= row) ? (col - row) : (row - col + size - 1);
    }

    [[nodiscard]] size_t linear_index(size_t third_toeplitz_index, size_t second_toeplitz_index,
                                      size_t first_toeplitz_index) const noexcept {
        assert(first_toeplitz_index < first_layer_toeplitz_size_);
        assert(second_toeplitz_index < second_layer_toeplitz_size_);
        assert(third_toeplitz_index < third_layer_toeplitz_size_);
        return first_toeplitz_index +
               first_layer_toeplitz_size_ *
                   (second_toeplitz_index + second_layer_toeplitz_size_ * third_toeplitz_index);
    }

  public:
    TripleToeplitz3x3Block() = default;

    TripleToeplitz3x3Block(size_t first_layer_size, size_t second_layer_size, size_t third_layer_size)
        : first_layer_size_(first_layer_size), second_layer_size_(second_layer_size),
          third_layer_size_(third_layer_size), first_layer_toeplitz_size_(toeplitz_storage_size(first_layer_size)),
          second_layer_toeplitz_size_(toeplitz_storage_size(second_layer_size)),
          third_layer_toeplitz_size_(toeplitz_storage_size(third_layer_size)),
          blocks_(first_layer_toeplitz_size_ * second_layer_toeplitz_size_ * third_layer_toeplitz_size_) {}

    [[nodiscard]] size_t first_layer_size() const noexcept { return first_layer_size_; }
    [[nodiscard]] size_t second_layer_size() const noexcept { return second_layer_size_; }
    [[nodiscard]] size_t third_layer_size() const noexcept { return third_layer_size_; }

    [[nodiscard]] size_t first_layer_toeplitz_size() const noexcept { return first_layer_toeplitz_size_; }
    [[nodiscard]] size_t second_layer_toeplitz_size() const noexcept { return second_layer_toeplitz_size_; }
    [[nodiscard]] size_t third_layer_toeplitz_size() const noexcept { return third_layer_toeplitz_size_; }

    [[nodiscard]] size_t rows_in_block() const noexcept { return 3; }
    [[nodiscard]] size_t cols_in_block() const noexcept { return 3; }
    [[nodiscard]] size_t rows() const noexcept { return 3 * first_layer_size_ * second_layer_size_ * third_layer_size_; }
    [[nodiscard]] size_t cols() const noexcept { return rows(); }
    [[nodiscard]] size_t stored_blocks_count() const noexcept { return blocks_.size(); }

    [[nodiscard]] block_type &get_toeplitz_block(size_t third_toeplitz_index, size_t second_toeplitz_index,
                                                 size_t first_toeplitz_index) noexcept {
        return blocks_[linear_index(third_toeplitz_index, second_toeplitz_index, first_toeplitz_index)].value;
    }

    [[nodiscard]] const block_type &get_toeplitz_block(size_t third_toeplitz_index, size_t second_toeplitz_index,
                                                       size_t first_toeplitz_index) const noexcept {
        return blocks_[linear_index(third_toeplitz_index, second_toeplitz_index, first_toeplitz_index)].value;
    }

    [[nodiscard]] block_type &get_block(size_t row3, size_t col3, size_t row2, size_t col2, size_t row1,
                                        size_t col1) noexcept {
        return get_toeplitz_block(toeplitz_index(row3, col3, third_layer_size_),
                                  toeplitz_index(row2, col2, second_layer_size_),
                                  toeplitz_index(row1, col1, first_layer_size_));
    }

    [[nodiscard]] const block_type &get_block(size_t row3, size_t col3, size_t row2, size_t col2, size_t row1,
                                              size_t col1) const noexcept {
        return get_toeplitz_block(toeplitz_index(row3, col3, third_layer_size_),
                                  toeplitz_index(row2, col2, second_layer_size_),
                                  toeplitz_index(row1, col1, first_layer_size_));
    }
};

template <typename T>
TripleToeplitz3x3Block<T> ZeroTripleToeplitz3x3Block(size_t first_layer_size, size_t second_layer_size,
                                                     size_t third_layer_size) {
    return TripleToeplitz3x3Block<T>(first_layer_size, second_layer_size, third_layer_size);
}

} // namespace EMW::Math::LinAgl::Matrix

#endif // TRIPLETOEPLITZ3X3BLOCK_HPP
