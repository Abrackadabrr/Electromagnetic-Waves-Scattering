#ifndef TRIPLE_TOEPLITZ_3X3_FOURIER_HPP
#define TRIPLE_TOEPLITZ_3X3_FOURIER_HPP

#include "math/matrix/Matrix.hpp"
#include "third_party/eigen/unsupported/Eigen/FFT"
#include "types/Types.hpp"

#include <array>
#include <cassert>
#include <stdexcept>
#include <type_traits>

namespace EMW::Math::Fourier {

/**
 * Compact storage for the unique coefficients of a 3-level Toeplitz matrix
 * with dense 3x3 internal blocks: levels_x * levels_y * levels_z * 3 * 3.
 */
template <typename scalar_t>
class TripleToeplitz3x3Tensor {
private:
    static constexpr Types::index block_size_ = 3;

    Types::index levels_x_ = 0;
    Types::index levels_y_ = 0;
    Types::index levels_z_ = 0;
    Containers::vector<scalar_t> data_;

    [[nodiscard]] Types::index flat_index(Types::index lx, Types::index ly, Types::index lz,
                                          Types::index row, Types::index col) const noexcept {
        return (((lz * levels_y_ + ly) * levels_x_ + lx) * block_size_ + row) * block_size_ + col;
    }

public:
    TripleToeplitz3x3Tensor() = default;

    TripleToeplitz3x3Tensor(Types::index levels_x, Types::index levels_y, Types::index levels_z,
                            scalar_t init = scalar_t{})
        : levels_x_(levels_x), levels_y_(levels_y), levels_z_(levels_z),
          data_(levels_x * levels_y * levels_z * block_size_ * block_size_, init) {
    }

    [[nodiscard]] Types::index levels_x() const noexcept { return levels_x_; }
    [[nodiscard]] Types::index levels_y() const noexcept { return levels_y_; }
    [[nodiscard]] Types::index levels_z() const noexcept { return levels_z_; }

    [[nodiscard]] const Containers::vector<scalar_t>& data() const noexcept { return data_; }

    [[nodiscard]] scalar_t& operator()(Types::index lx, Types::index ly, Types::index lz,
                                       Types::index row, Types::index col) noexcept {
        assert(lx < levels_x_);
        assert(ly < levels_y_);
        assert(lz < levels_z_);
        assert(row < block_size_);
        assert(col < block_size_);
        return data_[flat_index(lx, ly, lz, row, col)];
    }

    [[nodiscard]] const scalar_t& operator()(Types::index lx, Types::index ly, Types::index lz,
                                             Types::index row, Types::index col) const noexcept {
        assert(lx < levels_x_);
        assert(ly < levels_y_);
        assert(lz < levels_z_);
        assert(row < block_size_);
        assert(col < block_size_);
        return data_[flat_index(lx, ly, lz, row, col)];
    }
};

/**
 * Matrix-vector product for a triple Toeplitz matrix with 3x3 inner blocks.
 * The first 3 dimensions (Toeplitz levels) are multiplied via 3D FFT.
 */
template <typename scalar_t>
class TripleToeplitz3x3Fourier {
public:
    using scalar_type = scalar_t;
    using complex_type = Types::complex_d;
    using vector_type = Types::VectorX<scalar_t>;
    using tensor_type = TripleToeplitz3x3Tensor<scalar_t>;
    using triple_toeplitz_block_type = LinAgl::Matrix::TripleToeplitzBlock<scalar_t>;

    TripleToeplitz3x3Fourier() = default;

private:
    static constexpr Types::index block_size_ = 3;

    tensor_type levels_;

    Types::index nx_ = 0;
    Types::index ny_ = 0;
    Types::index nz_ = 0;

    Types::index shift_x_ = 0;
    Types::index shift_y_ = 0;
    Types::index shift_z_ = 0;

    Types::index fft_x_ = 0;
    Types::index fft_y_ = 0;
    Types::index fft_z_ = 0;

    std::array<Containers::vector<complex_type>, block_size_ * block_size_> kernel_spectrum_;

    // Cell numbering follows mesh convention:
    // x changes fastest, then y, then z:
    // idx = x + size_x * (y + size_y * z)
    [[nodiscard]] static Types::index flatten_xyz_x_major(Types::index x, Types::index y, Types::index z,
                                                          Types::index size_x, Types::index size_y) noexcept {
        return x + size_x * (y + size_y * z);
    }

    [[nodiscard]] Types::index flatten_fft(Types::index x, Types::index y, Types::index z) const noexcept {
        return flatten_xyz_x_major(x, y, z, fft_x_, fft_y_);
    }

    [[nodiscard]] Types::index flatten_cell(Types::index x, Types::index y, Types::index z) const noexcept {
        return flatten_xyz_x_major(x, y, z, nx_, ny_);
    }

    [[nodiscard]] Types::index flatten_cell_component(Types::index x, Types::index y, Types::index z,
                                                      Types::index component) const noexcept {
        return flatten_cell(x, y, z) * block_size_ + component;
    }

    [[nodiscard]] static complex_type to_complex(const scalar_t& value) noexcept {
        if constexpr (std::is_same_v<scalar_t, complex_type>) {
            return value;
        } else {
            return complex_type{static_cast<Types::scalar>(value), 0.0};
        }
    }

    [[nodiscard]] static scalar_t from_complex(const complex_type& value) noexcept {
        if constexpr (std::is_same_v<scalar_t, complex_type>) {
            return value;
        } else {
            return static_cast<scalar_t>(value.real());
        }
    }

    [[nodiscard]] Types::index fft_size() const noexcept {
        return fft_x_ * fft_y_ * fft_z_;
    }

    [[nodiscard]] static tensor_type tensor_from_triple_toeplitz_block(const triple_toeplitz_block_type& matrix) {
        const Types::index nz = matrix.rows_in_toeplitrz();
        if (nz == 0 || nz != matrix.cols_in_toeplitrz()) {
            throw std::invalid_argument("TripleToeplitzBlock must have square non-empty third Toeplitz level");
        }

        const auto& second_level_sample = matrix.get_block(0, 0);
        const Types::index ny = second_level_sample.rows_in_toeplitrz();
        if (ny == 0 || ny != second_level_sample.cols_in_toeplitrz()) {
            throw std::invalid_argument("TripleToeplitzBlock must have square non-empty second Toeplitz level");
        }

        const auto& first_level_sample = second_level_sample.get_block(0, 0);
        const Types::index nx = first_level_sample.rows_in_toeplitrz();
        if (nx == 0 || nx != first_level_sample.cols_in_toeplitrz()) {
            throw std::invalid_argument("TripleToeplitzBlock must have square non-empty first Toeplitz level");
        }

        const auto& internal_block_sample = first_level_sample.get_block(0, 0);
        if (internal_block_sample.rows() != block_size_ || internal_block_sample.cols() != block_size_) {
            throw std::invalid_argument("TripleToeplitzBlock internal blocks must be 3x3");
        }

        tensor_type levels(2 * nx - 1, 2 * ny - 1, 2 * nz - 1, scalar_t{});

        const auto row_col_from_level = [](Types::index level, Types::index base_size,
                                           Types::index& row, Types::index& col) {
            const Types::integer signed_shift =
                static_cast<Types::integer>(level) - static_cast<Types::integer>(base_size - 1);
            if (signed_shift >= 0) {
                row = static_cast<Types::index>(signed_shift);
                col = 0;
            } else {
                row = 0;
                col = static_cast<Types::index>(-signed_shift);
            }
        };

        for (Types::index lz = 0; lz < levels.levels_z(); ++lz) {
            Types::index i3 = 0;
            Types::index j3 = 0;
            row_col_from_level(lz, nz, i3, j3);
            const auto& second_level = matrix.get_block(i3, j3);

            if (second_level.rows_in_toeplitrz() != ny || second_level.cols_in_toeplitrz() != ny) {
                throw std::invalid_argument("TripleToeplitzBlock has inconsistent second-level Toeplitz sizes");
            }

            for (Types::index ly = 0; ly < levels.levels_y(); ++ly) {
                Types::index i2 = 0;
                Types::index j2 = 0;
                row_col_from_level(ly, ny, i2, j2);
                const auto& first_level = second_level.get_block(i2, j2);

                if (first_level.rows_in_toeplitrz() != nx || first_level.cols_in_toeplitrz() != nx) {
                    throw std::invalid_argument("TripleToeplitzBlock has inconsistent first-level Toeplitz sizes");
                }

                for (Types::index lx = 0; lx < levels.levels_x(); ++lx) {
                    Types::index i1 = 0;
                    Types::index j1 = 0;
                    row_col_from_level(lx, nx, i1, j1);
                    const auto& internal_block = first_level.get_block(i1, j1);

                    if (internal_block.rows() != block_size_ || internal_block.cols() != block_size_) {
                        throw std::invalid_argument("TripleToeplitzBlock has non-3x3 internal block");
                    }

                    for (Types::index out_comp = 0; out_comp < block_size_; ++out_comp) {
                        for (Types::index in_comp = 0; in_comp < block_size_; ++in_comp) {
                            levels(lx, ly, lz, out_comp, in_comp) = internal_block(out_comp, in_comp);
                        }
                    }
                }
            }
        }

        return levels;
    }

    void fft3_inplace(Containers::vector<complex_type>& data, bool inverse) const {
        Eigen::FFT<Types::scalar> fft;

        Containers::vector<complex_type> line_in;
        Containers::vector<complex_type> line_out;

        line_in.resize(fft_x_);
        line_out.resize(fft_x_);
        for (Types::index z = 0; z < fft_z_; ++z) {
            for (Types::index y = 0; y < fft_y_; ++y) {
                for (Types::index x = 0; x < fft_x_; ++x) {
                    line_in[x] = data[flatten_fft(x, y, z)];
                }

                if (inverse) {
                    fft.inv(line_out, line_in);
                } else {
                    fft.fwd(line_out, line_in);
                }

                for (Types::index x = 0; x < fft_x_; ++x) {
                    data[flatten_fft(x, y, z)] = line_out[x];
                }
            }
        }

        line_in.resize(fft_y_);
        line_out.resize(fft_y_);
        for (Types::index z = 0; z < fft_z_; ++z) {
            for (Types::index x = 0; x < fft_x_; ++x) {
                for (Types::index y = 0; y < fft_y_; ++y) {
                    line_in[y] = data[flatten_fft(x, y, z)];
                }

                if (inverse) {
                    fft.inv(line_out, line_in);
                } else {
                    fft.fwd(line_out, line_in);
                }

                for (Types::index y = 0; y < fft_y_; ++y) {
                    data[flatten_fft(x, y, z)] = line_out[y];
                }
            }
        }

        line_in.resize(fft_z_);
        line_out.resize(fft_z_);
        for (Types::index y = 0; y < fft_y_; ++y) {
            for (Types::index x = 0; x < fft_x_; ++x) {
                for (Types::index z = 0; z < fft_z_; ++z) {
                    line_in[z] = data[flatten_fft(x, y, z)];
                }

                if (inverse) {
                    fft.inv(line_out, line_in);
                } else {
                    fft.fwd(line_out, line_in);
                }

                for (Types::index z = 0; z < fft_z_; ++z) {
                    data[flatten_fft(x, y, z)] = line_out[z];
                }
            }
        }
    }

    void init_layout_and_kernel_spectrum() {
        if (levels_.levels_x() == 0 || levels_.levels_y() == 0 || levels_.levels_z() == 0) {
            throw std::invalid_argument("TripleToeplitz3x3Fourier requires non-zero level sizes");
        }

        if (levels_.levels_x() % 2 == 0 || levels_.levels_y() % 2 == 0 || levels_.levels_z() % 2 == 0) {
            throw std::invalid_argument(
                "Level sizes must be odd, because they represent Toeplitz shifts [-n+1, n-1]");
        }

        nx_ = (levels_.levels_x() + 1) / 2;
        ny_ = (levels_.levels_y() + 1) / 2;
        nz_ = (levels_.levels_z() + 1) / 2;

        shift_x_ = nx_ - 1;
        shift_y_ = ny_ - 1;
        shift_z_ = nz_ - 1;

        fft_x_ = levels_.levels_x() + nx_ - 1;
        fft_y_ = levels_.levels_y() + ny_ - 1;
        fft_z_ = levels_.levels_z() + nz_ - 1;

        const Types::index total_fft_size = fft_size();

        for (Types::index out_comp = 0; out_comp < block_size_; ++out_comp) {
            for (Types::index in_comp = 0; in_comp < block_size_; ++in_comp) {
                Containers::vector<complex_type> spatial(total_fft_size, complex_type{0.0, 0.0});

                for (Types::index lz = 0; lz < levels_.levels_z(); ++lz) {
                    for (Types::index ly = 0; ly < levels_.levels_y(); ++ly) {
                        for (Types::index lx = 0; lx < levels_.levels_x(); ++lx) {
                            spatial[flatten_fft(lx, ly, lz)] = to_complex(levels_(lx, ly, lz, out_comp, in_comp));
                        }
                    }
                }

                fft3_inplace(spatial, false);
                kernel_spectrum_[out_comp * block_size_ + in_comp] = std::move(spatial);
            }
        }
    }

public:
    explicit TripleToeplitz3x3Fourier(const tensor_type& levels)
        : levels_(levels) {
        init_layout_and_kernel_spectrum();
    }

    explicit TripleToeplitz3x3Fourier(tensor_type&& levels)
        : levels_(std::move(levels)) {
        init_layout_and_kernel_spectrum();
    }

    explicit TripleToeplitz3x3Fourier(const triple_toeplitz_block_type& matrix)
        : levels_(tensor_from_triple_toeplitz_block(matrix)) {
        init_layout_and_kernel_spectrum();
    }

    [[nodiscard]] const tensor_type& levels() const noexcept { return levels_; }

    [[nodiscard]] Types::index nx() const noexcept { return nx_; }
    [[nodiscard]] Types::index ny() const noexcept { return ny_; }
    [[nodiscard]] Types::index nz() const noexcept { return nz_; }

    [[nodiscard]] Types::index rows() const noexcept { return block_size_ * nx_ * ny_ * nz_; }
    [[nodiscard]] Types::index cols() const noexcept { return rows(); }

    [[nodiscard]] vector_type matvec(const vector_type& vec) const {
        if (static_cast<Types::index>(vec.size()) != cols()) {
            throw std::invalid_argument("Vector size does not match TripleToeplitz3x3Fourier::cols()");
        }

        const Types::index total_fft_size = fft_size();

        std::array<Containers::vector<complex_type>, block_size_> vec_spectrum;
        for (auto& component : vec_spectrum) {
            component.assign(total_fft_size, complex_type{0.0, 0.0});
        }

        for (Types::index z = 0; z < nz_; ++z) {
            for (Types::index y = 0; y < ny_; ++y) {
                for (Types::index x = 0; x < nx_; ++x) {
                    const Types::index grid_index = flatten_cell(x, y, z);
                    const Types::index fft_index = flatten_fft(x, y, z);
                    for (Types::index in_comp = 0; in_comp < block_size_; ++in_comp) {
                        vec_spectrum[in_comp][fft_index] = to_complex(vec(flatten_cell_component(x, y, z, in_comp)));
                    }
                }
            }
        }

        for (Types::index in_comp = 0; in_comp < block_size_; ++in_comp) {
            fft3_inplace(vec_spectrum[in_comp], false);
        }

        std::array<Containers::vector<complex_type>, block_size_> result_spectrum;
        for (auto& component : result_spectrum) {
            component.assign(total_fft_size, complex_type{0.0, 0.0});
        }

        for (Types::index out_comp = 0; out_comp < block_size_; ++out_comp) {
            for (Types::index in_comp = 0; in_comp < block_size_; ++in_comp) {
                const auto& kernel = kernel_spectrum_[out_comp * block_size_ + in_comp];
                const auto& source = vec_spectrum[in_comp];
                auto& destination = result_spectrum[out_comp];

                for (Types::index freq_index = 0; freq_index < total_fft_size; ++freq_index) {
                    destination[freq_index] += kernel[freq_index] * source[freq_index];
                }
            }
        }

        for (Types::index out_comp = 0; out_comp < block_size_; ++out_comp) {
            fft3_inplace(result_spectrum[out_comp], true);
        }

        vector_type result = vector_type::Zero(rows());

        for (Types::index z = 0; z < nz_; ++z) {
            for (Types::index y = 0; y < ny_; ++y) {
                for (Types::index x = 0; x < nx_; ++x) {
                    const Types::index fft_index = flatten_fft(x + shift_x_, y + shift_y_, z + shift_z_);

                    for (Types::index out_comp = 0; out_comp < block_size_; ++out_comp) {
                        result(flatten_cell_component(x, y, z, out_comp)) =
                            from_complex(result_spectrum[out_comp][fft_index]);
                    }
                }
            }
        }

        return result;
    }

    [[nodiscard]] vector_type operator*(const vector_type& vec) const {
        return matvec(vec);
    }
};

template <typename scalar_t>
[[nodiscard]] Types::VectorX<scalar_t> operator*(const TripleToeplitz3x3Fourier<scalar_t>& matrix,
                                                 const Types::VectorX<scalar_t>& vec) {
    return matrix.matvec(vec);
}

} // namespace EMW::Math::Fourier

#endif // TRIPLE_TOEPLITZ_3X3_FOURIER_HPP
