#include <gtest/gtest.h>

#include "math/fourier/TripleToeplitz3x3Fourier.hpp"
#include "math/matrix/Matrix.hpp"

#include <mesh/volume_mesh/CubeMeshWithData.hpp>
#include <operators/volume/OperatorK.hpp>

using namespace EMW;

template <typename scalar_t>
using tensor_t = Math::Fourier::TripleToeplitz3x3Tensor<scalar_t>;

template <typename scalar_t>
using matrix_t = Math::Fourier::TripleToeplitz3x3Fourier<scalar_t>;

template <typename scalar_t>
Types::VectorX<scalar_t> direct_matvec(const tensor_t<scalar_t> &levels, const Types::VectorX<scalar_t> &x) {
    const Types::index nx = (levels.levels_x() + 1) / 2;
    const Types::index ny = (levels.levels_y() + 1) / 2;
    const Types::index nz = (levels.levels_z() + 1) / 2;

    const Types::index shift_x = nx - 1;
    const Types::index shift_y = ny - 1;
    const Types::index shift_z = nz - 1;

    const auto flatten = [](Types::index ix, Types::index iy, Types::index iz, Types::index sx, Types::index sy) {
        return (iz * sy + iy) * sx + ix;
    };

    Types::VectorX<scalar_t> y = Types::VectorX<scalar_t>::Zero(3 * nx * ny * nz);

    for (Types::index oz = 0; oz < nz; ++oz) {
        for (Types::index oy = 0; oy < ny; ++oy) {
            for (Types::index ox = 0; ox < nx; ++ox) {
                const Types::index out_site = flatten(ox, oy, oz, nx, ny);

                for (Types::index out_comp = 0; out_comp < 3; ++out_comp) {
                    scalar_t acc = scalar_t{};

                    for (Types::index iz = 0; iz < nz; ++iz) {
                        for (Types::index iy = 0; iy < ny; ++iy) {
                            for (Types::index ix = 0; ix < nx; ++ix) {
                                const Types::index in_site = flatten(ix, iy, iz, nx, ny);

                                const Types::integer lx_signed = static_cast<Types::integer>(ox) -
                                                                 static_cast<Types::integer>(ix) + static_cast<
                                                                     Types::integer>(shift_x);
                                const Types::integer ly_signed = static_cast<Types::integer>(oy) -
                                                                 static_cast<Types::integer>(iy) + static_cast<
                                                                     Types::integer>(shift_y);
                                const Types::integer lz_signed = static_cast<Types::integer>(oz) -
                                                                 static_cast<Types::integer>(iz) + static_cast<
                                                                     Types::integer>(shift_z);

                                const Types::index lx = static_cast<Types::index>(lx_signed);
                                const Types::index ly = static_cast<Types::index>(ly_signed);
                                const Types::index lz = static_cast<Types::index>(lz_signed);

                                for (Types::index in_comp = 0; in_comp < 3; ++in_comp) {
                                    acc += levels(lx, ly, lz, out_comp, in_comp) * x(in_site * 3 + in_comp);
                                }
                            }
                        }
                    }

                    y(out_site * 3 + out_comp) = acc;
                }
            }
        }
    }

    return y;
}

TEST(TRIPLE_TOEPLITZ_3X3_FOURIER, MATCHES_DIRECT_REAL) {
    tensor_t<Types::scalar> levels(3, 3, 3, 0.0);

    for (Types::index lz = 0; lz < levels.levels_z(); ++lz) {
        for (Types::index ly = 0; ly < levels.levels_y(); ++ly) {
            for (Types::index lx = 0; lx < levels.levels_x(); ++lx) {
                for (Types::index out_comp = 0; out_comp < 3; ++out_comp) {
                    for (Types::index in_comp = 0; in_comp < 3; ++in_comp) {
                        levels(lx, ly, lz, out_comp, in_comp) =
                            0.7 * static_cast<Types::scalar>(lx + 1) +
                            0.3 * static_cast<Types::scalar>(ly + 1) -
                            0.2 * static_cast<Types::scalar>(lz + 1) +
                            0.05 * static_cast<Types::scalar>(out_comp + 1) -
                            0.04 * static_cast<Types::scalar>(in_comp + 1);
                    }
                }
            }
        }
    }

    matrix_t<Types::scalar> matrix(levels);
    Types::VectorX<Types::scalar> x = Types::VectorX<Types::scalar>::Zero(matrix.cols());

    for (Types::index i = 0; i < static_cast<Types::index>(x.size()); ++i) {
        x(i) = 0.5 + 0.1 * static_cast<Types::scalar>(i);
    }

    const auto y_fft = matrix * x;
    const auto y_direct = direct_matvec(levels, x);

    ASSERT_NEAR((y_fft - y_direct).norm(), 0.0, 1e-10);
}

TEST(TRIPLE_TOEPLITZ_3X3_FOURIER, MATCHES_DIRECT_COMPLEX_NON_CUBIC_GRID) {
    tensor_t<Types::complex_d> levels(3, 5, 3, Types::complex_d{0.0, 0.0});

    for (Types::index lz = 0; lz < levels.levels_z(); ++lz) {
        for (Types::index ly = 0; ly < levels.levels_y(); ++ly) {
            for (Types::index lx = 0; lx < levels.levels_x(); ++lx) {
                for (Types::index out_comp = 0; out_comp < 3; ++out_comp) {
                    for (Types::index in_comp = 0; in_comp < 3; ++in_comp) {
                        const Types::scalar real =
                            0.21 * static_cast<Types::scalar>(lx + 1) -
                            0.17 * static_cast<Types::scalar>(ly + 1) +
                            0.09 * static_cast<Types::scalar>(lz + 1) +
                            0.03 * static_cast<Types::scalar>(out_comp + 1) -
                            0.02 * static_cast<Types::scalar>(in_comp + 1);

                        const Types::scalar imag =
                            -0.11 * static_cast<Types::scalar>(lx + 1) +
                            0.07 * static_cast<Types::scalar>(ly + 1) -
                            0.05 * static_cast<Types::scalar>(lz + 1) +
                            0.04 * static_cast<Types::scalar>(out_comp + 1) +
                            0.06 * static_cast<Types::scalar>(in_comp + 1);

                        levels(lx, ly, lz, out_comp, in_comp) = Types::complex_d{real, imag};
                    }
                }
            }
        }
    }

    matrix_t<Types::complex_d> matrix(levels);
    Types::VectorX<Types::complex_d> x = Types::VectorX<Types::complex_d>::Zero(matrix.cols());

    for (Types::index i = 0; i < static_cast<Types::index>(x.size()); ++i) {
        x(i) = Types::complex_d{
            0.15 + 0.03 * static_cast<Types::scalar>(i),
            -0.2 + 0.01 * static_cast<Types::scalar>(i % 7)
        };
    }

    const auto y_fft = matrix * x;
    const auto y_direct = direct_matvec(levels, x);

    ASSERT_NEAR((y_fft - y_direct).norm(), 0.0, 1e-10 * y_direct.norm());
}

TEST(TRIPLE_TOEPLITZ_3X3_FOURIER, CENTRAL_LEVEL_REDUCES_TO_LOCAL_BLOCK_ACTION) {
    constexpr Types::index nx = 3;
    constexpr Types::index ny = 2;
    constexpr Types::index nz = 2;

    tensor_t<Types::scalar> levels(2 * nx - 1, 2 * ny - 1, 2 * nz - 1, 0.0);

    Types::Matrix3d local_block;
    local_block << 1.5, -0.3, 0.2,
        -0.7, 0.8, 0.1,
        0.4, 0.6, 1.2;

    const Types::index cx = nx - 1;
    const Types::index cy = ny - 1;
    const Types::index cz = nz - 1;

    for (Types::index out_comp = 0; out_comp < 3; ++out_comp) {
        for (Types::index in_comp = 0; in_comp < 3; ++in_comp) {
            levels(cx, cy, cz, out_comp, in_comp) = local_block(out_comp, in_comp);
        }
    }

    matrix_t<Types::scalar> matrix(levels);
    Types::VectorX<Types::scalar> x = Types::VectorX<Types::scalar>::Zero(matrix.cols());

    for (Types::index i = 0; i < static_cast<Types::index>(x.size()); ++i) {
        x(i) = -0.4 + 0.07 * static_cast<Types::scalar>(i);
    }

    Types::VectorX<Types::scalar> y_expected = Types::VectorX<Types::scalar>::Zero(matrix.rows());
    const Types::index total_sites = nx * ny * nz;

    for (Types::index site = 0; site < total_sites; ++site) {
        Types::Vector3d local_x;
        local_x(0) = x(site * 3 + 0);
        local_x(1) = x(site * 3 + 1);
        local_x(2) = x(site * 3 + 2);

        const Types::Vector3d local_y = local_block * local_x;

        y_expected(site * 3 + 0) = local_y(0);
        y_expected(site * 3 + 1) = local_y(1);
        y_expected(site * 3 + 2) = local_y(2);
    }

    const auto y_fft = matrix * x;

    ASSERT_NEAR((y_fft - y_expected).norm(), 0.0, 1e-10);
}

TEST(TRIPLE_TOEPLITZ_3X3_FOURIER, CONSTRUCTOR_FROM_TRIPLE_TOEPLITZ_BLOCK_MATCHES_DENSE_MATVEC) {
    constexpr Types::index nx = 30;
    constexpr Types::index ny = 30;
    constexpr Types::index nz = 30;

    auto toeplitz = Math::LinAgl::Matrix::ZeroTripleToeplitzBlock<Types::complex_d>(nx, ny, nz, 3);

    std::cout << "Initializing toeplitz matrix ..." << std::endl;

    for (Types::index i3 = 0; i3 < nz; ++i3) {
        for (Types::index j3 = 0; j3 < nz; ++j3) {
            for (Types::index i2 = 0; i2 < ny; ++i2) {
                for (Types::index j2 = 0; j2 < ny; ++j2) {
                    for (Types::index i1 = 0; i1 < nx; ++i1) {
                        for (Types::index j1 = 0; j1 < nx; ++j1) {
                            auto &block = toeplitz.get_block(i3, j3).get_block(i2, j2).get_block(i1, j1);
                            block = Types::Matrix3c::Random();
                        }
                    }
                }
            }
        }
    }
    std::cout << "Initializing tensor for FFT ..." << std::endl;
    const matrix_t<Types::complex_d> fft_matrix(toeplitz);

    Types::VectorXc x = Types::VectorXd::Zero(toeplitz.cols());
    for (Types::index i = 0; i < static_cast<Types::index>(x.size()); ++i) {
        x(i) = 3 - 2 * static_cast<Types::scalar>(i) + 10 * static_cast<Types::scalar>(i % 11);
    }

    std::cout << "dense start" << std::endl;
    const Types::VectorXc y_dense = toeplitz.matvec(x);
    std::cout << "dense done" << std::endl;
    const Types::VectorXc y_fft = fft_matrix * x;
    ASSERT_NEAR((y_dense - y_fft).norm(), 0.0, 1e-14 * y_dense.norm());
}

constexpr Types::scalar SPHERE_RADUIS = 0.5;

Types::scalar permittivity_distribution(const Types::point_t &x) {
    return x.norm() < SPHERE_RADUIS ? 2.56 : 1;
}

#include "Utils.hpp"

TEST(TRIPLE_TOEPLITZ_3X3_FOURIER, VIE_TEST) {
    std::cout << "Initializing toeplitz matrix ..." << std::endl;

    constexpr Types::scalar cube_length = 2.1 * SPHERE_RADUIS;
    constexpr Types::index Nx = 9;
    constexpr Types::index Ny = 9;
    constexpr Types::index Nz = 9;
    constexpr Types::scalar mesh_one_axis_size = cube_length / (Nx - 1);
    constexpr Types::scalar basis_fn_module = 1. / (mesh_one_axis_size * std::sqrt(mesh_one_axis_size));
    Mesh::VolumeMesh::CubeMeshWithData mesh{Types::point_t{-cube_length / 2, -cube_length / 2, -cube_length / 2},
                                            (Nx - 1) * mesh_one_axis_size,
                                            (Ny - 1) * mesh_one_axis_size, (Nz - 1) * mesh_one_axis_size, Nx, Ny, Nz};

    Operators::Volume::operator_K_over_cube_mesh operator_K{{1., 0.}, mesh};
    auto toeplitz = operator_K.compute_galerkin_matrix(basis_fn_module);

    auto [mat, perm] = operator_K.compute_galerkin_matrix_custom_blocksize(
        4, 4, 4, basis_fn_module);

    std::cout << "Initializing tensor for FFT ..." << std::endl;
    const matrix_t<Types::complex_d> fft_matrix(toeplitz);

    Types::VectorXc x = Types::VectorXc::Zero(toeplitz.cols());
    for (Types::index i = 0; i < static_cast<Types::index>(x.size()); ++i) {
        x(i) = 3 - 2 * static_cast<Types::scalar>(i) + 10 * static_cast<Types::scalar>(i % 11);
    }

    std::cout << "dense start" << std::endl;
    // Плотная умножалка
    const Types::VectorXc x_permuted = perm * x;
    auto start = std::chrono::steady_clock::now();
    const Types::VectorXc y_dense_permuted = mat.matvec(x_permuted);
    auto end = std::chrono::steady_clock::now();
    const Types::VectorXc y_dense = perm.transpose() * y_dense_permuted;
    std::cout << "dense done: " << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() <<
        std::endl;

    // Фурье-умножалка
    start = std::chrono::steady_clock::now();
    const Types::VectorXc y_fft = fft_matrix * x;
    end = std::chrono::steady_clock::now();
    std::cout << "FFT done: " << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() <<
        std::endl;

    // Сравнение для плотной умножалки через toeplitz
    const auto y_from_dense_toeplitz = toeplitz.matvec(x);
    ASSERT_NEAR((y_dense - y_from_dense_toeplitz).norm(), 0.0, 1e-12 * y_from_dense_toeplitz.norm());

    ASSERT_NEAR((y_dense - y_fft).norm(), 0.0, 1e-12 * y_dense.norm());
}

TEST(TRIPLE_TOEPLITZ_3X3_FOURIER, VIE_COMPARISON_WITH_SKELETON_FORMAT) {
    std::cout << "Initializing toeplitz matrix ..." << std::endl;

    constexpr Types::scalar cube_length = 2.1 * SPHERE_RADUIS;
    constexpr Types::index Nx = 41;
    constexpr Types::index Ny = 41;
    constexpr Types::index Nz = 41;
    constexpr Types::scalar mesh_one_axis_size = cube_length / (Nx - 1);
    constexpr Types::scalar basis_fn_norm = 1. / (mesh_one_axis_size * std::sqrt(mesh_one_axis_size));
    Mesh::VolumeMesh::CubeMeshWithData mesh{Types::point_t{-cube_length / 2, -cube_length / 2, -cube_length / 2},
                                            (Nx - 1) * mesh_one_axis_size,
                                            (Ny - 1) * mesh_one_axis_size, (Nz - 1) * mesh_one_axis_size, Nx, Ny, Nz};

    Operators::Volume::operator_K_over_cube_mesh operator_K{{1., 0.}, mesh};
    auto toeplitz = operator_K.compute_galerkin_matrix(basis_fn_norm);
    // Части, на которые режем матрицу
    size_t nx, ny, nz;
    nx = 4;
    ny = 4;
    nz = 4;
    const auto [mat, perm] = operator_K.compute_galerkin_matrix_custom_blocksize_compressed(
        nx, ny, nz, 1. / (mesh_one_axis_size * std::sqrt(mesh_one_axis_size)), 1e-4);

    std::cout << "Initializing tensor for FFT ..." << std::endl;

    const matrix_t<Types::complex_d> fft_matrix(toeplitz);

    Types::VectorXc x = Types::VectorXc::Random(toeplitz.cols());

    std::cout << "// --- PERFORMANCE REPORT --- //" << std::endl;

    // Плотная умножалка с замером времени
    Types::VectorXc x_permuted = perm * x;
    Types::VectorXc y_skeleton = Types::VectorXc::Zero(mat.rows());
    auto start = std::chrono::steady_clock::now();
    y_skeleton = mat.matvec(x_permuted);
    auto end = std::chrono::steady_clock::now();
    y_skeleton = perm.transpose() * y_skeleton;
    std::cout << "dense done in " << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() <<
        std::endl;

    // Фурье-умножалка с замером времени
    start = std::chrono::steady_clock::now();
    const Types::VectorXc y_fft = fft_matrix.matvec(x);
    end = std::chrono::steady_clock::now();
    std::cout << "FFT done in " << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() <<
        std::endl;

    // Смотрим на точность умножения на аппроксимацию матрицы

    // Для сравнение умножаем на плотную матрицу
    const auto y_dense = toeplitz.matvec(x);

    // Сравниваемся с умножением на тёплицеву матрицу
    std::cout << "DENSE vs SKELETON = " << (y_dense - y_skeleton).norm() / y_dense.norm() << std::endl;

    // Метод через Фурье должен совпадать с умножением на трижды тёплицеву матрицу
    std::cout << "FFT vs DENSE = " << (y_dense - y_fft).norm() / y_dense.norm() << std::endl;

    // Ошибка в умножении на аппроксимацию матрицы по сравнению с FFT
    std::cout << "FFT vs SKELETON = " << (y_fft - y_skeleton).norm() / y_dense.norm() << std::endl;

    // Параметры сжатой матрицы (анализ параметризации)
    std::cout << Utils::get_memory_usage(mat) << std::endl;
    std::cout << "Мозаичный ранг = " << Utils::get_elements_for_parametrization(mat) / mat.cols() << std::endl;
    std::cout << "log2(Nx * Ny * Nz) = " << std::log2(Nx * Ny * Nz) << std::endl;

#if 0
    // Ошибка в аппроксимации матрицы:
    const auto [mat_permuted_dense, _] = operator_K.compute_galerkin_matrix_custom_blocksize(
        nx, ny, nz, 1. / (mesh_one_axis_size * std::sqrt(mesh_one_axis_size)));

    std::cout << "DENSE_MAT vs SKELETON_MAT = " << Utils::relative_frobenius_error(mat_permuted_dense, mat) <<
        std::endl;
#endif
}
