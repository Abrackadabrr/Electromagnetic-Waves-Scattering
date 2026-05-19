//
// Created by evgen on 17.01.2026.
//

#ifndef OPERATORK_HPP
#define OPERATORK_HPP

#include "./Utils.hpp"

#include "mesh/volume_mesh/CubeMesh.hpp"

#include "math/matrix/Matrix.hpp"

namespace EMW::Operators::Volume {
/**
 *
 * Дискретизация интегрального оператора кусочно постоянными базисными функциями
 * с константными значениями, равными 1. Сетка подразумевается кубическая из одинаковых кубов,
 * базисные функции являются равномерно линейно независимыми (матрица там диагональная)
 *
 */
class operator_K_over_cube_mesh {
    const Mesh::VolumeMesh::CubeMesh &mesh;
    Types::complex_d wave_number;
    Types::complex_d wave_number_sqr;
    Types::index nearnes_tresholds = 3;
    Types::scalar rTol = 1e-6; // 1e-6
    Types::scalar aTol = 1e-20;
    size_t max_integration_level = 4;    // 4
    size_t max_6d_integration_level = 2; // 2

    struct Idx3d {
        size_t Nx, Ny, Nz;

        struct Int3d {
            Types::integer Nx, Ny, Nz;
        };

        Int3d operator-(const Idx3d &rhs) const {
            return Int3d{static_cast<Types::integer>(Nx) - static_cast<Types::integer>(rhs.Nx),
                         static_cast<Types::integer>(Ny) - static_cast<Types::integer>(rhs.Ny),
                         static_cast<Types::integer>(Nz) - static_cast<Types::integer>(rhs.Nz)};
        }

        bool operator==(const Idx3d &rhs) const { return Nx == rhs.Nx && Ny == rhs.Ny && Nz == rhs.Nz; }

        bool is_near(const Idx3d &rhs, size_t nx, size_t ny, size_t nz) const {
            const auto diff = *this - rhs;
            return std::abs(diff.Nx) <= nx && std::abs(diff.Ny) <= ny && std::abs(diff.Nz) <= nz;
        }
    };

    template <typename matrix_t> struct matrix_and_permutation {
        matrix_t matrix;
        Types::PermutationMatrix P;
    };

    [[nodiscard]] static constexpr decltype(auto) vector_stop_criterion(Types::scalar rTol, Types::scalar aTol) {
        return [rTol, aTol](const Types::Vector3c &v1, const Types::Vector3c &v2) {
            return (v1 - v2).norm() < rTol * v2.norm() + aTol;
        };
    }

    template <typename scalar_type>
    [[nodiscard]] static constexpr decltype(auto) scalar_stop_criterion(Types::scalar rTol, Types::scalar aTol) {
        return [rTol, aTol](scalar_type v1, scalar_type v2) { return std::abs(v1 - v2) < rTol * std::abs(v2) + aTol; };
    }

  public:
    /**
     * @brief Матрица поверхностной части интегрального оператора
     *
     * @param k линейный номер первого куба
     * @param p линейный номер второго куба
     *
     * @return матрица из коэффициентов для заданных кубов
     * @note Интегрирование идет с выделением особенности для близких кубов.
     * Если расстояние между кубами большое (в некотором смысле), то там идёт интегрирование без выделения особенности.
     */
    Types::Matrix3c matrix_2_coef(Types::index k, Types::index p) const;

    /**
     * Расчет объемного интеграла по двум кубам.
     * @param k линейный индекс первого куба
     * @param p линейный индекс второго куба
     * @return значение объемного интеграла, а не матрицу коэффициентов, потому что эта матрица шаровая
     * @note Интегрирование идет с выделением особенности для близких кубов.
     * Если расстояние между кубами большое (в некотором смысле), то там идёт интегрирование без выделения особенности.
     */
    Types::complex_d matrix_3_coef(Types::index k, Types::index p) const;

    /**
     * Расчет объемного члена оператора с выделением особенности
     *
     * @param k_corner, p_corner вершины кубов с минимальными значениями координат
     */
    Types::complex_d volume_part_singularity_extraction(const Types::point_t &k_corner, const Types::point_t &k_center,
                                                        const Types::point_t &p_corner) const;

    /**
     * Расчет матрицы взаимодействия двух кубов при достаточно большом расстоянии между ними
     * Достаточно большое расстояние <=> grad_x div_x F(x-y) достаточно гладкая функция
     *
     * @param k индекс куба с приемниками
     * @param p индекс куба с источниками
     *
     * @return блок взаимодействия между кубами
     */
    Types::Matrix3c far_zone_interaction(Types::index k, Types::index p) const;

    explicit operator_K_over_cube_mesh(Types::complex_d k, const Mesh::VolumeMesh::CubeMesh &mesh)
        : mesh(mesh), wave_number(k), wave_number_sqr(k * k){};

    /**
     * Аппроксимация объемного оператора методом Галёркина.
     * Считается полная матрица вне зависимости от тёплицевой структуры
     */
    [[nodiscard]] Types::MatrixXc compute_galerkin_matrix_dense(Types::scalar l1_basis_function_norm) const;

    [[nodiscard]] Types::Matrix3c galerkin_block_for_cubes(size_t k, size_t p) const;

    void compute_galerkin_matrix_dense_inplace(Types::MatrixXc *p_mat) const;

    /**
     * Аппроксимация методом Галеркина оператора К.
     * Используется специальная трижды тёплицева структура.
     */
    [[nodiscard]] Math::LinAgl::Matrix::TripleToeplitzBlock<Types::complex_d>
    compute_galerkin_matrix(Types::scalar basis_fucntion_module) const;

    /**
     * Расчет оператор K методом Галеркина для одинаковой коллекции кубов, которые отличаются начлаьными индексами
     *
     * @param start_i индексы нижнего левого куба для области приёмников
     * @param start_j индексы нижнего левого куба для области источников
     * @param sizes размеры областей в штуках кубов по трём осям
     * @param l1_basis_function_norm L1 норма кусочно-постоянных базисных функций
     *
     * @return Объект трижды-тёплицевой матрицы, соответствующий блоку взаимодействия
     * кубов источников и кубов приёмников
     *
     * @note Расчет матрицы, соответствующей всей кубической сетке, получается при
     * start_i = {0, 0, 0}; start_j = {0, 0, 0}; sizes = {mesh.nx() - 1, mesh.ny() - 1, mesh.nz() - 1}
     */
    [[nodiscard]] Math::LinAgl::Matrix::TripleToeplitzBlock<Types::complex_d>
    compute_galerkin_matrix(Idx3d start_i, Idx3d start_j, Idx3d sizes, Types::scalar l1_basis_function_norm = 1) const;

    /**
     * Аппроксимация оператора методом Галеркина. Используется трижды тёплицева структура с
     * внутренними блоками размера (3 * Nx) * (3 * Ny) * (3 * Nz). То есть теперь тёплицева структура будет считаться
     * не по отдельным кубам, а по их некоторому объединению в количестве (Nx, Ny, Nz).
     *
     * Такое представление матрицы влияет на нумерацию компонент в векторе неизвестных: матрица фактически бьется на
     * блоки и внутри блоков происходит отдельная нумерация.
     *
     * @param Nx,Ny,Nz размер блока кубов по каждой их координат
     * @param basis_fn_module модуль кусочно-постоянной базисной функции, используемой при дискретизации
     * @return Пару из трижды тёплицевой матрицы и матрицы перестановки, отвечающей
     * за правильную нумерацию неизвестных на сетке
     */
    [[nodiscard]] matrix_and_permutation<Math::LinAgl::Matrix::TripleToeplitzBlock<Types::complex_d>>
        compute_galerkin_matrix_custom_blocksize(size_t Nx, size_t Ny, size_t Nz, Types::scalar basis_fn_module) const;

        /**
         * Аппроксимация оператора методом Галеркина. Используется трижды тёплицева структура с
         * внутренними блоками размера (3 * Nx) * (3 * Ny) * (3 * Nz). То есть теперь тёплицева структура будет считаться
         * не по отдельным кубам, а по их некоторому объединению в количестве (Nx, Ny, Nz).
         *
         * Такое представление матрицы влияет на нумерацию компонент в векторе неизвестных: матрица фактически бьется на
         * блоки и внутри блоков происходит отдельная нумерация.
         *
         * Дополнительно к каждому блоку применяется сжатие через ACA с параметром epsilon, что позволяет
         * максимально сжать исходно плотную матрицу.
         *
         * @param Nx,Ny,Nz размер блока кубов по каждой их координат
         * @param basis_fn_module модуль кусочно-постоянной базисной функции, используемой при дискретизации
         * @return Пару из трижды тёплицевой матрицы и матрицы перестановки, отвечающей
         * за правильную нумерацию неизвестных на сетке
         */
        [[nodiscard]] matrix_and_permutation<Math::LinAgl::Matrix::TripleToeplitzFactoredBlock<Types::complex_d>>
        compute_galerkin_matrix_custom_blocksize_compressed(size_t Nx, size_t Ny, size_t Nz,
                                                            Types::scalar basis_fn_module,
                                                            Types::scalar epsilon) const;

        [[nodiscard]] matrix_and_permutation<Math::LinAgl::Matrix::TripleToeplitzFactoredBlock<Types::complex_d>>
        compute_galerkin_matrix_custom_blocksize_compressed(size_t Nx, size_t Ny, size_t Nz,
                                                    Types::scalar basis_fn_module,
                                                    Types::scalar epsilon, Math::LinAgl::Matrix::TripleToeplitzBlock<Types::complex_d>& dense_mat) const;
    };
} // namespace EMW::Operators::Volume

#endif // OPERATORK_HPP
