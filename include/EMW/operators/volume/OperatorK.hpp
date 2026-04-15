//
// Created by evgen on 17.01.2026.
//

#ifndef OPERATORK_HPP
#define OPERATORK_HPP

#include "./Utils.hpp"

#include "mesh/volume_mesh/CubeMesh.hpp"

#include "math/matrix/Matrix.hpp"

namespace EMW::Operators::Volume
{
    /**
     *
     * Дискретизация интегрального оператора кусочно постоянными базисными функциями
     * с константными значениями, равными 1. Сетка подразумевается кубическая из одинаковых кубов,
     * базисные функции являются равномерно линейно независимыми (матрица там диагональная)
     *
     */
    class operator_K_over_cube_mesh
    {
        const Mesh::VolumeMesh::CubeMesh& mesh;
        Types::complex_d wave_number;
        Types::index nearnes_tresholds = 2;

        struct Idx3d
        {
            size_t Nx, Ny, Nz;

            struct Int3d
            {
                Types::integer Nx, Ny, Nz;
            };

            Int3d operator-(const Idx3d& rhs) const
            {
                return Int3d{
                    static_cast<Types::integer>(Nx) - static_cast<Types::integer>(rhs.Nx),
                    static_cast<Types::integer>(Ny) - static_cast<Types::integer>(rhs.Ny),
                    static_cast<Types::integer>(Nz) - static_cast<Types::integer>(rhs.Nz)
                };
            }

            bool operator==(const Idx3d& rhs) const
            {
                return Nx == rhs.Nx && Ny == rhs.Ny && Nz == rhs.Nz;
            }

            bool is_near(const Idx3d& rhs, size_t nx, size_t ny, size_t nz) const
            {
                const auto diff = *this - rhs;
                return std::abs(diff.Nx) <= nx && std::abs(diff.Ny) <= ny && std::abs(diff.Nz) <= nz;
            }
        };

        template <typename matrix_t>
        struct matrix_and_permutation
        {
            matrix_t matrix;
            Types::PermutationMatrix P;
        };

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
         * @brief Ньютонов потенциал куба с линейным номером k в сетке mesh
         *
         * @param k линейный индекс куба в сетке mesh
         * @param r точка, в которой считаем значение Ньютонова потенциала
         *
         * @note: эта функция должна быть заменена аналитическим интегрированием по кубу после того, как оно будет отлажено
         * и перевыведено из статей с выводом этого дела, а пока что там какая-то ошибпбочка.
         */
        [[nodiscard]] Types::scalar newton_potential_of_cube(Types::index k, Types::point_t r) const;

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

        explicit operator_K_over_cube_mesh(Types::complex_d k,
                                           const Mesh::VolumeMesh::CubeMesh& mesh) : mesh(mesh), wave_number(k)
        {
        };

        /**
         * Аппроксимация объемного оператора методом Галёркина.
         * Считается полная матрица вне зависимости от тёплицевой структуры
         */
        [[nodiscard]] Types::MatrixXc compute_galerkin_matrix_dense(Types::scalar l1_basis_function_norm) const;

        [[nodiscard]] Types::Matrix3c galerkin_block_for_cubes(size_t k, size_t p) const;

        void compute_galerkin_matrix_dense_inplace(Types::MatrixXc* p_mat) const;

        /**
         * Аппроксимация методом Галеркина оператора К.
         * Используется специальная трижды тёплицева структура.
         */
        [[nodiscard]] Math::LinAgl::Matrix::TripleToeplitzBlock<Types::complex_d> compute_galerkin_matrix(
            Types::scalar basis_fucntion_module) const;

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
        [[nodiscard]] Math::LinAgl::Matrix::TripleToeplitzBlock<Types::complex_d> compute_galerkin_matrix(
            Idx3d start_i, Idx3d start_j, Idx3d sizes, Types::scalar l1_basis_function_norm = 1) const;

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
