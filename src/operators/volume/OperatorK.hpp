//
// Created by evgen on 17.01.2026.
//

#ifndef OPERATORK_HPP
#define OPERATORK_HPP

#include "./Utils.hpp"

#include "mesh/volume_mesh/CubeMesh.hpp"

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

  public:
    explicit operator_K_over_cube_mesh(Types::complex_d k, const Mesh::VolumeMesh::CubeMesh &mesh) : mesh(mesh), wave_number(k) {};

    [[nodiscard]] Types::MatrixXc get_galerkin_matrix() const;

    void get_galerkin_matrix_inplace(Types::MatrixXc * p_mat) const;
};

} // namespace EMW::Operators::Volume

#endif // OPERATORK_HPP
