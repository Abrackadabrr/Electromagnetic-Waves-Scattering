//
// Created by evgen on 31.01.24.
//

#ifndef ELECTROMAGNETIC_WAVES_SCATTERING_MATRIXGENERATIONFUNCTION_HPP
#define ELECTROMAGNETIC_WAVES_SCATTERING_MATRIXGENERATIONFUNCTION_HPP

#include "mesh/MeshTypes.hpp"
#include "mesh/SurfaceMesh.hpp"
#include "types/Types.hpp"

namespace EMW::Matrix {

enum class operator_t { R = 0, K = 1 };

/**
 * Contain parts of matrix coefficients
 * first figure means number of tau_i (i - number of an equation)
 * second figure means number of tau_j (j - number of a part in sum)
 */
struct MatrixCoefs {
    Types::complex_d a11;
    Types::complex_d a12;
    Types::complex_d a21;
    Types::complex_d a22;
};

namespace DiscreteK {
struct ContourIntegralParts {
    Types::Vector3c ab;
    Types::Vector3c bc;
    Types::Vector3c cd;
    Types::Vector3c da;
};

/***
 * Returning the first part of matrix coefficient (surface integral) divided by k^2
 * (aka dot product of K_1{ tau[m]_j, \sigma_j} and tau[m]_i)
 * @param i - index of an equation
 * @param j - index inside the equation
 * @param k - wave number
 * @param cells - cells in the mesh (with collocation nodes)
 * @return
 */
Types::complex_d getFirstPartIntegral(Types::index i, Types::index j, Types::complex_d k,
                                      const Containers::vector<Mesh::IndexedCell> &cells);

/**
 * Returning the zero part of matrix coefficients (contour integral)
 * @param i - index of an equation
 * @param j - index inside the equation
 * @param k - wave number
 * @param cells - cells in the mesh (with collocation nodes)
 * @return
 */
Types::Matrix3c getZeroPartIntegral(Types::index i, Types::index j, Types::complex_d k,
                                    const Containers::vector<Mesh::IndexedCell> &cells);

MatrixCoefs getMatrixCoefs(Types::index i, Types::index j, Types::complex_d k,
                           const Containers::vector<Mesh::IndexedCell> &cells);
} // namespace DiscreteK

namespace DiscreteR {
MatrixCoefs getMatrixCoefs(Types::index i, Types::index j, Types::complex_d k,
                           const Containers::vector<Mesh::IndexedCell> &cells);

/**
 * Расчет коэффициентов матрицы от ячейки cell_j в точке cell_i.collPoint_
 * @param cell_i ячейка, в которой есть точка коллокации
 * @param cell_j ячейка, по которой производится расчет интеграла
 * @param k волновое число
 * @return четыре коэффициента в матрицу
 *
 * @brief Функция, созданная для расчета матриц, где поверхность интегрирования не совпадает с поверьностью, где
 * лежат точки коллокации
 */
MatrixCoefs getMatrixCoefs(const Mesh::IndexedCell &cell_i, const Mesh::IndexedCell &cell_j, Types::complex_d k);

} // namespace DiscreteR

/**
* Расчет матрицы оператора K для токов и значений образа на одной и той же поверхности
*/
Types::MatrixXc getMatrixK(Types::complex_d k, const Mesh::SurfaceMesh &surface_mesh);

/**
* Расчет матрицы оператора K для токов и значений образа на разных поверхностях
*/
Types::MatrixXc getMatrixK(Types::complex_d k, const Mesh::SurfaceMesh &integration_mesh,
                                               const Mesh::SurfaceMesh &mesh_with_coll_points);

/**
* Расчет матрицы оператора R для токов и значений образа на одной и той же поверхности
*/
Types::MatrixXc getMatrixR(Types::complex_d k, const Mesh::SurfaceMesh &surface_mesh);

/**
* Расчет матрицы оператора R для токов и значений образа на разных поверхностях
*/
Types::MatrixXc getMatrixR(Types::complex_d k, const Mesh::SurfaceMesh &integration_mesh,
                                               const Mesh::SurfaceMesh &mesh_with_coll_points);

/**
* Расчет матрицы оператора векторного умножения на нормаль справа
*/
Types::MatrixXc getMatrixCrossNormal(Types::complex_d k, const Mesh::SurfaceMesh &surface_mesh);

/**
 * Расчет матрицы единичного оператора на поверхности
 */
Types::MatrixXc getMatrixIdentity(Types::complex_d k, const Mesh::SurfaceMesh &surface_mesh);

} // namespace EMW::Matrix

#endif // ELECTROMAGNETIC_WAVES_SCATTERING_MATRIXGENERATIONFUNCTION_HPP
