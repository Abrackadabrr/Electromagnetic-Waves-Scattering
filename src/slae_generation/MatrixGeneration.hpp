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
 * Контейнер для частей коэффициентов в матрице
 * первая цифра означает номер tau_i (i - номер уравнения, почти что номер строки в матрице)
 * вторая цифра означает номер tau_j (j - номер внутри суммы, почти что номер столбца)
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

/**
 * Рассчитывает первую часть коэффициентов в матрице (поверхностный интеграл) делённый на k^2
 * (скалярное произведение K_1{\tau[m]_j,\sigma_j} and tau[m]_i)
 * @param i - индекс уравнения
 * @param j - индекс внутри уравнения (индекс переменной)
 * @param k - волновое число
 * @param cells - ячейки сетки (в порядке как в SurfaceMesh)
 * @return
 */
Types::complex_d getFirstPartIntegral(Types::index i, Types::index j, Types::complex_d k,
                                      const Containers::vector<Mesh::IndexedCell> &cells);

Types::complex_d getFirstPartIntegral(const Mesh::IndexedCell &cell_i, const Mesh::IndexedCell &cell_j,
                                      Types::complex_d k);

/**
 * Рассчитывает нулевую часть коэффициентов в матрицу (контурный интеграл)
 * @param i - индекс уравнения
 * @param j - индекс внутри уравнения (интекс переменной)
 * @param k - волновое число
 * @param cells - ячейки сетки (в порядке, котором они лежат в SurfaceMesh)
 * @return
 */
Types::Matrix3c getZeroPartIntegral(Types::index i, Types::index j, Types::complex_d k,
                                    const Containers::vector<Mesh::IndexedCell> &cells);

Types::Matrix3c getZeroPartIntegral(const Mesh::IndexedCell &cell_i, const Mesh::IndexedCell &cell_j,
                                    Types::complex_d k);

MatrixCoefs getMatrixCoefs(Types::index i, Types::index j, Types::complex_d k,
                           const Containers::vector<Mesh::IndexedCell> &cells);

MatrixCoefs getMatrixCoefs(const Mesh::IndexedCell &cell_i, const Mesh::IndexedCell &cell_j, Types::complex_d k);

Containers::array<Types::complex_d, 4> getMatrixCoefsInArray(const Mesh::IndexedCell &cell_i, const Mesh::IndexedCell &cell_j, Types::complex_d k);
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

Containers::array<Types::complex_d, 4> getMatrixCoefsInArray(const Mesh::IndexedCell &cell_i, const Mesh::IndexedCell &cell_j, Types::complex_d k);

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
