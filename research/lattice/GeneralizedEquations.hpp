//
// Created by evgen on 18.01.2025.
//

#ifndef GENERALIZEDEQUATIONS_HPP
#define GENERALIZEDEQUATIONS_HPP

#include "mesh/SurfaceMesh.hpp"

#include "slae_generation/MatrixGeneration.hpp"

#include "types/Types.hpp"

#include "math/MathConstants.hpp"

#include "geometry/PeriodicStructure.hpp"

using namespace EMW;

namespace GeneralizedEquations {

namespace detail {
/**
 * Вычисляем матрицу для "самодействия"
 */
inline EMW::Types::MatrixXc diagonal(const Mesh::SurfaceMesh &mesh_all, const Mesh::SurfaceMesh &mesh_zero,
                                     Types::scalar a, Types::complex_d k) {
    // описываем размеры итоговой матрицы
    const Types::index N = mesh_all.getCells().size();
    const Types::index K = mesh_zero.getCells().size();
    const Types::index size = 2 * (N + K);
    Types::MatrixXc A(size, size);

    const Types::complex_d epsilon{1., 0};
    const Types::complex_d mu{1., 0};

    // расчет коэффициента импеданса и вспомогательных констант
    const Types::complex_d beta = std::sqrt(k * k - EMW::Math::Constants::PI_square<Types::scalar>() / (a * a));
#ifdef OLD
    const Types::complex_d z = -k * Math::Constants::one_div_c / beta;
    const Types::complex_d with_e = Math::Constants::i * Math::Constants::mu_0_c / k;
    const Types::complex_d with_mu = Math::Constants::i * Math::Constants::e_0_c / k;
#else
    const Types::complex_d z = -k * mu / beta;
    const Types::complex_d with_e = Math::Constants::i / (k * epsilon);
    const Types::complex_d with_mu = Math::Constants::i / (k * mu);
#endif

    // генерируем необходимые подматрицы
    auto K_all = Matrix::getMatrixK(k, mesh_all);
    auto R_0 = Matrix::getMatrixR(k, mesh_zero, mesh_all);
    auto K_0 = Matrix::getMatrixK(k, mesh_zero);
    auto R_all = Matrix::getMatrixR(k, mesh_all, mesh_zero);
    // тривиальные операторы нужных размеров из уравнений
    // оператор векторного умножения на нормаль справа
    const Types::MatrixXc C = 0.5 * Matrix::getMatrixCrossNormal(k, mesh_zero);
    // единичный оператор
    const auto I_m = Matrix::getMatrixIdentity(k, mesh_zero);

    // 1) Дописываем слагаемые в уравнение на поверхности sigma_0
    // ВНИМАНИЕ: тут идёт существенный упор на то, что мы знаем, что sigma_0 лежит последняя в sigma_all
    R_0.block(N - K, 0, K, 2 * K) -= C.block(0, 0, K, 2 * K);
    R_0.block(2 * N - K, 0, K, 2 * K) -= C.block(K, 0, K, 2 * K);

    R_all.block(0, N - K, 2 * K, K) += C.block(0, 0, 2 * K, K);
    R_all.block(0, 2 * N - K, 2 * K, K) += C.block(0, K, 2 * K, K);

    // 2) Собираем матрицу
    A.block(0, 0, 2 * N, 2 * N) = with_e * K_all;
    A.block(0, 2 * N, 2 * N, 2 * K) = -R_0;
    A.block(2 * N, 0, 2 * K, 2 * N) = z * R_all;
    A.block(2 * N, 2 * N, 2 * K, 2 * K) = z * with_mu * K_0 + I_m;

    // 3) Готово
    return A;
}

/**
 * Вспомогательная функция для вычисления взаимодействия двух антенн
 * Сетки mesh_to_integrate и mesh_with_collocation_points не обязаны совпадать,
 * функция будет работать для произвольных сеток
 *
 * Результат этой функции при передаче индентичных сеток не эквивалентен расчету матрицы "самодействия",
 * так как в этом расчете не учитываются неаддитивные орепаторы, зависящие только от сетки
 * с точками коллокации (единичный оператор и оператор векторного умножения на нормаль справа)
 */
inline Types::MatrixXc submatrix(const Mesh::SurfaceMesh &mest_to_integrate,
                                 const Mesh::SurfaceMesh &mest_with_collocation_points, const Types::scalar a,
                                 const Types::complex_d k) {
    // описываем специфические части сеток
    // части с электрическими токами
    const Mesh::SurfaceMesh &mesh_to_integrate_all = mest_to_integrate;
    const Mesh::SurfaceMesh mesh_to_integrate_zero =
        mest_to_integrate.getSubmesh(Mesh::IndexedCell::Tag::WAVEGUIDE_CROSS_SECTION);
    const Mesh::SurfaceMesh &mesh_collocation_all = mest_with_collocation_points;
    const Mesh::SurfaceMesh mesh_collocation_zero =
        mest_with_collocation_points.getSubmesh(Mesh::IndexedCell::Tag::WAVEGUIDE_CROSS_SECTION);

    // Надо понять какие размеры будем иметь итоговая матрица
    const Types::index N_cols = mesh_to_integrate_all.getCells().size();
    const Types::index K_cols = mesh_to_integrate_zero.getCells().size();
    const Types::index cols = 2 * (N_cols + K_cols);
    const Types::index N_rows = mesh_collocation_all.getCells().size();
    const Types::index K_rows = mesh_collocation_zero.getCells().size();
    const Types::index rows = 2 * (N_rows + K_rows);
    Types::MatrixXc A(rows, cols);

    const Types::complex_d epsilon{1., 0};
    const Types::complex_d mu{1., 0};

    // расчет коэффициента импеданса и вспоногательных констант
    const Types::complex_d beta = std::sqrt(k * k - EMW::Math::Constants::PI_square<Types::scalar>() / (a * a));
#ifdef OLD
    const Types::complex_d z = -k * Math::Constants::one_div_c / beta;
    const Types::complex_d with_e = Math::Constants::i * Math::Constants::mu_0_c / k;
    const Types::complex_d with_mu = Math::Constants::i * Math::Constants::e_0_c / k;
#else
    const Types::complex_d z = -k * mu / beta;
    const Types::complex_d with_e = Math::Constants::i / (k * epsilon);
    const Types::complex_d with_mu = Math::Constants::i / (k * mu);
#endif

    // 1) генерируем необходимые подматрицы
    // матрица оператора K размерами 2 * N_rows x 2 * N_cols
    auto K_all = Matrix::getMatrixK(k, mesh_to_integrate_all, mesh_collocation_all);
    // матрица оператора R размерами 2 * N_rows x 2 * K_cols
    auto R_0 = Matrix::getMatrixR(k, mesh_to_integrate_zero, mesh_collocation_all);
    // матрица оператора K размерами 2 * K_rows x 2 * K_cols
    auto K_0 = Matrix::getMatrixK(k, mesh_to_integrate_zero, mesh_collocation_zero);
    // матрица оператора K размерами 2 * K_rows x 2 * N_cols
    auto R_all = Matrix::getMatrixR(k, mesh_to_integrate_all, mesh_collocation_zero);

    // 2) Неаддитивные по поверхности операторы отваливаются при расчете недиагональных подматриц
    // поэтому их тут не будет. Далее все идентично как и в случае одной антенны, блоки распихиваем точно так же

    // 3) Собираем матрицу
    A.block(0, 0, 2 * N_rows, 2 * N_cols) = with_e * K_all;
    A.block(0, 2 * N_cols, 2 * N_rows, 2 * K_cols) = -R_0;
    A.block(2 * N_rows, 0, 2 * K_rows, 2 * N_cols) = z * R_all;
    A.block(2 * N_rows, 2 * N_cols, 2 * K_rows, 2 * K_cols) = z * with_mu * K_0;

    // Готово
    return A;
}
} // namespace detail

/**
 * Вычисляем матрицу для системы уравнений для периодической структуры втупую, без хранения и аппроксимации
 */
template <Types::index N1, Types::index N2>
EMW::Types::MatrixXc getMatrix(const Geometry::PeriodicStructure<N1, N2> &geometry, Types::scalar a,
                               Types::complex_d k) {
    // Создаем итоговый результат
    // Размер итоговой матрицы
    const Types::index size_of_block = 2 * (geometry.get_mesh_matrix()[0][0].getCells().size() +
        geometry.get_mesh_matrix()[0][0].getSubmesh(Mesh::IndexedCell::Tag::WAVEGUIDE_CROSS_SECTION).getCells().size());
    const Types::index size = size_of_block * geometry.size();
    std::cout << size << std::endl;
    Types::MatrixXc A(size, size);

    // 2) Поэтапный расчет внедиагональных блоков без учета того, что какие-то из них одинаковые
    for (int i = 0; i < geometry.size(); i++) {
        for (int j = 0; j < geometry.size(); j++)
            if (i != j) {
                // если у нас недиагональный блок, то считаем матрицу
                // A_ij блок показывает как k сетка влияет на поле в точках коллокации на j сетке
                const auto A_ij = detail::submatrix(geometry.get(j), geometry.get(i), a, k);
                // рассчитывем место, где этот блок должен находится
                const Types::index first_row = i * size_of_block;
                const Types::index first_col = j * size_of_block;
                A.block(first_row, first_col, size_of_block, size_of_block) = A_ij;
            }
        A.block(i * size_of_block, i * size_of_block, size_of_block, size_of_block) = detail::diagonal(geometry.get(i),
        geometry.get(i).getSubmesh(Mesh::IndexedCell::Tag::WAVEGUIDE_CROSS_SECTION), a, k);
    }
    return A;
}
}
#endif //GENERALIZEDEQUATIONS_HPP
