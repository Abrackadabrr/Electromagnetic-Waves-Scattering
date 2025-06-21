//
// Created by evgen on 18.01.2025.
//

#ifndef GENERALIZEDEQUATIONS_HPP
#define GENERALIZEDEQUATIONS_HPP

#include "mesh/SurfaceMesh.hpp"

#include "slae_generation/MatrixGeneration.hpp"

#include "types/Types.hpp"

#include "math/MathConstants.hpp"
#include "math/fields/SurfaceVectorField.hpp"

namespace WaveGuideWithActiveSection {

using namespace EMW;

/**
 * Вычисляем матрицу для "самодействия"
 */
inline EMW::Types::MatrixXc diagonal(const Mesh::SurfaceMesh &mesh_all, Types::scalar a, Types::complex_d k) {
    // вытаскиваем сетку для магнитных токов
    const auto mesh_zero = mesh_all.getSubmesh(Mesh::IndexedCell::Tag::WAVEGUIDE_CROSS_SECTION);
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
    // TODO: переписать код, чтобы матрицы генерировалась in-place
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
    const Mesh::SurfaceMesh &mesh_to_integrate_zero =
        mest_to_integrate.getSubmesh(Mesh::IndexedCell::Tag::WAVEGUIDE_CROSS_SECTION);
    const Mesh::SurfaceMesh &mesh_collocation_all = mest_with_collocation_points;
    const Mesh::SurfaceMesh &mesh_collocation_zero =
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

/**
 * Не прямо-таки элемент в матрице, а его "брат-близнец" в верхней части соотвествующего блока в матрице
 */
inline Types::complex_d
element_of_submatrix(Types::index i, Types::index j, const Mesh::SurfaceMesh &mesh_to_integrate_all,
                     const Mesh::SurfaceMesh &mesh_to_integrate_zero, const Mesh::SurfaceMesh &mesh_collocation_all,
                     const Mesh::SurfaceMesh &mesh_collocation_zero, const Types::scalar a, const Types::complex_d k) {
    // Надо понять какие размеры будем иметь итоговая матрица
    const Types::index N_cols = mesh_to_integrate_all.getCells().size();
    const Types::index K_cols = mesh_to_integrate_zero.getCells().size();
    const Types::index cols = 2 * (N_cols + K_cols);
    const Types::index N_rows = mesh_collocation_all.getCells().size();
    const Types::index K_rows = mesh_collocation_zero.getCells().size();
    const Types::index rows = 2 * (N_rows + K_rows);
    // Надо понять какие константы будут фигурировать в расчете
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
    // Расчет значения элемента непосредственно

    const auto local_array_index = [](Types::index i, Types::index j, Types::index rows, Types::index cols) {
        // есть двоичное число i / rows, j / cols
        // десятичное будет в виде (i / rows * 2) + (j / cols)
        return (i / rows) * 2 + (j / cols);
    };

    if (i < 2 * N_rows) {
        if (j < 2 * N_cols) {
            return with_e * Matrix::DiscreteK::getMatrixCoefsInArray(mesh_to_integrate_all.getCells()[j % N_cols],
                                                              mesh_collocation_all.getCells()[i % N_rows], k)
                                                                [local_array_index(i, j, N_rows, N_cols)];
            // а давайте попробуем приближать матрицу, используя только один набор коэффициентов
            // таким образом, вместо настоящего коэффициента,
            // мы всегда будем возвращать соответствующий элемент в блоке (0, N_rows) x (0, N_cols)
        } else {
            Types::index local_j = j - 2 * N_cols;
            return -Matrix::DiscreteR::getMatrixCoefsInArray(mesh_to_integrate_zero.getCells()[local_j % K_cols],
                                                      mesh_collocation_all.getCells()[i % N_rows], k)
                                                        [local_array_index(i, local_j, N_rows, K_cols)];
        }
    } else {
        Types::index local_i = i - 2 * N_rows;
        if (j < 2 * N_cols) {
            return z * Matrix::DiscreteR::getMatrixCoefsInArray(mesh_to_integrate_all.getCells()[j % N_cols],
                                                         mesh_collocation_zero.getCells()[local_i % K_rows], k)
                                                         [local_array_index(local_i, j, K_rows, N_cols)];
        } else {
            Types::index local_j = j - 2 * N_cols;
            return z * with_mu *
                   Matrix::DiscreteK::getMatrixCoefsInArray(mesh_to_integrate_zero.getCells()[local_j % K_cols],
                                                     mesh_collocation_zero.getCells()[local_i % K_rows], k)
                                                         [local_array_index(local_i, local_j, K_rows, K_cols)];
        }
    }
}

inline Types::VectorXc getRhsBase(const Mesh::SurfaceMesh &mesh_all, Types::scalar a, Types::complex_d k) {

    // Собираем сабсетку с активным сечением
    const auto &mesh_zero = mesh_all.getSubmesh(Mesh::IndexedCell::Tag::WAVEGUIDE_CROSS_SECTION);
    // Собираем функцию поля на активном сечении волновода
    const Types::complex_d beta = std::sqrt(k * k - (EMW::Math::Constants::PI_square<Types::scalar>() / (a * a)));
    const auto get_e_h10_mode = [&k, &a, &beta](const Mesh::point_t &x) {
        const auto pi = Math::Constants::PI<Types::scalar>();
        const Types::complex_d mult = Math::Constants::i * (a / pi) * k / Math::Constants::e_0_c;
        const Types::complex_d exp = std::exp(Math::Constants::i * beta * (x.z() + 0.17));
        const Types::scalar cos = std::cos(pi * x.x() / a);
        return Types::Vector3c{Types::complex_d{0, 0}, mult * cos * exp, Types::complex_d{0, 0}};
    };

    const Types::index N = mesh_all.getCells().size();
    const Types::index K = mesh_zero.getCells().size();

    const Math::SurfaceVectorField direct_e_field{mesh_zero, get_e_h10_mode};
    const Types::VectorXc non_zer0_rhs = -2 * direct_e_field.normalCrossField().asVector();

    Types::VectorXc res = Types::VectorXc::Zero(2 * N + 2 * K);
    res.block(2 * N, 0, 2 * K, 1) = non_zer0_rhs;
    return res;
}

template<Types::index N>
Types::VectorXc getRhs(std::array<Types::complex_d, N> phases,
    const Mesh::SurfaceMesh &mesh_all, Types::scalar a, Types::complex_d k) {
    const auto base_rhs = getRhsBase(mesh_all, a, k);

    Types::VectorXc res{N * base_rhs.size()};

    for (int i = 0; i < N; ++i) {
        res.block(i * base_rhs.rows(), 0, base_rhs.rows(), 1) = phases[i] * base_rhs;
    }
    return res;
}

using diagonal_t = decltype(diagonal);
using submatrix_t = decltype(submatrix);
}

#endif //GENERALIZEDEQUATIONS_HPP
