//
// Created by evgen on 24.01.24.
//

#ifndef ELECTROMAGNETIC_WAVES_SCATTERING_TYPES_H
#define ELECTROMAGNETIC_WAVES_SCATTERING_TYPES_H

#include "third_party/eigen/Eigen/Core"
#include "third_party/eigen/Eigen/Geometry"
#include <array>
#include <vector>
#include <set>

namespace EMW::Types {
    using scalar = double;
    using index = std::size_t;
    using integer = long int;
    using complex_d = std::complex<scalar>;

    template<typename t>
    using Vector3 = Eigen::Matrix<t, 3, 1>;

    using Vector3d = Vector3<scalar>;
    using Vector2d = Eigen::Matrix<scalar, 2, 1>;
    using Vector3c = Vector3<complex_d>;
    template<typename t> using VectorX = Eigen::Matrix<t, Eigen::Dynamic, 1>;
    using VectorXc = Eigen::VectorXcd;
    using VectorXd = Eigen::VectorXd;

    template<typename t> using MatrixX = Eigen::Matrix<t, Eigen::Dynamic, Eigen::Dynamic>;
    using MatrixXc = Eigen::MatrixX<Types::complex_d>;
    using MatrixXd = Eigen::MatrixX<Types::scalar>;
    using Matrix3d = Eigen::Matrix3<scalar>;
    using Matrix3c = Eigen::Matrix3<complex_d>;
    template<typename t> using DiagonalMatrixX = Eigen::DiagonalMatrix<t, Eigen::Dynamic>;
    using DiagonalMatrixXd = DiagonalMatrixX<Types::scalar>;
    using DiagonalMatrixXc = DiagonalMatrixX<Types::complex_d>;
};

namespace EMW::Containers {
    template<Types::index n>
    using array_d = std::array<Types::scalar, n>;
    template<typename type, Types::index n>
    using array = std::array<type, n>;
    using vector_d = std::vector<Types::scalar>;
    template<typename type>
    using vector = std::vector<type>;
    template<typename ... types>
    using tuple = std::tuple<types...>;
    template<typename Type>
    using set = std::set<Type>;
};

#endif //ELECTROMAGNETIC_WAVES_SCATTERING_TYPES_H
