//
// Created by evgen on 24.01.24.
//

#ifndef ELECTROMAGNETIC_WAVES_SCATTERING_TYPES_H
#define ELECTROMAGNETIC_WAVES_SCATTERING_TYPES_H

#include "third_party/eigen/Eigen/Core"
#include "third_party/eigen/Eigen/Geometry"
#include <array>
#include <vector>

namespace EMW::Types {
    using scalar = double;
    using index = std::size_t;
    using integer = std::int_fast32_t;
    using complex_d = std::complex<scalar>;

    template<typename t>
    using Vector3 = Eigen::Matrix<t, 3, 1>;

    using Vector3d = Vector3<scalar>;
    using Vector2d = Eigen::Matrix<scalar, 2, 1>;
    using Vector3c = Vector3<complex_d>;
    using VectorXc = Eigen::VectorXcd;
    using VectorXd = Eigen::VectorXd;

    using MatrixXc = Eigen::MatrixXcd;
    using MatrixXd = Eigen::MatrixXd;
    using Matrix3d = Eigen::Matrix3<scalar>;
    using Matrix3c = Eigen::Matrix3<complex_d>;

}

namespace EMW::Containers {
    template<Types::index n>
    using array_d = std::array<Types::scalar, n>;
    template<typename type, Types::index n>
    using array = std::array<type, n>;
    using vector_d = std::vector<Types::scalar>;
    template<typename type>
    using vector = std::vector<type>;
}

#endif //ELECTROMAGNETIC_WAVES_SCATTERING_TYPES_H
