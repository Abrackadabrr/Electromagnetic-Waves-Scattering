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
    using Vector3d = Eigen::Matrix<scalar, 3, 1>;
    using Vector2d = Eigen::Matrix<scalar, 2, 1>;
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
