//
// Created by evgen on 14.01.2025.
//

#include "meshes/plate/PlateGrid.hpp"
#include "gtest/gtest.h"
#include "math/Productions.hpp"
#include "math/fields/SurfaceVectorField.hpp"

using namespace EMW;
using namespace EMW::Types;

class PRODUCTIONS_TESTS : public testing::Test {
};

Vector3c crossCheck(const Types::Vector3c& a, const Types::Vector3c& b) {
    return Vector3c{ a.y() * b.z() - a.z() * b.y(),
                     a.z() * b.x() - a.x() * b.z(),
                     a.x() * b.y() - a.y() * b.x() };
}

Vector3c z_unit() {
    return Vector3c{std::complex<scalar>{1., 0}, std::complex<scalar>{0., 0}, std::complex<scalar>{0., 0}};
}

Vector3c y_unit() {
    return Vector3c{std::complex<scalar>{0., 0}, std::complex<scalar>{1., 0}, std::complex<scalar>{0., 0}};
}

TEST_F(PRODUCTIONS_TESTS, EIGEN_COMPARISON) {
    const Types::Vector3c z{complex_d{1, 2}, complex_d{2, 3}, complex_d{4, 5}};
    const Types::Vector3c y{complex_d{-14, 7}, complex_d{-2, 56}, complex_d{2, 13}};
    const Types::Vector3c eigen_result = z.cross(y);
    const Types::Vector3c check_result = crossCheck(z, y);
    const Types::Vector3c my_res = Math::cross(z, y);
    std::cout << eigen_result << std::endl;
    std::cout << check_result << std::endl;
    std::cout << my_res << std::endl;
    ASSERT_EQ(check_result, my_res);
}
#if 0
TEST_F(PRODUCTIONS_TESTS, EIGEN_COMPARISON_CD) {
    const Types::Vector3c z{complex_d{1, 2}, complex_d{2, 3}, complex_d{4, 5}};
    const Types::Vector3d y{1,1,1};

    const Types::Vector3c eigen_result = z.cross(y);
    const Types::Vector3c my_res = Math::crossCD(z, y);
    std::cout << eigen_result << std::endl;
    std::cout << my_res << std::endl;
    ASSERT_EQ(eigen_result, my_res);
}

TEST_F(PRODUCTIONS_TESTS, EIGEN_COMPARISON_CC) {
    const Types::Vector3c z{complex_d{1, 2}, complex_d{2, 3}, complex_d{4, 5}};
    const Types::Vector3c y{complex_d{-14, 7}, complex_d{-2, 56}, complex_d{2, 13}};

    const Types::Vector3c eigen_result = z.cross(y);
    const Types::Vector3c my_res = Math::crossCC(z, y);
    std::cout << eigen_result << std::endl;
    std::cout << my_res << std::endl;
    ASSERT_EQ(eigen_result, my_res);
}
#endif