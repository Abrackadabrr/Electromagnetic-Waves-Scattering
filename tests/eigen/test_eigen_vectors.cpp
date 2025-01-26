//
// Created by evgen on 04.02.24.
//

#include "gtest/gtest.h"
#include "types/Types.hpp"


using namespace EMW::Types;

template<std::size_t ... indexes>
Vector3c
sum_fold(const std::array<Vector3c, 3> &vec, const std::array<scalar, 3> &w,
         std::index_sequence<indexes...>) {
    return ((vec[indexes] * w[indexes]) + ...);
}

TEST(EigenVectors, CastTest) {
    Vector3c a{std::complex<scalar>{1, 0}, std::complex<scalar>{1, 0}, std::complex<scalar>{0, 1}};
    scalar b = 1;
    const Vector3c c = a * b;
    std::cout << c.norm() << std::endl;
}

TEST(EigenVectors, FoldTest) {
    std::array<Vector3c, 3> vectors{
            Vector3c{std::complex<scalar>{1, 0}, std::complex<scalar>{1, 0}, std::complex<scalar>{0, 1}},
            Vector3c{std::complex<scalar>{1, 0}, std::complex<scalar>{1, 0}, std::complex<scalar>{0, 1}},
            Vector3c{std::complex<scalar>{1, 0}, std::complex<scalar>{1, 0}, std::complex<scalar>{0, 1}}};
    std::array<scalar, 3> weights{1, 2, 3};
    const Vector3c sum = sum_fold(vectors, weights, std::make_index_sequence<3>{});
//    const Vector3c sum_for = sum(vectors, weights, std::make_index_sequence<3>{});
}

TEST(EigenVectors, StingTranspose) {
    const Vector3c a = {complex_d{1, 0}, complex_d{1, 0}, complex_d{1, 0}};
    const Vector3d b = {1, 1, 1};
    const auto bT = b.transpose();
    const complex_d bTa = bT * a;
    const Eigen::Matrix3<complex_d> abT = a * bT;
}

TEST(EigenVectors, DotProduction) {
    const Vector3c a = {complex_d{0, 1}, complex_d{0, 1}, complex_d{0, 1}};
    std::cout << a << std::endl;
    std::cout << a.sum() << std::endl;
    const Vector3d b = {1, 1, 1};
    std::cout << a.dot(b) << std::endl;
    std::cout << b.dot(a) << std::endl;
}


TEST(EigenVectors, INverseEigenVector) {
    const Vector3c a = {complex_d{10, 0}, complex_d{10, 0}, complex_d{10, 0}};
    const Vector3c a_inv = a.cwiseInverse();
    std::cout << a_inv << std::endl;
    std::cout << a_inv.cwiseProduct(a) << std::endl;
}
