//
// Created by evgen on 18.05.2026.
//

#include <gtest/gtest.h>

#include <cmath>

#include "EMW/types/Types.hpp"

#include "EMW/math/integration/simplicial/SimplicialRange.hpp"

using namespace EMW;

template <typename point_t>
using TriangleRange = Math::Integration::Numerical::Simplicial::QuadratureUtils::TriangleRange<point_t>;

namespace {

Types::scalar triangle_area(const Types::point_t &a, const Types::point_t &b, const Types::point_t &c) {
    const Types::scalar ux = b.x() - a.x();
    const Types::scalar uy = b.y() - a.y();

    const Types::scalar vx = c.x() - a.x();
    const Types::scalar vy = c.y() - a.y();

    return std::abs(ux * vy - uy * vx) * Types::scalar(0.5);
}

bool inside_canonical_simplex(const Types::point_t &p) {
    constexpr Types::scalar eps = Types::scalar(1e-12);

    return p.x() >= -eps && p.y() >= -eps && p.x() + p.y() <= Types::scalar(1) + eps;
}

} // namespace

TEST(TriangleRange, LevelZero) {
    const Containers::array root{Types::point_t{0.0, 0.0, 0}, Types::point_t{1.0, 0.0, 0}, Types::point_t{0.0, 1.0, 0}};

    size_t level = 0;
    const std::size_t N = std::size_t{1} << level;
    const std::size_t expected_count = N * N;

    const Types::scalar expected_area = 1. / (2 * static_cast<Types::scalar>(expected_count));

    std::size_t count = 0;
    auto total_area = static_cast<Types::scalar>(0);

    for (const auto small_tri : TriangleRange(root, level)) {
        ASSERT_EQ(root[2], small_tri.a);
        ASSERT_EQ(root[0], small_tri.b);
        ASSERT_EQ(root[1], small_tri.c);

        const Types::scalar area = triangle_area(small_tri.a, small_tri.b, small_tri.c);

        EXPECT_NEAR(area, 0.5, 1e-12);

        total_area += area;
        ++count;
    }

    EXPECT_EQ(count, expected_count);

    EXPECT_NEAR(total_area, 0.5, 1e-12);
}

TEST(TriangleRange, CanonicalSimplexAreaAndCount) {
    const Containers::array root{Types::point_t{0.0, 0.0, 0}, Types::point_t{1.0, 0.0, 0}, Types::point_t{0.0, 1.0, 0}};

    for (std::size_t level = 0; level <= 8; ++level) {
        const std::size_t N = std::size_t{1} << level;
        const std::size_t expected_count = N * N;

        const Types::scalar expected_small_area = 1. / (2 * static_cast<Types::scalar>(expected_count));

        std::size_t count = 0;
        auto total_area = static_cast<Types::scalar>(0);

        for (const auto small_tri : TriangleRange(root, level)) {
            EXPECT_TRUE(inside_canonical_simplex(small_tri.a));
            EXPECT_TRUE(inside_canonical_simplex(small_tri.b));
            EXPECT_TRUE(inside_canonical_simplex(small_tri.c));

            const Types::scalar area = triangle_area(small_tri.a, small_tri.b, small_tri.c);

            EXPECT_NEAR(area, expected_small_area, 1e-12);

            total_area += area;
            ++count;
        }

        EXPECT_EQ(count, expected_count);

        EXPECT_NEAR(total_area, 0.5, 1e-12);
    }
}

Types::scalar signed_area_2d(const Types::point_t &a, const Types::point_t &b, const Types::point_t &c) {
    const Types::scalar ux = b.x() - a.x();
    const Types::scalar uy = b.y() - a.y();

    const Types::scalar vx = c.x() - a.x();
    const Types::scalar vy = c.y() - a.y();

    return Types::scalar(0.5) * (ux * vy - uy * vx);
}

TEST(TriangleRange, SmallTrianglesHaveConsistentOrientation) {
    const Containers::array root{Types::point_t{0.0, 0.0, 0.}, Types::point_t{1.0, 0.0, 0},
                                 Types::point_t{0.0, 1.0, 0}};

    const Types::scalar root_signed_area = signed_area_2d(root[0], root[1], root[2]);

    ASSERT_GT(root_signed_area, Types::scalar(0));

    for (std::size_t level = 0; level <= 8; ++level) {
        for (const auto small_tri : TriangleRange(root, level)) {
            const Types::scalar small_signed_area = signed_area_2d(small_tri.a, small_tri.b, small_tri.c);

            EXPECT_GT(small_signed_area, Types::scalar(0)) << "level = " << level;
        }
    }
}

TEST(TriangleRange, SmallTrianglesPreserveNegativeOrientation) {
    const Containers::array root{Types::point_t{0.0, 0.0, 0}, Types::point_t{0.0, 1.0, 0}, Types::point_t{1.0, 0.0, 0}};

    const Types::scalar root_signed_area = signed_area_2d(root[0], root[1], root[2]);

    ASSERT_LT(root_signed_area, Types::scalar(0));

    for (std::size_t level = 0; level <= 8; ++level) {
        for (const auto small_tri : TriangleRange(root, level)) {
            const Types::scalar small_signed_area = signed_area_2d(small_tri.a, small_tri.b, small_tri.c);

            EXPECT_LT(small_signed_area, Types::scalar(0)) << "level = " << level;
        }
    }
}