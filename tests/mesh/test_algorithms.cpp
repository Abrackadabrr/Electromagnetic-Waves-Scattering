//
// Created by evgen on 19.08.24.
//
#include "gtest/gtest.h"
#include "mesh/Algorithms.hpp"

using namespace EMW;

class ALGORITHMS_TESTS : public testing::Test {
};

TEST_F(ALGORITHMS_TESTS, POINT_IN_TRIANGLE) {
    const Containers::vector<Mesh::point_t> triangle = {
            Mesh::point_t{1, 1, 0},
            Mesh::point_t{4, 1, 0},
            Mesh::point_t{6, 4, 0}
    };

    const Containers::vector<Mesh::point_t> points = {
            Mesh::point_t{2, 2, 0},
            Mesh::point_t{1, 4, 0},
            Mesh::point_t{6, 2, 0},
            Mesh::point_t{3, 3, 0},
            Mesh::point_t{4, 3, 0},
            Mesh::point_t{2, 1, 0},
    };

    const Containers::vector<bool> answers = {false, false, false, false, false, true};

    for (int i = 0; i != points.size(); ++i) {
        const bool res = Mesh::Algorithm::PointInTriangle(points[i], Mesh::Cell{triangle[0], triangle[1], triangle[2],
                                                                                triangle[2]});
        ASSERT_TRUE(Mesh::Algorithm::logical_nor(res, answers[i]) || (res && answers[i])) << i << std::endl;
    }
}
