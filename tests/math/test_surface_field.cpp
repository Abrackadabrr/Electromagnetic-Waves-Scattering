//
// Created by evgen on 02.10.2024.
//

#include "meshes/plate/PlateGrid.hpp"
#include "gtest/gtest.h"

#include "math/fields/SurfaceField.hpp"

using namespace EMW;
using namespace EMW::Types;

class SURFACE_FILEDS_TESTS : public testing::Test {
};

Vector3c z_unit(const Vector3d&) {
    return Vector3c{{0., 0}, {0., 0}, {1., 0}};
}

Vector3c y_unit(const Vector3d&) {
    return Vector3c{{0., 0}, {1., 0}, {0., 0}};
}

TEST_F(SURFACE_FILEDS_TESTS, CROSS_WITH_NORMAL) {
    int N1 = 43;
    int N2 = 31;

    const auto surfaceMesh = Mesh::SurfaceMesh{EMW::Examples::Plate::generateRectangularMesh(N1, N2, 1./(N1 -1), 1./(N2 -1))};

    const auto field_z = Math::SurfaceField(surfaceMesh, [](Vector3d x)->Vector3c{return z_unit(x);});
    const auto field_y = Math::SurfaceField(surfaceMesh, [](Vector3d x)->Vector3c{return y_unit(x);});
    const auto field_diff = field_z.crossWithNormalField() - field_y;
    ASSERT_NEAR(field_diff.supNorm(), 0, 1e-14);
}