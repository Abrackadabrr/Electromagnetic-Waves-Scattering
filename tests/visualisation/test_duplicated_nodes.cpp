//
// Created by evgen on 24.12.2024.
//

#include "VTKFunctions.hpp"

#include "gtest/gtest.h"
#include "types/Types.hpp"
#include "mesh/SurfaceMesh.hpp"
#include "mesh/Utils.hpp"

#include "mesh/Parser.hpp"
#include <ranges>

/**
 * Хочу проверить тезис о том, что в сетку можно загрузить совпадающие узлы и она корректно отобразится
 */
TEST(MESH_TESTING_WITH_DUPLICATED_NODES, SIMPLE_TEST) {
    // Генерация прямоугольной сетки с совпадающими узлами
    std::vector<EMW::Mesh::point_t> meshgrid;
    std::vector<EMW::Containers::array<EMW::Types::index, 4>> cells;
    const EMW::Types::integer N = 10;
    const EMW::Types::scalar h = 1. / (N - 1);
    meshgrid.reserve(4 * N * N);
    cells.reserve(4 * N * N);

    // Заполнение сетки через одновременное создание узлов и ячеек
    for (EMW::Types::scalar i = 0; i != N - 1; i += 1) {
        for (EMW::Types::scalar j = 0; j != N - 1; j += 1) {
            EMW::Types::index current = meshgrid.size();
            // нужно сгенерировать 4 точки
            const auto ij_point = EMW::Types::Vector3d{i, j, 0};
            const auto ij_p_point = EMW::Types::Vector3d{i, j + 1, 0};
            const auto i_p_j_p_point = EMW::Types::Vector3d{i + 1, j + 1, 0};
            const auto i_p_j_point = EMW::Types::Vector3d{i + 1, j, 0};

            cells.push_back(EMW::Containers::array<EMW::Types::index, 4>{current, current + 1, current + 2, current + 3});

            meshgrid.push_back(ij_point);
            meshgrid.push_back(ij_p_point);
            meshgrid.push_back(i_p_j_p_point);
            meshgrid.push_back(i_p_j_point);
        }
    }

    // Создать объект сетки
    EMW::Mesh::SurfaceMesh mesh{meshgrid, cells};

    EMW::Math::SurfaceVectorField normal = EMW::Math::SurfaceVectorField::NormalField(mesh);

    VTK::united_snapshot({normal}, {},mesh, "/home/evgen/Education/MasterDegree/thesis/"
                                "Electromagnetic-Waves-Scattering/tests/visualisation/test.vtk");
}

/**
 * Хочу проверить тезис о том, что в сетку можно загрузить совпадающие узлы и она корректно отобразится
 */
TEST(MESH_TESTING_WITH_DUPLICATED_NODES, RUPOR_MESH_TEST) {
    const std::string nodesFile = "/home/evgen/Education/MasterDegree/thesis/Electromagnetic-Waves-Scattering/meshes/"
                                  "rupor/15200_nodes.csv";
    const std::string cellsFile = "/home/evgen/Education/MasterDegree/thesis/Electromagnetic-Waves-Scattering/meshes/"
                                  "rupor/3800_cells.csv";
    const EMW::Types::index nNodes = 15200;
    const EMW::Types::index nCells = 3800;

    const auto parser_out = EMW::Parser::parseMesh(nodesFile, cellsFile, nNodes, nCells);
    auto surfaceMesh = EMW::Mesh::SurfaceMesh{parser_out.first, parser_out.second};
    // surfaceMesh = surfaceMesh.getSubmesh(EMW::Mesh::IndexedCell::Tag::WAVEGUIDE_CROSS_SECTION);

    EMW::Math::SurfaceScalarField tag_field = EMW::Math::SurfaceScalarField::tagField(surfaceMesh);
    tag_field.setName("tag");
    EMW::Math::SurfaceScalarField number_field = EMW::Math::SurfaceScalarField::sequenceNumberField(surfaceMesh);
    number_field.setName("number");
    VTK::united_snapshot({}, {tag_field, number_field}, surfaceMesh, "/home/evgen/Education/MasterDegree/thesis/"
                               "Electromagnetic-Waves-Scattering/tests/visualisation/");
}
