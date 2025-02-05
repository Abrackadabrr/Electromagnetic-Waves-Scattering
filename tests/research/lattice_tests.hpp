//
// Created by evgen on 04.02.2025.
//

#ifndef LATTICE_TESTS_HPP
#define LATTICE_TESTS_HPP

#include <gtest/gtest.h>

#include "mesh/SurfaceMesh.hpp"

#include "experiment/PhysicalCondition.hpp"

#include "mesh/Parser.hpp"

using namespace EMW;

class AssemblingTests : public ::testing::Test {
protected:
    EMW::Mesh::SurfaceMesh mesh_base;
    const Types::scalar a = 0.07;
    const Types::scalar freq = Math::Constants::c / 1e8;
    const Types::complex_d k{Physics::get_k_on_frquency(freq), 0};
public:
    void SetUp() override {
        const std::string nodesFile = "/home/evgen/Education/MasterDegree/thesis/Electromagnetic-Waves-Scattering/meshes/"
                                          "lattice/8000_nodes.csv";
        const std::string cellsFile = "/home/evgen/Education/MasterDegree/thesis/Electromagnetic-Waves-Scattering/meshes/"
                                      "lattice/2000_cells.csv";
        const EMW::Types::index nNodes = 8000;
        const EMW::Types::index nCells = 2000;

        // собираем сетки
        const auto parser_out = EMW::Parser::parseMesh(nodesFile, cellsFile, nNodes, nCells);
        mesh_base = EMW::Mesh::SurfaceMesh{parser_out.first, parser_out.second};
    };
};

#endif //LATTICE_TESTS_HPP
