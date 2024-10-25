//
// Created by evgen on 25.01.24.
//

#ifndef ELECTROMAGNETIC_WAVES_SCATTERING_OBJPARSER_H
#define ELECTROMAGNETIC_WAVES_SCATTERING_OBJPARSER_H

#include "MeshTypes.hpp"
#include "SurfaceMesh.hpp"
#include "../types/Types.hpp"
#include "third_party/csv/csv.h"
#include <string>

namespace EMW::Parser {
    Mesh::SurfaceMesh parseMesh(const std::string & csvNodes, const std::string & csvCells, int nNodes, int nCells);

}

#endif //ELECTROMAGNETIC_WAVES_SCATTERING_OBJPARSER_H
