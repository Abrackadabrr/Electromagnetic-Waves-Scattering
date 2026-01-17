//
// Created by evgen on 17.01.2026.
//

#ifndef CUBEMESH_HPP
#define CUBEMESH_HPP

#include "../MeshBase.hpp"
#include "VolumeCells.hpp"

namespace EMW::Mesh::VolumeMesh {

class CubeMesh : public MeshBase<VolumeCells::IndexedCube> {

};

};

#endif //CUBEMESH_HPP
