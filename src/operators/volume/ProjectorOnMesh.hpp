//
// Created by evgen on 27.01.2026.
//

#ifndef PROJECTORONMESH_HPP
#define PROJECTORONMESH_HPP

#include "math/"

namespace EMW::Operators::Volume {

class ProjectorOnMesh {
    using mesh_t = Mesh::VolumeMesh::CubeMesh
    const mesh_t &mesh_;

    public:
    ProjectorOnMesh(const mesh_t &mesh): mesh_(mesh) {}
};

}
#endif //PROJECTORONMESH_HPP

