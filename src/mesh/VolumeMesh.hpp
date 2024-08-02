//
// Created by evgen on 08.02.24.
//

#ifndef ELECTROMAGNETIC_WAVES_SCATTERING_VOLUMEMESH_HPP
#define ELECTROMAGNETIC_WAVES_SCATTERING_VOLUMEMESH_HPP

#include "types/Types.hpp"
#include "mesh/MeshTypes.hpp"
#include "math/SurfaceField.hpp"

namespace EMW::Mesh {
    class VolumeMesh {
        Containers::vector<Node> nodes_;
        std::string name_ = "default_volume_mesh_name";

    public:
        VolumeMesh(const Mesh::SurfaceMesh &surfaceMesh, const Containers::vector<Node> &nodes) :
                nodes_(nodes) {};

        void evaluateOperator(Types::scalar k, const Math::SurfaceField &field);

        void addInitialField(const std::function<Types::Vector3c(Mesh::point_t)> &initial);

        void calculateFullField(Types::scalar k, const Math::SurfaceField &field,
                                const std::function<Types::Vector3c(Mesh::point_t)> &initial);

        [[nodiscard]] const Containers::vector<Node> &getNodes() const noexcept {
            return nodes_;
        }

        [[nodiscard]] std::string getName() const noexcept {
            return name_;
        }

        void setName(const std::string &n) noexcept {
            name_ = n;
        }
    };
}

#endif //ELECTROMAGNETIC_WAVES_SCATTERING_VOLUMEMESH_HPP
