//
// Created by evgen on 08.02.24.
//

#ifndef ELECTROMAGNETIC_WAVES_SCATTERING_VOLUMEMESH_HPP
#define ELECTROMAGNETIC_WAVES_SCATTERING_VOLUMEMESH_HPP

#include "types/Types.hpp"
#include "mesh/MeshTypes.hpp"
#include "mesh/Mesh.hpp"

namespace EMW::Mesh {
    class VolumeMesh {
        Containers::vector<Node> nodes_;
        const Mesh::SurfaceMesh &surfaceMesh_;
        std::string name_ = "volume_mesh";

        Types::Vector3c
        getFirstPartIntegral(const Point &point, Types::complex_d k);

        Types::Vector3c
        getZeroPartIntegral(const Point &point, Types::complex_d k);

        Types::Vector3c
        getRIntegral(const Point &point, Types::complex_d k);

    public:
        VolumeMesh() = delete;

        VolumeMesh(const Mesh::SurfaceMesh &surfaceMesh, const Containers::vector<Node> &nodes) :
                surfaceMesh_(surfaceMesh), nodes_(nodes) {};

        void calculateAll(const Types::Vector3d &polarization, const Types::Vector3d &k_vec,
                          Types::complex_d k);

        void calculateESS(const Types::Vector3d &polarization, const Types::Vector3d &k_vec,
                          Types::complex_d k);

        [[nodiscard]] const Containers::vector<Node> &getNodes() const noexcept {
            return nodes_;
        }

        [[nodiscard]] const Mesh::SurfaceMesh &getSurface() const noexcept {
            return surfaceMesh_;
        }

        [[nodiscard]] std::string getName() const noexcept {
            return name_;
        }


    };
}

#endif //ELECTROMAGNETIC_WAVES_SCATTERING_VOLUMEMESH_HPP
