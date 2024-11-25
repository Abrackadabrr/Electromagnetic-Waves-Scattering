//
// Created by evgen on 25.01.24.
//

#ifndef ELECTROMAGNETIC_WAVES_SCATTERING_SURFACEMESH_HPP
#define ELECTROMAGNETIC_WAVES_SCATTERING_SURFACEMESH_HPP

#include <cmath>
#include <utility>
#include "types/Types.hpp"
#include "MeshTypes.hpp"

namespace EMW::Mesh {
    /**
     * Класс расчетной поверхностной сетки, дискретный аналог поверхности
     */
    class SurfaceMesh {
        Containers::vector<point_t> nodes_;
        Containers::vector<IndexedCell> cells_;
        std::string name = "default_mesh_name";

    public:

        SurfaceMesh(Containers::vector<point_t> nodes, Containers::vector<Containers::array<Types::index, 4>> cells);

        SurfaceMesh(Containers::vector<point_t> nodes, Containers::vector<Containers::array<Types::index, 4>> cells,
        const std::function<point_t(const Containers::array<Types::index, 4> &,
                                                 const Containers::vector<point_t> &)> &getPoint);

        void customLocalBasis(const std::function<std::array<Types::Vector3d, 3>(const Mesh::IndexedCell&)>& func);

        void customCollocationPpoints(const std::function<Types::Vector3d(const Mesh::IndexedCell &)>& func);

        [[nodiscard]] constexpr const Containers::vector<IndexedCell> &getCells() const { return cells_; }

        [[nodiscard]] constexpr const Containers::vector<point_t> &getNodes() const { return nodes_; }

        [[nodiscard]] std::string getName() const { return name; }

        void setName(const std::string &n) {name = n;}
    };
}
#endif //ELECTROMAGNETIC_WAVES_SCATTERING_SURFACEMESH_HPP
