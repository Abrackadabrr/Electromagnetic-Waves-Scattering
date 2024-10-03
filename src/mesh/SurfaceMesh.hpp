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
        bool jFille = false;

    public:

        SurfaceMesh(Containers::vector<point_t> nodes, Containers::vector<Containers::array<Types::index, 4>> cells);

        void fillJ(const Types::VectorXc &j);

        void customLocalBasis(std::function<std::array<Types::Vector3d, 3>(const Mesh::IndexedCell&)> func);

        [[nodiscard]] bool jFilled() const {return jFille;}

        [[nodiscard]] constexpr const Containers::vector<IndexedCell> &getCells() const { return cells_; }

        [[nodiscard]] constexpr const Containers::vector<point_t> &getNodes() const { return nodes_; }

        [[nodiscard]] std::string getName() const { return name; }

        void setName(const std::string &n) {name = n;}

        // deprecated
        //        SurfaceMesh(Containers::vector<Point> nodes,
        //                    Containers::vector<Containers::array<Types::index, 4>> cells,
        //                    Containers::vector<Types::Vector3c> E_field, Containers::vector<Types::Vector3c> H_field);
    };
}
#endif //ELECTROMAGNETIC_WAVES_SCATTERING_SURFACEMESH_HPP
