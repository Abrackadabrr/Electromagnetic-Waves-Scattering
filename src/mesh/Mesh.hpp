//
// Created by evgen on 25.01.24.
//

#ifndef ELECTROMAGNETIC_WAVES_SCATTERING_MESH_HPP
#define ELECTROMAGNETIC_WAVES_SCATTERING_MESH_HPP

#include <cmath>
#include <utility>
#include "types/Types.hpp"
#include "MeshTypes.hpp"

namespace EMW::Mesh {
    /**
     * Класс расчетной поверхностной сетки
     */
    class SurfaceMesh {
        Containers::vector<Point> nodes_;
        Containers::vector<IndexedCell> cells_;
        std::string name = "default_mesh_name";
        bool jFilled_;
    public:
        SurfaceMesh(Containers::vector<Point> nodes, Containers::vector<Containers::array<Types::index, 4>> cells);

        void fillJ(const Types::VectorXc &j);

        [[nodiscard]] constexpr const Containers::vector<IndexedCell> &getCells() const { return cells_; }

        [[nodiscard]] constexpr const Containers::vector<Point> &getNodes() const { return nodes_; }

        [[nodiscard]] std::string getName() const { return name; }

        void setName(const std::string &n) {name = n;}

        [[nodiscard]] bool jFilled() const {return jFilled_;}

        void basisHack();

        // deprecated
        //        SurfaceMesh(Containers::vector<Point> nodes,
        //                    Containers::vector<Containers::array<Types::index, 4>> cells,
        //                    Containers::vector<Types::Vector3c> E_field, Containers::vector<Types::Vector3c> H_field);
    };
}
#endif //ELECTROMAGNETIC_WAVES_SCATTERING_MESH_HPP
