//
// Created by evgen on 02.08.24.
//

#ifndef ELECTROMAGNETIC_WAVES_SCATTERING_SURFACEFIELD_HPP
#define ELECTROMAGNETIC_WAVES_SCATTERING_SURFACEFIELD_HPP

#include "types/Types.hpp"
#include "mesh/SurfaceMesh.hpp"
#include "mesh/MeshTypes.hpp"

namespace EMW::Math {
    class SurfaceField {
        using manifold_t = Mesh::SurfaceMesh;
        using field_t = Types::Vector3c;

        const manifold_t &manifold_;
        Containers::vector<field_t> field_data_;
        bool initialized;
        std::string name_ = "default_field_name";

    public:
        /**
         * Консткруктор нулевого поля на поверхности
         * @param manifold_ref
         */
        explicit SurfaceField(const manifold_t &manifold_ref) : manifold_(manifold_ref),
                                                                field_data_(Containers::vector<field_t>{}),
                                                                initialized(false) {}

        /**
         * Конструктор поля по форматированным данным
         * @param man_ref
         * @param field_data
         */
        SurfaceField(const manifold_t &man_ref, const Containers::vector<field_t> &field_data) : manifold_(man_ref),
                                                                                                 field_data_(
                                                                                                         field_data),
                                                                                                 initialized(true) {}

        /**
         * Конструктор поля по аналитическому заданию в точках поверхности
         * @param man_ref
         * @param function
         */
        SurfaceField(const manifold_t &man_ref, const std::function<field_t(const Mesh::point_t &)> &function);

        // --- Methods --- //
        [[nodiscard]] SurfaceField surfaceProjection() const;

        [[nodiscard]] Types::VectorXc asSLAERHS() const;

        // --- Getters & Setters --- //
        [[nodiscard]] const manifold_t &getManifold() const { return manifold_; };

        [[nodiscard]] const Containers::vector<field_t> &getField() const { return field_data_; };

        [[nodiscard]] std::string getName() const {return name_;};

        void setName(const std::string &name) { name_ = name; };

        // --- Fabric --- //
        static SurfaceField ZeroField(const manifold_t &manifold);

        static SurfaceField TangentField(const manifold_t &manifold, const Types::VectorXc &fieldProjections);

        static SurfaceField NormalField(const manifold_t &manifold);
    };
}

#endif //ELECTROMAGNETIC_WAVES_SCATTERING_SURFACEFIELD_HPP
