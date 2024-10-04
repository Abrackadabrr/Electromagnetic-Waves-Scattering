//
// Created by evgen on 02.08.24.
//

#ifndef ELECTROMAGNETIC_WAVES_SCATTERING_SurfaceVectorField_HPP
#define ELECTROMAGNETIC_WAVES_SCATTERING_SurfaceVectorField_HPP

#include "SurfaceFieldBase.h"
#include "mesh/MeshTypes.hpp"
#include "mesh/SurfaceMesh.hpp"
#include "types/Types.hpp"

namespace EMW::Math {
class SurfaceVectorField : public SurfaceFieldBase<Types::Vector3c> {
    using Base = SurfaceFieldBase<Types::Vector3c>;
    using manifold_t = SurfaceFieldBase<Types::Vector3c>::manifold_t;
    using field_t = SurfaceFieldBase<Types::Vector3c>::field_t;

  public:
    SurfaceVectorField(const manifold_t &manifold_ref) : Base(manifold_){};
    SurfaceVectorField(const manifold_t &man_ref, const Containers::vector<field_t> &field_data)
        : Base(man_ref, field_data){};
    SurfaceVectorField(const manifold_t &manifold_ref, const std::function<field_t(const Mesh::point_t &)> &function)
        : Base(manifold_ref, function){};
    SurfaceVectorField(const manifold_t &manifold_ref,
                       const std::function<field_t(const Mesh::IndexedCell &)> &function)
        : Base(manifold_ref, function){};

    // --- Methods --- //
    [[nodiscard]] SurfaceVectorField surfaceProjection() const;

    [[nodiscard]] SurfaceVectorField crossWithNormalField() const;

    [[nodiscard]] Types::VectorXc asSLAERHS() const;

    [[nodiscard]] Types::scalar supNorm() const;

    // --- Fabric --- //
    static SurfaceVectorField ZeroField(const manifold_t &manifold);

    static SurfaceVectorField TangentField(const manifold_t &manifold, const Types::VectorXc &fieldProjections);

    static SurfaceVectorField NormalField(const manifold_t &manifold);
};

// --- Operators --- //
/*
 * Пока что эти функции не безопасные для использования, их можно успользовать
 * только если поля заданы на одном и том же многообразии
 */
SurfaceVectorField operator-(const SurfaceVectorField &lhs, const SurfaceVectorField &rhs);

SurfaceVectorField operator+(const SurfaceVectorField &lhs, const SurfaceVectorField& rhs);
}

#endif //ELECTROMAGNETIC_WAVES_SCATTERING_SurfaceVectorField_HPP
