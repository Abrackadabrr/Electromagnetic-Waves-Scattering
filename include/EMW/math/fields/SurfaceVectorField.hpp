//
// Created by evgen on 02.08.24.
//

#ifndef ELECTROMAGNETIC_WAVES_SCATTERING_SurfaceVectorField_HPP
#define ELECTROMAGNETIC_WAVES_SCATTERING_SurfaceVectorField_HPP

#include "SurfaceFieldBase.h"
#include "SurfaceScalarField.h"
#include "mesh/MeshTypes.hpp"
#include "mesh/SurfaceMesh.hpp"
#include "types/Types.hpp"

namespace EMW::Math {
class SurfaceVectorField : public SurfaceFieldBase<Types::Vector3c> {
    using Base = SurfaceFieldBase<Types::Vector3c>;
    using manifold_t = SurfaceFieldBase<Types::Vector3c>::manifold_t;
    using field_t = SurfaceFieldBase<Types::Vector3c>::field_t;

  public:
    SurfaceVectorField(const manifold_t &manifold_ref) : Base(manifold_ref){};
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

    [[nodiscard]] SurfaceVectorField normalCrossField() const;

    [[nodiscard]] Types::complex_d supNorm() const;

    /**
     * Представление касательной проекции поля в виде вектора из компонент разложения по базису в касательной плоскости
     * Сначала идут все компоненты, соответствующие "первым" базисным векторам, потом всем "вторым" базисным векторам
     *
     * Такое представление в виде вектора очевидно зависит от нумерации ячеек многообразия manifold,
     * ссылка на которое хранится в классе поля
     */
    [[nodiscard]] Types::VectorXc asVector() const;

    [[nodiscard]] SurfaceScalarField<Types::complex_d> fieldNorm(const std::string name) const;

    [[nodiscard]] SurfaceVectorField
    pointwiseMultiplication(const std::function<Types::scalar(const Types::Vector3d &)> &function) const;

    [[nodiscard]] SurfaceVectorField
    pointwiseMultiplication(const std::function<Types::scalar(const Mesh::IndexedCell &)> &function) const;

    void multiply(const std::function<Types::scalar(const Types::Vector3d &)> &function);

    void multiply(const std::function<Types::scalar(const Mesh::IndexedCell &)> &function);

    // --- Fabric --- //
    static SurfaceVectorField ZeroField(const manifold_t &manifold);

    static SurfaceVectorField TangentField(const manifold_t &manifold, const Types::VectorXc &fieldProjections,
                                           std::string name = "tangent_field");

    static SurfaceVectorField NormalField(const manifold_t &manifold);

    // --- Служебные methods --- //
    static constexpr int getNumberOfComponents() { return 3; };
};

// --- Operators --- //
/*
 * Пока что эти функции не безопасные для использования, их можно использовать
 * только если поля заданы на одном и том же многообразии
 */
SurfaceVectorField operator-(const SurfaceVectorField &lhs, const SurfaceVectorField &rhs);

SurfaceVectorField operator+(const SurfaceVectorField &lhs, const SurfaceVectorField &rhs);

SurfaceVectorField operator-(const SurfaceVectorField &lhs);

SurfaceVectorField operator*(const Types::complex_d &lhs, const SurfaceVectorField &rhs);

SurfaceVectorField operator*(const Types::scalar &lhs, const SurfaceVectorField &rhs);
}

#endif //ELECTROMAGNETIC_WAVES_SCATTERING_SurfaceVectorField_HPP
