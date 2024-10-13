//
// Created by evgen on 04.10.2024.
//

#ifndef SURFACESCALARFIELD_H
#define SURFACESCALARFIELD_H

#include "SurfaceFieldBase.h"

namespace EMW::Math {
class SurfaceScalarField : public SurfaceFieldBase<Types::complex_d> {
    using Base = SurfaceFieldBase<Types::complex_d>;
    using manifold_t = Base::manifold_t;
    using field_t = Base::field_t;

  public:
    SurfaceScalarField(const manifold_t &manifold_ref) : Base(manifold_){};
    SurfaceScalarField(const manifold_t &man_ref, const Containers::vector<field_t> &field_data)
        : Base(man_ref, field_data){}
    SurfaceScalarField(const manifold_t &manifold_ref, const std::function<field_t(const Mesh::point_t &)> &function)
        : Base(manifold_ref, function){}
    SurfaceScalarField(const manifold_t &manifold_ref,
                       const std::function<field_t(const Mesh::IndexedCell &)> &function)
        : Base(manifold_ref, function){}

    static constexpr int getNumberOfComponents() {return 1;}
};
} // namespace EMW::Math

#endif // SURFACESCALARFIELD_H
