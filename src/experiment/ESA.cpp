//
// Created by evgen on 15.02.24.
//

#include "ESA.hpp"
#include "math/MathConstants.hpp"
#include "math/Productions.hpp"
#include "math/fields/SurfaceVectorField.hpp"
#include "math/integration/Quadrature.hpp"
#include "math/integration/gauss_quadrature/GaussLegenderPoints.hpp"
#include "math/integration/newton_cotess/Rectangular.hpp"
#include "mesh/MeshTypes.hpp"
#include "mesh/SurfaceMesh.hpp"
#include "operators/Functions.hpp"
#include "types/Types.hpp"

namespace EMW::ESA {

// Для корректной работы редукции в цикле
#pragma omp declare reduction(+ : Types::Vector3c : omp_out += omp_in) initializer(omp_priv = decltype(omp_orig)::Zero())

Types::Vector3c sigmaOverCell(Types::complex_d k, const Types::Vector3d &tau, const Mesh::IndexedCell &cell,
                              const Types::Vector3c &j_e, const Types::Vector3c &j_m) {
    const auto phi = [&](Types::scalar p, Types::scalar q) -> Types::Vector3c {
        const Mesh::point_t y = cell.parametrization(p, q);
        const Types::scalar mul = cell.multiplier(p, q);
        return Helmholtz::sigmaKernel(k, tau, y, j_e, j_m) * mul;
    };
    return DefiniteIntegrals::integrate<DefiniteIntegrals::GaussLegendre::Quadrature<2, 2>>(phi, {0, 0}, {1, 1});
}

Types::Vector3c sigmaOverCube(Types::complex_d k, const Types::Vector3d &tau,
                              const Types::Vector3d &left_down_corner, Types::scalar dx, Types::scalar dy,
                              Types::scalar dz,
                              const Types::Vector3c &j_e, const Types::Vector3c &j_m, Types::complex_d eps) {
    const auto phi = [&](Types::scalar x, Types::scalar y, Types::scalar z) -> Types::Vector3c {
        return (eps - 1.) * Helmholtz::sigmaKernel(k, tau, {x, y, z}, j_e, j_m, eps);
    };

    return DefiniteIntegrals::integrate<DefiniteIntegrals::NewtonCotess::Quadrature<4, 4, 4>>(
        phi, {left_down_corner.x(), left_down_corner.y(), left_down_corner.z()}, {dx, dy, dz});
}

Types::scalar calculateESA(const Types::Vector3d &tau, Types::complex_d k, const Math::SurfaceVectorField &j_e,
                           const Math::SurfaceVectorField &j_m) {
    const auto &cells_e = j_e.getManifold().getCells();
    const auto &field_e = j_e.getField();
    const auto &cells_m = j_m.getManifold().getCells();
    const auto &field_m = j_m.getField();
    Types::Vector3c result = Types::Vector3c::Zero();
#pragma omp parallel for reduction(+ : result) num_threads(14)
    for (int i = 0; i != field_e.size(); i++) {
        result += sigmaOverCell(k, tau, cells_e[i], field_e[i], Types::Vector3c::Zero());
    }

#pragma omp parallel for reduction(+ : result) num_threads(14)
    for (int i = 0; i != field_m.size(); i++) {
        result += sigmaOverCell(k, tau, cells_m[i], Types::Vector3c::Zero(), field_m[i]);
    }
    return Math::Constants::inverse_4PI<Types::scalar>() * result.squaredNorm();
}

Types::scalar calculateESA(const Types::Vector3d &tau, Types::complex_d k,
                           const Containers::vector<Math::SurfaceVectorField> &j_e_v,
                           const Containers::vector<Math::SurfaceVectorField> &j_m_v) {
    Types::Vector3c result = Types::Vector3c::Zero();

    for (const auto &j_e : j_e_v) {
        const auto &cells_e = j_e.getManifold().getCells();
        const auto &field_e = j_e.getField();
#pragma omp parallel for reduction(+:result) num_threads(14)
        for (int i = 0; i != field_e.size(); i++) {
            result += sigmaOverCell(k, tau, cells_e[i], field_e[i], Types::Vector3c::Zero());
        }
    }

    for (const auto &j_m : j_m_v) {
        const auto &cells_m = j_m.getManifold().getCells();
        const auto &field_m = j_m.getField();
#pragma omp parallel for reduction(+:result) num_threads(14)
        for (int i = 0; i != field_m.size(); i++) {
            result += sigmaOverCell(k, tau, cells_m[i], Types::Vector3c::Zero(), field_m[i]);
        }
    }
    return Math::Constants::inverse_4PI<Types::scalar>() * result.squaredNorm();
}


Types::scalar calculateESA(const Types::Vector3d &tau, Types::complex_d k, const Math::SurfaceVectorField &j_e) {
    const auto &cells_e = j_e.getManifold().getCells();
    const auto &field_e = j_e.getField();
    Types::Vector3c result = Types::Vector3c::Zero();
#pragma omp parallel for reduction(+ : result) num_threads(14)
    for (int i = 0; i != field_e.size(); i++) {
        result += sigmaOverCell(k, tau, cells_e[i], field_e[i], Types::Vector3c::Zero());
    }
    return Math::Constants::inverse_4PI<Types::scalar>() * result.squaredNorm();
}

Types::scalar calculateRSP(const Types::Vector3d &tau, Types::complex_d k, const std::string &field_name,
                           const Mesh::VolumeMesh::CubeMeshWithData &cube_mesh) {
    auto&& j_data = cube_mesh.getVectorData(field_name);
    auto&& eps_data = cube_mesh.getScalarData("eps");
    Types::Vector3c result = Types::Vector3c::Zero();
#pragma omp parallel for reduction(+ : result) num_threads(14)
    for (size_t i = 0; i < cube_mesh.getCells().size(); i++) {
        result += sigmaOverCube(k, tau, cube_mesh.leftDownCorner(i), cube_mesh.dx(), cube_mesh.dy(), cube_mesh.dz(),
                                j_data[i], Types::Vector3c::Zero(), eps_data[i]);
    }
    return Math::Constants::inverse_4PI<Types::scalar>() * result.squaredNorm();
}

} // namespace EMW::ESA
