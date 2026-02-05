//
// Created by evgen on 27.01.2026.
//

#ifndef PROJECTORONMESH_HPP
#define PROJECTORONMESH_HPP

#include "math/integration/gauss_quadrature/GaussLegenderPoints.hpp"
#include "math/integration/Quadrature.hpp"

#include "types/Types.hpp"


namespace EMW::Operators::Volume
{
    /**
    * @brief Оператор проекции на кубическую сетку в пространстве через Галеркинскую проекцию
    * с кусочно-постоянными базисными функциями
     */
    class ProjectorOnMesh
    {
        using mesh_t = Mesh::VolumeMesh::CubeMesh;
        const mesh_t& mesh_;
        using vector_t = Types::VectorXc;

    public:
        ProjectorOnMesh(const mesh_t& mesh): mesh_(mesh)
        {
        }

        template <typename Callable>
        vector_t operator ()(Callable&& function) const;
    };

    template <typename Callable>
    auto ProjectorOnMesh::operator()(Callable&& function) const -> vector_t
    {
        // Проиходимся по всем кубам и оформляем себе скаларные произведения с кусочно-постоянными базисными функциями
        const Types::index n_cubes = mesh_.getCells().size();
        vector_t result = vector_t::Zero(3 * n_cubes);
        for (Types::index i = 0; i < n_cubes; ++i)
        {
            const auto cube_i = mesh_.leftDownCorner(i);
            const auto integrand = [function](
                Types::scalar x, Types::scalar y, Types::scalar z)
            {
                return function({x, y, z});
            };
            const Types::Vector3c value = DefiniteIntegrals::integrate<DefiniteIntegrals::GaussLegendre::Quadrature<
                4, 4, 4>>(
                integrand, {cube_i[0], cube_i[1], cube_i[2]},
                {mesh_.dx(), mesh_.dy(), mesh_.dz()});
            result.block(3 * i, 0, 3, 1) = value;
        }
        return result;
    }
}
#endif //PROJECTORONMESH_HPP
