//
// Created by evgen on 10.02.2026.
//

#ifndef CUBE_MESH_WITH_DATA_HPP
#define CUBE_MESH_WITH_DATA_HPP

#include <ranges>

#include "types/Types.hpp"
#include "mesh/volume_mesh/CubeMesh.hpp"
#include "math/integration/decart/Quadrature.hpp"

namespace EMW::Mesh::VolumeMesh
{
    class CubeMeshWithData : public CubeMesh
    {
        using Base = CubeMesh;

        using scalarDataContainer = std::vector<Types::complex_d>;
        using vectorDataContainer = std::vector<Types::Vector3c>;

        Containers::map<std::string, scalarDataContainer> scalar_data;
        Containers::map<std::string, vectorDataContainer> vector_data;
        Types::index nodesCount;
        Types::index cellsCount;

        CubeMeshWithData(NodesContainer_t&& nodes, CellsContainer_t&& cells): Base(std::move(nodes), std::move(cells)),
                                                                              cellsCount(cells_.size()),
                                                                              nodesCount(nodes_.size())
        {
        }

        template <typename Container>
        void loadScalarData(Container&& specific_scalar_data, const std::string& name)
        {
            static_assert(std::is_convertible_v<typename Container::value_type, Types::complex_d>);
            if (specific_scalar_data.size() != cells_.size())
                throw std::runtime_error("Scalar data size mismatch with mesh's cells amount");
            // load
            scalar_data[name] = std::forward<Container>(specific_scalar_data);
        }

    public:
        CubeMeshWithData() = default;

        CubeMeshWithData(const Types::Vector3d& minCorner, Types::scalar xs, Types::scalar ys, Types::scalar zs,
                         std::size_t nx,
                         std::size_t ny, std::size_t nz) : Base(minCorner, xs, ys, zs, nx, ny, nz),
                                                           nodesCount(nodes_.size()), cellsCount(cells_.size())
        {
        }

        CubeMeshWithData(const Types::Vector3d& minCorner, Types::scalar xs, Types::index ns) :
            CubeMeshWithData(minCorner, xs, xs, xs, ns, ns, ns)
        {
        };

        // --- Setters ----- //
        template <typename Container>
        void setScalarData(const std::string& name, Container&& data);
        template <typename Container>
        void setVectorData(const std::string& name, Container&& data);
        // --- Invokers --- //
        template <typename Callable>
        void invokeScalarData(const std::string& name, Callable&& cell_function) requires std::is_invocable_v<
            Callable, CellsType>;
        template <typename Callable>
        void invokeScalarData(const std::string& name, Callable&& point_function) requires std::is_invocable_v<
            Callable, Types::point_t>;
        template <typename Quadrature, typename Callable>
        void smoothScalarData(const std::string& name, Callable&& point_function) requires std::is_invocable_v<
            Callable, Types::point_t>;

        template <typename Callable>
        void invokeVectorData(const std::string& name, Callable&& point_function) requires std::is_invocable_v<
            Callable, CellsType>;
        template <typename Callable>
        void invokeVectorData(const std::string& name, Callable&& point_function) requires std::is_invocable_v<
            Callable, Types::point_t>;

        // ---- Getters ---- //
        [[nodiscard]] const auto& getScalarData(const std::string& name) const
        {
            return scalar_data.find(name)->second;
        };

        [[nodiscard]] const auto& getVectorData(const std::string& name) const
        {
            return vector_data.find(name)->second;
        };
        [[nodiscard]] const auto& getScalarData() const { return scalar_data; }
        [[nodiscard]] const auto& getVectorData() const { return vector_data; };

        /**
        * Возвращает cell-data на сетке в виде вектора с нумерацией, согласно нумерации ячеек в сетке
        *
        * @param name - имя поля данных, которые хотим получить
        * @return комплексный вектор размерности количества ячеек в сетке
         */
        [[nodiscard]] Types::VectorXc getScalarDataAsVector(const std::string& name) const;

        template <typename Selector>
        [[nodiscard]] CubeMeshWithData getSubmeshBasedOnScalarData(Selector&& select_cell, const std::string& data_name)
            requires std::is_invocable_v<Selector, Types::complex_d>;
    };

    template <typename Container>
    void CubeMeshWithData::setScalarData(const std::string& name, Container&& data)
    {
        static_assert(std::is_same_v<typename Container::value_type, Types::complex_d>);
        if (scalar_data.contains(name))
            throw std::runtime_error("CubeMeshWithData::setVectorData: scalar has already been added. Rename");
        this->scalar_data[name] = std::forward<Container>(data);
    }

    template <typename Container>
    void CubeMeshWithData::setVectorData(const std::string& name, Container&& data)
    {
        static_assert(std::is_same_v<typename Container::value_type, Types::Vector3c>);
        if (vector_data.contains(name))
            throw std::runtime_error("CubeMeshWithData::setVectorData: vector has already been added. Rename");
        this->vector_data[name] = std::forward<Container>(data);
    }

    template <typename Callable>
    void CubeMeshWithData::invokeScalarData(const std::string& name, Callable&& cell_function) requires std::
        is_invocable_v<Callable, VolumeCells::IndexedCube>
    {
        // generate a Containers::vector<Types::complex_d> from cell function
        Containers::vector<Types::complex_d> data;
        data.reserve(cellsCount);
        for (std::size_t i = 0; i < cellsCount; ++i)
            data.emplace_back(std::forward<Callable>(cell_function)(cells_[i]));
        setScalarData(name, std::move(data));
    }

    template <typename Callable>
    void CubeMeshWithData::invokeScalarData(const std::string& name, Callable&& point_function) requires std::
        is_invocable_v<Callable, Types::point_t>
    {
        // generate a Containers::vector<Types::complex_d> from point function
        Containers::vector<Types::complex_d> data;
        data.reserve(cellsCount);
        for (std::size_t i = 0; i < cellsCount; ++i)
            data.emplace_back(std::forward<Callable>(point_function)(cells_[i].center_));
        setScalarData(name, std::move(data));
    }

    template <typename Quadrature, typename Callable>
    void CubeMeshWithData::smoothScalarData(const std::string& name, Callable&& point_function) requires
        std::is_invocable_v<Callable, Types::point_t>
    {
        // интегрально сгладить функцию по кубу
        const auto function_to_integrate = [point_function](Types::scalar x, Types::scalar y, Types::scalar z)
        {
            return point_function({x, y, z});
        };

        const auto integral_over_cube = [&function_to_integrate, this](const CellsType& cell)
        {
            const auto ldc = this->nodes_[cell.nodes_[0]];
            return Math::Integration::Numerical::Decart::integrate<Quadrature>(function_to_integrate, {ldc.x(), ldc.y(), ldc.z()},
                                                            {this->dx(), this->dy(), this->dz()}) /
                (this->dx() * this->dy() * this->dz());
        };

        invokeScalarData(name, integral_over_cube);
    }

    template <typename Callable>
    void CubeMeshWithData::invokeVectorData(const std::string& name, Callable&& point_function) requires std::
        is_invocable_v<Callable, VolumeCells::IndexedCube>
    {
        throw std::runtime_error("CubeMeshWithData::invokeVectorData: callable over cell is not implemented yet");
    }

    template <typename Callable>
    void CubeMeshWithData::invokeVectorData(const std::string& name, Callable&& point_function) requires std::
        is_invocable_v<Callable, Eigen::Matrix<double, 3, 1>>
    {
        // generate a Containers::vector<Types::complex_d> from cube function
        Containers::vector<Types::complex_d> data;
        data.reserve(cellsCount);
        for (std::size_t i = 0; i < cellsCount; ++i)
            data.emplace_back(std::forward<Callable>(point_function)(cells_[i].center_));
        setVectorData(name, std::move(data));
    }

    template <typename Selector>
    CubeMeshWithData CubeMeshWithData::getSubmeshBasedOnScalarData(Selector&& select_cell, const std::string& data_name)
        requires std::is_invocable_v<
            Selector, Types::complex_d>
    {
        decltype(cells_) new_cells;
        new_cells.reserve(cells_.size());
        Containers::vector<Types::complex_d> new_data;
        new_data.reserve(cells_.size());
        auto&& selective_data = getScalarData(data_name);
        size_t total_new_cells = 0;
        for (size_t idx = 0; idx < cellsCount; ++idx)
        {
            if (select_cell(selective_data[idx]))
            {
                total_new_cells++;
                // save specified data
                new_data.push_back(selective_data[idx]);
                new_cells.push_back(cells_[idx]);
            }
        }
        auto new_nodes = nodes_;
        CubeMeshWithData submesh{std::move(new_nodes), std::move(new_cells)};
        submesh.loadScalarData(std::move(new_data), data_name);
        return submesh;
    }
}

#endif // CUBE_MESH_WITH_DATA_HPP
