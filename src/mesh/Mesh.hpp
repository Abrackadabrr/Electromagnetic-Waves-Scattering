//
// Created by evgen on 25.01.24.
//

#ifndef ELECTROMAGNETIC_WAVES_SCATTERING_MESH_HPP
#define ELECTROMAGNETIC_WAVES_SCATTERING_MESH_HPP

#include <cmath>
#include <utility>
#include "Types.hpp"
#include "MeshTypes.hpp"

#include <vtkDoubleArray.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkXMLStructuredGridWriter.h>
#include <vtkStructuredGrid.h>
#include <vtkSmartPointer.h>

namespace EMW::Mesh {
    /**
     * Класс расчетной поверхностной сетки
     */
    class SurfaceMesh {
        Containers::vector<Point> nodes_;
        Containers::vector<IndexedCell> cells_;
        Containers::vector<Node> collocationPoints_;
    public:
        SurfaceMesh() = default;

        SurfaceMesh(Containers::vector<Point> nodes,
                    Containers::vector<IndexedCell> cells,
                    Containers::vector<Node> collocationPoints);
    };
}
#endif //ELECTROMAGNETIC_WAVES_SCATTERING_MESH_HPP
