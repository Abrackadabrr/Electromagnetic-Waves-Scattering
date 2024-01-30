//
// Created by evgen on 30.01.24.
//

#include "VTKFunctions.hpp"
#include <vtkDoubleArray.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkPolygon.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkUnstructuredGrid.h>
#include <vtkSmartPointer.h>

namespace EMW::Mesh {
    void Mesh::test_snapshot(Types::index snap_number, const SurfaceMesh &mesh) {
        // VTK grid
        vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
        // VTK points
        vtkSmartPointer<vtkPoints> dumpPoints = vtkSmartPointer<vtkPoints>::New();

        auto E = vtkSmartPointer<vtkDoubleArray>::New();
        E->SetName("E");
        E->SetNumberOfComponents(3);
//        auto H = vtkSmartPointer<vtkDoubleArray>::New();
//        H->SetNumberOfComponents(3);
//        H->SetName("H");
//        auto J = vtkSmartPointer<vtkDoubleArray>::New();
//        J->SetNumberOfComponents(3);
//        J->SetName("J");

        const Containers::vector<Point> &nodes = mesh.getNodes();
        const Containers::vector<IndexedCell> &cells = mesh.getCells();

        // Обходим точки коллакации нашей сетки
        for (auto &cell : cells) {
            dumpPoints->InsertNextPoint(cell.collPoint_.x, cell.collPoint_.y, cell.collPoint_.z);
            E->InsertNextValue(Types::Vector3d{0, 0, 1});
        }

        // Обходим все точки нашей расчётной сетки
        for (auto &node : nodes) {
            // Вставляем новую точку в сетку VTK-снапшота
            dumpPoints->InsertNextPoint(node.x(), node.y(), node.z());
        }

        // Грузим точки в сетку
        unstructuredGrid->SetPoints(dumpPoints);

        // А теперь пишем, как наши точки объединены в четырехугольники (поверхностныее полигоны)
        for (auto &cell : cells) {
            auto poly = vtkSmartPointer<vtkPolygon>::New();
            poly->GetPointsIds()->SetNumberOfIds(4);
            poly->GetPointIds()->SetId(0, cell.points_[0] + cells.size());
            poly->GetPointIds()->SetId(1, cell.points_[1] + cells.size());
            poly->GetPointIds()->SetId(2, cell.points_[2] + cells.size());
            poly->GetPointIds()->SetId(3, cell.points_[3] + cells.size());
//        tetra->GetPoints()->SetPoint(0, _nodes[cell.nodes[0]].x, _nodes[cell.nodes[0]].y, _nodes[cell.nodes[0]].z);
//        tetra->GetPoints()->SetPoint(1, _nodes[cell.nodes[1]].x, _nodes[cell.nodes[1]].y, _nodes[cell.nodes[1]].z);
//        tetra->GetPoints()->SetPoint(2, _nodes[cell.nodes[2]].x, _nodes[cell.nodes[2]].y, _nodes[cell.nodes[2]].z);
//        tetra->GetPoints()->SetPoint(3, _nodes[cell.nodes[3]].x, _nodes[cell.nodes[3]].y, _nodes[cell.nodes[3]].z);
            unstructuredGrid->InsertNextCell(tetra->GetCellType(), tetra->GetPointIds());
        }

        unstructuredGrid->GetPointData()->AddArray(E);

        // Создаём снапшот в файле с заданным именем
        std::string fileName = mesh.getName() + "-step-" + std::to_string(snap_number) + ".vtu";
        vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
        writer->SetFileName(fileName.c_str());
        writer->SetInputData(unstructuredGrid);
        writer->Write();
    }
}