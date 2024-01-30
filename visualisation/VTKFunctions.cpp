//
// Created by evgen on 30.01.24.
//

#include "VTKFunctions.hpp"
#include <vtkDoubleArray.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkPolygon.h>
#include <vtkTetra.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkUnstructuredGrid.h>
#include <vtkSmartPointer.h>

namespace VTK {
    void
    test_snapshot(EMW::Types::index snap_number, const EMW::Mesh::SurfaceMesh &mesh, const std::string &path_to_file) {
        // VTK grid
        vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
        unstructuredGrid->Allocate(100);
        // VTK points
        vtkSmartPointer<vtkPoints> dumpPoints = vtkSmartPointer<vtkPoints>::New();

        auto E = vtkSmartPointer<vtkDoubleArray>::New();
        E->SetName("E");
        E->SetNumberOfComponents(3);
        auto H = vtkSmartPointer<vtkDoubleArray>::New();
        H->SetNumberOfComponents(3);
        H->SetName("H");
        auto J = vtkSmartPointer<vtkDoubleArray>::New();
        J->SetNumberOfComponents(3);
        J->SetName("J");

        const EMW::Containers::vector<EMW::Mesh::Point> &nodes = mesh.getNodes();
        const EMW::Containers::vector<EMW::Mesh::IndexedCell> &cells = mesh.getCells();

        // Обходим точки коллакации нашей сетки
        for (auto &cell: cells) {
            dumpPoints->InsertNextPoint(cell.collPoint_.point_.x(), cell.collPoint_.point_.y(),
                                        cell.collPoint_.point_.z());
            double e[3] = {cell.collPoint_.E_.x(), cell.collPoint_.E_.y(), cell.collPoint_.E_.z()};
            double h[3] = {cell.collPoint_.H_.x(), cell.collPoint_.H_.y(), cell.collPoint_.H_.z()};
            double j[3] = {cell.collPoint_.J_.x(), cell.collPoint_.J_.y(), cell.collPoint_.J_.z()};
            E->InsertNextTuple(e);
            H->InsertNextTuple(h);
            J->InsertNextTuple(j);
        }

        // Обходим все точки нашей расчётной сетки
        for (auto &node: nodes) {
            // Вставляем новую точку в сетку VTK-снапшота
            dumpPoints->InsertNextPoint(node.x(), node.y(), node.z());
            double zero_field[3] = {0, 0, 0};
            E->InsertNextTuple(zero_field);
            H->InsertNextTuple(zero_field);
            J->InsertNextTuple(zero_field);
        }

        // Грузим точки в сетку
        unstructuredGrid->SetPoints(dumpPoints);

        // А теперь пишем, как наши точки объединены в четырехугольники (поверхностныее полигоны)
        for (auto &cell: cells) {
            auto poly = vtkSmartPointer<vtkPolygon>::New();
            poly->GetPointIds()->SetNumberOfIds(4);
            poly->GetPointIds()->SetId(0, cell.points_[0] + cells.size());
            poly->GetPointIds()->SetId(1, cell.points_[1] + cells.size());
            poly->GetPointIds()->SetId(2, cell.points_[2] + cells.size());
            poly->GetPointIds()->SetId(3, cell.points_[3] + cells.size());
            unstructuredGrid->InsertNextCell(poly->GetCellType(), poly->GetPointIds());
        }

        unstructuredGrid->GetPointData()->AddArray(E);
        unstructuredGrid->GetPointData()->AddArray(H);
        unstructuredGrid->GetPointData()->AddArray(J);

        // Создаём снапшот в файле с заданным именем
        std::string fileName = mesh.getName() + "-step-" + std::to_string(snap_number) + ".vtu";
        vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
        writer->SetFileName((path_to_file + fileName).c_str());
        writer->SetInputData(unstructuredGrid);
        writer->Write();
    }
}
