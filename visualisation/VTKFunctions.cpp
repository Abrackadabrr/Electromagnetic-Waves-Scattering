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

        // Электромагнитное поле
        auto E = vtkSmartPointer<vtkDoubleArray>::New();
        E->SetName("E");
        E->SetNumberOfComponents(3);
        auto H = vtkSmartPointer<vtkDoubleArray>::New();
        H->SetNumberOfComponents(3);
        H->SetName("H");
        // Токи
        auto J = vtkSmartPointer<vtkDoubleArray>::New();
        J->SetNumberOfComponents(3);
        J->SetName("J");
        // Локальные базисы
        auto tau1 = vtkSmartPointer<vtkDoubleArray>::New();
        tau1->SetNumberOfComponents(3);
        tau1->SetName("tau1");
        auto tau2 = vtkSmartPointer<vtkDoubleArray>::New();
        tau2->SetNumberOfComponents(3);
        tau2->SetName("tau2");
        auto n = vtkSmartPointer<vtkDoubleArray>::New();
        n->SetNumberOfComponents(3);
        n->SetName("n");

        const EMW::Containers::vector<EMW::Mesh::Point> &nodes = mesh.getNodes();
        const EMW::Containers::vector<EMW::Mesh::IndexedCell> &cells = mesh.getCells();

        // Обходим точки коллакации нашей сетки
        for (auto &cell: cells) {
            dumpPoints->InsertNextPoint(cell.collPoint_.point_.x(), cell.collPoint_.point_.y(),
                                        cell.collPoint_.point_.z());
            E->InsertNextTuple(cell.collPoint_.E_.data());
            H->InsertNextTuple(cell.collPoint_.H_.data());
            J->InsertNextTuple(cell.collPoint_.J_.data());
            tau1->InsertNextTuple(cell.tau1.data());
            tau2->InsertNextTuple(cell.tau2.data());
            n->InsertNextTuple(cell.normal.data());
            assert(std::abs(cell.normal.norm() - 1) < 1e-10);
            assert(std::abs(cell.tau1.norm() - 1) < 1e-10);
            assert(std::abs(cell.tau2.norm() - 1) < 1e-10);
        }

        // Обходим все точки нашей расчётной сетки
        for (auto &node: nodes) {
            // Вставляем новую точку в сетку VTK-снапшота
            dumpPoints->InsertNextPoint(node.x(), node.y(), node.z());
            double zero[3] = {0, 0, 0};
            E->InsertNextTuple(zero);
            H->InsertNextTuple(zero);
            J->InsertNextTuple(zero);
            tau1->InsertNextTuple(zero);
            tau2->InsertNextTuple(zero);
            n->InsertNextTuple(zero);
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
        unstructuredGrid->GetPointData()->AddArray(tau1);
        unstructuredGrid->GetPointData()->AddArray(tau2);
        unstructuredGrid->GetPointData()->AddArray(n);

        // Создаём снапшот в файле с заданным именем
        std::string fileName = mesh.getName() + "-step-" + std::to_string(snap_number) + ".vtu";
        vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
        writer->SetFileName((path_to_file + fileName).c_str());
        writer->SetInputData(unstructuredGrid);
        writer->Write();
    }
}