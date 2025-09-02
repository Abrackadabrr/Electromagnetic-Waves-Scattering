//
// Created by evgen on 30.01.24.
//

#include "VTKFunctions.hpp"
#include <ranges>
#include <vtkPoints.h>
#include <vtkPolygon.h>
#include <vtkUnstructuredGrid.h>

namespace VTK {

vtkSmartPointer<vtkUnstructuredGrid> detail::formUnstructuredGrid(const EMW::Mesh::SurfaceMesh &mesh) {
    // VTK grid
    vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
    unstructuredGrid->Allocate(mesh.getNodes().size());
    // VTK points
    vtkSmartPointer<vtkPoints> dumpPoints = vtkSmartPointer<vtkPoints>::New();
    // Поля, специфичные для поверхностной сетки: локальные базисы
    auto tau1 = vtkSmartPointer<vtkDoubleArray>::New();
    tau1->SetNumberOfComponents(3);
    tau1->SetName("tau1");
    auto tau2 = vtkSmartPointer<vtkDoubleArray>::New();
    tau2->SetNumberOfComponents(3);
    tau2->SetName("tau2");
    auto n = vtkSmartPointer<vtkDoubleArray>::New();
    n->SetNumberOfComponents(3);
    n->SetName("n");

    const EMW::Containers::vector<EMW::Mesh::point_t> &nodes = mesh.getNodes();
    const EMW::Containers::vector<EMW::Mesh::IndexedCell> &cells = mesh.getCells();

    // Обходим точки коллакации нашей сетки
    for (const auto &cell : cells) {
        dumpPoints->InsertNextPoint(cell.collPoint_.x(), cell.collPoint_.y(), cell.collPoint_.z());
        tau1->InsertNextTuple(cell.tau[0].data());
        tau2->InsertNextTuple(cell.tau[1].data());
        n->InsertNextTuple(cell.normal.data());
        assert(std::abs(cell.normal.norm() - 1) < 1e-10);
        assert(std::abs(cell.tau[0].norm() - 1) < 1e-10);
        assert(std::abs(cell.tau[1].norm() - 1) < 1e-10);
    }

    // Обходим все точки нашей расчётной сетки
    for (auto &node : nodes) {
        // Вставляем новую точку в сетку VTK-снапшота
        dumpPoints->InsertNextPoint(node.x(), node.y(), node.z());
    }

    // Грузим точки в сетку
    unstructuredGrid->SetPoints(dumpPoints);

    // А теперь пишем, как наши точки объединены в четырехугольники (поверхностныее полигоны)
    for (const auto &cell : cells) {
        auto poly = vtkSmartPointer<vtkPolygon>::New();
        poly->GetPointIds()->SetNumberOfIds(4);
        poly->GetPointIds()->SetId(0, cell.points_[0] + cells.size());
        poly->GetPointIds()->SetId(1, cell.points_[1] + cells.size());
        poly->GetPointIds()->SetId(2, cell.points_[2] + cells.size());
        poly->GetPointIds()->SetId(3, cell.points_[3] + cells.size());
        unstructuredGrid->InsertNextCell(poly->GetCellType(), poly->GetPointIds());
    }

    unstructuredGrid->GetCellData()->AddArray(tau1);
    unstructuredGrid->GetCellData()->AddArray(tau2);
    unstructuredGrid->GetCellData()->AddArray(n);

    return unstructuredGrid;
}

void surface_snapshot(const EMW::Mesh::SurfaceMesh &mesh, const std::string &path_to_file) {
    // VTK grid
    vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid = detail::formUnstructuredGrid(mesh);

    // Создаём снапшот в файле с заданным именем
    std::string fileName = mesh.getName() + ".vtu";
    vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
    writer->SetFileName((path_to_file + fileName).c_str());
    writer->SetInputData(unstructuredGrid);
    writer->Write();
}

void field_snapshot(const EMW::Math::SurfaceVectorField &field, const std::string &path_to_file) {
    // VTK grid
    vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
    unstructuredGrid->Allocate(field.getManifold().getCells().size());
    // VTK points
    vtkSmartPointer<vtkPoints> dumpPoints = vtkSmartPointer<vtkPoints>::New();

    auto real_field = vtkSmartPointer<vtkDoubleArray>::New();
    real_field->SetName((field.getName() + "_real").c_str());
    real_field->SetNumberOfComponents(3);
    auto imag_field = vtkSmartPointer<vtkDoubleArray>::New();
    imag_field->SetNumberOfComponents(3);
    imag_field->SetName((field.getName() + "_imag").c_str());

    // Обходим все точки пространственной окружающей расчётной сетки
    const auto &cells = field.getManifold().getCells();
    for (auto [i, cell] : cells | std::views::enumerate) {
        const auto node = cell.collPoint_;
        // Вставляем новую точку в сетку VTK-снапшота
        dumpPoints->InsertNextPoint(node.x(), node.y(), node.z());
        const auto f = field.getField()[i];
        double f_real[3] = {f(0).real(), f(1).real(), f(2).real()};
        double f_imag[3] = {f(0).imag(), f(1).imag(), f(2).imag()};
        real_field->InsertNextTuple(f_real);
        imag_field->InsertNextTuple(f_imag);
    }

    // Грузим точки в сетку
    unstructuredGrid->SetPoints(dumpPoints);

    unstructuredGrid->GetPointData()->AddArray(real_field);
    unstructuredGrid->GetPointData()->AddArray(imag_field);

    // Создаём снапшот в файле с заданным именем
    std::string fileName = field.getName() + ".vtu";
    vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
    writer->SetFileName((path_to_file + fileName).c_str());
    writer->SetInputData(unstructuredGrid);
    writer->Write();
}

void field_in_points_snapshot(const std::vector<std::vector<EMW::Types::Vector3c>> &fields,
                              const std::vector<std::vector<EMW::Types::scalar>> &scalar_fields,
                              const std::vector<std::string> &vector_names,
                              const std::vector<std::string> &scalar_names,
                              const std::vector<EMW::Types::Vector3d> &points, const std::string &mesh_name,
                              const std::string &path_to_file) {
    // 1) формируем сетку из никак не соединенных точек
    vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
    unstructuredGrid->Allocate(points.size());
    // VTK points
    vtkSmartPointer<vtkPoints> dumpPoints = vtkSmartPointer<vtkPoints>::New();

    // Обходим точки коллакации нашей сетки
    for (const auto &point : points) {
        dumpPoints->InsertNextPoint(point.x(), point.y(), point.z());
    }

    // Грузим точки в сетку и готово
    unstructuredGrid->SetPoints(dumpPoints);

    // 2) Записываем поля в эту сетку
    for (const auto [field, name] : std::views::zip(fields, vector_names)) {
        detail::VTKFieldsInitialiser<std::vector<EMW::Types::Vector3c>>(field, name).dumpField(unstructuredGrid);
    }

    for (const auto [field, name] : std::views::zip(scalar_fields, scalar_names)) {
        detail::VTKFieldsInitialiser<std::vector<EMW::Types::scalar>>(field, name).dumpField(unstructuredGrid);
    }

    // 3) Делаем запись
    std::string fileName = mesh_name + ".vtu";
    vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
    writer->SetFileName((path_to_file + fileName).c_str());
    writer->SetInputData(unstructuredGrid);
    writer->Write();
}


}


