//
// Created by evgen on 30.01.24.
//

#include "visualisation/include/VTKFunctions.hpp"
#include <algorithm>
#include <cmath>
#include <memory>
#include <numeric>
#include <ranges>
#include <stdexcept>
#include <string_view>
#include <unordered_map>
#include <vtkHexahedron.h>
#include <vtkPoints.h>
#include <vtkPolygon.h>
#include <vtkTriangle.h>
#include <vtkCell.h>
#include <vtkDataArray.h>
#include <vtkDataSetAttributes.h>
#include <vtkCellType.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <filesystem>


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
        // TODO: вернуть vtkPolygon назад
        auto poly = vtkSmartPointer<vtkTriangle>::New();
        poly->GetPointIds()->SetNumberOfIds(3);
        poly->GetPointIds()->SetId(0, cell.points_[0] + cells.size());
        poly->GetPointIds()->SetId(1, cell.points_[1] + cells.size());
        poly->GetPointIds()->SetId(2, cell.points_[2] + cells.size());
        //        poly->GetPointIds()->SetId(3, cell.points_[3] + cells.size());
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

void volume_mesh_snapshot(const EMW::Mesh::VolumeMesh::CubeMesh &mesh, const std::string &path_to_file) {
    // VTK grid
    vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
    unstructuredGrid->Allocate(mesh.getNodes().size());
    // VTK points
    vtkSmartPointer<vtkPoints> dumpPoints = vtkSmartPointer<vtkPoints>::New();

    const EMW::Containers::vector<EMW::Mesh::point_t> &nodes = mesh.getNodes();
    const auto &cells = mesh.getCells();

    // Обходим все точки нашей расчётной сетки
    for (auto &node : nodes) {
        // Вставляем новую точку в сетку VTK-снапшота
        dumpPoints->InsertNextPoint(node.x(), node.y(), node.z());
    }

    // Грузим точки в сетку
    unstructuredGrid->SetPoints(dumpPoints);

    // А теперь пишем, как наши точки объединены в четырехугольники (поверхностныее полигоны)
    for (const auto &cell : cells) {
        auto poly = vtkSmartPointer<vtkHexahedron>::New();
        vtkIdType ids[8] = {
            cell.nodes_[0], cell.nodes_[1], cell.nodes_[3], cell.nodes_[2],
            cell.nodes_[4], cell.nodes_[5], cell.nodes_[7], cell.nodes_[6],
        };
        unstructuredGrid->InsertNextCell(VTK_HEXAHEDRON, 8, ids);
    }

    // Создаём снапшот в файле с заданным именем
    std::string fileName = mesh.getName() + ".vtu";
    vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
    writer->SetFileName((path_to_file + fileName).c_str());
    writer->SetInputData(unstructuredGrid);
    writer->Write();
}

void volume_mesh_withdata_snapshot(const EMW::Mesh::VolumeMesh::CubeMeshWithData &mesh,
                                   const std::string &path_to_file) {
    // VTK grid
    vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
    unstructuredGrid->Allocate(mesh.getNodes().size());
    // VTK points
    vtkSmartPointer<vtkPoints> dumpPoints = vtkSmartPointer<vtkPoints>::New();

    const EMW::Containers::vector<EMW::Mesh::point_t> &nodes = mesh.getNodes();
    const EMW::Containers::vector<EMW::Mesh::VolumeCells::IndexedCube> &cells = mesh.getCells();

    // Обходим все точки нашей расчётной сетки
    for (auto &node : nodes) {
        // Вставляем новую точку в сетку VTK-снапшота
        dumpPoints->InsertNextPoint(node.x(), node.y(), node.z());
    }

    // Грузим точки в сетку
    unstructuredGrid->SetPoints(dumpPoints);

    // А теперь пишем, как наши точки объединены в четырехугольники (поверхностныее полигоны)
    for (const auto &cell : cells) {
        auto poly = vtkSmartPointer<vtkHexahedron>::New();
        vtkIdType ids[8] = {
            cell.nodes_[0], cell.nodes_[1], cell.nodes_[3], cell.nodes_[2],
            cell.nodes_[4], cell.nodes_[5], cell.nodes_[7], cell.nodes_[6],
        };
        unstructuredGrid->InsertNextCell(VTK_HEXAHEDRON, 8, ids);
    }

    // Записываем данные из сетки в cellData
    for (auto &&[name, scalar_data] : mesh.getScalarData()) {
        auto real_field = vtkSmartPointer<vtkDoubleArray>::New();
        real_field->SetName((name + "_real").c_str());
        real_field->SetNumberOfComponents(1);
        auto imag_field = vtkSmartPointer<vtkDoubleArray>::New();
        imag_field->SetNumberOfComponents(1);
        imag_field->SetName((name + "_imag").c_str());

        for (auto &&value : scalar_data) {
            const double real = value.real();
            const double imag = value.imag();
            real_field->InsertNextValue(real);
            imag_field->InsertNextValue(imag);
        }
        unstructuredGrid->GetCellData()->AddArray(real_field);
        unstructuredGrid->GetCellData()->AddArray(imag_field);
    }

    // Добавляем нумерацию ячеек в сетке
    auto cell_idx_field = vtkSmartPointer<vtkIntArray>::New();
    cell_idx_field->SetName("cell_idx");
    for (size_t idx = 0; idx < mesh.getCells().size(); ++idx) {
        cell_idx_field->InsertNextValue(idx);
    }
    unstructuredGrid->GetCellData()->AddArray(cell_idx_field);

    for (auto &&[name, vector_data] : mesh.getVectorData()) {
        auto real_field = vtkSmartPointer<vtkDoubleArray>::New();
        real_field->SetName((name + "_real").c_str());
        real_field->SetNumberOfComponents(3);
        auto imag_field = vtkSmartPointer<vtkDoubleArray>::New();
        imag_field->SetNumberOfComponents(3);
        imag_field->SetName((name + "_imag").c_str());

        for (auto &&value : vector_data) {
            const double real[3] = {value.real()[0], value.real()[1], value.real()[2]};
            const double imag[3] = {value.imag()[0], value.imag()[1], value.imag()[2]};
            real_field->InsertNextTuple(real);
            imag_field->InsertNextTuple(imag);
        }
        unstructuredGrid->GetCellData()->AddArray(real_field);
        unstructuredGrid->GetCellData()->AddArray(imag_field);
    }

    // Создаём снапшот в файле с заданным именем
    const std::string fileName = mesh.getName() + ".vtu";
    vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
    writer->SetFileName((path_to_file + fileName).c_str());
    writer->SetInputData(unstructuredGrid);
    writer->Write();
}



namespace {

constexpr std::string_view kRealSuffix = "_real";
constexpr std::string_view kImagSuffix = "_imag";

struct ComplexFieldArrays {
    vtkDataArray *real = nullptr;
    vtkDataArray *imag = nullptr;
};

struct StructuredGridParams {
    EMW::Types::Vector3d minCorner;
    EMW::Types::scalar xs;
    EMW::Types::scalar ys;
    EMW::Types::scalar zs;
    std::size_t nx;
    std::size_t ny;
    std::size_t nz;
};

EMW::Types::scalar axisTolerance(const EMW::Types::scalar value) {
    return std::max<EMW::Types::scalar>(1e-10, 1e-6 * std::max<EMW::Types::scalar>(1.0, std::abs(value)));
}

std::vector<EMW::Types::scalar> uniqueAxisValues(std::vector<EMW::Types::scalar> axisValues) {
    std::sort(axisValues.begin(), axisValues.end());
    std::vector<EMW::Types::scalar> uniqueValues;
    uniqueValues.reserve(axisValues.size());

    for (const auto value : axisValues) {
        if (uniqueValues.empty()) {
            uniqueValues.push_back(value);
            continue;
        }

        if (std::abs(value - uniqueValues.back()) > axisTolerance(uniqueValues.back())) {
            uniqueValues.push_back(value);
        }
    }

    return uniqueValues;
}

void validateUniformAxis(const std::vector<EMW::Types::scalar> &axisValues, const std::string &axisName) {
    if (axisValues.size() < 2) {
        throw std::runtime_error("volume_mesh_withdata_from_vtu: axis " + axisName + " has less than 2 points");
    }

    const auto fullSize = axisValues.back() - axisValues.front();
    const auto expectedStep = fullSize / static_cast<EMW::Types::scalar>(axisValues.size() - 1);
    if (expectedStep <= 0) {
        throw std::runtime_error("volume_mesh_withdata_from_vtu: axis " + axisName + " is degenerate");
    }

    const auto tol = axisTolerance(expectedStep);
    for (std::size_t idx = 1; idx < axisValues.size(); ++idx) {
        const auto step = axisValues[idx] - axisValues[idx - 1];
        if (std::abs(step - expectedStep) > tol) {
            throw std::runtime_error("volume_mesh_withdata_from_vtu: axis " + axisName + " is not uniform");
        }
    }
}

StructuredGridParams inferStructuredGridParams(vtkUnstructuredGrid *unstructuredGrid) {
    const auto nPoints = unstructuredGrid->GetNumberOfPoints();
    if (nPoints == 0) {
        throw std::runtime_error("volume_mesh_withdata_from_vtu: VTU contains no points");
    }

    std::vector<EMW::Types::scalar> xAxis;
    std::vector<EMW::Types::scalar> yAxis;
    std::vector<EMW::Types::scalar> zAxis;
    xAxis.reserve(static_cast<std::size_t>(nPoints));
    yAxis.reserve(static_cast<std::size_t>(nPoints));
    zAxis.reserve(static_cast<std::size_t>(nPoints));

    for (vtkIdType pointIdx = 0; pointIdx < nPoints; ++pointIdx) {
        double point[3]{};
        unstructuredGrid->GetPoint(pointIdx, point);
        xAxis.push_back(point[0]);
        yAxis.push_back(point[1]);
        zAxis.push_back(point[2]);
    }

    xAxis = uniqueAxisValues(std::move(xAxis));
    yAxis = uniqueAxisValues(std::move(yAxis));
    zAxis = uniqueAxisValues(std::move(zAxis));

    validateUniformAxis(xAxis, "X");
    validateUniformAxis(yAxis, "Y");
    validateUniformAxis(zAxis, "Z");

    const auto nx = xAxis.size();
    const auto ny = yAxis.size();
    const auto nz = zAxis.size();
    const auto expectedPoints = static_cast<vtkIdType>(nx * ny * nz);
    if (nPoints != expectedPoints) {
        throw std::runtime_error("volume_mesh_withdata_from_vtu: point count does not match structured grid");
    }

    StructuredGridParams params{};
    params.minCorner = EMW::Types::Vector3d{xAxis.front(), yAxis.front(), zAxis.front()};
    params.xs = xAxis.back() - xAxis.front();
    params.ys = yAxis.back() - yAxis.front();
    params.zs = zAxis.back() - zAxis.front();
    params.nx = nx;
    params.ny = ny;
    params.nz = nz;
    return params;
}

std::runtime_error surfaceReadError(const std::string &message) {
    return std::runtime_error("surface_mesh_with_vector_fields_from_vtu: " + message);
}

EMW::Types::Vector3d vtkPoint(vtkUnstructuredGrid *unstructuredGrid, const vtkIdType pointIdx) {
    double point[3]{};
    unstructuredGrid->GetPoint(pointIdx, point);
    return {point[0], point[1], point[2]};
}

vtkSmartPointer<vtkUnstructuredGrid> readUnstructuredGrid(const std::string &path_to_file) {
    auto unstructuredGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
    const auto extension = std::filesystem::path(path_to_file).extension().string();

    if (extension == ".vtk") {
        auto reader = vtkSmartPointer<vtkUnstructuredGridReader>::New();
        reader->SetFileName(path_to_file.c_str());
        reader->ReadAllVectorsOn();
        reader->Update();
        if (!reader->GetOutput()) {
            throw surfaceReadError("cannot read legacy VTK file: " + path_to_file);
        }
        unstructuredGrid->ShallowCopy(reader->GetOutput());
        return unstructuredGrid;
    }

    auto reader = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
    reader->SetFileName(path_to_file.c_str());
    reader->Update();
    if (!reader->GetOutput()) {
        throw surfaceReadError("cannot read VTU file: " + path_to_file);
    }
    unstructuredGrid->ShallowCopy(reader->GetOutput());
    return unstructuredGrid;
}

std::vector<vtkIdType> inferLeadingCollocationPointIds(vtkUnstructuredGrid *unstructuredGrid) {
    const auto nCells = unstructuredGrid->GetNumberOfCells();
    const auto nPoints = unstructuredGrid->GetNumberOfPoints();
    if (nCells == 0 || nPoints < nCells) {
        return {};
    }

    for (vtkIdType cellIdx = 0; cellIdx < nCells; ++cellIdx) {
        auto *cell = unstructuredGrid->GetCell(cellIdx);
        if (!cell) {
            throw surfaceReadError("cannot read cell " + std::to_string(cellIdx));
        }

        const auto nCellPoints = cell->GetNumberOfPoints();
        for (vtkIdType localPointIdx = 0; localPointIdx < nCellPoints; ++localPointIdx) {
            if (cell->GetPointId(localPointIdx) < nCells) {
                return {};
            }
        }
    }

    std::vector<vtkIdType> collocationPointIds(static_cast<std::size_t>(nCells));
    std::iota(collocationPointIds.begin(), collocationPointIds.end(), vtkIdType{0});
    return collocationPointIds;
}

void validateVectorArray(vtkDataArray *array, const std::string &arrayName, const vtkIdType expectedTuples,
                         const std::string &context) {
    if (!array) {
        throw surfaceReadError("missing " + context + " array " + arrayName);
    }
    if (array->GetNumberOfComponents() != 3) {
        throw surfaceReadError(context + " array " + arrayName + " must have 3 components");
    }
    if (array->GetNumberOfTuples() != expectedTuples) {
        throw surfaceReadError(context + " array " + arrayName + " has invalid tuple count");
    }
}

void applySurfaceCellData(vtkUnstructuredGrid *unstructuredGrid, EMW::Mesh::SurfaceMesh &mesh,
                          const std::vector<vtkIdType> &collocationPointIds) {
    if (!collocationPointIds.empty()) {
        for (vtkIdType cellIdx = 0; cellIdx < unstructuredGrid->GetNumberOfCells(); ++cellIdx) {
            mesh.getCells()[static_cast<std::size_t>(cellIdx)].collPoint_ =
                vtkPoint(unstructuredGrid, collocationPointIds[static_cast<std::size_t>(cellIdx)]);
        }
    }

    auto *cellData = unstructuredGrid->GetCellData();
    if (!cellData) {
        return;
    }

    auto *tau1 = cellData->GetArray("tau1");
    auto *tau2 = cellData->GetArray("tau2");
    auto *normal = cellData->GetArray("n");
    if (!tau1 && !tau2 && !normal) {
        return;
    }

    const auto nCells = unstructuredGrid->GetNumberOfCells();
    validateVectorArray(tau1, "tau1", nCells, "basis");
    validateVectorArray(tau2, "tau2", nCells, "basis");
    validateVectorArray(normal, "n", nCells, "basis");

    for (vtkIdType cellIdx = 0; cellIdx < nCells; ++cellIdx) {
        auto &cell = mesh.getCells()[static_cast<std::size_t>(cellIdx)];
        cell.tau[0] = {tau1->GetComponent(cellIdx, 0), tau1->GetComponent(cellIdx, 1),
                       tau1->GetComponent(cellIdx, 2)};
        cell.tau[1] = {tau2->GetComponent(cellIdx, 0), tau2->GetComponent(cellIdx, 1),
                       tau2->GetComponent(cellIdx, 2)};
        cell.normal = {normal->GetComponent(cellIdx, 0), normal->GetComponent(cellIdx, 1),
                       normal->GetComponent(cellIdx, 2)};
    }
}

std::shared_ptr<EMW::Mesh::SurfaceMesh>
readSurfaceMesh(vtkUnstructuredGrid *unstructuredGrid, const std::string &path_to_file,
                std::vector<vtkIdType> &collocationPointIds) {
    if (unstructuredGrid->GetNumberOfCells() == 0) {
        throw surfaceReadError("mesh file contains no cells");
    }

    EMW::Containers::vector<EMW::Mesh::point_t> nodes;
    EMW::Containers::vector<EMW::Mesh::IndexedCell::nodes_t> cells;
    nodes.reserve(static_cast<std::size_t>(unstructuredGrid->GetNumberOfPoints()));
    cells.reserve(static_cast<std::size_t>(unstructuredGrid->GetNumberOfCells()));

    std::unordered_map<vtkIdType, EMW::Types::index> pointIdToNodeIdx;

    const auto getNodeIdx = [&](const vtkIdType vtkPointId) -> EMW::Types::index {
        const auto [it, inserted] =
            pointIdToNodeIdx.emplace(vtkPointId, static_cast<EMW::Types::index>(nodes.size()));
        if (inserted) {
            nodes.push_back(vtkPoint(unstructuredGrid, vtkPointId));
        }
        return it->second;
    };

    for (vtkIdType cellIdx = 0; cellIdx < unstructuredGrid->GetNumberOfCells(); ++cellIdx) {
        auto *cell = unstructuredGrid->GetCell(cellIdx);
        if (!cell) {
            throw surfaceReadError("cannot read cell " + std::to_string(cellIdx));
        }

        const auto nCellPoints = cell->GetNumberOfPoints();
        if (nCellPoints != 3 && nCellPoints != 4) {
            throw surfaceReadError("only triangle and quadrilateral surface cells are supported");
        }

        EMW::Mesh::IndexedCell::nodes_t cellNodes{};
        for (vtkIdType localPointIdx = 0; localPointIdx < nCellPoints; ++localPointIdx) {
            cellNodes[static_cast<std::size_t>(localPointIdx)] = getNodeIdx(cell->GetPointId(localPointIdx));
        }
        if (nCellPoints == 3) {
            cellNodes[3] = cellNodes[2];
        }

        cells.push_back(cellNodes);
    }

    collocationPointIds = inferLeadingCollocationPointIds(unstructuredGrid);

    auto mesh = std::make_shared<EMW::Mesh::SurfaceMesh>(std::move(nodes), std::move(cells));
    mesh->setName(std::filesystem::path(path_to_file).stem().string());
    applySurfaceCellData(unstructuredGrid, *mesh, collocationPointIds);
    return mesh;
}

std::pair<vtkDataArray *, vtkDataArray *> findComplexVectorArrays(vtkDataSetAttributes *attributes,
                                                                 const std::string &fieldName,
                                                                 const std::string &dataLocation) {
    if (!attributes) {
        return {nullptr, nullptr};
    }

    const auto realName = fieldName + std::string(kRealSuffix);
    const auto imagName = fieldName + std::string(kImagSuffix);
    auto *real = attributes->GetArray(realName.c_str());
    auto *imag = attributes->GetArray(imagName.c_str());
    if (!real && !imag) {
        return {nullptr, nullptr};
    }
    if (!real || !imag) {
        throw surfaceReadError("field " + fieldName + " requires arrays " + realName + " and " + imagName +
                               " in " + dataLocation + " data");
    }
    if (real->GetNumberOfComponents() != 3 || imag->GetNumberOfComponents() != 3) {
        throw surfaceReadError("field " + fieldName + " must have 3 components in " + dataLocation + " data");
    }
    if (real->GetNumberOfTuples() != imag->GetNumberOfTuples()) {
        throw surfaceReadError("real and imaginary tuple count mismatch for field " + fieldName);
    }
    return {real, imag};
}

EMW::Math::SurfaceVectorField readSurfaceVectorField(const EMW::Mesh::SurfaceMesh &mesh,
                                                     vtkDataArray *realArray, vtkDataArray *imagArray,
                                                     const std::vector<vtkIdType> &tupleIds,
                                                     const std::string &fieldName) {
    EMW::Containers::vector<EMW::Types::Vector3c> fieldData(tupleIds.size());

    for (std::size_t cellIdx = 0; cellIdx < tupleIds.size(); ++cellIdx) {
        const auto tupleIdx = tupleIds[cellIdx];
        if (tupleIdx < 0 || tupleIdx >= realArray->GetNumberOfTuples()) {
            throw surfaceReadError("tuple index is out of range for field " + fieldName);
        }

        EMW::Types::Vector3c value;
        for (int component = 0; component < 3; ++component) {
            value[component] = EMW::Types::complex_d(realArray->GetComponent(tupleIdx, component),
                                                     imagArray->GetComponent(tupleIdx, component));
        }
        fieldData[cellIdx] = value;
    }

    EMW::Math::SurfaceVectorField field(mesh, fieldData);
    field.setName(fieldName);
    return field;
}

std::vector<EMW::Math::SurfaceVectorField>
readSurfaceVectorFields(const EMW::Mesh::SurfaceMesh &mesh, vtkUnstructuredGrid *unstructuredGrid,
                        const std::vector<std::string> &fieldNames,
                        const std::vector<vtkIdType> &collocationPointIds) {
    std::vector<EMW::Math::SurfaceVectorField> fields;
    fields.reserve(fieldNames.size());

    std::vector<vtkIdType> tupleIds(mesh.getCells().size());
    std::iota(tupleIds.begin(), tupleIds.end(), vtkIdType{0});

    for (const auto &fieldName : fieldNames) {
        if (auto [real, imag] = findComplexVectorArrays(unstructuredGrid->GetCellData(), fieldName, "cell");
            real && imag) {
            if (real->GetNumberOfTuples() != static_cast<vtkIdType>(mesh.getCells().size())) {
                throw surfaceReadError("cell field " + fieldName + " tuple count does not match mesh cells");
            }
            fields.push_back(readSurfaceVectorField(mesh, real, imag, tupleIds, fieldName));
            continue;
        }

        if (auto [real, imag] = findComplexVectorArrays(unstructuredGrid->GetPointData(), fieldName, "point");
            real && imag) {
            const auto nTuples = real->GetNumberOfTuples();
            if (nTuples == static_cast<vtkIdType>(mesh.getCells().size())) {
                fields.push_back(readSurfaceVectorField(mesh, real, imag, tupleIds, fieldName));
                continue;
            }

            if (!collocationPointIds.empty() && nTuples == unstructuredGrid->GetNumberOfPoints()) {
                fields.push_back(readSurfaceVectorField(mesh, real, imag, collocationPointIds, fieldName));
                continue;
            }

            throw surfaceReadError("point field " + fieldName + " cannot be mapped to mesh cells");
        }

        const auto realName = fieldName + std::string(kRealSuffix);
        const auto imagName = fieldName + std::string(kImagSuffix);
        throw surfaceReadError("field " + fieldName + " requires arrays " + realName + " and " + imagName);
    }

    return fields;
}

} // namespace

EMW::Mesh::VolumeMesh::CubeMeshWithData volume_mesh_withdata_from_vtu(const std::string &path_to_file) {
    auto reader = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
    reader->SetFileName(path_to_file.c_str());
    reader->Update();

    auto unstructuredGrid = reader->GetOutput();
    if (!unstructuredGrid) {
        throw std::runtime_error("volume_mesh_withdata_from_vtu: cannot read VTU file: " + path_to_file);
    }

    const auto params = inferStructuredGridParams(unstructuredGrid);
    const auto nCells = unstructuredGrid->GetNumberOfCells();

    const auto expectedCells = static_cast<vtkIdType>((params.nx - 1) * (params.ny - 1) * (params.nz - 1));
    if (nCells != expectedCells) {
        throw std::runtime_error("volume_mesh_withdata_from_vtu: cell count does not match structured grid");
    }

    for (vtkIdType cellIdx = 0; cellIdx < nCells; ++cellIdx) {
        if (unstructuredGrid->GetCellType(cellIdx) != VTK_HEXAHEDRON) {
            throw std::runtime_error("volume_mesh_withdata_from_vtu: only VTK_HEXAHEDRON cells are supported");
        }
    }

    auto mesh = EMW::Mesh::VolumeMesh::CubeMeshWithData(params.minCorner, params.xs, params.ys, params.zs, params.nx,
                                                        params.ny, params.nz);
    mesh.setName(std::filesystem::path(path_to_file).stem().string());

    std::unordered_map<std::string, ComplexFieldArrays> fieldsByName;
    auto cellData = unstructuredGrid->GetCellData();
    if (!cellData) {
        return mesh;
    }

    std::vector<std::size_t> dataCellToMeshCell(static_cast<std::size_t>(nCells));
    for (vtkIdType idx = 0; idx < nCells; ++idx) {
        dataCellToMeshCell[static_cast<std::size_t>(idx)] = static_cast<std::size_t>(idx);
    }

    if (auto *cellIdxArray = cellData->GetArray("cell_idx")) {
        if (cellIdxArray->GetNumberOfComponents() != 1 || cellIdxArray->GetNumberOfTuples() != nCells) {
            throw std::runtime_error("volume_mesh_withdata_from_vtu: invalid cell_idx array shape");
        }

        std::vector<bool> used(static_cast<std::size_t>(nCells), false);
        for (vtkIdType dataIdx = 0; dataIdx < nCells; ++dataIdx) {
            const auto meshCellIdx = static_cast<long long>(std::llround(cellIdxArray->GetComponent(dataIdx, 0)));
            if (meshCellIdx < 0 || meshCellIdx >= nCells) {
                throw std::runtime_error("volume_mesh_withdata_from_vtu: cell_idx value is out of range");
            }
            if (used[static_cast<std::size_t>(meshCellIdx)]) {
                throw std::runtime_error("volume_mesh_withdata_from_vtu: duplicated value in cell_idx");
            }
            used[static_cast<std::size_t>(meshCellIdx)] = true;
            dataCellToMeshCell[static_cast<std::size_t>(dataIdx)] = static_cast<std::size_t>(meshCellIdx);
        }
    }

    for (int arrIdx = 0; arrIdx < cellData->GetNumberOfArrays(); ++arrIdx) {
        auto *array = cellData->GetArray(arrIdx);
        if (!array || !array->GetName()) {
            continue;
        }

        const std::string arrayName = array->GetName();
        if (arrayName == "cell_idx") {
            continue;
        }

        const std::string_view nameView = arrayName;
        if (nameView.ends_with(kRealSuffix)) {
            const auto fieldName = arrayName.substr(0, arrayName.size() - kRealSuffix.size());
            if (fieldName.empty()) {
                throw std::runtime_error("volume_mesh_withdata_from_vtu: invalid real field name");
            }
            auto &field = fieldsByName[fieldName];
            if (field.real != nullptr) {
                throw std::runtime_error("volume_mesh_withdata_from_vtu: duplicated real part for field " + fieldName);
            }
            field.real = array;
        } else if (nameView.ends_with(kImagSuffix)) {
            const auto fieldName = arrayName.substr(0, arrayName.size() - kImagSuffix.size());
            if (fieldName.empty()) {
                throw std::runtime_error("volume_mesh_withdata_from_vtu: invalid imaginary field name");
            }
            auto &field = fieldsByName[fieldName];
            if (field.imag != nullptr) {
                throw std::runtime_error("volume_mesh_withdata_from_vtu: duplicated imaginary part for field " +
                                         fieldName);
            }
            field.imag = array;
        }
    }

    for (auto &&[name, arrays] : fieldsByName) {
        if (!arrays.real || !arrays.imag) {
            throw std::runtime_error("volume_mesh_withdata_from_vtu: complex field is incomplete: " + name);
        }

        if (arrays.real->GetNumberOfTuples() != nCells || arrays.imag->GetNumberOfTuples() != nCells) {
            throw std::runtime_error("volume_mesh_withdata_from_vtu: tuple count mismatch for field " + name);
        }

        if (arrays.real->GetNumberOfComponents() != arrays.imag->GetNumberOfComponents()) {
            throw std::runtime_error("volume_mesh_withdata_from_vtu: component mismatch for field " + name);
        }

        const auto nComponents = arrays.real->GetNumberOfComponents();
        if (nComponents == 1) {
            EMW::Containers::vector<EMW::Types::complex_d> scalarData(static_cast<std::size_t>(nCells));
            for (vtkIdType dataCellIdx = 0; dataCellIdx < nCells; ++dataCellIdx) {
                const auto meshCellIdx = dataCellToMeshCell[static_cast<std::size_t>(dataCellIdx)];
                scalarData[meshCellIdx] = {arrays.real->GetComponent(dataCellIdx, 0),
                                           arrays.imag->GetComponent(dataCellIdx, 0)};
            }
            mesh.setScalarData(name, std::move(scalarData));
            continue;
        }

        if (nComponents == 3) {
            EMW::Containers::vector<EMW::Types::Vector3c> vectorData(static_cast<std::size_t>(nCells));
            for (vtkIdType dataCellIdx = 0; dataCellIdx < nCells; ++dataCellIdx) {
                EMW::Types::Vector3c value;
                for (int comp = 0; comp < 3; ++comp) {
                    value[comp] = EMW::Types::complex_d(
                        arrays.real->GetComponent(dataCellIdx, comp), arrays.imag->GetComponent(dataCellIdx, comp));
                }
                const auto meshCellIdx = dataCellToMeshCell[static_cast<std::size_t>(dataCellIdx)];
                vectorData[meshCellIdx] = value;
            }
            mesh.setVectorData(name, std::move(vectorData));
            continue;
        }

        throw std::runtime_error(
            "volume_mesh_withdata_from_vtu: unsupported number of components in field " + name);
    }

    return mesh;
}

SurfaceMeshVTUData surface_mesh_with_vector_fields_from_vtu(const std::string &path_to_file,
                                                            const std::vector<std::string> &first_vector_field_names) {
    auto unstructuredGrid = readUnstructuredGrid(path_to_file);

    SurfaceMeshVTUData result{};
    std::vector<vtkIdType> collocationPointIds;
    result.mesh = readSurfaceMesh(unstructuredGrid, path_to_file, collocationPointIds);
    result.first_vector_fields =
        readSurfaceVectorFields(*result.mesh, unstructuredGrid, first_vector_field_names, collocationPointIds);
    return result;
}

}
