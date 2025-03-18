//
// Created by evgen on 30.01.24.
//

#ifndef ELECTROMAGNETIC_WAVES_SCATTERING_VTKFUNCTIONS_HPP
#define ELECTROMAGNETIC_WAVES_SCATTERING_VTKFUNCTIONS_HPP

#include "mesh/SurfaceMesh.hpp"
#include "math/fields/SurfaceVectorField.hpp"
#include "types/Types.hpp"

#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkXMLMultiBlockDataWriter.h>

namespace VTK {
namespace detail {
vtkSmartPointer<vtkUnstructuredGrid> formUnstructuredGrid(const EMW::Mesh::SurfaceMesh &mesh);

template <typename field_t>
void dumpField(const vtkSmartPointer<vtkUnstructuredGrid> &unstructuredGrid, const field_t &field) {
    auto real_field = vtkSmartPointer<vtkDoubleArray>::New();
    real_field->SetName((field.getName() + "_real").c_str());
    real_field->SetNumberOfComponents(field.getNumberOfComponents());
    auto imag_field = vtkSmartPointer<vtkDoubleArray>::New();
    imag_field->SetNumberOfComponents(field.getNumberOfComponents());
    imag_field->SetName((field.getName() + "_imag").c_str());

    // Обходим все точки пространственной окружающей расчётной сетки
    if constexpr (field.getNumberOfComponents() >= 3) {
        for (const auto &f : field.getField()) {
            double f_real[3] = {f(0).real(), f(1).real(), f(2).real()};
            double f_imag[3] = {f(0).imag(), f(1).imag(), f(2).imag()};
            real_field->InsertNextTuple(f_real);
            imag_field->InsertNextTuple(f_imag);
        }
    } else {
        for (const auto &f : field.getField()) {
            double f_real[1] = {f.real()};
            double f_imag[1] = {f.imag()};
            real_field->InsertNextTuple(f_real);
            imag_field->InsertNextTuple(f_imag);
        }
    }
    unstructuredGrid->GetCellData()->AddArray(real_field);
    unstructuredGrid->GetCellData()->AddArray(imag_field);
    std::cout << "Field " << field.getName() << " dumped" << std::endl;
};

} // namespace detail

void surface_snapshot(const EMW::Mesh::SurfaceMesh &mesh, const std::string &part_to_file);

void field_snapshot(const EMW::Math::SurfaceVectorField &field, const std::string &part_to_file);

void united_snapshot(const std::vector<EMW::Math::SurfaceVectorField> &vectorFields,
                     const std::vector<EMW::Math::SurfaceScalarField> &scalarFiends, const EMW::Mesh::SurfaceMesh &mesh,
                     const std::string &path_to_file, int number = 0);

void field_in_points_snapshot(const std::vector<std::vector<EMW::Types::Vector3c>> &fields,
                              const std::vector<std::string> &names,
                              const std::vector<EMW::Types::Vector3d> &points,
                              const std::string& mesh_name,
                              const std::string &path_to_file);

template<typename TopologicalStructure>
void geometry_snapshot(const TopologicalStructure& geometry, const std::string &part_to_file) {
    vtkSmartPointer<vtkMultiBlockDataSet> multiBlockDataSet = vtkSmartPointer<vtkMultiBlockDataSet>::New();

    // add each data set
    for (int i = 0; i != geometry.size(); ++i) {
        const auto pointer_to_mesh = detail::formUnstructuredGrid(geometry.get(i));
        multiBlockDataSet->SetBlock(i, pointer_to_mesh);
    }
    // write the result
    vtkSmartPointer<vtkXMLMultiBlockDataWriter> writer = vtkSmartPointer<vtkXMLMultiBlockDataWriter>::New();
    writer->SetFileName(part_to_file.c_str());
    writer->SetInputData(multiBlockDataSet);
    writer->Write();
}

template<typename field_set_t>
void set_of_fields_snapshot(const field_set_t& fields, const std::string &part_to_file) {
    vtkSmartPointer<vtkMultiBlockDataSet> multiBlockDataSet = vtkSmartPointer<vtkMultiBlockDataSet>::New();
    // дампим сначала все сетки и соответвующие поля
    // электрические
    for (int i = 0; i != fields.get_electric_fields().size(); ++i) {
        const auto pointer_to_mesh = detail::formUnstructuredGrid(fields.get_electric_fields()[i].getManifold());
        detail::dumpField(pointer_to_mesh, fields.get_electric_fields()[i]);
        multiBlockDataSet->SetBlock(i, pointer_to_mesh);
    }
    // магнитные
    for (int i = 0; i != fields.get_magnetic_fields().size(); ++i) {
        const auto pointer_to_mesh = detail::formUnstructuredGrid(fields.get_magnetic_fields()[i].getManifold());
        detail::dumpField(pointer_to_mesh,fields.get_magnetic_fields()[i]);
        multiBlockDataSet->SetBlock(i + fields.get_electric_fields().size(), pointer_to_mesh);
    }
    // write the result
    vtkSmartPointer<vtkXMLMultiBlockDataWriter> writer = vtkSmartPointer<vtkXMLMultiBlockDataWriter>::New();
    writer->SetFileName(part_to_file.c_str());
    writer->SetInputData(multiBlockDataSet);
    writer->Write();
}
} // namespace VTK

#endif //ELECTROMAGNETIC_WAVES_SCATTERING_VTKFUNCTIONS_HPP
