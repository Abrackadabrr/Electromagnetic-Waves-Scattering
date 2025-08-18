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
#include <vtkXMLUnstructuredGridWriter.h>

namespace VTK {
namespace detail {

template<typename Field_type>
struct VTKFieldsInitialiser {};

template<>
struct VTKFieldsInitialiser<EMW::Math::SurfaceScalarField<EMW::Types::scalar>> {
    using Field_Type = EMW::Math::SurfaceScalarField<EMW::Types::scalar>;
    using vtkPointer = vtkSmartPointer<vtkDoubleArray>;
protected:
    static constexpr EMW::Types::index n_fields = 1;
    static constexpr EMW::Types::index n_components = 1;
    EMW::Containers::array<vtkPointer, n_fields> fields;
    const std::string name;

public:
    explicit VTKFieldsInitialiser(const Field_Type& field_to_display): name(field_to_display.getName()) {
        fields[0] = vtkSmartPointer<vtkDoubleArray>::New();
        fields[0] ->SetName((field_to_display.getName()).c_str());
        fields[0] ->SetNumberOfComponents(field_to_display.getNumberOfComponents());

        for (const auto &f : field_to_display.getField()) {
            double f_real[1] = {f};
            fields[0] ->InsertNextTuple(f_real);
        }
    }

    [[nodiscard]] const EMW::Containers::array<vtkPointer, n_fields>& getFields() const { return fields; }
    void dumpField(const vtkSmartPointer<vtkUnstructuredGrid> &unstructuredGrid) const {
        for (auto&& field : fields)
            unstructuredGrid->GetCellData()->AddArray(field);
        std::cout << "Field " << name << " dumped" << std::endl;
    }
};

template<>
struct VTKFieldsInitialiser<EMW::Math::SurfaceScalarField<EMW::Types::complex_d>> {
    using Field_Type = EMW::Math::SurfaceScalarField<EMW::Types::complex_d>;
    using vtkPointer = vtkSmartPointer<vtkDoubleArray>;
protected:
    static constexpr EMW::Types::index n_fields = 2;
    static constexpr EMW::Types::index n_components = 1;
    EMW::Containers::array<vtkPointer, n_fields> fields;
    const std::string name;

public:
    explicit VTKFieldsInitialiser(const Field_Type& field_to_display): name(field_to_display.getName()) {
        fields[0] = vtkSmartPointer<vtkDoubleArray>::New();
        fields[0] ->SetName((field_to_display.getName() + "_real").c_str());
        fields[0] ->SetNumberOfComponents(field_to_display.getNumberOfComponents());

        fields[1] = vtkSmartPointer<vtkDoubleArray>::New();
        fields[1] ->SetName((field_to_display.getName() + "_imag").c_str());
        fields[1] ->SetNumberOfComponents(field_to_display.getNumberOfComponents());

        for (const auto &f : field_to_display.getField()) {
            double f_real[1] = {f.real()};
            double f_imag[3] = {f.imag()};
            fields[0] ->InsertNextTuple(f_real);
            fields[1] ->InsertNextTuple(f_imag);
        }
    }

    [[nodiscard]] const EMW::Containers::array<vtkPointer, n_fields>& getFields() const { return fields; }
    void dumpField(const vtkSmartPointer<vtkUnstructuredGrid> &unstructuredGrid) const {
        for (auto&& field : fields)
            unstructuredGrid->GetCellData()->AddArray(field);
        std::cout << "Field " << name << " dumped" << std::endl;
    }
};

template<>
struct VTKFieldsInitialiser<EMW::Math::SurfaceVectorField> {
    using Field_Type = EMW::Math::SurfaceVectorField;
    using vtkPointer = vtkSmartPointer<vtkDoubleArray>;
protected:
    static constexpr EMW::Types::index n_fields = 2;
    static constexpr EMW::Types::index n_components = 3;
    EMW::Containers::array<vtkPointer, n_fields> fields;
    const std::string name;

public:
    explicit VTKFieldsInitialiser(const Field_Type& field_to_display): name(field_to_display.getName()) {
        fields[0] = vtkSmartPointer<vtkDoubleArray>::New();
        fields[0] ->SetName((field_to_display.getName() + "_real").c_str());
        fields[0] ->SetNumberOfComponents(field_to_display.getNumberOfComponents());
        fields[1] = vtkSmartPointer<vtkDoubleArray>::New();
        fields[1] ->SetName((field_to_display.getName() + "_imag").c_str());
        fields[1] ->SetNumberOfComponents(field_to_display.getNumberOfComponents());

        for (const auto &f : field_to_display.getField()) {
            double f_real[3] = {f[0].real(), f[1].real(), f[2].real()};
            double f_imag[3] = {f[0].imag(), f[1].imag(), f[2].imag()};
            fields[0] ->InsertNextTuple(f_real);
            fields[1] ->InsertNextTuple(f_imag);
        }
    }

    [[nodiscard]] const EMW::Containers::array<vtkPointer, n_fields>& getFields() const { return fields; }
    void dumpField(const vtkSmartPointer<vtkUnstructuredGrid> &unstructuredGrid) const {
        for (auto&& field : fields)
            unstructuredGrid->GetCellData()->AddArray(field);
        std::cout << "Field " << name << " dumped" << std::endl;
    }
};

vtkSmartPointer<vtkUnstructuredGrid> formUnstructuredGrid(const EMW::Mesh::SurfaceMesh &mesh);

};

void surface_snapshot(const EMW::Mesh::SurfaceMesh &mesh, const std::string &part_to_file);

void field_snapshot(const EMW::Math::SurfaceVectorField &field, const std::string &part_to_file);

template<typename ScalarField>
void united_snapshot(const std::vector<ScalarField> &scalarFields,
                     const std::vector<EMW::Math::SurfaceVectorField> &vectorFields,
                     const EMW::Mesh::SurfaceMesh &mesh,
                     const std::string &path_to_file, int number = 0) {
    // Создаем поверхность
    vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid = detail::formUnstructuredGrid(mesh);

    // Создаем поля
    for (auto&& field : vectorFields) {
        const auto dumper = detail::VTKFieldsInitialiser<EMW::Math::SurfaceVectorField>(field);
        dumper.dumpField(unstructuredGrid);
    }

    for (const auto &field : scalarFields) {
        const auto dumper = detail::VTKFieldsInitialiser<ScalarField>(field);
        dumper.dumpField(unstructuredGrid);
    }

    // Создаём снапшот в файле с заданным именем
    std::string fileName = mesh.getName() + (number ? std::to_string(number) : std::string("")) + ".vtu";
    vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
    writer->SetFileName((path_to_file + fileName).c_str());
    writer->SetInputData(unstructuredGrid);
    writer->Write();
}

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
        detail::VTKFieldsInitialiser(fields.get_electric_fields()[i]).dumpField(pointer_to_mesh);
        multiBlockDataSet->SetBlock(i, pointer_to_mesh);
    }
    // магнитные
    for (int i = 0; i != fields.get_magnetic_fields().size(); ++i) {
        const auto pointer_to_mesh = detail::formUnstructuredGrid(fields.get_magnetic_fields()[i].getManifold());
        detail::VTKFieldsInitialiser(fields.get_magnetic_fields()[i]).dumpField(pointer_to_mesh);
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
