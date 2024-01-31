//
// Created by evgen on 24.01.24.
//

#include "types/Types.hpp"
#include "PlateGrid.hpp"
#include "mesh/Mesh.hpp"
#include "visualisation/VTKFunctions.hpp"

int main() {
    const auto mesh = EMW::Examples::Plate::generatePlatePrimaryMesh(10, 0.1);
    VTK::test_snapshot(0, mesh,
                       "/media/evgen/SecondLinuxDisk/4_level/Electromagnetic-Waves-Scattering/vtk_files/examples/plate/");
}
