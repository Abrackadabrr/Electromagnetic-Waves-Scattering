//
// Created by evgen on 15.02.24.
//

#include "slae_generation/MatrixGeneration.hpp"
#include "meshes/plate/PlateGrid.hpp"
#include "experiment/ESA.hpp"
#include "math/MathConstants.hpp"
#include "experiment/PhysicalCondition.hpp"

#include <Eigen/Core>
#include <Eigen/IterativeLinearSolvers>
#include <unsupported/Eigen/IterativeSolvers>
#include <iostream>
#include <fstream>

using namespace EMW;
using namespace EMW::Types;

template<typename Container>
void to_csv(const Container &cont1, const Container &cont2, std::ofstream &str) {
    str << "sigma,alpha\n";
    for (int i = 0; i < cont1.size(); i++) {
        str << cont1[i] << ',' << cont2[i] << '\n';
    }
}

template<typename Range1, typename Range2, typename OutputIterator>
void cartesian_productXY(Range1 const &r1, Range2 const &r2, OutputIterator out, Types::index N,
                         Types::scalar h) {
    using std::begin;
    using std::end;

    for (auto i = begin(r1); i != end(r1); ++i) {
        for (auto j = begin(r2); j != end(r2); ++j) {

            *out++ = Types::Vector3d{static_cast<Types::scalar>(*i) - static_cast<Types::scalar>(N / 2.),
                                     static_cast<Types::scalar>(*j) - static_cast<Types::scalar>(N / 2.),
                                     0} * h;

        }
    }
}

int main() {
    int N = 51;
    scalar h = 1. / (N-1);
    // сетка
    auto surfaceMesh = Mesh::SurfaceMesh{EMW::Examples::Plate::generatePlatePrimaryMesh(N, h)};

    surfaceMesh.setName("surface_mesh_" + std::to_string(N));

    // физика
    EMW::Physics::planeWaveCase physics{
            Vector3d{0, 0, 1}.normalized(),
            complex_d{4 * Math::Constants::PI<scalar>(), 0},
            Vector3d{-1, 0, 0}.normalized()
    };
    const auto initial_field_function = [physics](const Mesh::point_t &point) -> Vector3c {
        return physics.value(point);
    };
    // след падающего поля на расчетной поверхности
    const Math::SurfaceVectorField incidentField(surfaceMesh, initial_field_function);

    // на мелкой
    const VectorXc b = -incidentField.asVector();
    const MatrixXc A = Matrix::getMatrixK(physics.k, surfaceMesh);
//    std::cout << A3.operatorNorm() << std::endl << std::endl;
//    std::cout << newMatrix << std::endl << std::endl;
//    std::cout << newMatrix.inverse().norm() << std::endl;
    std::cout << "Matrix has been built" << std::endl;
//    std::cout << "Matrix conditional value: " << (A3.operatorNorm()) * (A3.inverse().operatorNorm()) << std::endl;
//    std::cout << "Matrix result value: " << newMatrix.operatorNorm() * (newMatrix.inverse().operatorNorm()) << std::endl;

    auto method = Eigen::GMRES<MatrixXc>{};
    method.setMaxIterations(500);
    std::cout << method.maxIterations() << std::endl;
    method.setTolerance(1e-2);
    method.set_restart(500);
    method.compute(A);
    const auto j = VectorXc{method.solveWithGuess(b, b)};
    std::cout << "total iterations: " << method.iterations() << std::endl;
    std::cout << "total error: " << method.error() << std::endl;
    std::cout << "Info: " << static_cast<int>(method.info()) << std::endl;

    auto j_field = Math::SurfaceVectorField::TangentField(surfaceMesh, j);
    j_field.setName("j");

    int samples = 360;
    Containers::vector_d esas;
    esas.reserve(samples);
    Containers::vector_d angles;
    angles.reserve(samples);

    for (int i = 0; i < samples; i++) {
        Types::scalar angle = i * Math::Constants::PI<scalar>() * 2 / samples;
        Types::Vector3d tau = {std::cos(angle), std::sin(angle), 0};
        esas.push_back(ESA::calculateESA(tau, physics.k, j_field));
        angles.push_back(angle);
    }
    std::ofstream sigma(
            "/home/evgen/Education/MasterDegree/thesis/Electromagnetic-Waves-Scattering/vtk_files/examples/plane/sigmas/sigma_" +
            std::to_string(N) + ".csv");

    to_csv(esas, angles, sigma);
}
