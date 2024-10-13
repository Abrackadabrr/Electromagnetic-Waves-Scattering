//
// Created by evgen on 02.08.24.
//

#include "SurfaceVectorField.hpp"

#include "math/MathConstants.hpp"
#include "math/Productions.hpp"

#include <ranges>

EMW::Types::VectorXc EMW::Math::SurfaceVectorField::asSLAERHS() const {
    const auto &cells = manifold_.getCells();
    const long N = static_cast<long>(cells.size());
    Types::VectorXc result = Types::VectorXc::Zero(2 * N);
    for (int i = 0; i < cells.size(); i++) {
        result(i) = -Math::quasiDot(cells[i].tau[0], field_data_[i]);
        result(i + N) = -Math::quasiDot(cells[i].tau[1], field_data_[i]);
    }
    return result;
}

EMW::Types::scalar EMW::Math::SurfaceVectorField::supNorm() const {
    const auto element = std::max_element(field_data_.begin(), field_data_.end(),
                                          [](const field_t &v1, const field_t &v2) { return v1.norm() < v2.norm(); });
    if (element == field_data_.end()) {
        return -1;
    }
    return element->norm();
}
EMW::Math::SurfaceScalarField
EMW::Math::SurfaceVectorField::fieldNorm(const std::string name = "delfault_field_name") const {
    const auto field_view = field_data_ | std::views::transform([](const field_t &v) -> Types::complex_d {
                                return Types::complex_d{v.real().norm(), v.imag().norm()};
                            });
    SurfaceScalarField result{manifold_, std::vector<Types::complex_d>{field_view.begin(), field_view.end()}};
    result.setName(name);
    return result;
}

EMW::Math::SurfaceVectorField
EMW::Math::SurfaceVectorField::TangentField(const EMW::Math::SurfaceVectorField::manifold_t &manifold,
                                            const EMW::Types::VectorXc &fieldProjections) {
    SurfaceVectorField result(manifold);
    const long N = static_cast<long>(manifold.getCells().size());
    for (auto [i, cell] : manifold.getCells() | std::views::enumerate) {
        result.field_data_.emplace_back(fieldProjections(i) * cell.tau[0] + fieldProjections(i + N) * cell.tau[1]);
    }
    result.initialized = true;
    return result;
}

EMW::Math::SurfaceVectorField EMW::Math::SurfaceVectorField::pointwiseMultiplication(
    const std::function<Types::scalar(const Mesh::IndexedCell &)> &function) const {
    const auto field_data = std::views::iota(0, static_cast<int>(manifold_.getCells().size())) |
                            std::views::transform([&](int i) -> Types::Vector3c {
                                const auto cell = manifold_.getCells()[i];
                                return field_data_[i] * function(cell);
                            });
    return {manifold_, {field_data.begin(), field_data.end()}};
}

EMW::Math::SurfaceVectorField EMW::Math::SurfaceVectorField::pointwiseMultiplication(
    const std::function<Types::scalar(const EMW::Types::Vector3d &)> &function) const {
    const auto field_data = std::views::iota(0, static_cast<int>(manifold_.getCells().size())) |
                            std::views::transform([&](int i) -> Types::Vector3c {
                                const auto cell = manifold_.getCells()[i];
                                return field_data_[i] * function(cell.collPoint_.point_);
                            });
    return {manifold_, {field_data.begin(), field_data.end()}};
}

void EMW::Math::SurfaceVectorField::multiply(
    const std::function<Types::scalar(const EMW::Types::Vector3d &)> &function) {
    for (int i = 0; i < static_cast<int>(manifold_.getCells().size()); i++) {
        const auto cell = manifold_.getCells()[i];
        field_data_[i] = field_data_[i] * function(cell.collPoint_.point_);
    }
}

void EMW::Math::SurfaceVectorField::multiply(
    const std::function<Types::scalar(const EMW::Mesh::IndexedCell &)> &function) {
    for (int i = 0; i < static_cast<int>(manifold_.getCells().size()); i++) {
        const auto cell = manifold_.getCells()[i];
        field_data_[i] = field_data_[i] * function(cell);
    }
}

EMW::Math::SurfaceVectorField
EMW::Math::SurfaceVectorField::ZeroField(const EMW::Math::SurfaceVectorField::manifold_t &manifold) {
    SurfaceVectorField result(manifold);
    std::fill(result.field_data_.begin(), result.field_data_.end(), Types::Vector3c::Zero());
    result.initialized = true;
    return result;
}

EMW::Math::SurfaceVectorField EMW::Math::SurfaceVectorField::surfaceProjection() const {
    const auto field_data =
        std::views::iota(0, static_cast<int>(manifold_.getCells().size())) | std::views::transform([&](int i) {
            const auto cell = manifold_.getCells()[i];
            const Types::complex_d c0 = Math::quasiDot(field_data_[i], cell.tau[0]);
            const Types::complex_d c1 = Math::quasiDot(field_data_[i], cell.tau[1]);
            return c0 * cell.tau[0] + c1 * cell.tau[1];
        });
    return {manifold_, {field_data.begin(), field_data.end()}};
}

EMW::Math::SurfaceVectorField EMW::Math::SurfaceVectorField::crossWithNormalField() const {
    const auto field_data = std::views::iota(0, static_cast<int>(manifold_.getCells().size())) |
                            std::views::transform([&](int k) -> Types::Vector3c {
                                const auto cell = manifold_.getCells()[k];
                                const Types::Vector3d real = field_data_[k].real().cross(cell.normal);
                                const Types::Vector3d imag = field_data_[k].imag().cross(cell.normal);
                                return (real + Math::Constants::i * imag);
                            });
    return {manifold_, {field_data.begin(), field_data.end()}};
}

EMW::Math::SurfaceVectorField EMW::Math::SurfaceVectorField::normalCrossField() const {
    const auto field_data = std::views::iota(0, static_cast<int>(manifold_.getCells().size())) |
                            std::views::transform([&](int k) -> Types::Vector3c {
                                const auto cell = manifold_.getCells()[k];
                                const Types::Vector3d real = cell.normal.cross(field_data_[k].real());
                                const Types::Vector3d imag = cell.normal.cross(field_data_[k].imag());;
                                return (real + Math::Constants::i * imag);
                            });
    return {manifold_, {field_data.begin(), field_data.end()}};
}

EMW::Math::SurfaceVectorField
EMW::Math::SurfaceVectorField::NormalField(const EMW::Math::SurfaceVectorField::manifold_t &manifold) {
    const auto normals_view = manifold.getCells() | std::views::transform([](const Mesh::IndexedCell &cell) -> field_t {
                                  return cell.normal;
                              });
    return {manifold, {normals_view.begin(), normals_view.end()}};
}

EMW::Math::SurfaceVectorField EMW::Math::operator-(const EMW::Math::SurfaceVectorField &lhs,
                                                   const EMW::Math::SurfaceVectorField &rhs) {
    Containers::vector<Types::Vector3c> field_data{};
    field_data.reserve(lhs.getField().size());
    for (int i = 0; i != lhs.getField().size(); i++) {
        field_data.emplace_back(lhs.getField()[i] - rhs.getField()[i]);
    }
    return {lhs.getManifold(), field_data};
}

EMW::Math::SurfaceVectorField EMW::Math::operator+(const EMW::Math::SurfaceVectorField &lhs,
                                                   const EMW::Math::SurfaceVectorField &rhs) {
    Containers::vector<Types::Vector3c> field_data{};
    field_data.reserve(lhs.getField().size());
    for (int i = 0; i != lhs.getField().size(); i++) {
        field_data.emplace_back(lhs.getField()[i] + rhs.getField()[i]);
    }
    return {lhs.getManifold(), field_data};
}
