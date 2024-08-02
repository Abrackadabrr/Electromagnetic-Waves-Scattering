//
// Created by evgen on 02.08.24.
//

#include "SurfaceField.hpp"
#include "MathConstants.hpp"
#include "Productions.hpp"

#include <ranges>

EMW::Math::SurfaceField::SurfaceField(const EMW::Math::SurfaceField::manifold_t &man_ref,
                                      const std::function<field_t(const Mesh::point_t &)> &function) : manifold_(
        man_ref), initialized(true) {
    field_data_.reserve(man_ref.getCells().size());
    const auto coll_points_data_view = man_ref.getCells() | std::views::transform(
            [](const Mesh::IndexedCell &cell) -> Mesh::point_t {
                return cell.collPoint_.point_;
            });
    // странно, поскольку не могу применить подряд несколько трансформов
    auto points = Containers::vector<Mesh::point_t>{coll_points_data_view.begin(),
                                                    coll_points_data_view.end()};

    const auto filed_data_view = points | std::views::transform(function);
    field_data_ = std::vector<field_t>{filed_data_view.begin(), filed_data_view.end()};
}

EMW::Types::VectorXc EMW::Math::SurfaceField::asSLAERHS() const {
    const auto &cells = manifold_.getCells();
    const long N = static_cast<long>(cells.size());
    Types::VectorXc result = Types::VectorXc::Zero(2 * N);
    for (int i = 0; i < cells.size(); i++) {
        result(i) = -Math::quasiDot(cells[i].tau[0], field_data_[i]);
        result(i + N) = -Math::quasiDot(cells[i].tau[1], field_data_[i]);
    }
    return result;
}

EMW::Math::SurfaceField
EMW::Math::SurfaceField::TangentField(const EMW::Math::SurfaceField::manifold_t &manifold,
                                      const EMW::Types::VectorXc &fieldProjections) {
    SurfaceField result(manifold);
    const long N = static_cast<long>(manifold.getCells().size());
    for (auto [i, cell]: manifold.getCells() | std::views::enumerate) {
        result.field_data_.emplace_back(fieldProjections(i) * cell.tau[0] + fieldProjections(i + N) * cell.tau[1]);
    }
    result.initialized = true;
    return result;
}

EMW::Math::SurfaceField EMW::Math::SurfaceField::ZeroField(const EMW::Math::SurfaceField::manifold_t &manifold) {
    SurfaceField result(manifold);
    std::fill(result.field_data_.begin(), result.field_data_.end(), Types::Vector3c::Zero());
    result.initialized = true;
    return result;
}

EMW::Math::SurfaceField EMW::Math::SurfaceField::surfaceProjection() const {
    const auto field_data =
            std::views::iota(0, static_cast<int>(manifold_.getCells().size())) | std::views::transform([&](int i){
                const auto cell = manifold_.getCells()[i];
                const Types::complex_d c0 = Math::quasiDot(field_data_[i], cell.tau[0]);
                const Types::complex_d c1 = Math::quasiDot(field_data_[i], cell.tau[1]);
                return c0 * cell.tau[0] + c1 * cell.tau[1];
            });
    return {manifold_, {field_data.begin(), field_data.end()}};
}

EMW::Math::SurfaceField EMW::Math::SurfaceField::NormalField(const EMW::Math::SurfaceField::manifold_t &manifold) {
    const auto normals_view = manifold.getCells() | std::views::transform([](const Mesh::IndexedCell &cell) -> field_t {
        return cell.normal;
    });
    return {manifold, {normals_view.begin(), normals_view.end()}};
};

