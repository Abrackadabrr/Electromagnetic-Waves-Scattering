#include <gtest/gtest.h>

#include "math/matrix/ToeplitzMultilevelNew.hpp"

#include <Eigen/Core>

#include <array>
#include <cmath>
#include <span>

namespace {

using Index = Eigen::Index;
using Dense = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;
using Vector = Eigen::VectorXd;
using Matrix = EMW::Math::LinAgl::Matrix::MultilevelToeplitz<Dense>;

constexpr Index grid_extent = 4;
constexpr Index leaf_extent = 200;

std::pair<double, double> coefficients(std::span<const Index> row, std::span<const Index> col) {
    const double first = static_cast<double>(col[0] - row[0]);
    const double second = static_cast<double>(col[1] - row[1]);
    const double third = static_cast<double>(col[2] - row[2]);
    return {1.0 + 0.031 * first - 0.017 * second + 0.023 * third,
            0.0002 * (2.0 * first + 3.0 * second - third)};
}

Dense make_leaf(std::span<const Index> row, std::span<const Index> col) {
    const auto [diagonal, constant] = coefficients(row, col);
    Dense result = Dense::Constant(leaf_extent, leaf_extent, constant);
    result.diagonal().array() += diagonal;
    return result;
}

Index logical_index(Index first, Index second, Index third) {
    return (first * grid_extent + second) * grid_extent + third;
}

} // namespace

TEST(MultilevelToeplitzNew, Triple4x4x4With200x200Leaves) {
    Eigen::setNbThreads(1);
    Matrix matrix(3, grid_extent, grid_extent, make_leaf);

    ASSERT_EQ(matrix.levels(), 3);
    ASSERT_EQ(matrix.unique_leaf_blocks(), 7 * 7 * 7);
    ASSERT_EQ(matrix.rows(), grid_extent * grid_extent * grid_extent * leaf_extent);
    ASSERT_EQ(matrix.cols(), matrix.rows());
    ASSERT_EQ(matrix.get_storage().packed_coefficients().size(),
              static_cast<std::size_t>(7 * 7 * 7 * leaf_extent * leaf_extent));

    Vector rhs(matrix.cols());
    for (Index index = 0; index < rhs.size(); ++index)
        rhs[index] = std::sin(0.013 * static_cast<double>(index)) +
                     0.25 * std::cos(0.007 * static_cast<double>(index));

    Vector expected = Vector::Zero(matrix.rows());
    for (Index row_first = 0; row_first < grid_extent; ++row_first)
        for (Index row_second = 0; row_second < grid_extent; ++row_second)
            for (Index row_third = 0; row_third < grid_extent; ++row_third) {
                const Index output = logical_index(row_first, row_second, row_third) * leaf_extent;
                auto output_segment = expected.segment(output, leaf_extent);
                for (Index col_first = 0; col_first < grid_extent; ++col_first)
                    for (Index col_second = 0; col_second < grid_extent; ++col_second)
                        for (Index col_third = 0; col_third < grid_extent; ++col_third) {
                            const std::array<Index, 3> row{row_first, row_second, row_third};
                            const std::array<Index, 3> col{col_first, col_second, col_third};
                            const auto [diagonal, constant] = coefficients(row, col);
                            const Index input = logical_index(col_first, col_second, col_third) * leaf_extent;
                            const auto input_segment = rhs.segment(input, leaf_extent);
                            output_segment.array() += diagonal * input_segment.array() +
                                                      constant * input_segment.sum();
                        }
            }

    Matrix::workspace_type workspace;
    Vector result(matrix.rows());
    matrix.matvec(rhs, result, workspace);
    const double relative_error = (result - expected).norm() / expected.norm();
    EXPECT_LE(relative_error, 2e-13);
}
