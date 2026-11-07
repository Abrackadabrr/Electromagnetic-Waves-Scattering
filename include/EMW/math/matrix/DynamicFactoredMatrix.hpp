//
// Created by evgen on 20.03.2025.
//

#ifndef DYMANICFACTOREDMATRIX_HPP
#define DYMANICFACTOREDMATRIX_HPP

#include "MatrixTraits.hpp"
#include "types/Types.hpp"

#include <algorithm>
#include <cassert>
#include <iostream>
#include <limits>

namespace EMW::Math::LinAgl::Matrix
{
    template <typename factor_t>
    class DynamicFactoredMatrix
    {
        using factor_traits = MatrixTraits<factor_t>;

        Containers::vector<factor_t> factors_;
        Containers::vector<bool> transposed_;
        size_t rows_, cols_;

    public:
        DynamicFactoredMatrix() = default;

        DynamicFactoredMatrix(size_t rows, size_t cols): factors_{}, transposed_{}, rows_(rows), cols_(cols)
        {
        };

        DynamicFactoredMatrix(Containers::vector<factor_t>&& factors)
            : factors_(std::move(factors)), transposed_(factors_.size(), false), rows_(factors_[0].rows()),
              cols_(factors_.back().cols())
        {
        };

        DynamicFactoredMatrix(Containers::vector<factor_t>&& factors, Containers::vector<bool>&& transposed)
            : factors_(std::move(factors)), transposed_(std::move(transposed)),
              rows_(transposed_.front() ? factors_[0].cols() : factors_[0].rows()),
              cols_(transposed_.back() ? factors_.back().rows() : factors_.back().cols())
        {
        };

        DynamicFactoredMatrix(const Containers::vector<factor_t>& factors)
            : factors_(factors), transposed_(factors_.size(), false), rows_(factors_[0].rows()),
              cols_(factors_.back().cols())
        {
            std::cout << "FactoredMatrix constructor with copying" << std::endl;
        };

        template <typename vector_t>
        vector_t matvec(const vector_t& vector) const;

        void matvec(Eigen::Ref<const typename factor_traits::vector_t> vector,
                    Eigen::Ref<typename factor_traits::vector_t> dest) const;

        void matvec(typename factor_traits::element_t* vector, size_t vector_size,
                    Eigen::Ref<typename factor_traits::vector_t> dest) const;

        void matvec(typename factor_traits::element_t* vector, size_t vector_size,
                    typename factor_traits::element_t* dest, size_t dest_size) const;

        void matmul_wise(typename factor_traits::element_t* rhs, size_t rhs_rows, size_t rhs_cols,
                         typename factor_traits::element_t* dest, size_t dest_rows, size_t dest_cols) const;
        void matmul_wise(Eigen::Ref<const typename factor_traits::matrix_t> rhs,
                         Eigen::Ref<typename factor_traits::matrix_t> dest) const;

        /** Number of scalar entries required in each external scratch buffer by
         * matmul_wise_noalloc(). The buffers allow a nested Toeplitz owner to
         * allocate all ACA intermediates when that matrix is constructed.
         */
        [[nodiscard]] size_t matmul_wise_workspace_size(size_t rhs_cols) const noexcept;

        /** dest += this * rhs without allocating intermediate factor products.
         * workspace_first and workspace_second must each contain at least
         * matmul_wise_workspace_size(rhs.cols()) scalar entries.
         */
        void matmul_wise_noalloc(Eigen::Ref<const typename factor_traits::matrix_t> rhs,
                                 Eigen::Ref<typename factor_traits::matrix_t> dest,
                                 typename factor_traits::element_t* workspace_first,
                                 typename factor_traits::element_t* workspace_second,
                                 size_t workspace_size) const noexcept;

        void matmull(Eigen::Map<factor_t> mat, Eigen::Map<factor_t> dest,
                     typename factor_traits::element_t* additional_space) const noexcept;

        typename factor_traits::production_t compute() const noexcept;

        template <Types::index I>
        const factor_t& get() const;

        typename factor_traits::vector_t diagonal() const;

        // ---- Selectors ---- //
        inline Types::index factor_number() const { return factors_.size(); };
        inline Types::scalar memory_usage() const;

        inline Types::index rows() const { return rows_; }

        inline Types::index cols() const { return cols_; }

        // --- Aux methods --- //
        decltype(auto) to_dense() const noexcept { return compute(); };

        // --- Доступ к элементам матрицы --- //
    };

    template <typename factor_t>
    template <typename vector_t>
    vector_t DynamicFactoredMatrix<factor_t>::matvec(const vector_t& vector) const
    {
        if (factor_number() == 0)
        {
            return vector_t::Zero(rows_);
        }

        vector_t result = vector;
        for (Types::integer i = static_cast<Types::integer>(factor_number()) - 1; i >= 0; --i)
        {
            if (transposed_[i])
                result = factors_[i].transpose() * result;
            else
                result = factors_[i] * result;
        }
        return result;
    }

    template <typename factor_t>
    void DynamicFactoredMatrix<factor_t>::matvec(
        Eigen::Ref<const typename factor_traits::vector_t> vector,
        Eigen::Ref<typename factor_traits::vector_t> dest) const
    {
        assert(static_cast<size_t>(vector.size()) == cols_);
        if (factor_number() == 0)
        {
            return;
        }

        typename factor_traits::vector_t result;
        typename factor_traits::vector_t additional_workspace;
        if (transposed_.back())
            result.noalias() = factors_.back().transpose() * vector;
        else
            result.noalias() = factors_.back() * vector;

        for (Types::integer i = static_cast<Types::integer>(factor_number()) - 2; i >= 0; --i)
        {
            if (transposed_[i])
                additional_workspace.noalias() = factors_[i].transpose() * result;
            else
                additional_workspace.noalias() = factors_[i] * result;
            result.noalias() = additional_workspace;
        }
        dest.noalias() += result;
    }

    template <typename factor_t>
    void DynamicFactoredMatrix<factor_t>::matvec(
        typename factor_traits::element_t* vector, size_t vector_size,
        Eigen::Ref<typename factor_traits::vector_t> dest) const
    {
        assert(vector != nullptr);
        assert(vector_size == cols_);
        Eigen::Map<const typename factor_traits::vector_t> mapped_vector(vector,
                                                                         static_cast<Eigen::Index>(vector_size));
        matvec(mapped_vector, dest);
    }

    template <typename factor_t>
    void DynamicFactoredMatrix<factor_t>::matvec(typename factor_traits::element_t* vector, size_t vector_size,
                                                 typename factor_traits::element_t* dest, size_t dest_size) const
    {
        assert(vector != nullptr);
        assert(vector_size == cols_);
        Eigen::Map<typename factor_traits::vector_t> mapped_vector(vector,
                                                                   static_cast<Eigen::Index>(vector_size));
        Eigen::Map<typename factor_traits::vector_t> dest_vector(dest,
                                                                 static_cast<Eigen::Index>(dest_size));
        matvec(mapped_vector, dest_vector);
    }

    template <typename factor_t>
    void DynamicFactoredMatrix<factor_t>::matmul_wise(typename factor_traits::element_t* rhs, size_t rhs_rows,
                                                      size_t rhs_cols, typename factor_traits::element_t* dest,
                                                      size_t dest_rows, size_t dest_cols) const
    {
        assert(rhs != nullptr);
        assert(dest != nullptr);
        assert(rhs_rows == cols_);
        assert(dest_rows == rows_);
        assert(rhs_cols == dest_cols);

        Eigen::Map<const typename factor_traits::matrix_t> rhs_matrix(rhs,
                                                                      static_cast<Eigen::Index>(rhs_rows),
                                                                      static_cast<Eigen::Index>(rhs_cols));
        Eigen::Map<typename factor_traits::matrix_t> dest_matrix(dest,
                                                                 static_cast<Eigen::Index>(dest_rows),
                                                                 static_cast<Eigen::Index>(dest_cols));
        matmul_wise(rhs_matrix, dest_matrix);
    }

    template <typename factor_t>
    void DynamicFactoredMatrix<factor_t>::matmul_wise(
        Eigen::Ref<const typename factor_traits::matrix_t> rhs,
        Eigen::Ref<typename factor_traits::matrix_t> dest) const
    {
        assert(rhs.rows() == static_cast<Eigen::Index>(cols_));
        assert(dest.rows() == static_cast<Eigen::Index>(rows_));
        assert(rhs.cols() == dest.cols());

        if (factor_number() == 0)
        {
            return;
        }

        typename factor_traits::matrix_t result;
        typename factor_traits::matrix_t additional;
        if (transposed_.back())
            result.noalias() = factors_.back().transpose() * rhs;
        else
            result.noalias() = factors_.back() * rhs;

        for (Types::integer i = static_cast<Types::integer>(factor_number()) - 2; i >= 0; --i)
        {
            if (transposed_[i])
                additional.noalias() = factors_[i].transpose() * result;
            else
                additional.noalias() = factors_[i] * result;
            result.noalias() = additional;
        }
        dest.noalias() += result;
    }

    template <typename factor_t>
    size_t DynamicFactoredMatrix<factor_t>::matmul_wise_workspace_size(size_t rhs_cols) const noexcept
    {
        if (factor_number() < 2 || rhs_cols == 0)
        {
            return 0;
        }

        size_t maximum_rows = 0;
        // Factor zero writes directly to dest. Only products of the remaining
        // factors have to live in scratch storage.
        for (size_t i = 1; i < factors_.size(); ++i)
        {
            const size_t output_rows = static_cast<size_t>(
                transposed_[i] ? factors_[i].cols() : factors_[i].rows());
            maximum_rows = std::max(maximum_rows, output_rows);
        }
        if (maximum_rows > std::numeric_limits<size_t>::max() / rhs_cols)
        {
            std::terminate();
        }
        return maximum_rows * rhs_cols;
    }

    template <typename factor_t>
    void DynamicFactoredMatrix<factor_t>::matmul_wise_noalloc(
        Eigen::Ref<const typename factor_traits::matrix_t> rhs,
        Eigen::Ref<typename factor_traits::matrix_t> dest,
        typename factor_traits::element_t* workspace_first,
        typename factor_traits::element_t* workspace_second,
        size_t workspace_size) const noexcept
    {
        assert(rhs.rows() == static_cast<Eigen::Index>(cols_));
        assert(dest.rows() == static_cast<Eigen::Index>(rows_));
        assert(rhs.cols() == dest.cols());

        if (factor_number() == 0)
        {
            return;
        }
        if (factor_number() == 1)
        {
            if (transposed_[0])
                dest.noalias() += factors_[0].transpose() * rhs;
            else
                dest.noalias() += factors_[0] * rhs;
            return;
        }

        const size_t required = matmul_wise_workspace_size(static_cast<size_t>(rhs.cols()));
        assert(workspace_first != nullptr);
        assert(workspace_second != nullptr || factor_number() == 2);
        assert(workspace_size >= required);

        using matrix_t = typename factor_traits::matrix_t;
        const auto output_rows = [&](size_t i) -> Eigen::Index {
            return transposed_[i] ? factors_[i].cols() : factors_[i].rows();
        };

        const size_t last = factors_.size() - 1;
        Eigen::Map<matrix_t> first(workspace_first, output_rows(last), rhs.cols());
        if (transposed_[last])
            first.noalias() = factors_[last].transpose() * rhs;
        else
            first.noalias() = factors_[last] * rhs;

        bool result_in_first = true;
        for (size_t i = last; i-- > 1;)
        {
            const Eigen::Index previous_rows = output_rows(i + 1);
            const Eigen::Index current_rows = output_rows(i);
            Eigen::Map<const matrix_t> previous(result_in_first ? workspace_first : workspace_second,
                                                 previous_rows, rhs.cols());
            Eigen::Map<matrix_t> current(result_in_first ? workspace_second : workspace_first,
                                         current_rows, rhs.cols());
            if (transposed_[i])
                current.noalias() = factors_[i].transpose() * previous;
            else
                current.noalias() = factors_[i] * previous;
            result_in_first = !result_in_first;
        }

        const Eigen::Index intermediate_rows = output_rows(1);
        Eigen::Map<const matrix_t> intermediate(result_in_first ? workspace_first : workspace_second,
                                                 intermediate_rows, rhs.cols());
        if (transposed_[0])
            dest.noalias() += factors_[0].transpose() * intermediate;
        else
            dest.noalias() += factors_[0] * intermediate;
    }


    template <typename factor_t>
    void DynamicFactoredMatrix<factor_t>::matmull(Eigen::Map<factor_t> mat, Eigen::Map<factor_t> dest,
                                                  typename factor_traits::element_t* additional_space) const noexcept {
        // Учимся умножать на плотную матрицу
        Eigen::Map<factor_t> matmul_result{additional_space, dest.rows(), dest.cols()};

        if (transposed_.back())
            matmul_result.noalias() = factors_.back().transpose() * mat;
        else
            matmul_result.noalias() = factors_.back() * mat;

        for (Types::integer i = static_cast<Types::integer>(factor_number()) - 2; i >= 0; --i)
        {
            if (transposed_[i])
                matmul_result = factors_[i].transpose() * matmul_result;
            else
                matmul_result = factors_[i] * matmul_result;
        }
        dest.noalias() += matmul_result;
    }


    template <typename factor_t>
    typename DynamicFactoredMatrix<factor_t>::factor_traits::vector_t DynamicFactoredMatrix<factor_t>::diagonal() const
    {
        if (factor_number() == 1)
            return factors_.front().diagonal();
        throw std::runtime_error(
            "diagonal() is time-consuming operation while running with more that 1 factor numbers");
    }

    template <typename factor_t>
    typename DynamicFactoredMatrix<factor_t>::factor_traits::production_t
    DynamicFactoredMatrix<factor_t>::compute() const noexcept
    {
        factor_t result{};
        if (transposed_[0])
            result = factors_[0].transpose();
        else
            result = factors_[0];
        for (Types::index i = 1; i < factor_number(); ++i)
        {
            if (transposed_[i])
                result *= factors_[i].transpose();
            else
                result *= factors_[i];
        }
        return result;
    }

    template <typename factor_t>
    template <Types::index I>
    const factor_t& DynamicFactoredMatrix<factor_t>::get() const
    {
        return factors_[I];
    }

    // --- Selectors --- //
    template <typename factor_t>
    Types::scalar DynamicFactoredMatrix<factor_t>::memory_usage() const
    {
        Types::index result = 0;
        for (auto& factor : factors_)
        {
            result += factor.size();
        }
        return result;
    }

    template <typename factor_t, typename vector_t>
    vector_t operator*(const DynamicFactoredMatrix<factor_t>& mat, const vector_t& vector)
    {
        return mat.matvec(vector);
    }
} // namespace EMW::Math::LinAgl::Matrix

#endif // DYMANICFACTOREDMATRIX_HPP
