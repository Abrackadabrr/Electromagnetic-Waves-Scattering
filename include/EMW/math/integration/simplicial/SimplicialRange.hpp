//
// Created by evgen on 06.05.2026.
//

#ifndef ELECTROMAGNETIC_WAVES_SCATTERING_SIMPLICIALRANGE_HPP
#define ELECTROMAGNETIC_WAVES_SCATTERING_SIMPLICIALRANGE_HPP

namespace EMW::Math::Integration::Numerical::Simplicial::QuadratureUtils {

#include <cstddef>
#include <iterator>

template <typename point_t> class TriangleRange {
  public:
    struct Triangle {
        point_t a, b, c;

        operator Containers::array<point_t, 3>() const { return Containers::array<point_t, 3>{a, b, c}; }
    };

    TriangleRange(const Triangle &root, size_t level) : root_(root), level_(level) {}

    TriangleRange(const Containers::array<point_t, 3> &root, size_t level)
        : TriangleRange(Triangle{root[0], root[1], root[2]}, level) {}

    TriangleRange(const point_t &a, const point_t &b, const point_t &c, size_t level)
        : TriangleRange(Triangle{a, b, c}, level) {}

    auto begin() const { return Iterator(root_, level_, false); }

    auto end() const { return Iterator(root_, level_, true); }

  private:
    class Iterator {
      public:
        using iterator_category = std::forward_iterator_tag;
        using value_type = Triangle;
        using difference_type = std::ptrdiff_t;
        using pointer = void;
        using reference = Triangle;

        Iterator(const Triangle &root, size_t level, bool end)
            : root_(root), N_(size_t{1} << level), phase_(Phase::Upper), i_(0), j_(0), k_(N_), index_(0),
              total_count_(N_ * N_) {
            if (end) {
                // end -- итератор на треугольник "за последним"
                index_ = total_count_;
            }
        }

        Triangle operator*() const { return phase_ == Phase::Upper ? upper_triangle() : lower_triangle(); }

        Iterator &operator++() {
            ++index_;

            if (index_ == total_count_) {
                // после end operator++ работает некорректно
                return *this;
            }

            if (phase_ == Phase::Upper) {
                advance_upper();
            } else {
                advance_lower();
            }

            return *this;
        }

        Iterator operator++(int) {
            Iterator tmp = *this;
            ++(*this);
            return tmp;
        }

        bool operator==(const Iterator &other) const { return index_ == other.index_; }

        bool operator!=(const Iterator &other) const { return !(*this == other); }

      private:
        enum class Phase { Upper, Lower };

        [[nodiscard]] point_t barycentric(size_t l0, size_t l1, size_t l2) const {
            const double invN = 1.0 / static_cast<double>(N_);

            return root_.a * (static_cast<double>(l0) * invN) + root_.b * (static_cast<double>(l1) * invN) +
                   root_.c * (static_cast<double>(l2) * invN);
        }

        [[nodiscard]] Triangle upper_triangle() const {
            return Triangle{barycentric(i_, j_, k_), barycentric(i_ + 1, j_, k_ - 1), barycentric(i_, j_ + 1, k_ - 1)};
        }

        [[nodiscard]] Triangle lower_triangle() const {
            return Triangle{barycentric(i_ + 1, j_, k_), barycentric(i_ + 1, j_ + 1, k_ - 1),
                            barycentric(i_, j_ + 1, k_)};
        }

        void advance_upper() {
            if (j_ < N_ - i_ - 1) {
                ++j_;
                --k_;
            } else {
                ++i_;
                j_ = 0;

                if (i_ < N_) {
                    k_ = N_ - i_;
                } else {
                    phase_ = Phase::Lower;
                    i_ = 0;
                    j_ = 0;
                    k_ = N_ - 1;
                }
            }
        }

        void advance_lower() {
            if (j_ < N_ - i_ - 2) {
                ++j_;
                --k_;
            } else {
                ++i_;
                j_ = 0;
                k_ = N_ - 1 - i_;
            }
        }

        Triangle root_;
        size_t N_{};

        Phase phase_{Phase::Upper};

        size_t i_{};
        size_t j_{};
        size_t k_{};

        size_t index_{};
        size_t total_count_{};
    };

    Triangle root_;
    size_t level_;
};

// deduction hint
template <class point_t>
TriangleRange(typename TriangleRange<point_t>::Triangle, size_t) -> TriangleRange<point_t>;

} // namespace EMW::Math::Integration::Numerical::Simplicial::QuadratureUtils

#endif // ELECTROMAGNETIC_WAVES_SCATTERING_SIMPLICIALRANGE_HPP
