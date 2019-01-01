#pragma once

/**
 * Kahan summation implementation.
 *
 * The algorithm is taken from
 * https://en.wikipedia.org/wiki/Kahan_summation_algorithm
 *
 * @tparam Numeric A numeric type which supports addition and subtraction. The
 * most sensible values are `double` and `std::complex<double>`.
 */
template <typename Numeric>
class KahanAccumulator {
 public:
  Numeric value() const { return sum_; }

  KahanAccumulator &operator+=(Numeric const right) {
    Numeric const y = right - c_;
    Numeric const t = sum_ + y;
    c_ = (t - sum_) - y;
    sum_ = t;

    return *this;
  }

 private:
  Numeric sum_ = 0.0;
  Numeric c_ = 0.0;
};
