#pragma once

/**
 * Kahan summation implementation.
 *
 * The algorithm is taken from
 * https://en.wikipedia.org/wiki/Kahan_summation_algorithm
 */
class KahanAccumulator {
 public:
  double value() const { return sum_; }

  KahanAccumulator &operator+=(double const right) {
    double const y = right - c_;
    double const t = sum_ + y;
    c_ = (t - sum_) - y;
    sum_ = t;

    return *this;
  }

 private:
  double sum_ = 0.0;
  double c_ = 0.0;
};
