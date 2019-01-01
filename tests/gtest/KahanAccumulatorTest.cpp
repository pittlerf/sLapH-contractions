#include "KahanAccumulator.hpp"

#include <gtest/gtest.h>

#include <cmath>
#include <iostream>

// A simple test where the normal double precision is sufficient.
TEST(KahanAccumulator, simple) {
  KahanAccumulator acc;
  double sum = 0.0;

  for (int i = 1; i != 10; ++i) {
    acc += i;
    sum += i;
  }

  EXPECT_NEAR(acc.value(), sum, 1e-20);
}

// Add a fraction of the double epsilon to 1.0 many times. With the
// conventional addition the result will still say 1.0, with the Kahan
// summation we expect an exact result.
TEST(KahanAccumulator, hard) {
  KahanAccumulator acc;
  double sum = 0.0;

  acc += 1.0;
  sum += 1.0;

  int const iters = 10000;

  double const eps = std::numeric_limits<double>::epsilon();
  double const summand = eps / 100;
  double const exact = 1.0 + iters * summand;

  for (int i = 0; i < iters; ++i) {
    acc += summand;
    sum += summand;
  }

  // One has to be very careful with exact equality on floating point numbers.
  // In this particular case we do want an exact equality, though.
  EXPECT_EQ(acc.value(), exact);
  EXPECT_NE(sum, exact);
}
