#include "DilutedProductFactory.hpp"

void DilutedProductFactoryQ0Q2::build(Key const &time_key,
                                      std::array<ssize_t, 2> const &key) {
  // Extract time keys.
  int constexpr nt1 = 1;
  int constexpr nt2 = 3;
  std::array<int, nt1> time_key1;
  std::array<int, nt2> time_key2;
  std::copy_n(std::begin(time_key) + 0, nt1, std::begin(time_key1));
  std::copy_n(std::begin(time_key) + nt1, nt2, std::begin(time_key2));

  multiply<1, 1>(
      Q0Q2_[time_key], key, factory_q0_[time_key1], factory_q2_[time_key2]);
}
