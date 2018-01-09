#pragma once

#include "DilutedFactor.h"
#include "Gamma.h"
#include "OperatorsForMesons.h"
#include "Perambulator.h"
#include "dilution-iterator.h"
#include "typedefs.h"

#include "Eigen/Dense"
#include "boost/circular_buffer.hpp"
#include "boost/multi_array.hpp"

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

template <QuarkLineType qlt>
class QuarkLineBlock2 {
 public:
  using Key = std::array<int, QuarkLineIndices<qlt>::num_times>;
  using Value = OperatorToFactorMap<1, 0>;

  QuarkLineBlock2(RandomVector const &random_vector,
                  Perambulator const &perambulator,
                  OperatorsForMesons const &_meson_operator,
                  size_t const dilT,
                  size_t const dilE,
                  size_t const nev,
                  typename QuarkLineIndices<qlt>::type const &quarkline_indices);

  Value const &operator[](Key const &key) const {
    if (Ql.count(key) == 0) {
      std::cout << "Tried to access an element which does not exist\n";
      std::cout << "Size of the map: " << Ql.size() << std::endl;
      abort();
    }

    return Ql.at(key);
  }

  void clear() { Ql.clear(); }

  void build_Q1_one_t(const int t_source, const int t_sink);

  void build_block_pair(DilutionIterator const &block_pair);

 private:
  std::map<Key, Value> Ql;

  RandomVector const &rnd_vec;
  Perambulator const &peram;
  OperatorsForMesons const &meson_operator;
  const size_t dilT, dilE, nev;
  typename QuarkLineIndices<qlt>::type const &quarkline_indices;

  static int constexpr dilD = 4;
};
