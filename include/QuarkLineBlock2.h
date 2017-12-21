#pragma once

#include "Gamma.h"
#include "OperatorsForMesons.h"
#include "Perambulator.h"
#include "dilution-iterator.h"
#include "typedefs.h"
#include "DilutedFactor.h"
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
  QuarkLineBlock2(const size_t dilT,
                  const size_t dilE,
                  const size_t nev,
                  const typename QuarkLineIndices<qlt>::type &quarkline_indices,
                  const std::vector<RandomIndexCombinationsQ2> &ric_lookup);

  std::vector<DilutedFactor<0>> const &operator()(const int t,
                                               const int b,
                                               const size_t op_id) const {
    /*! @todo catch when t,b is an invalid index */
    auto const time_key = std::make_pair(t, b);
    typename OperatorToFactorMap<1, 0>::key_type const key{op_id};

    auto it = Ql.find(time_key);
    if (it == Ql.end()) {
      std::cout << "Tried to access " << t << " " << b << std::endl;
      std::cout << "Size of the map: " << Ql.size() << std::endl;
      abort();
    }

    auto const &Ql_elem = Ql.at(time_key);

    auto it2 = Ql_elem.find(key);
    if (it2 == Ql_elem.end()) {
      std::cout << "Tried to access " << to_string<1, 0>(key) << std::endl;
      std::cout << "Size of the map: " << Ql_elem.size() << std::endl;
      print(Ql_elem);
      abort();
    }

    return Ql_elem.at(key);
  }

  OperatorToFactorMap<1, 0> const &operator()(const int t, const int b) const {
    /*! @todo catch when t,b is an invalid index */
    auto const time_key = std::make_pair(t, b);

    auto it = Ql.find(time_key);
    if (it == Ql.end()) {
      std::cout << "Tried to access " << t << " " << b << std::endl;
      std::cout << "Size of the map: " << Ql.size() << std::endl;
      abort();
    }

    return Ql.at(time_key);
  }

  void clear() { Ql.clear(); }

  // ----------------- INTERFACE FOR BUILDING QUARKLINES -----------------------
  // ---------------------------------------------------------------------------
  void build_Q1_one_t(const Perambulator &peram,
                      const OperatorsForMesons &meson_operator,
                      const int t_source,
                      const int t_sink,
                      const typename QuarkLineIndices<qlt>::type &ql_lookup,
                      const std::vector<RandomIndexCombinationsQ2> &ric_lookup);

  void build(const Perambulator &peram,
             const OperatorsForMesons &meson_operator,
             const int t_source,
             const int t_sink,
             const typename QuarkLineIndices<qlt>::type &ql_lookup,
             const std::vector<RandomIndexCombinationsQ2> &ric_lookup);

  void build_block_pair(Perambulator const &peram,
                        OperatorsForMesons const &meson_operator,
                        DilutionIterator const &block_pair,
                        typename QuarkLineIndices<qlt>::type const &ql_lookup,
                        std::vector<RandomIndexCombinationsQ2> const &ric_lookup);
  // Overload for Q0
  void build_block_pair(RandomVector const &rnd_vec,
                        OperatorsForMesons const &meson_operator,
                        DilutionIterator const &block_pair,
                        typename QuarkLineIndices<qlt>::type const &ql_lookup,
                        std::vector<RandomIndexCombinationsQ2> const &ric_lookup);

 private:
  /*!
    Containers for the three types of quark lines.

    Indices:

    1. Time slice
    */
  std::map<std::pair<int, int>, OperatorToFactorMap<1, 0>> Ql;

  const size_t dilT, dilE, nev;

  static int constexpr dilD = 4;
};
