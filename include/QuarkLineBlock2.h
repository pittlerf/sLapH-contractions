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

namespace LapH {


template <QuarkLineType qlt>
class QuarkLineBlock2 {
 public:
  QuarkLineBlock2(const size_t dilT,
                  const size_t dilE,
                  const size_t nev,
                  const typename QuarkLineIndices<qlt>::type &quarkline_indices,
                  const std::vector<RandomIndexCombinationsQ2> &ric_lookup);

  std::vector<DilutedFactor> const &operator()(const int t,
                                                  const int b,
                                                  const int op_id) const {
    auto const id =
        std::find(Ql_id.begin(), Ql_id.end(), std::pair<int, int>(t, b)) - Ql_id.begin();
    /*! @todo catch when t,b is an invalid index */
    typename OperatorToFactorMap<1>::key_type const key{op_id};
    return Ql[id].at(key);
  }

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
  std::vector<OperatorToFactorMap<1>> Ql;

  boost::circular_buffer<std::pair<int, int>> Ql_id;

  const size_t dilT, dilE, nev;

  static int constexpr dilD = 4;
};

}  // end of namespace
