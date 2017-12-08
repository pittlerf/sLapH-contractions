#pragma once

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

namespace LapH {

/*! typetrait class which allows to use QuarklineQ1Indices for Q1 and
 *  QuarklineQ2Indices for Q2L and Q2V
 */
template <QuarkLineType qlt>
struct QuarkLineIndices {};

/*! @todo QuarkLineType is a bad name in this case. That's a proxy for
 *        CorrInfo.lookup
 */
template <>
struct QuarkLineIndices<QuarkLineType::Q0> {
  typedef std::vector<VdaggerVRandomLookup> type;
};

template <>
struct QuarkLineIndices<QuarkLineType::Q1> {
  typedef std::vector<QuarklineQ1Indices> type;
};

template <>
struct QuarkLineIndices<QuarkLineType::Q2L> {
  typedef std::vector<QuarklineQ2Indices> type;
};

template <>
struct QuarkLineIndices<QuarkLineType::Q2V> {
  typedef std::vector<QuarklineQ2Indices> type;
};

template <QuarkLineType qlt>
class QuarkLineBlock2 {
 public:
  QuarkLineBlock2(const size_t dilT,
                  const size_t dilE,
                  const size_t nev,
                  const typename QuarkLineIndices<qlt>::type &quarkline_indices,
                  const std::vector<RandomIndexCombinationsQ2> &ric_lookup);

  Eigen::MatrixXcd const &operator()(const int t,
                                     const int b,
                                     const int op_id,
                                     const int rnd) const {
    auto const id =
        std::find(Ql_id.begin(), Ql_id.end(), std::pair<int, int>(t, b)) - Ql_id.begin();
    /*! @todo catch when t,b is an invalid index */
    return Ql[id][op_id].at(rnd);
  }

  std::vector<Eigen::MatrixXcd> const &operator()(const int t,
                                                  const int b,
                                                  const int op_id) const {
    auto const id =
        std::find(Ql_id.begin(), Ql_id.end(), std::pair<int, int>(t, b)) - Ql_id.begin();
    /*! @todo catch when t,b is an invalid index */
    return Ql[id][op_id];
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

    `vector<vector>` rather than `multiarray`, because `nb_rnd` depends on the operator.

    Indices:

    1. Time slice
    2. Operator
    3. Random Vector
    */
  std::vector<std::vector<std::vector<Eigen::MatrixXcd>>> Ql;

  // using MyQuarkLine = Eigen::MatrixXcd;
  // using MyQuarkLineRandomVec = std::vector<MyQuarkLine>;
  // using MyQuarkLineRandomVecOperator = std::vector<MyQuarkLineRandomVec>;

  // std::vector<std::pair<std::pair<int, int>, MyQuarkLineRandomVecOperator>> Ql;

  boost::circular_buffer<std::pair<int, int>> Ql_id;

  const size_t dilT, dilE, nev;

  static int constexpr dilD = 4;
};

}  // end of namespace
