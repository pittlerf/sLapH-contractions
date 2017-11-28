#pragma once

#include "OperatorsForMesons.h"
#include "Perambulator.h"
#include "Quarklines.h"
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

template <QuarkLineType qlt>
class QuarkLineBlock {
 public:
  QuarkLineBlock(const size_t dilT,
                 const size_t dilE,
                 const size_t nev,
                 const typename QuarkLineIndices<qlt>::type& quarkline_indices,
                 const std::vector<RandomIndexCombinationsQ2>& ric_lookup);

  ~QuarkLineBlock(){};

  const cmplx& return_gamma_val(const size_t gamma_id, const size_t row) const {
    return gamma[gamma_id].value[row];
  }

  int const& return_gamma_row(const size_t gamma_id, const size_t row) const {
    return gamma[gamma_id].row[row];
  }

  Eigen::MatrixXcd const &operator()(const int t,
                                     const int b,
                                     const int op_id,
                                     const int rnd) const {
    auto const id =
        std::find(Ql_id.begin(), Ql_id.end(), std::pair<int, int>(t, b)) - Ql_id.begin();
    /*! @todo catch when t,b is an invalid index */
    return Ql[id][op_id].at(rnd);
  }

  // ----------------- INTERFACE FOR BUILDING QUARKLINES -----------------------
  // ---------------------------------------------------------------------------
  void build_Q1_one_t(const Perambulator& peram,
                      const OperatorsForMesons& meson_operator,
                      const int t_source,
                      const int t_sink,
                      const typename QuarkLineIndices<qlt>::type& ql_lookup,
                      const std::vector<RandomIndexCombinationsQ2>& ric_lookup);

  void build(const Perambulator& peram,
             const OperatorsForMesons& meson_operator,
             const int t_source,
             const int t_sink,
             const typename QuarkLineIndices<qlt>::type& ql_lookup,
             const std::vector<RandomIndexCombinationsQ2>& ric_lookup);

  void build_block_pair(const Perambulator &peram,
                        const OperatorsForMesons &meson_operator,
                        DilutionIterator const &block_pair,
                        const typename QuarkLineIndices<qlt>::type &ql_lookup,
                        const std::vector<RandomIndexCombinationsQ2> &ric_lookup);

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
  std::vector<LapH::gamma_lookup> gamma;

  static int constexpr dilD = 4;
};

}  // end of namespace
