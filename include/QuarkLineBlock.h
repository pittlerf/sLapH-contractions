#pragma once

#include "OperatorsForMesons.h"
#include "Perambulator.h"
#include "Quarklines.h"
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

  /*! @todo Rewrite that with bracket overload */
  Eigen::MatrixXcd const& return_Ql(const size_t t,
                                    const size_t op_id,
                                    const size_t rnd) const {
    return Ql[t][op_id].at(rnd);
  }

  // ----------------- INTERFACE FOR BUILDING QUARKLINES -----------------------
  // ---------------------------------------------------------------------------
  void build_Q1_one_t(const Perambulator& peram,
                      const OperatorsForMesons& meson_operator,
                      size_t pos,
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

 private:
  // containers for the three types of quark lines
  // vector<vector> rather than multiarray, because nb_rnd depends on the operator
  boost::circular_buffer< std::vector< std::vector<Eigen::MatrixXcd> > > Ql;
  boost::circular_buffer< std::pair<int,int> > Ql_id;

  const size_t dilT, dilE, nev;
  std::vector<LapH::gamma_lookup> gamma;
};

}  // end of namespace
