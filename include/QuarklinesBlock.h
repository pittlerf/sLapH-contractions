#ifndef QUARKLINES_H_
#define OUARKLINES_H_

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

#include "boost/multi_array.hpp"
#include "Eigen/Dense"

#include "OperatorsForMesons.h"
#include "Perambulator.h"
#include "typedefs.h"

namespace LapH {

struct gamma_lookup {
  std::array<int, 4> row;
  std::array<cmplx, 4> value;
};

class Quarklines {

private:
  // containers for the three types of quark lines
  array_quarkline Q1;
  array_quarkline Q2V;
  array_quarkline Q2L;
  const size_t Lt, dilT, dilE, nev;
  std::vector<LapH::gamma_lookup>  gamma;

  void build_Q1(const Perambulator& peram,
                const OperatorsForMesons& meson_operator,
                const std::vector<QuarklineQ1Indices>& ql_lookup,
                const std::vector<RandomIndexCombinationsQ2>& ric_lookup);
  void build_Q2V(const Perambulator& peram,
                 const OperatorsForMesons& meson_operator,
                 const std::vector<QuarklineQ2Indices>& ql_lookup,
                 const std::vector<RandomIndexCombinationsQ2>& ric_lookup);
  void build_Q2L(const Perambulator& peram,
                 const OperatorsForMesons& meson_operator,
                 const std::vector<QuarklineQ2Indices>& ql_lookup,
                 const std::vector<RandomIndexCombinationsQ2>& ric_lookup);

public:

  Quarklines (const size_t Lt, const size_t dilT, const size_t dilE, 
              const size_t nev, const QuarklineLookup& quarkline_lookuptable,
              const std::vector<RandomIndexCombinationsQ2>& ric_lookup);
  ~Quarklines () {}; // dtor

  inline const cmplx& return_gamma_val(const size_t gamma_id, 
                                       const size_t row) const {
    return gamma[gamma_id].value[row];
  }
  inline const int& return_gamma_row(const size_t gamma_id, 
                                     const size_t row) const{
    return gamma[gamma_id].row[row];
  }
  inline const Eigen::MatrixXcd& return_Q1(const size_t t1, const size_t t2,
                                   const size_t op_id, const size_t rnd) const {
    return Q1[t1][t2][op_id].at(rnd);
  }
  inline const Eigen::MatrixXcd& return_Q2V(const size_t t1, const size_t t2,
                                   const size_t op_id, const size_t rnd) const {
    return Q2V[t1][t2][op_id].at(rnd);
  }
  inline const Eigen::MatrixXcd& return_Q2L(const size_t t1, const size_t t2,
                                   const size_t op_id, const size_t rnd) const {
    return Q2L[t1][t2][op_id].at(rnd);
  }

  void create_quarklines(const Perambulator& peram, 
                      const OperatorsForMesons& meson_operator,
                      const QuarklineLookup& quarkline_lookuptable,
                      const std::vector<RandomIndexCombinationsQ2>& ric_lookup);
};


// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

//class Quarklines_one_t {
//
//private:
//  // containers for the three types of quark lines
//  boost::multi_array<std::vector<Eigen::MatrixXcd>, 2> Q1;
//  boost::multi_array<std::vector<Eigen::MatrixXcd>, 2> Q2V;
//  boost::multi_array<std::vector<Eigen::MatrixXcd>, 2> Q2L;
//
//  const size_t dilT, dilE, nev;
//  std::vector<LapH::gamma_lookup>  gamma;
//
//public:
//
//  Quarklines_one_t (const size_t dilT, const size_t dilE, 
//              const size_t nev, const QuarklineLookup& quarkline_lookuptable,
//              const std::vector<RandomIndexCombinationsQ2>& ric_lookup);
//  ~Quarklines_one_t () {}; // dtor
//
//  inline const cmplx& return_gamma_val(const size_t gamma_id, 
//                                       const size_t row) const {
//    return gamma[gamma_id].value[row];
//  }
//  inline const int& return_gamma_row(const size_t gamma_id, 
//                                     const size_t row) const{
//    return gamma[gamma_id].row[row];
//  }
//  inline const Eigen::MatrixXcd& return_Q1(const size_t t,
//                                   const size_t op_id, const size_t rnd) const {
//    return Q1[t][op_id].at(rnd);
//  }
//  inline const Eigen::MatrixXcd& return_Q2V(const size_t t,
//                                   const size_t op_id, const size_t rnd) const {
//    return Q2V[t][op_id].at(rnd);
//  }
//  inline const Eigen::MatrixXcd& return_Q2L(const size_t t,
//                                   const size_t op_id, const size_t rnd) const {
//    return Q2L[t][op_id].at(rnd);
//  }
//
//  // ----------------- INTERFACE FOR BUILDING QUARKLINES -----------------------
//  // ---------------------------------------------------------------------------
//  void build_Q1_one_t(const Perambulator& peram,
//                      const OperatorsForMesons& meson_operator, size_t pos,
//                      const int t_source, const int t_sink,
//                      const std::vector<QuarklineQ1Indices>& ql_lookup,
//                      const std::vector<RandomIndexCombinationsQ2>& ric_lookup);
//  void build_Q1_mult_t(const Perambulator& peram,
//                       const OperatorsForMesons& meson_operator,
//                       const int t_source, const int t_sink,
//                       const std::vector<QuarklineQ1Indices>& ql_lookup,
//                       const std::vector<RandomIndexCombinationsQ2>& ric_lookup);
//  void build_Q2V_one_t(const Perambulator& peram,
//                       const OperatorsForMesons& meson_operator,
//                       const int t1_block, const int t2_block,
//                       const std::vector<QuarklineQ2Indices>& ql_lookup,
//                       const std::vector<RandomIndexCombinationsQ2>& ric_lookup);
//  void build_Q2L_one_t(const Perambulator& peram,
//                       const OperatorsForMesons& meson_operator,
//                       const int t1_block, const int t2_block,
//                       const std::vector<QuarklineQ2Indices>& ql_lookup,
//                       const std::vector<RandomIndexCombinationsQ2>& ric_lookup);
//};
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------


/*! typetrait class which allows to use QuarklineQ1Indices for Q1 and 
 *  QuarklineQ2Indices for Q2L and Q2V
 */
template <QuarkLineType qlt>
struct QuarkLineIndices {};

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
class QuarkLine_one_t {

private:

  // containers for the three types of quark lines
  boost::multi_array<std::vector<Eigen::MatrixXcd>, 2> Ql;

  const size_t dilT, dilE, nev;
  std::vector<LapH::gamma_lookup>  gamma;

public:

  QuarkLine_one_t (const size_t dilT, const size_t dilE, 
              const size_t nev, const typename QuarkLineIndices<qlt>::type& quarkline_indices,
              const std::vector<RandomIndexCombinationsQ2>& ric_lookup);
  ~QuarkLine_one_t () {}; // dtor

  inline const cmplx& return_gamma_val(const size_t gamma_id, 
                                       const size_t row) const {
    return gamma[gamma_id].value[row];
  }
  inline const int& return_gamma_row(const size_t gamma_id, 
                                     const size_t row) const{
    return gamma[gamma_id].row[row];
  }
  /*! @todo Rewrite that with bracket overload */
  inline const Eigen::MatrixXcd& return_Ql(const size_t t,
                                   const size_t op_id, const size_t rnd) const {
    return Ql[t][op_id].at(rnd);
  }

  // ----------------- INTERFACE FOR BUILDING QUARKLINES -----------------------
  // ---------------------------------------------------------------------------
  void build_Q1_one_t(const Perambulator& peram,
                      const OperatorsForMesons& meson_operator, size_t pos,
                      const int t_source, const int t_sink,
                      const typename QuarkLineIndices<qlt>::type& ql_lookup,
                      const std::vector<RandomIndexCombinationsQ2>& ric_lookup);
  void build(const Perambulator& peram,
                       const OperatorsForMesons& meson_operator,
                       const int t_source, const int t_sink,
                       const typename QuarkLineIndices<qlt>::type& ql_lookup,
                       const std::vector<RandomIndexCombinationsQ2>& ric_lookup);
};
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

} // end of namespace

#endif

