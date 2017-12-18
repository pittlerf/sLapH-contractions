/*! @file Operator.h
 *  Declaration of Operator multiplication
 *
 *  @author Markus Werner
 *
 */

#pragma once

#include <iosfwd>
#include <set>
#include <sstream>
#include <vector>

#include "Eigen/Dense"

#include "OperatorsForMesons.h"
#include "QuarkLineBlock.h"
#include "typedefs.h"

namespace LapH {

struct DilutedFactor {
  using Data = Eigen::MatrixXcd;
  using RndId = int8_t;

  Data data;
  std::pair<RndId, RndId> ric;
  std::vector<RndId> used_rnd_ids;
};

struct DilutedTrace {
  using Data = cmplx;
  using RndId = int8_t;

  Data data;
  std::vector<RndId> used_rnd_ids;
};

/*! Product yielding the off-diagonal elements.

  From the sets of DilutedFactor elements, the product set of DilutedFactor is build such
  that it only contains elements with _unequal_ left and right random vector index. This
  set is intended to be used as an intermediate result.
  */
std::vector<DilutedFactor> operator*(std::vector<DilutedFactor> const &left_vec,
                                     std::vector<DilutedFactor> const &right_vec);

template <int n>
using OperatorToFactorMap = std::map<std::array<size_t, n>, std::vector<DilutedFactor>>;

#if 0
template <int n1, int n2>
OperatorToFactorMap<n1 + n2> operator*(OperatorToFactorMap<n1> const &left_map,
                                       OperatorToFactorMap<n2> const &right_map) {
  OperatorToFactorMap<n1 + n2> result;

  for (auto const &left : left_map) {
    for (auto const &right : right_map) {
      // Concatenate the two keys from the left and right element into the new key.
      typename OperatorToFactorMap<n1 + n2>::key_type key;
      auto out_it =
          std::copy(std::begin(left.first), std::end(left.first), std::begin(key));
      std::copy(std::begin(right.first), std::end(right.first), out_it);

      // Do the actual multiplication.
      result[key] = left.second * right.second;
    }
  }
}
#endif

// Proposed:
// template <int n>
// using OperatorToFactorMap = std::map<std::array<QuantumNumbers, n>,
// std::vector<DilutedFactor>>;

template <int n>
std::string to_string(typename OperatorToFactorMap<n>::key_type const &array) {
  std::ostringstream oss;
  oss << "{";
  for (int i = 0; i < n; ++i) {
    if (i != 0) {
      oss << ", ";
    }
    oss << array[i];
  }
  oss << "}";

  return oss.str();
}

template <int n>
void print(OperatorToFactorMap<n> const &otfm) {
  std::cout << "OperatorToFactorMap, size = " << otfm.size() << "\n";
  for (auto const &elem : otfm) {
    std::cout << "  " << to_string<n>(elem.first) << " -> "
              << "std::vector(size = " << elem.second.size() << ")\n";
  }
}

/*! @todo Be more restrictive with lookup tables. .Q2V etc. is enough */
template <QuarkLineType qlt>
void check_random_combinations(std::string const &diagram,
                               std::vector<size_t> const &lookup,
                               std::vector<RandomIndexCombinationsQ2> const &ric_lookup,
                               std::vector<VdaggerVRandomLookup> const &rvdaggervr_lookup,
                               std::vector<QuarklineQ2Indices> const &Q2_lookup);

/*! Multiply (Q2V*rVdaggerV) and take trace
 *  - corrC
 */
template <QuarkLineType qlt1, QuarkLineType qlt2>
std::vector<cmplx> trace(QuarkLineBlock<qlt1> const &quarkline1,
                         QuarkLineBlock<qlt2> const &quarkline2,
                         int const t1,
                         int const b2,
                         int const t2,
                         std::vector<size_t> const &lookup,
                         std::vector<RandomIndexCombinationsQ2> const &ric_lookup,
                         std::vector<size_t> const &ric_ids,
                         int const gamma,
                         size_t const dilE,
                         size_t const dilD);

/*! Multiply (Q1*Q1) and take trace
 *  - corr0
 */
std::vector<cmplx> trace(std::vector<Eigen::MatrixXcd> const &quarkline1,
                         std::vector<Eigen::MatrixXcd> const &quarkline2,
                         std::vector<RandomIndexCombinationsQ2> const &ric_lookup,
                         std::vector<size_t> const &ric_ids);

/*! Multiply two traces of two Quarklines each: tr(QQ) * tr(QQ)
 */
compcomp_t trtr(std::vector<cmplx> const &factor1,
                std::vector<cmplx> const &factor2,
                std::vector<RandomIndexCombinationsQ2> const &ric_lookup,
                std::vector<size_t> const &ric_ids);

/******************************************************************************/

/*! Create vector<MatrixXcd> with Q1 for all rnd_vecs
 *  - C1
 *  - C3c
 *  - C30
 */
void Q1(std::vector<Eigen::MatrixXcd> &result,
        std::vector<Eigen::MatrixXcd> const &quarklines,
        std::vector<RandomIndexCombinationsQ2> const &ric_lookup,
        std::vector<size_t> const &ric_ids,
        size_t const dilE,
        size_t const dilD);

void Q1xQ1(std::vector<DilutedFactor> &result,
           std::vector<DilutedFactor> const &quarkline1,
           std::vector<DilutedFactor> const &quarkline2,
           std::vector<RandomIndexCombinationsQ2> const &ric_lookup,
           std::vector<size_t> const ric_ids,
           size_t const dilE,
           size_t const dilD);

/*! Create vector<MatrixXcd> with Q1*Q1 for all rnd vecs not equal
 *  - C30
 *  - C40B
 *  - C40C
 */
void Q1xQ1(std::vector<Eigen::MatrixXcd> &result,
           std::vector<Eigen::MatrixXcd> const &quarkline1,
           std::vector<Eigen::MatrixXcd> const &quarkline2,
           std::vector<RandomIndexCombinationsQ2> const &ric_lookup,
           std::vector<size_t> const ric_ids,
           size_t const dilE,
           size_t const dilD);

/*! Create vector<MatrixXcd> with Q0*Q2 for all rnd vecs not equal
 *  - (corrC)
 *  - C4cB
 *  - C4cC
 *  - C3c
 */
void rVdaggerVrxQ2(std::vector<Eigen::MatrixXcd> &result,
                   std::vector<Eigen::MatrixXcd> const &quarkline1,
                   std::vector<Eigen::MatrixXcd> const &quarkline2,
                   std::vector<RandomIndexCombinationsQ2> const &ric_lookup,
                   std::vector<size_t> const &ric_ids,
                   size_t const dilE,
                   size_t const dilD);

/*! Multiply (QQ)*(Q) and take trace
 *  - C3c
 *  - C30
 *  @calls M1xM2 for Optimization
 */
cmplx trace_3pt(std::vector<Eigen::MatrixXcd> const &M1,
                std::vector<Eigen::MatrixXcd> const &M2,
                std::vector<RandomIndexCombinationsQ2> const &ric_lookup,
                std::vector<size_t> const &ric_ids,
                size_t const dilE,
                size_t const dilD);

/*! Multiply (QQ)*(QQ). or
 *  and take trace
 *  - C4cC
 *  - C4cB
 *  - C40C
 *  - C40B
 *  @calls M1xM2 for Optimization
 */
cmplx trace(std::vector<Eigen::MatrixXcd> const &M1,
            std::vector<Eigen::MatrixXcd> const &M2,
            std::vector<RandomIndexCombinationsQ2> const &ric_lookup,
            std::vector<size_t> const &ric_ids,
            size_t const dilE,
            size_t const dilD);

cmplx trace(std::vector<DilutedFactor> const &left_vec,
            std::vector<DilutedFactor> const &right_vec);

compcomp_t inner_product(std::vector<DilutedTrace> const &left_vec,
                         std::vector<DilutedTrace> const &right_vec);

}  // end of namespace
