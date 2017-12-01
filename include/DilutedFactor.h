/*! @file Operator.h
 *  Declaration of Operator multiplication
 *
 *  @author Markus Werner
 *
 */

#pragma once

#include <vector>
#include "Eigen/Dense"

#include "OperatorsForMesons.h"
#include "QuarkLineBlock.h"
#include "typedefs.h"

namespace LapH{

/*! @todo Be more restrictive with lookup tables. .Q2V etc. is enough */
template <QuarkLineType qlt>
void check_random_combinations(std::string const &diagram,
                               std::vector<size_t> const &lookup,
                               std::vector<RandomIndexCombinationsQ2> const &ricQ2_lookup,
                               std::vector<VdaggerVRandomLookup> const &rvdaggervr_lookup,
                               std::vector<QuarklineQ2Indices> const &Q2_lookup);

/*! Create vector<MatrixXcd> with Q1 for all rnd_vecs
 *  - C1
 *  - C3c
 *  - C30
 */
template <QuarkLineType qlt>
void Q1(std::vector<Eigen::MatrixXcd> &result, 
                    QuarkLineBlock<qlt> const &quarklines,
                    int const t1,
                    int const b2,
                    std::array<size_t, 2> const look,
                    std::vector<RandomIndexCombinationsQ2> const &ricQ2_lookup,
                    std::vector<QuarklineQ1Indices> const &Q1_lookup,
                    size_t const dilE,
                    size_t const dilD);

/*! Create vector<MatrixXcd> with Q1*Q1 for all rnd vecs not equal
 *  - C30
 *  - C40B
 *  - C40C
 */
template <QuarkLineType qlt>
void Q1xQ1(std::vector<Eigen::MatrixXcd> &result, 
           QuarkLineBlock<qlt> const &quarklines,
           int const t1,
           int const b1,
           int const t2,
           int const b2,
           std::array<size_t, 3> const look,
           std::vector<RandomIndexCombinationsQ2> const &ricQ2_lookup,
           std::vector<QuarklineQ1Indices> const &Q1_lookup,
           size_t const dilE,
           size_t const dilD);

/*! Create vector<MatrixXcd> with Q1*Q1 for all rnd vecs not equal
 *  - (corr0)
 *  - C3c
 *  - C40B
 *  - C40C
 */
template <QuarkLineType qlt>
void Q2xrVdaggerVr(std::vector<Eigen::MatrixXcd> &result, 
                    QuarkLineBlock<qlt> const &quarklines,
                    OperatorsForMesons const &meson_operator,
                    int const b2,
                    int const t1,
                    int const t2,
                    std::array<size_t, 3> const look,
                    std::vector<RandomIndexCombinationsQ2> const &ricQ2_lookup,
                    std::vector<VdaggerVRandomLookup> const &rvdaggervr_lookup,
                    std::vector<QuarklineQ2Indices> const &Q2V_lookup,
                    size_t const dilE,
                    size_t const dilD);

/*! Create vector<MatrixXcd> with Q1*Q1 for all rnd vecs not equal
 *  - (corrC)
 *  - C4cB
 *  - C4cC
 */
template <QuarkLineType qlt>
void rVdaggerVrxQ2(std::vector<Eigen::MatrixXcd> &result, 
                    QuarkLineBlock<qlt> const &quarklines,
                    OperatorsForMesons const &meson_operator,
                    int const t1,
                    int const b2,
                    std::array<size_t, 3> const look,
                    std::vector<RandomIndexCombinationsQ2> const &ricQ2_lookup,
                    std::vector<VdaggerVRandomLookup> const &rvdaggervr_lookup,
                    std::vector<QuarklineQ2Indices> const &Q2V_lookup,
                    size_t const dilE,
                    size_t const dilD);

/*! Multiply (Q2V*rVdaggerV) and take trace
 *  - corrC
 */
template <QuarkLineType qlt1, QuarkLineType qlt2>
std::vector<cmplx> trace(
    QuarkLineBlock<qlt1> const &quarklines,
    OperatorsForMesons const &meson_operator,
    int const t1,
    int const b2,
    int const t2,
    std::vector<size_t> const &lookup,
    std::vector<RandomIndexCombinationsQ2> const &ricQ2_lookup,
    std::vector<VdaggerVRandomLookup> const &rvdaggervr_lookup,
    std::vector<QuarklineQ2Indices> const &Q2_lookup,
    int const gamma,
    size_t const dilE,
    size_t const dilD);

/*! Multiply (Q1*Q1) and take trace
 *  - corr0
 */
template <QuarkLineType qlt1, QuarkLineType qlt2>
std::vector<cmplx> trace(
    QuarkLineBlock<qlt1> const &quarklines,
    int const t1,
    int const b2,
    int const t2,
    int const b1,
    std::vector<size_t> const &lookup,
    std::vector<RandomIndexCombinationsQ2> const &ricQ2_lookup,
    std::vector<QuarklineQ1Indices> const &Q1_lookup);

/*! Multiply (Q2V*rVdaggerV)*(Q2V*rVdaggerVr). or 
 *           (rVdaggerVr*Q2L)*(rVdaggerVr*Q2L) (two implementations) 
 *  and take trace
 *  - C4cC
 *  - C4cB
 *  @calls M1xM2 for Optimization
 */
template <QuarkLineType qlt1, QuarkLineType qlt2>
cmplx trace(std::vector<Eigen::MatrixXcd> const &M1, 
           std::vector<Eigen::MatrixXcd> const &M2, 
           std::vector<size_t> const &lookup,
           std::vector<RandomIndexCombinationsQ2> const &ricQ2_lookup,
           std::vector<VdaggerVRandomLookup> const &rvdaggervr_lookup,
           std::vector<QuarklineQ2Indices> const &Q2_lookup,
           size_t const dilE,
           size_t const dilD);

/*! Multiply (Q2V*rVdaggerV)*(Q1) and take trace
 *  - C3c
 */
template <QuarkLineType qlt1, QuarkLineType qlt2>
cmplx trace(std::vector<Eigen::MatrixXcd> const &M1, 
           std::vector<Eigen::MatrixXcd> const &M2, 
           std::vector<size_t> const &lookup,
           std::vector<RandomIndexCombinationsQ2> const &ricQ2_lookup,
           std::vector<RandomIndexCombinationsQ1> const &ricQ1_lookup,
           std::vector<VdaggerVRandomLookup> const &rvdaggervr_lookup,
           std::vector<QuarklineQ1Indices> const &Q1_lookup,
           std::vector<QuarklineQ2Indices> const &Q2_lookup,
           size_t const dilE,
           size_t const dilD);

/*! Multiply (Q1*Q1)*(Q1) and take trace
 *  - C40B
 */
template <QuarkLineType qlt1, QuarkLineType qlt2>
cmplx trace(std::vector<Eigen::MatrixXcd> const &L1, 
            std::vector<Eigen::MatrixXcd> const &L2, 
            std::vector<size_t> const &lookup,
            std::vector<RandomIndexCombinationsQ2> const &ricQ2_lookup,
            std::vector<QuarklineQ1Indices> const &Q1_lookup,
            size_t const dilE,
            size_t const dilD);

/*! Multiply (Q1*Q1)*(Q1*Q1) and take trace
 *  - C30
 */
template <QuarkLineType qlt1, QuarkLineType qlt2>
cmplx trace(std::vector<Eigen::MatrixXcd> const &L1, 
            std::vector<Eigen::MatrixXcd> const &L2, 
            std::vector<size_t> const &lookup,
            std::vector<RandomIndexCombinationsQ2> const &ricQ2_lookup,
            std::vector<QuarklineQ1Indices> const &Q1_lookup);

/*! Multiply tr(Q1*Q1) * tr(Q1*Q1) 
 *  - C40D
 *  - C40V
 */
template <QuarkLineType qlt1, QuarkLineType qlt2>
void trtr(compcomp_t &result,
          std::vector<cmplx> const &factor1,
          std::vector<cmplx> const &factor2,
          std::vector<size_t> const &lookup,
          std::vector<RandomIndexCombinationsQ2> const &ricQ2_lookup,
          std::vector<QuarklineQ1Indices> const &Q1_lookup);

t/*! Multiply tr(Q2V*rVdaggerVr) * tr(Q2V*rVdaggerVr) 
 *  - C40D
 *  - C40V
 *  rVdaggerVr_lookup is not necessary as it contains information redundant
 *  to Q2V_lookup
 */
emplate<QuarkLineType qlt1, QuarkLineType qlt2>
void trtr(compcomp_t &result,
          std::vector<cmplx> const &factor1,
          std::vector<cmplx> const &factor2,
          std::vector<size_t> const &lookup,
          std::vector<RandomIndexCombinationsQ2> const &ricQ2_lookup,
          std::vector<QuarklineQ2Indices> const &Q2_lookup);

}  // end of namespace
