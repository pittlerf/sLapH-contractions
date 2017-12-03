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
                               std::vector<RandomIndexCombinationsQ2> const &ric_lookup,
                               std::vector<VdaggerVRandomLookup> const &rvdaggervr_lookup,
                               std::vector<QuarklineQ2Indices> const &Q2_lookup);

/*! Multiply (Q2V*rVdaggerV) and take trace
 *  - corrC
 */
template <QuarkLineType qlt1, QuarkLineType qlt2>
std::vector<cmplx> trace(
    QuarkLineBlock<qlt1> const &quarkline1,
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
std::vector<cmplx> trace(
    std::vector<Eigen::MatrixXcd> const &quarkline1,
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



}  // end of namespace
