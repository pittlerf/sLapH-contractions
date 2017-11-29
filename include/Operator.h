/*! @file Operator.h
 *  Declaration of Operator multiplication
 *
 *  @author Martin Ueding
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
void check_random_combinations(std::string const &diagram,
                               std::vector<size_t> const &lookup,
                               std::vector<RandomIndexCombinationsQ2> const &ricQ2_lookup,
                               std::vector<VdaggerVRandomLookup> const &rvdaggervr_lookup,
                               std::vector<QuarklineQ2Indices> const &Q2V_lookup);

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

cmplx trace(std::vector<Eigen::MatrixXcd> const &M1, 
           std::vector<Eigen::MatrixXcd> const &M2, 
           std::vector<size_t> const &lookup,
           std::vector<RandomIndexCombinationsQ2> const &ricQ2_lookup,
           std::vector<VdaggerVRandomLookup> const &rvdaggervr_lookup,
           std::vector<QuarklineQ2Indices> const &Q2V_lookup,
           size_t const dilE,
           size_t const dilD);

}  // end of namespace
