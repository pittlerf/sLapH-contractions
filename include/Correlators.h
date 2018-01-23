/*! @file Correlators.h
 *  Class declaration of Correlators
 *
 *  @author Bastian Knippschild
 *  @author Martin Ueding
 *  @author Markus Werner
 *
 */

#pragma once

#include <complex>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include <sys/types.h>  
#include <sys/stat.h>   

#include "boost/multi_array.hpp"
#include "Eigen/Dense"

#include "OperatorsForMesons.h"
#include "Perambulator.h"
#include "typedefs.h"
#include "DiagramForward.h"
#include "QuarkLineBlock2.h"


/*! Locally replaces QuarklineLookup extended by lookuptable for rVdaggerVr */
struct DilutedFactorLookup{
  std::vector<QuarklineIndices> const Q0;  
  std::vector<QuarklineIndices> const Q1;
  std::vector<QuarklineIndices> const Q2V;
  std::vector<QuarklineIndices> const Q2L;
};

/******************************************************************************/

/*! Calculates correlation functions according to the stochastic Laplacian 
 *  Heaviside (sLapH) method.
 *
 *  Within the sLapH framework every hadronic correlation function can be 
 *  rewritten as trace of a product of perambulators and operators. In this 
 *  class said traces are calculated.
 *
 *  The central function is performed in contract(). Here, all the single 
 *  diagrams are calculated successively. To this end, lookup tables for 
 *  Correlators, Quarklines and Operators are needed that specify which 
 *  combinations of physical quantum numbers are to be evaluated. These lists 
 *  are built from the infile in GlobalData::init_lookup_tables().
 *  
 *  Additionally the necessary data is passed in the form of 
 *  instances of Quarklines, OperatorsForMesons and 
 *  Perambulators 
 *
 *  The diagrams corr0, corrC (and thus C20, C2+) as well as C3c and C4cB are 
 *  memory optimized calling quarklines within an outer loop over time and thus
 *  only need 1/Lt the memory.
 *
 *  @todo check whether other correlators are still functional
 *  @todo make other correlators call one_t quarklines as well
 *
 */

void contract(const size_t Lt,
              const size_t dilT,
              const size_t dilE,
              const size_t nev,
              OperatorsForMesons const &meson_operator,
              RandomVector const &randomvectors,
              Perambulator const &perambulators,
              OperatorLookup const &operator_lookup,
              CorrelatorLookup const &corr_lookup,
              DilutedFactorIndicesCollection const &quark_lookup,
              std::string const output_path,
              std::string const output_filename);
