/*! @file Correlators.h
 *  Class declaration of LapH::Correlators
 *
 *  @author Bastian Knippschild
 *  @author Markus Werner
 *
 */

#ifndef CORRELATORS_H_
#define CORRELATORS_H_

#include <complex>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include <sys/types.h>  
#include <sys/stat.h>   

#include "boost/multi_array.hpp"
#include "boost/filesystem.hpp"
#include "Eigen/Dense"

#include "OperatorsForMesons.h"
#include "Quarklines.h"
#include "typedefs.h"

#include "H5Cpp.h"

namespace LapH {

/******************************************************************************/ 
// This is just a workaround for complex numbers to get it running for hdf5
/*! @TODO: Why is this in a named namespace in the include file and not in an
 *         unnamed namespace in the source file?
 */
typedef struct { 
  double re; 
  double im; 
} complex_t; 

// This is the datatype to write 4pt functions and that alike directly
struct compcomp_t { 
  double rere;   
  double reim;
  double imre;
  double imim;   
  compcomp_t(const double rere, const double reim, 
             const double imre, const double imim) : 
                              rere(rere), reim(reim), imre(imre), imim(imim) {};
}; 
/******************************************************************************/

/*! Calculates correlation functions according to the stochastic Laplacian 
 *  Heaviside method
 *
 *  Given lookup tables for Correlators, Quarklines and Operators and the data 
 *  specified within this is the central computation.
 *
 *  The diagrams corr0, corrC (and thus C20, C2+) as well as C3c and C4cB are 
 *  memory optimized calling quarklines within an outer loop over time and thus
 *  only need 1/Lt the memory.
 *
 *  @TODO check whether other correlators are still functional
 *  @TODO make other correlators call one_t quarklines as well
 *
 */
class Correlators {

private:
  array_corr corrC, corr0; // this is only needed intermideately
  const size_t Lt, dilT, dilE, nev;

  // this is just for intermediate steps
  void build_corr0(const OperatorsForMesons& meson_operator, 
                   const Perambulator& perambulators,
                   const std::vector<CorrInfo>& corr_lookup,
                   const QuarklineLookup& quark_lookup,
                   const OperatorLookup& operator_lookup);
  void build_corrC(const Perambulator& perambulators,
                   const OperatorsForMesons& meson_operator,
                   const OperatorLookup& OperatorLookup,
                   const std::vector<CorrInfo>& corr_lookup,
                   const QuarklineLookup& quark_lookup);

  // functions to build correlation functions
  void build_C1(const Quarklines& quarklines, 
                const std::vector<CorrInfo>& corr_lookup,
                const QuarklineLookup& quark_lookup,
                const std::vector<RandomIndexCombinationsQ2>& ric_lookup);
  void build_C20(const std::vector<CorrInfo>& corr_lookup);
  void build_C30(const Quarklines& quarklines, 
                 const std::vector<CorrInfo>& corr_lookup,
                 const QuarklineLookup& quark_lookup,
                 const std::vector<RandomIndexCombinationsQ2>& ric_lookup);
  void build_C40D(const OperatorLookup& operator_lookup, 
                  const CorrelatorLookup& corr_lookup,
                  const QuarklineLookup& quark_lookup);
  void build_C40V(const OperatorLookup& operator_lookup, 
                  const CorrelatorLookup& corr_lookup,
                  const QuarklineLookup& quark_lookup);
  void build_C40C(const Quarklines& quarklines, 
                  const std::vector<CorrInfo>& corr_lookup,
                  const QuarklineLookup& quark_lookup,
                  const std::vector<RandomIndexCombinationsQ2>& ric_lookup);
  void build_C40B(const Quarklines& quarklines, 
                  const std::vector<CorrInfo>& corr_lookup,
                  const QuarklineLookup& quark_lookup,
                  const std::vector<RandomIndexCombinationsQ2>& ric_lookup);
  void build_C2c(const std::vector<CorrInfo>& corr_lookup);
  void build_C3c(const OperatorsForMesons& meson_operator,
                 const Perambulator& perambulators,
                 const OperatorLookup& operator_lookup,
                 const std::vector<CorrInfo>& corr_lookup,
                 const QuarklineLookup& quark_lookup);
  void build_C4cD(const OperatorLookup& operator_lookup, 
                  const CorrelatorLookup& corr_lookup,
                  const QuarklineLookup& quark_lookup);
  void build_C4cV(const OperatorLookup& operator_lookup, 
                  const CorrelatorLookup& corr_lookup,
                  const QuarklineLookup& quark_lookup);
  void build_C4cC(const Quarklines& quarklines, 
                  const OperatorsForMesons& meson_operator,
                  const OperatorLookup& operator_lookup,
                  const std::vector<CorrInfo>& corr_lookup,
                  const QuarklineLookup& quark_lookup);
  void build_C4cB(const OperatorsForMesons& meson_operator,
                  const Perambulator& perambulators,
                  const OperatorLookup& operator_lookup,
                  const std::vector<CorrInfo>& corr_lookup,
                  const QuarklineLookup& quark_lookup);

public:
  // Constructor
  Correlators (const size_t Lt, const size_t dilT, const size_t dilE, 
               const size_t nev, const CorrelatorLookup& corr_lookup) :
               Lt(Lt), dilT(dilT), dilE(dilE), nev(nev) {};
  // Standard Destructor
  ~Correlators () {};

  /*! Call all functions building a correlator */
  void contract(Quarklines& quarklines, 
                const OperatorsForMesons& meson_operator,
                const Perambulator& perambulators,
                const OperatorLookup& operator_lookup,
                const CorrelatorLookup& corr_lookup, 
                const QuarklineLookup& quark_lookup);
};
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

} // end of namespace

#endif // CORRELATORS_H_ 
