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


/*! Locally replaces QuarklineLookup extended by lookuptable for rVdaggerVr */
struct DilutedFactorLookup{
  std::vector<QuarklineIndices> const Q0;  
  std::vector<QuarklineIndices> const Q1;
  std::vector<QuarklineQ2Indices> const Q2V;
  std::vector<QuarklineQ2Indices> const Q2L;
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
class Correlators {

public:
  // Constructor
  Correlators (const size_t Lt, const size_t dilT, const size_t dilE, 
               const size_t nev, const CorrelatorLookup& corr_lookup,
               OperatorLookup const &operator_lookup,
               QuarklineLookup const &quark_lookup);
  // Standard Destructor
  ~Correlators () {};

  /*! Call all functions building a correlator */
  void contract(OperatorsForMesons const &meson_operator,
                RandomVector const &randomvectors,
                Perambulator const &perambulators,
                OperatorLookup const &operator_lookup,
                CorrelatorLookup const &corr_lookup, 
                QuarklineLookup const &quark_lookup,
                std::string const output_path,
                std::string const output_filename);

private:
  /*! @todo that should not be here but taken from Globaldata */
  const size_t Lt, dilT, dilE, nev;

  DilutedFactorLookup const dil_fac_lookup;

  /*! Temporal memory for tr(rVdaggerV*Q1*rVdaggerV*Q1) */
  DilutedTraceCollection<2> corr0; 
  /*! Calculate tr(rVdaggerV*Q1*rVdaggerV*Q1) */
  void build_corr0(RandomVector const &randomvectors,
                   OperatorsForMesons const &meson_operator, 
                   Perambulator const &perambulators,
                   std::vector<CorrInfo> const &corr_lookup);
  /*! Temporal memory for tr(Q2V*rVdaggerVr) */
  DilutedTraceCollection<2> corrC; 
  /*! Calculate tr(Q2V*rVdaggerVr) */
  void build_corrC(RandomVector const &randomvectors,
                   Perambulator const &perambulators,
                   OperatorsForMesons const &meson_operator,
                   std::vector<CorrInfo> const &corr_lookup);

  
  /*! Temporal memory for tr(Q1) */
  DilutedTraceCollection2<1> corr_part_trQ1;
  /*! Build 1pt loops */
  void build_part_trQ1(RandomVector const &randomvectors,
                   OperatorsForMesons const &meson_operator, 
                Perambulator const &perambulators,
                std::vector<CorrInfo> const &corr_lookup,
                std::string const output_path,
                std::string const output_filename);

  // Functions to build correlation functions

  /*! Build neutral 2pt correlation function 
   *  @f{align}{
   *    C = \langle D_\mathtt{Q0}^{-1}(t'|t) \Gamma_\mathtt{Op0} 
   *                D_\mathtt{Q1}^{-1}(t|t') \Gamma_\mathtt{Op1} \rangle
   *  @f}
   */
  void build_C20(std::vector<CorrInfo> const &corr_lookup, 
                 std::string const output_path, 
                 std::string const filename);

  /*! Build neutral 2pt correlation function 
   *  @f{align}{
   *    C = \langle D_\mathtt{Q0}^{-1}(t|t) \Gamma_\mathtt{Op0} \rangle \cdot
   *        \langle D_\mathtt{Q1}^{-1}(t'|t') \Gamma_\mathtt{Op1} \rangle
   *  @f}
   */
  void build_C20V(std::vector<CorrInfo> const &corr_lookup,
                  std::string const output_path,
                  std::string const output_filename);


  /*! Build neutral 4pt correlation function: Direct diagram
   *  @f{align}{
   *    C = \langle D_\mathtt{Q0}^{-1}(t'|t) \Gamma_\mathtt{Op0} 
   *                D_\mathtt{Q1}^{-1}(t|t') \Gamma_\mathtt{Op1} \rangle \cdot
   *        \langle D_\mathtt{Q2}^{-1}(t'|t) \Gamma_\mathtt{Op2} 
   *                D_\mathtt{Q3}^{-1}(t|t') \Gamma_\mathtt{Op3} \rangle
   *  @f}
   */
  void build_C40D(std::vector<CorrInfo> const &corr_lookup,
                  std::string const output_path,
                  std::string output_filename);
  /*! Build neutral 4pt correlation function: Vacuum diagram
   *  @f{align}{
   *    C = \langle D_\mathtt{Q0}^{-1}(t|t) \Gamma_\mathtt{Op0} 
   *                D_\mathtt{Q1}^{-1}(t|t) \Gamma_\mathtt{Op1} \rangle \cdot
   *        \langle D_\mathtt{Q2}^{-1}(t'|t') \Gamma_\mathtt{Op2} 
   *                D_\mathtt{Q3}^{-1}(t'|t') \Gamma_\mathtt{Op3} \rangle
   *  @f}
   */
  void build_C40V(std::vector<CorrInfo> const &corr_lookup,
                  std::string const output_path,
                  std::string const output_filename);

  /*! Build charged 2pt correlation function 
   *  @f{align}{
   *    C = \langle \gamma_5 D_\mathtt{Q0}^{-1}(t|t')^\dagger \gamma_5  \Gamma_\mathtt{Op0} 
   *                D_\mathtt{Q1}^{-1}(t|t') \Gamma_\mathtt{Op1} \rangle
   *  @f}
   */
  void build_C2c(std::vector<CorrInfo> const &corr_lookup,
                 std::string const output_path,
                 std::string const output_filename);

  /*! Build charged 4pt correlation function: Direct diagram
   *  @f{align}{
   *    C = \langle \gamma_5 D_\mathtt{Q0}^{-1}(t|t')^\dagger \gamma_5 \Gamma_\mathtt{Op0} 
   *                D_\mathtt{Q1}^{-1}(t|t') \Gamma_\mathtt{Op1} \rangle \cdot
   *        \langle \gamma_5 D_\mathtt{Q2}^{-1}(t|t')^\dagger \gamma_5 \Gamma_\mathtt{Op2} 
   *                D_\mathtt{Q3}^{-1}(t|t') \Gamma_\mathtt{Op3} \rangle
   *  @f}
   */
  void build_C4cD(CorrelatorLookup const &corr_lookup,
                  std::string const output_path,
                  std::string const output_filename);
  /*! Build charged 4pt correlation function: Vacuum diagram
   *  @f{align}{
   *    C = \langle \gamma_5 D_\mathtt{Q0}^{-1}(t|t)^\dagger \gamma_5 \Gamma_\mathtt{Op0} 
   *                D_\mathtt{Q1}^{-1}(t|t) \Gamma_\mathtt{Op1} \rangle \cdot
   *        \langle \gamma_5 D_\mathtt{Q2}^{-1}(t'|t')^\dagger \gamma_5 \Gamma_\mathtt{Op2} 
   *                D_\mathtt{Q3}^{-1}(t'|t') \Gamma_\mathtt{Op3} \rangle
   *  @f}
   */
  void build_C4cV(CorrelatorLookup const &corr_lookup,
                  std::string const output_path,
                  std::string const output_filename);
};

template <typename Numeric>
void build_diagram(typename DiagramTraits<Numeric>::Diagram &diagram,
                   RandomVector const &randomvectors,
                   OperatorsForMesons const &meson_operator,
                   Perambulator const &perambulators,
                   std::string const output_path,
                   std::string const output_filename,
                   const size_t Lt,
                   const size_t dilT,
                   const size_t dilE,
                   const size_t nev,
                   DilutedFactorLookup const &dil_fac_lookup,
                   DilutedTraceCollection<2> &corr0,
                   DilutedTraceCollection<2> &corrC,
                   DilutedTraceCollection2<1> &corr_part_trQ1);
