/*! @file Correlators.h
 *  Class declaration of LapH::Correlators
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
 *  instances of LapH::Quarklines, LapH::OperatorsForMesons and 
 *  LapH::Perambulators 
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

private:
  /*! @todo that should not be here but taken from Globaldata */
  const size_t Lt, dilT, dilE, nev;

  /*! Temporal memory for Q2V*rVdaggerVr (without trace!) */
  array_corr corrC;
  /*! Calculate Q2V*rVdaggerVr (without trace!) */
  void build_corr0(const OperatorsForMesons& meson_operator, 
                   const Perambulator& perambulators,
                   const std::vector<CorrInfo>& corr_lookup,
                   const QuarklineLookup& quark_lookup,
                   const OperatorLookup& operator_lookup);
  /*! Temporal memory for Q1*VdaggerVr*Q1*VdaggerVr (without trace!) */
  array_corr corr0; 
  /*! Calculate Q2V*rVdaggerVr (without trace!) */
  void build_corrC(const Perambulator& perambulators,
                   const OperatorsForMesons& meson_operator,
                   const OperatorLookup& OperatorLookup,
                   const std::vector<CorrInfo>& corr_lookup,
                   const QuarklineLookup& quark_lookup);

  // Functions to build correlation functions
  
  /*! Build 1pt correlation function 
   *  @f{align}{
   *    C = \langle D_\mathtt{Q0}^{-1}(t|t') \Gamma_\mathtt{Op0} \rangle
   *  @f}
   */
  void build_C1(const Quarklines& quarklines, 
                const std::vector<CorrInfo>& corr_lookup,
                const QuarklineLookup& quark_lookup,
                const std::vector<RandomIndexCombinationsQ2>& ric_lookup);
  /*! Build neutral 2pt correlation function 
   *  @f{align}{
   *    C = \langle D_\mathtt{Q0}^{-1}(t'|t) \Gamma_\mathtt{Op0} 
   *                D_\mathtt{Q1}^{-1}(t|t') \Gamma_\mathtt{Op1} \rangle
   *  @f}
   */
  void build_C20(const std::vector<CorrInfo>& corr_lookup, 
                     const std::string output_path, const std::string filename);

  /*! Build neutral 3pt correlation function 
   *  @f{align}{
   *    C = \langle D_\mathtt{Q0}^{-1}(t|t) \Gamma_\mathtt{Op0} 
   *                D_\mathtt{Q1}^{-1}(t|t') \Gamma_\mathtt{Op1} 
   *                D_\mathtt{Q2}^{-1}(t'|t) \Gamma_\mathtt{Op2} \rangle
   *  @f}
   */
  void build_C30(OperatorsForMesons const &meson_operator,
                  Perambulator const &perambulators,
                  std::vector<CorrInfo> const &corr_lookup,
                  QuarklineLookup const &quark_lookup,
                  OperatorLookup const &operator_lookup,
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
  void build_C40D(const OperatorLookup& operator_lookup, 
                  const CorrelatorLookup& corr_lookup,
                  const QuarklineLookup& quark_lookup,
                  const std::string output_path,
                  const std::string output_filename);
  /*! Build neutral 4pt correlation function: Vacuum diagram
   *  @f{align}{
   *    C = \langle D_\mathtt{Q0}^{-1}(t|t) \Gamma_\mathtt{Op0} 
   *                D_\mathtt{Q1}^{-1}(t|t) \Gamma_\mathtt{Op1} \rangle \cdot
   *        \langle D_\mathtt{Q2}^{-1}(t'|t') \Gamma_\mathtt{Op2} 
   *                D_\mathtt{Q3}^{-1}(t'|t') \Gamma_\mathtt{Op3} \rangle
   *  @f}
   */
  void build_C40V(OperatorLookup const &operator_lookup, 
                  CorrelatorLookup const &corr_lookup,
                  QuarklineLookup const &quark_lookup,
                  std::string const output_path,
                  std::string const output_filename);
  /*! Build neutral 4pt correlation function: Cross diagram
   *  @f{align}{
   *    C = \langle D_\mathtt{Q0}^{-1}(t'|t) \Gamma_\mathtt{Op0} 
   *                D_\mathtt{Q1}^{-1}(t|t') \Gamma_\mathtt{Op1}
   *                D_\mathtt{Q2}^{-1}(t'|t) \Gamma_\mathtt{Op2} 
   *                D_\mathtt{Q3}^{-1}(t|t') \Gamma_\mathtt{Op3} \rangle
   *  @f}
   */
  void build_C40C(OperatorsForMesons const &meson_operator,
                  Perambulator const &perambulators,
                  std::vector<CorrInfo> const &corr_lookup,
                  QuarklineLookup const &quark_lookup,
                  OperatorLookup const &operator_lookup,
                  std::string const output_path,
                  std::string const output_filename);
  /*! Build neutral 4pt correlation function: Box diagram
   *  @f{align}{
   *    C = \langle D_\mathtt{Q0}^{-1}(t|t) \Gamma_\mathtt{Op0} 
   *                D_\mathtt{Q1}^{-1}(t|t') \Gamma_\mathtt{Op1}
   *                D_\mathtt{Q2}^{-1}(t'|t') \Gamma_\mathtt{Op2} 
   *                D_\mathtt{Q3}^{-1}(t'|t) \Gamma_\mathtt{Op3} \rangle
   *  @f}
   */
  void build_C40B(OperatorsForMesons const &meson_operator,
                  Perambulator const &perambulators,
                  std::vector<CorrInfo> const &corr_lookup,
                  QuarklineLookup const &quark_lookup,
                  OperatorLookup const &operator_lookup,
                  std::string const output_path,
                  std::string const output_filename);
  /*! Build charged 2pt correlation function 
   *  @f{align}{
   *    C = \langle \gamma_5 D_\mathtt{Q0}^{-1}(t|t')^\dagger \gamma_5  \Gamma_\mathtt{Op0} 
   *                D_\mathtt{Q1}^{-1}(t|t') \Gamma_\mathtt{Op1} \rangle
   *  @f}
   */
  void build_C2c(const std::vector<CorrInfo>& corr_lookup,
                 const std::string output_path,
                 const std::string output_filename);
  /*! Build neutral 3pt correlation function 
   *  @f{align}{
   *    C = \langle \gamma_5 D_\mathtt{Q0}^{-1}(t|t)^\dagger \gamma_5 \Gamma_\mathtt{Op0} 
   *                D_\mathtt{Q1}^{-1}(t|t') \Gamma_\mathtt{Op1} 
   *                D_\mathtt{Q2}^{-1}(t'|t) \Gamma_\mathtt{Op2} \rangle
   *  @f}
   */
  void build_C3c(const OperatorsForMesons& meson_operator,
                 const Perambulator& perambulators,
                 const OperatorLookup& operator_lookup,
                 const std::vector<CorrInfo>& corr_lookup,
                 const QuarklineLookup& quark_lookup,
                 const std::string output_path,
                 const std::string output_filename);
  /*! Build charged 4pt correlation function: Direct diagram
   *  @f{align}{
   *    C = \langle \gamma_5 D_\mathtt{Q0}^{-1}(t|t')^\dagger \gamma_5 \Gamma_\mathtt{Op0} 
   *                D_\mathtt{Q1}^{-1}(t|t') \Gamma_\mathtt{Op1} \rangle \cdot
   *        \langle \gamma_5 D_\mathtt{Q2}^{-1}(t|t')^\dagger \gamma_5 \Gamma_\mathtt{Op2} 
   *                D_\mathtt{Q3}^{-1}(t|t') \Gamma_\mathtt{Op3} \rangle
   *  @f}
   */
  void build_C4cD(const OperatorLookup& operator_lookup, 
                  const CorrelatorLookup& corr_lookup,
                  const QuarklineLookup& quark_lookup,
                  const std::string output_path,
                  const std::string output_filename);
  /*! Build charged 4pt correlation function: Vacuum diagram
   *  @f{align}{
   *    C = \langle \gamma_5 D_\mathtt{Q0}^{-1}(t|t)^\dagger \gamma_5 \Gamma_\mathtt{Op0} 
   *                D_\mathtt{Q1}^{-1}(t|t) \Gamma_\mathtt{Op1} \rangle \cdot
   *        \langle \gamma_5 D_\mathtt{Q2}^{-1}(t'|t')^\dagger \gamma_5 \Gamma_\mathtt{Op2} 
   *                D_\mathtt{Q3}^{-1}(t'|t') \Gamma_\mathtt{Op3} \rangle
   *  @f}
   */
  void build_C4cV(const OperatorLookup& operator_lookup, 
                  const CorrelatorLookup& corr_lookup,
                  const QuarklineLookup& quark_lookup,
                  const std::string output_path,
                  const std::string output_filename);
  /*! Build charged 4pt correlation function: Cross diagram
   *  @f{align}{
   *    C = \langle \gamma_5 D_\mathtt{Q0}^{-1}(t|t')^\dagger \gamma_5 \Gamma_\mathtt{Op0} 
   *                D_\mathtt{Q1}^{-1}(t|t') \Gamma_\mathtt{Op1}
   *                \gamma_5 D_\mathtt{Q2}^{-1}(t|t')^\dagger \gamma_5 \Gamma_\mathtt{Op2} 
   *                D_\mathtt{Q3}^{-1}(t|t') \Gamma_\mathtt{Op3} \rangle
   *  @f}
   */
  void build_C4cC(const OperatorsForMesons& meson_operator,
                  const Perambulator& perambulators,
                  const OperatorLookup& operator_lookup,
                  const std::vector<CorrInfo>& corr_lookup,
                  const QuarklineLookup& quark_lookup,
                  const std::string output_path,
                  const std::string output_filename);
  /*! Build charged 4pt correlation function: Box diagram
   *  @f{align}{
   *    C = \langle \gamma_5 D_\mathtt{Q0}^{-1}(t|t)^\dagger \gamma_5 \Gamma_\mathtt{Op0} 
   *                D_\mathtt{Q1}^{-1}(t|t') \Gamma_\mathtt{Op1}
   *                \gamma_5 D_\mathtt{Q2}^{-1}(t'|t')^\dagger \gamma_5 \Gamma_\mathtt{Op2} 
   *                D_\mathtt{Q3}^{-1}(t'|t) \Gamma_\mathtt{Op3} \rangle
   *  @f}
   */
  void build_C4cB(const OperatorsForMesons& meson_operator,
                  const Perambulator& perambulators,
                  const OperatorLookup& operator_lookup,
                  const std::vector<CorrInfo>& corr_lookup,
                  const QuarklineLookup& quark_lookup,
                  const std::string output_path,
                  const std::string output_filename);

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
                const QuarklineLookup& quark_lookup,
                const std::string output_path,
                const std::string output_filename);
};

} // end of namespace
