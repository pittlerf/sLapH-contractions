#ifndef CORRELATORS_H_
#define CORRELATORS_H_

#include <complex>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "boost/multi_array.hpp"
#include "Eigen/Dense"

#include "OperatorsForMesons.h"
#include "Quarklines.h"
#include "typedefs.h"

namespace LapH {

class Correlators {

private:
  array_cd_d2 C2c, C3c, C4cD, C4cV, C4cC, C4cB, C20, C30, C40D, C40V, C40C, C40B;
  array_C1 C1, C1T;
  array_corr corrC, corr0; // this is only needed intermideately
  const size_t Lt, dilT, dilE, nev;

  // this is just for intermediate steps
  void build_corr0(const Quarklines& quarklines, 
                   const std::vector<CorrInfo>& corr_lookup,
                   const QuarklineLookup& quark_lookup,
                   const std::vector<RandomIndexCombinationsQ2>& ric_lookup);
  void build_corrC(const Quarklines& quarklines, 
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
  void build_C4cB(const Quarklines& quarklines, 
                  const OperatorsForMesons& meson_operator,
                  const OperatorLookup& operator_lookup,
                  const std::vector<CorrInfo>& corr_lookup,
                  const QuarklineLookup& quark_lookup);

public:
  Correlators (const size_t Lt, const size_t dilT, const size_t dilE, 
               const size_t nev, const CorrelatorLookup& corr_lookup);
  ~Correlators () {}; // dtor

  void contract(const Quarklines& quarklines, 
                const OperatorsForMesons& meson_operator,
                const OperatorLookup& operator_lookup,
                const CorrelatorLookup& corr_lookup, 
                const QuarklineLookup& quark_lookup);
};
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

} // end of namespace

#endif // CORRELATORS_H_ 
