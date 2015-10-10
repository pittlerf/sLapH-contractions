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
  array_cd_d2 C2c, C20;
  array_corr corrC, corr0; // this is only needed intermideately
  const size_t Lt, dilT, dilE, nev;

  // this is just for intermediate steps
  void build_corr0(const Quarklines& quarklines, 
                   const std::vector<CorrInfo>& corr_lookup,
                   const QuarklineLookup& quark_lookup,
                   const std::vector<RandomIndexCombinationsQ2>& ric_lookup);
  void build_corrC();

  // functions to build correlation functions
  void build_C20();
  void build_C2c();

public:
  Correlators (const size_t Lt, const size_t dilT, const size_t dilE, 
               const size_t nev, const CorrelatorLookup& corr_lookup);
  ~Correlators () {}; // dtor

  void contract(const Quarklines& quarklines, 
                const OperatorsForMesons& meson_operator,
                const CorrelatorLookup& corr_lookup, 
                const QuarklineLookup& quark_lookup,
                const std::vector<RandomIndexCombinationsQ2>& ric_lookup);
};
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

} // end of namespace

#endif // CORRELATORS_H_ 
