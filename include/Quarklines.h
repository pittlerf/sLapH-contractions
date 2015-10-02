#ifndef QUARKLINES_H_
#define OUARKLINES_H_

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

#include "boost/multi_array.hpp"
#include "Eigen/Dense"

#include "OperatorsForMesons.h"
#include "Perambulator.h"
#include "typedefs.h"

namespace LapH {

class Quarklines {

private:
  // containers for the three types of quark lines
  array_quarkline Q1;
  array_quarkline Q2V;
  array_quarkline Q2L;
  const size_t Lt, dilT, dilE;

  void build_Q1(const Perambulator& peram,
                const OperatorsForMesons& meson_operator,
                const std::vector<QuarklineQ1Indices>& lookup);
  void build_Q2V(const Perambulator& peram,
                 const OperatorsForMesons& meson_operator,
                 const std::vector<QuarklineQ2Indices>& lookup);
  void build_Q2L(const Perambulator& peram,
                 const OperatorsForMesons& meson_operator,
                 const std::vector<QuarklineQ2Indices>& lookup);

public:

  Quarklines (const size_t Lt, const size_t dilT, const size_t dilE, 
              const QuarklineLookup& quarkline_lookuptable,
              const std::vector<RandomIndexCombinationsQ2>& ric_lookup);
  ~Quarklines () {}; // dtor

  // ----------------- INTERFACE FOR BUILDING QUARKLINES -----------------------
  // ---------------------------------------------------------------------------
  void create_quarklines(const Perambulator& peram, 
                         const OperatorsForMesons& meson_operator,
                         const QuarklineLookup& quarkline_lookuptable);

};
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

} // end of namespace

#endif

