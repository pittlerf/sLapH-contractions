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

struct gamma_lookup {
  std::array<int, 4> row;
  std::array<cmplx, 4> value;
};

class Quarklines {

private:
  // containers for the three types of quark lines
  array_quarkline Q1;
  array_quarkline Q2V;
  array_quarkline Q2L;
  const size_t Lt, dilT, dilE, nev;
  std::vector<LapH::gamma_lookup>  gamma;

  void build_Q1(const Perambulator& peram,
                const OperatorsForMesons& meson_operator,
                const std::vector<QuarklineQ1Indices>& ql_lookup,
                const std::vector<RandomIndexCombinationsQ2>& ric_lookup);
  void build_Q2V(const Perambulator& peram,
                 const OperatorsForMesons& meson_operator,
                 const std::vector<QuarklineQ2Indices>& ql_lookup,
                 const std::vector<RandomIndexCombinationsQ2>& ric_lookup);
  void build_Q2L(const Perambulator& peram,
                 const OperatorsForMesons& meson_operator,
                 const std::vector<QuarklineQ2Indices>& ql_lookup,
                 const std::vector<RandomIndexCombinationsQ2>& ric_lookup);

public:

  Quarklines (const size_t Lt, const size_t dilT, const size_t dilE, 
              const size_t nev, const QuarklineLookup& quarkline_lookuptable,
              const std::vector<RandomIndexCombinationsQ2>& ric_lookup);
  ~Quarklines () {}; // dtor

  inline const Eigen::MatrixXcd& return_Q1(const size_t t1, const size_t t2,
                                   const size_t op_id, const size_t rnd) const {
    return Q1[t1][t2][op_id].at(rnd);
  }
  // ----------------- INTERFACE FOR BUILDING QUARKLINES -----------------------
  // ---------------------------------------------------------------------------
  void create_quarklines(const Perambulator& peram, 
                      const OperatorsForMesons& meson_operator,
                      const QuarklineLookup& quarkline_lookuptable,
                      const std::vector<RandomIndexCombinationsQ2>& ric_lookup);

};
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

} // end of namespace

#endif

