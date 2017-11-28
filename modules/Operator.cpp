#include <iostream>
#include <stdexcept>
#include <string>

#include "Operator.h"

namespace LapH {

void check_random_combinations(std::string const &diagram,
                               std::vector<size_t> const &lookup,
                               OperatorLookup const &operator_lookup,
                               QuarklineLookup const &quark_lookup){
   const auto &ric0 =
       operator_lookup.ricQ2_lookup[quark_lookup.Q2V[lookup[0]].id_ric_lookup]
           .rnd_vec_ids;
   const auto &ric1 =
       operator_lookup
           .ricQ2_lookup[  
               operator_lookup.rvdaggervr_lookuptable[lookup[1]].id_ricQ_lookup].rnd_vec_ids;
   const auto &ric2 =
       operator_lookup.ricQ2_lookup[quark_lookup.Q2V[lookup[2]].id_ric_lookup]
           .rnd_vec_ids;
   const auto &ric3 =
       operator_lookup.ricQ2_lookup[ 
               operator_lookup.rvdaggervr_lookuptable[lookup[3]].id_ricQ_lookup].rnd_vec_ids;

   if (ric0.size() != ric1.size() || ric0.size() != ric2.size() ||
       ric0.size() != ric3.size()) {
     std::string error_message = 
       std::string("rnd combinations are not the same in ") + diagram;
     throw std::length_error(error_message);
   }
}

/*! 
 *  Multiply Operator with 2 Quarks and Operator with zero quarks
 *
 *  Dependency inversion principle: This interface will be a policy that must be 
 *  fullfilled by any linear algebra library to be used. 
 *
 *  @todo Layer of abstraction for Eigen in call parameter result
 */
template <>
void Q2xrVdaggerVr<QuarkLineType::Q2V>(std::vector<Eigen::MatrixXcd> &result, 
                    QuarkLineBlock<QuarkLineType::Q2V> const &quarklines,
                    OperatorsForMesons const &meson_operator,
                    int const b2,
                    int const t1,
                    int const t2,
                    std::array<size_t, 3> const look,
                    OperatorLookup const &operator_lookup,
                    QuarklineLookup const &quark_lookup,
                    size_t const dilE,
                    size_t const dilD){


  const auto &ric0 =
      operator_lookup.ricQ2_lookup[quark_lookup.Q2V[look[1]].id_ric_lookup]
          .rnd_vec_ids;
  const auto &ric1 =
      operator_lookup
          .ricQ2_lookup[  // needed only for checking
              operator_lookup.rvdaggervr_lookuptable[look[2]].id_ricQ_lookup]
          .rnd_vec_ids;
  size_t result_rnd_counter = 0;
  for (const auto &rnd0 : ric0) {
    for (const auto &rnd1 : ric1) {
      if (rnd0.first != rnd1.first && rnd0.second == rnd1.second) {
        const size_t idr0 = &rnd0 - &ric0[0];
        const size_t idr1 = &rnd1 - &ric1[0];
        for (size_t d = 0; d < 4; d++) {
          const cmplx value =
              quarklines.return_gamma_val(5, d);  // TODO: gamma hardcoded
          const size_t gamma_index =
              quarklines.return_gamma_row(5, d);  // TODO: gamma hardcoded
          //          const cmplx value =
          //          quarklines.return_gamma_val(c_look.gamma[0], d);
          //          const size_t gamma_index = quarklines.return_gamma_row(
          //                                                          c_look.gamma[0],
          //                                                          d);
          result[result_rnd_counter].block(0, d * dilE, dilD * dilE, dilE) =
              value *
              quarklines(
                  t1, b2, look[1], idr0)
                  .block(0, gamma_index * dilE, dilD * dilE, dilE) *
              meson_operator.return_rvdaggervr(look[2], t2, idr1)
                  .block(gamma_index * dilE, d * dilE, dilE, dilE);
        }
        result_rnd_counter++;
      }
    }
  }
}

}  // end of LapH namespace
