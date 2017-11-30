#include <cassert>
#include <iostream>
#include <stdexcept>
#include <string>

#include "Operator.h"

namespace {

template <QuarkLineType qlt>
void M1xM2(Eigen::MatrixXcd &result, 
           Eigen::MatrixXcd const &M1, 
           std::vector<Eigen::MatrixXcd> const &M2, 
           std::vector<size_t> const &lookup,
           std::vector<RandomIndexCombinationsQ2> const &ricQ2_lookup,
           std::vector<VdaggerVRandomLookup> const &rvdaggervr_lookup,
           std::vector<QuarklineQ2Indices> const &Q2_lookup,
           std::pair<size_t, size_t> const & rnd0,
           std::pair<size_t, size_t> const & rnd1,
           size_t const dilE,
           size_t const dilD);

template <>
void M1xM2<QuarkLineType::Q2V>(Eigen::MatrixXcd &result, 
           Eigen::MatrixXcd const &M1, 
           std::vector<Eigen::MatrixXcd> const &M2, 
           std::vector<size_t> const &lookup,
           std::vector<RandomIndexCombinationsQ2> const &ricQ2_lookup,
           std::vector<VdaggerVRandomLookup> const &rvdaggervr_lookup,
           std::vector<QuarklineQ2Indices> const &Q2_lookup,
           std::pair<size_t, size_t> const & rnd0,
           std::pair<size_t, size_t> const & rnd1,
           size_t const dilE,
           size_t const dilD){

   const auto &ric2 = ricQ2_lookup[Q2_lookup[lookup[2]].id_ric_lookup].rnd_vec_ids;
   const auto &ric3 = ricQ2_lookup[rvdaggervr_lookup[lookup[3]].id_ricQ_lookup].rnd_vec_ids;

   size_t M2_rnd_counter = 0;
   for (const auto &rnd2 : ric2) {
     for (const auto &rnd3 : ric3) {
       if (rnd2.first != rnd3.first && rnd2.second == rnd3.second) {

         // Check that no random vector is used in M1 and M2 at the same 
         // time
         if (rnd1.first == rnd2.first && rnd0.first == rnd3.first &&
             rnd0.second != rnd3.second) {
           result += M2[M2_rnd_counter];
         }
         M2_rnd_counter++;
       }
     }
   }

   result = (M1 * result);
}

template <>
void M1xM2<QuarkLineType::Q2L>(Eigen::MatrixXcd &result, 
           Eigen::MatrixXcd const &M1, 
           std::vector<Eigen::MatrixXcd> const &M2, 
           std::vector<size_t> const &lookup,
           std::vector<RandomIndexCombinationsQ2> const &ricQ2_lookup,
           std::vector<VdaggerVRandomLookup> const &rvdaggervr_lookup,
           std::vector<QuarklineQ2Indices> const &Q2_lookup,
           std::pair<size_t, size_t> const & rnd0,
           std::pair<size_t, size_t> const & rnd1,
           size_t const dilE,
           size_t const dilD){

  const auto &ric2 = ricQ2_lookup[Q2_lookup[lookup[2]].id_ric_lookup].rnd_vec_ids;
  const auto &ric3 = ricQ2_lookup[rvdaggervr_lookup[lookup[1]].id_ricQ_lookup].rnd_vec_ids;

  size_t M2_rnd_counter = 0;
  for (const auto &rnd2 : ric2) {
    for (const auto &rnd3 : ric3) {
      if (rnd2.first == rnd3.first && rnd2.second != rnd3.second) {

        if (rnd0.second == rnd3.second && rnd1.second == rnd2.second &&
            rnd0.first != rnd2.first) {
          result += M2[M2_rnd_counter];
        }
        M2_rnd_counter++;
      }
    }
  }

  result = (M1 * result);
}

} // end of anonymous namespace 

namespace LapH {

template<>
void check_random_combinations<QuarkLineType::Q2V>(std::string const &diagram,
                               std::vector<size_t> const &lookup,
                               std::vector<RandomIndexCombinationsQ2> const &ricQ2_lookup,
                               std::vector<VdaggerVRandomLookup> const &rvdaggervr_lookup,
                               std::vector<QuarklineQ2Indices> const &Q2_lookup){
   const auto &ric0 =
       ricQ2_lookup[Q2_lookup[lookup[0]].id_ric_lookup]
           .rnd_vec_ids;
   const auto &ric1 =
      ricQ2_lookup[rvdaggervr_lookup[lookup[1]].id_ricQ_lookup].rnd_vec_ids;
   const auto &ric2 =
       ricQ2_lookup[Q2_lookup[lookup[2]].id_ric_lookup].rnd_vec_ids;
   const auto &ric3 =
       ricQ2_lookup[rvdaggervr_lookup[lookup[3]].id_ricQ_lookup].rnd_vec_ids;

   if (ric0.size() != ric1.size() || ric0.size() != ric2.size() ||
       ric0.size() != ric3.size()) {
     std::string error_message = 
       std::string("rnd combinations are not the same in ") + diagram;
     throw std::length_error(error_message);
   }
}

template<>
void check_random_combinations<QuarkLineType::Q2L>(std::string const &diagram,
                               std::vector<size_t> const &lookup,
                               std::vector<RandomIndexCombinationsQ2> const &ricQ2_lookup,
                               std::vector<VdaggerVRandomLookup> const &rvdaggervr_lookup,
                               std::vector<QuarklineQ2Indices> const &Q2_lookup){
   const auto &ric0 =
       ricQ2_lookup[Q2_lookup[lookup[0]].id_ric_lookup]
           .rnd_vec_ids;
   const auto &ric1 =
      ricQ2_lookup[rvdaggervr_lookup[lookup[3]].id_ricQ_lookup].rnd_vec_ids;
   const auto &ric2 =
       ricQ2_lookup[Q2_lookup[lookup[2]].id_ric_lookup].rnd_vec_ids;
   const auto &ric3 =
       ricQ2_lookup[rvdaggervr_lookup[lookup[1]].id_ricQ_lookup].rnd_vec_ids;

   if (ric0.size() != ric1.size() || ric0.size() != ric2.size() ||
       ric0.size() != ric3.size()) {
     std::string error_message = 
       std::string("rnd combinations are not the same in ") + diagram;
     throw std::length_error(error_message);
   }
}

template <>
void Q1<QuarkLineType::Q1>(std::vector<Eigen::MatrixXcd> &result, 
                    QuarkLineBlock<QuarkLineType::Q1> const &quarklines,
                    int const t1,
                    int const b2,
                    std::array<size_t, 2> const look,
                    std::vector<RandomIndexCombinationsQ2> const &ricQ2_lookup,
                    std::vector<QuarklineQ1Indices> const &Q1_lookup,
                    size_t const dilE,
                    size_t const dilD){

  /*! @todo Why ricQ2 and not ricQ1? */ 
  const auto &ric1 = ricQ2_lookup[Q1_lookup[look[1]].id_ric_lookup].rnd_vec_ids;

  size_t M2_rnd_counter = 0;
  for (const auto &rnd1 : ric1){
    const auto idr1 = &rnd1 - &ric1[0];
    /*! @Note Allocation should be refactored */
    result.emplace_back(Eigen::MatrixXcd::Zero(dilE * dilD, dilE * dilD));

    result[M2_rnd_counter] = quarklines(t1, b2, look[1], idr1);
    ++M2_rnd_counter;
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
                    std::vector<RandomIndexCombinationsQ2> const &ricQ2_lookup,
                    std::vector<VdaggerVRandomLookup> const &rvdaggervr_lookup,
                    std::vector<QuarklineQ2Indices> const &Q2V_lookup,
                    size_t const dilE,
                    size_t const dilD){

  /*! Assume full dilution in Dirac space in the Loop over d */
  assert(dilD = 4);

  const auto &ric0 =
      ricQ2_lookup[Q2V_lookup[look[1]].id_ric_lookup].rnd_vec_ids;
  const auto &ric1 =
      ricQ2_lookup[rvdaggervr_lookup[look[2]].id_ricQ_lookup].rnd_vec_ids;

  size_t result_rnd_counter = 0;
  for (const auto &rnd0 : ric0) {
    for (const auto &rnd1 : ric1) {
      if (rnd0.first != rnd1.first && rnd0.second == rnd1.second) {

        /*! @Note Allocation should be refactored */
        result.emplace_back(Eigen::MatrixXcd::Zero(dilE * dilD, dilE * dilD));

        const size_t idr0 = &rnd0 - &ric0[0];
        const size_t idr1 = &rnd1 - &ric1[0];
        for (size_t d = 0; d < 4; d++) {
          // TODO: gamma hardcoded
          const cmplx value = quarklines.return_gamma_val(5, d);  
          const size_t gamma_index = quarklines.return_gamma_row(5, d);  
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

template <>
void rVdaggerVrxQ2<QuarkLineType::Q2L>(std::vector<Eigen::MatrixXcd> &result, 
                    QuarkLineBlock<QuarkLineType::Q2L> const &quarklines,
                    OperatorsForMesons const &meson_operator,
                    int const t1,
                    int const b2,
                    std::array<size_t, 3> const look,
                    std::vector<RandomIndexCombinationsQ2> const &ricQ2_lookup,
                    std::vector<VdaggerVRandomLookup> const &rvdaggervr_lookup,
                    std::vector<QuarklineQ2Indices> const &Q2_lookup,
                    size_t const dilE,
                    size_t const dilD){

  /*! Assume full dilution in Dirac space in the Loop over d */
  assert(dilD = 4);

  const auto &ric0 = ricQ2_lookup[Q2_lookup[look[2]].id_ric_lookup].rnd_vec_ids;
  const auto &ric1 = ricQ2_lookup[rvdaggervr_lookup[look[1]].id_ricQ_lookup].rnd_vec_ids;

  size_t result_rnd_counter = 0;
  for (const auto &rnd0 : ric0) {
    for (const auto &rnd1 : ric1) {
      if (rnd0.first == rnd1.first && rnd0.second != rnd1.second) {

        /*! @Note Allocation should be refactored */
        result.emplace_back(Eigen::MatrixXcd::Zero(dilE * dilD, dilE * dilD));

        const size_t idr0 = &rnd0 - &ric0[0];
        const size_t idr1 = &rnd1 - &ric1[0];

        for (size_t d = 0; d < 4; d++) {
          // TODO: gamma hardcoded
          const cmplx value = quarklines.return_gamma_val(5, d);
          const size_t gamma_index = quarklines.return_gamma_row(5, d);

          result[result_rnd_counter].block(d * dilE, 0, dilE, dilD * dilE) =
              value *
              meson_operator.return_rvdaggervr(look[1], t1, idr1)
                  .block(d * dilE, gamma_index * dilE, dilE, dilE) *
              quarklines(t1, b2, look[2], idr0)
                  .block(gamma_index * dilE, 0, dilE, dilD * dilE);
        }
        ++result_rnd_counter;
      }
    }
  }
}



template<>
cmplx trace<QuarkLineType::Q2V, QuarkLineType::Q2V>(std::vector<Eigen::MatrixXcd> const &M1, 
           std::vector<Eigen::MatrixXcd> const &M2, 
           std::vector<size_t> const &lookup,
           std::vector<RandomIndexCombinationsQ2> const &ricQ2_lookup,
           std::vector<VdaggerVRandomLookup> const &rvdaggervr_lookup,
           std::vector<QuarklineQ2Indices> const &Q2_lookup,
           size_t const dilE,
           size_t const dilD){

  /*! @todo unnessary allocation */
   Eigen::MatrixXcd M3 = Eigen::MatrixXcd::Zero(dilE * dilD, dilE * dilD);
   cmplx result = cmplx(.0,.0);

   const auto &ric0 = ricQ2_lookup[Q2_lookup[lookup[0]].id_ric_lookup].rnd_vec_ids;
   const auto &ric1 = ricQ2_lookup[rvdaggervr_lookup[lookup[1]].id_ricQ_lookup].rnd_vec_ids;
   size_t M1_rnd_counter = 0;
   for (const auto &rnd0 : ric0) {
     for (const auto &rnd1 : ric1) {
       if (rnd0.first != rnd1.first && rnd0.second == rnd1.second) {
         // setting matrix values to zero
          M3.setZero(dilE * 4, dilE * 4);

          M1xM2<QuarkLineType::Q2V>(M3, M1[M1_rnd_counter], M2, lookup, 
                                   ricQ2_lookup, rvdaggervr_lookup, 
                                   Q2_lookup, rnd0, rnd1, dilE, 4);

          result += M3.trace();
          ++M1_rnd_counter;
       }
     }
   }

  return result;
}

template<>
cmplx trace<QuarkLineType::Q2L, QuarkLineType::Q2L>(
           std::vector<Eigen::MatrixXcd> const &M1, 
           std::vector<Eigen::MatrixXcd> const &M2, 
           std::vector<size_t> const &lookup,
           std::vector<RandomIndexCombinationsQ2> const &ricQ2_lookup,
           std::vector<VdaggerVRandomLookup> const &rvdaggervr_lookup,
           std::vector<QuarklineQ2Indices> const &Q2_lookup,
           size_t const dilE,
           size_t const dilD){

  /*! @todo unnessary allocation */
  Eigen::MatrixXcd M3 = Eigen::MatrixXcd::Zero(dilE * dilD, dilE * dilD);
  cmplx result = cmplx(.0,.0);

  const auto &ric0 = ricQ2_lookup[Q2_lookup[lookup[0]].id_ric_lookup].rnd_vec_ids;
  const auto &ric1 = ricQ2_lookup[rvdaggervr_lookup[lookup[3]].id_ricQ_lookup].rnd_vec_ids;

  size_t M1_rnd_counter = 0;
  for (const auto &rnd0 : ric0) {
    for (const auto &rnd1 : ric1) {
      if (rnd0.first == rnd1.first && rnd0.second != rnd1.second) {

        // setting matrix values to zero
        M3.setZero(4 * dilE, 4 * dilE);  

        M1xM2<QuarkLineType::Q2L>(M3, M1[M1_rnd_counter], M2, lookup, 
                                 ricQ2_lookup, rvdaggervr_lookup, 
                                 Q2_lookup, rnd0, rnd1, dilE, 4);

        result += M3.trace();
        ++M1_rnd_counter;

        }
    }
  }

  return result;
}

template<>
std::vector<cmplx> trace<QuarkLineType::Q2V, QuarkLineType::Q0>(
    QuarkLineBlock<QuarkLineType::Q2V> const &quarklines,
    OperatorsForMesons const &meson_operator,
    int const t1,
    int const b2,
    int const t2,
    std::vector<size_t> const &lookup,
    std::vector<RandomIndexCombinationsQ2> const &ricQ2_lookup,
    std::vector<VdaggerVRandomLookup> const &rvdaggervr_lookup,
    std::vector<QuarklineQ2Indices> const &Q2_lookup,
    int const gamma,
    size_t const dilE,
    size_t const dilD){

  assert(dilD == 4);

  std::vector<cmplx> result;

  const auto &ric0 = 
    ricQ2_lookup[Q2_lookup[lookup[0]].id_ric_lookup].rnd_vec_ids;
  const auto &ric1 =
    ricQ2_lookup[rvdaggervr_lookup[lookup[1]].id_ricQ_lookup].rnd_vec_ids;

  for (const auto &rnd : ric0) {
    const auto idr0 = &rnd - &ric0[0];
    result.emplace_back(cmplx(0.0, 0.0));

    for (size_t d = 0; d < 4; d++) {
      const auto gamma_index = quarklines.return_gamma_row(gamma, d);
      result[idr0] +=
          quarklines.return_gamma_val(gamma, d) *
          ( quarklines(t1, b2, lookup[0], idr0)
               .block(d * dilE, gamma_index * dilE, dilE, dilE) *
            meson_operator.return_rvdaggervr(lookup[1], t2, idr0)
               .block(gamma_index * dilE, d * dilE, dilE, dilE))
          .trace();
    }
  }

  return result;
}

template<>
std::vector<cmplx> trace<QuarkLineType::Q1, QuarkLineType::Q1>(
    QuarkLineBlock<QuarkLineType::Q1> const &quarklines,
    int const t1,
    int const b2,
    int const t2,
    int const b1,
    std::vector<size_t> const &lookup,
    std::vector<RandomIndexCombinationsQ2> const &ricQ2_lookup,
    std::vector<QuarklineQ1Indices> const &Q1_lookup){

  assert(dilD == 4);

  std::vector<cmplx> result;

  const auto& ric0 =
    ricQ2_lookup[Q1_lookup[lookup[0]].id_ric_lookup].rnd_vec_ids;
  const auto& ric1 = 
    ricQ2_lookup[Q1_lookup[lookup[1]].id_ric_lookup].rnd_vec_ids;

  for (const auto& rnd : ric0) {
    const auto idr0 = &rnd - &ric0[0];
    result.emplace_back(cmplx(0.0, 0.0));

    // check that ric1 and ric0 are indeed a full trace 
    // i.e. ric1.first == ric2.second and ric1.second == ric2.first
    const auto it1 = std::find_if(
        ric1.begin(), ric1.end(), [&](std::pair<size_t, size_t> pair) {
          return (pair == std::make_pair(rnd.second, rnd.first));
        });
    if (it1 == ric1.end()) {
      std::cout << "something wrong with random vectors in build_corr0"
                << std::endl;
      exit(1);
    }
    const auto idr1 = it1 - ric1.begin();

    /*! @todo How do I properly get the block indices for sink? */
    result[idr0] += (quarklines(t1, b2, lookup[0], idr0) * 
                     quarklines(t2, b1, lookup[1], idr1))
                    .trace();
  }

  return result;
}


template<>
cmplx trace<QuarkLineType::Q2L, QuarkLineType::Q1>(
           std::vector<Eigen::MatrixXcd> const &M1, 
           std::vector<Eigen::MatrixXcd> const &M2, 
           std::vector<size_t> const &lookup,
           std::vector<RandomIndexCombinationsQ2> const &ricQ2_lookup,
           std::vector<RandomIndexCombinationsQ1> const &ricQ1_lookup,
           std::vector<VdaggerVRandomLookup> const &rvdaggervr_lookup,
           std::vector<QuarklineQ1Indices> const &Q1_lookup,
           std::vector<QuarklineQ2Indices> const &Q2_lookup,
           size_t const dilE,
           size_t const dilD){

  cmplx result = cmplx(.0,.0);

  const auto &ric0 = 
    ricQ2_lookup[Q2_lookup[lookup[0]].id_ric_lookup].rnd_vec_ids;
  const auto &ric1 = 
    ricQ2_lookup[Q1_lookup[lookup[1]].id_ric_lookup].rnd_vec_ids;
  const auto &ric2 = 
    ricQ2_lookup[rvdaggervr_lookup[lookup[2]].id_ricQ_lookup].rnd_vec_ids;

  size_t M1_rnd_counter = 0;
  for (const auto &rnd0 : ric0) {
    for (const auto &rnd2 : ric2) {
      if (rnd0.first == rnd2.first && rnd0.second != rnd2.second) {
        size_t M2_rnd_counter = 0;
        for (const auto &rnd1 : ric1) {
          const auto idr1 = &rnd1 - &ric1[0];
          if (rnd1.first != rnd2.first && rnd1.second == rnd2.second &&
              rnd1.first == rnd0.second && rnd1.second != rnd0.first) {
            result += (M1[M1_rnd_counter] *
                                M2[M2_rnd_counter]).trace();
          }
          ++M2_rnd_counter;
        }
      ++M1_rnd_counter;
      }
    }
  }

  return result;
}

}  // end of LapH namespace
