#include <cassert>
#include <iostream>
#include <stdexcept>
#include <string>

#include "DilutedFactor.h"

namespace {

void M1xM2(Eigen::MatrixXcd &result, 
           Eigen::MatrixXcd const &M1, 
           std::vector<Eigen::MatrixXcd> const &M2, 
           std::vector<RandomIndexCombinationsQ2> const &ric_lookup,
           std::vector<size_t> const &ric_ids,
           std::pair<size_t, size_t> const & rnd0,
           std::pair<size_t, size_t> const & rnd1,
           size_t const dilE,
           size_t const dilD){

  const auto &ric2 = ric_lookup[ric_ids[2]].rnd_vec_ids;
  const auto &ric3 = ric_lookup[ric_ids[3]].rnd_vec_ids;

  size_t M2_rnd_counter = 0;
  for (const auto &rnd2 : ric2) {
    for (const auto &rnd3 : ric3) {
      if (rnd2.second == rnd3.first && rnd2.first != rnd3.second) {

        if (rnd0.first == rnd3.second && rnd1.second == rnd2.first &&
            rnd0.second != rnd3.first) {
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
                               std::vector<RandomIndexCombinationsQ2> const &ric_lookup,
                               std::vector<VdaggerVRandomLookup> const &rvdaggervr_lookup,
                               std::vector<QuarklineQ2Indices> const &Q2_lookup){
   const auto &ric0 =
       ric_lookup[Q2_lookup[lookup[0]].id_ric_lookup]
           .rnd_vec_ids;
   const auto &ric1 =
      ric_lookup[rvdaggervr_lookup[lookup[1]].id_ric_lookup].rnd_vec_ids;
   const auto &ric2 =
       ric_lookup[Q2_lookup[lookup[2]].id_ric_lookup].rnd_vec_ids;
   const auto &ric3 =
       ric_lookup[rvdaggervr_lookup[lookup[3]].id_ric_lookup].rnd_vec_ids;

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
                               std::vector<RandomIndexCombinationsQ2> const &ric_lookup,
                               std::vector<VdaggerVRandomLookup> const &rvdaggervr_lookup,
                               std::vector<QuarklineQ2Indices> const &Q2_lookup){
   const auto &ric0 =
       ric_lookup[Q2_lookup[lookup[0]].id_ric_lookup]
           .rnd_vec_ids;
   const auto &ric1 =
      ric_lookup[rvdaggervr_lookup[lookup[3]].id_ric_lookup].rnd_vec_ids;
   const auto &ric2 =
       ric_lookup[Q2_lookup[lookup[2]].id_ric_lookup].rnd_vec_ids;
   const auto &ric3 =
       ric_lookup[rvdaggervr_lookup[lookup[1]].id_ric_lookup].rnd_vec_ids;

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
                    std::vector<RandomIndexCombinationsQ2> const &ric_lookup,
                    std::vector<QuarklineQ1Indices> const &Q1_lookup,
                    size_t const dilE,
                    size_t const dilD){

  /*! @todo Why ricQ2 and not ricQ1? */ 
  const auto &ric1 = ric_lookup[Q1_lookup[look[1]].id_ric_lookup].rnd_vec_ids;

  size_t M2_rnd_counter = 0;
  for (const auto &rnd1 : ric1){
    const auto idr1 = &rnd1 - &ric1[0];
    /*! @Note Allocation should be refactored */
    result.emplace_back(Eigen::MatrixXcd::Zero(dilE * dilD, dilE * dilD));

    result[M2_rnd_counter] = quarklines(t1, b2, look[1], idr1);
    ++M2_rnd_counter;
  }
}

template <>
void Q1xQ1<QuarkLineType::Q1>(
    std::vector<Eigen::MatrixXcd> &result, 
    QuarkLineBlock<QuarkLineType::Q1> const &quarklines,
    int const t1,
    int const b1,
    int const t2,
    int const b2,
    std::array<size_t, 3> const look,
    std::vector<RandomIndexCombinationsQ2> const &ric_lookup,
    std::vector<QuarklineQ1Indices> const &Q1_lookup,
    size_t const dilE,
    size_t const dilD){

  const auto &ric0 =
      ric_lookup[Q1_lookup[look[1]].id_ric_lookup].rnd_vec_ids;
  const auto &ric1 =
      ric_lookup[Q1_lookup[look[2]].id_ric_lookup].rnd_vec_ids;

  size_t result_rnd_counter = 0;
  for (const auto &rnd0 : ric0) {
    for (const auto &rnd1 : ric1) {

      auto const idr0 = &rnd0 - &ric0[0];
      auto const idr1 = &rnd1 - &ric1[0];

      if (rnd0.second == rnd1.first && rnd0.first != rnd1.second) {

        result.emplace_back(Eigen::MatrixXcd::Zero(dilE * dilD, dilE * dilD));

        result[result_rnd_counter] = quarklines(t1, b1, look[1], idr0) *
                        quarklines(t2, b2, look[2], idr1);
        ++result_rnd_counter;
      }
    }
  }

}

/*! 
 *  Multiply Operator with 0 quarks and Operator with 2 quarks
 *
 *  Dependency inversion principle: This interface will be a policy that must be 
 *  fullfilled by any linear algebra library to be used. 
 *
 *  @todo Layer of abstraction for Eigen in call parameter result
 */
void rVdaggerVrxQ2(std::vector<Eigen::MatrixXcd> &result, 
                   std::vector<Eigen::MatrixXcd> const &quarkline1,
                   std::vector<Eigen::MatrixXcd> const &quarkline2,
                   std::vector<RandomIndexCombinationsQ2> const &ric_lookup,
                   std::vector<size_t> const &ric_ids,
                   size_t const dilE,
                   size_t const dilD){

  /*! Assume full dilution in Dirac space in the Loop over d */
  assert(dilD == 4);

  gamma_lookup gamma_5{};
  gamma_5.row[0] = 0;
  gamma_5.value[0] = 1;
  gamma_5.row[1] = 1;
  gamma_5.value[1] = 1;
  gamma_5.row[2] = 2;
  gamma_5.value[2] = -1;
  gamma_5.row[3] = 3;
  gamma_5.value[3] = -1;


  const auto &ric0 = ric_lookup[ric_ids[0]].rnd_vec_ids;
  const auto &ric1 = ric_lookup[ric_ids[1]].rnd_vec_ids;

  size_t result_rnd_counter = 0;
  for (const auto &rnd0 : ric0) {
    for (const auto &rnd1 : ric1) {
      if (rnd0.second == rnd1.first && rnd0.first != rnd1.second) {

        /*! @Note Allocation should be refactored */
        result.emplace_back(Eigen::MatrixXcd::Zero(dilE * dilD, dilE * dilD));

        const size_t idr0 = &rnd0 - &ric0[0];
        const size_t idr1 = &rnd1 - &ric1[0];

        for (size_t d = 0; d < 4; d++) {
          // TODO: gamma hardcoded
          const cmplx value = gamma_5.value[d];
          const size_t gamma_index = gamma_5.row[d];

          result[result_rnd_counter].block(d * dilE, 0, dilE, dilD * dilE) =
              value *
              quarkline1[idr0].block(d * dilE, gamma_index * dilE, dilE, dilE) *
              quarkline2[idr1].block(gamma_index * dilE, 0, dilE, dilD * dilE);
        }
        ++result_rnd_counter;
      }
    }
  }

}

cmplx trace(std::vector<Eigen::MatrixXcd> const &M1, 
            std::vector<Eigen::MatrixXcd> const &M2, 
            std::vector<RandomIndexCombinationsQ2> const &ric_lookup,
            std::vector<size_t> const &ric_ids,
            size_t const dilE,
            size_t const dilD){

  /*! @todo unnessary allocation */
   Eigen::MatrixXcd M3 = Eigen::MatrixXcd::Zero(dilE * dilD, dilE * dilD);
   cmplx result = cmplx(.0,.0);

   const auto &ric0 = ric_lookup[ric_ids[0]].rnd_vec_ids;
   const auto &ric1 = ric_lookup[ric_ids[1]].rnd_vec_ids;

   size_t M1_rnd_counter = 0;
   for (const auto &rnd0 : ric0) {
     for (const auto &rnd1 : ric1) {
       if (rnd0.second == rnd1.first && rnd0.first != rnd1.second) {
         // setting matrix values to zero
          M3.setZero(dilE * 4, dilE * 4);

          M1xM2(M3, M1[M1_rnd_counter], M2, ric_lookup, ric_ids, rnd0, rnd1, dilE, 4);

          result += M3.trace();
          ++M1_rnd_counter;
       }
     }
   }

  return result;
}

/*! corrC */
template<>
std::vector<cmplx> trace<QuarkLineType::Q2V, QuarkLineType::Q0>(
    QuarkLineBlock<QuarkLineType::Q2V> const &quarkline1,
    QuarkLineBlock<QuarkLineType::Q0> const &quarkline2,
    int const t1,
    int const b2,
    int const t2,
    std::vector<size_t> const &lookup,
    std::vector<RandomIndexCombinationsQ2> const &ric_lookup,
    std::vector<size_t> const &ric_ids,
    int const gamma,
    size_t const dilE,
    size_t const dilD){

  assert(dilD == 4);

  std::vector<cmplx> result;

  const auto &ric0 = ric_lookup[ric_ids[0]].rnd_vec_ids;
  const auto &ric1 = ric_lookup[ric_ids[1]].rnd_vec_ids;

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
      std::cout << "something wrong with random vectors in build_corrC"
                << std::endl;
      exit(1);
    }
    const auto idr1 = it1 - ric1.begin();

    for (size_t d = 0; d < 4; d++) {
      const auto gamma_index = quarkline1.return_gamma_row(gamma, d);
      result[idr0] +=
          quarkline1.return_gamma_val(gamma, d) *
          ( quarkline1(t1, b2, lookup[0], idr0)
               .block(d * dilE, gamma_index * dilE, dilE, dilE) *
            /*! @warning idr0 instead of idr1 because rvdaggervr are interchanged */
            quarkline2(t2, -1, lookup[1], idr1)
               .block(gamma_index * dilE, d * dilE, dilE, dilE))
          .trace();
    }
  }

  return result;
}

/*! corr0 */
std::vector<cmplx> trace(
    std::vector<Eigen::MatrixXcd> const &quarkline1,
    std::vector<Eigen::MatrixXcd> const &quarkline2,
    std::vector<RandomIndexCombinationsQ2> const &ric_lookup,
    std::vector<size_t> const &ric_ids){

  std::vector<cmplx> result;

  const auto& ric0 =
    ric_lookup[ric_ids[0]].rnd_vec_ids;
  const auto& ric1 = 
    ric_lookup[ric_ids[1]].rnd_vec_ids;

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
    result[idr0] += (quarkline1[idr0] * 
                     quarkline2[idr1])
                    .trace();
  }

  return result;
}


/*! C3c */
template<>
cmplx trace<QuarkLineType::Q2L, QuarkLineType::Q1>(
           std::vector<Eigen::MatrixXcd> const &M1, 
           std::vector<Eigen::MatrixXcd> const &M2, 
           std::vector<size_t> const &lookup,
           std::vector<RandomIndexCombinationsQ2> const &ric_lookup,
           std::vector<VdaggerVRandomLookup> const &rvdaggervr_lookup,
           std::vector<QuarklineQ1Indices> const &Q1_lookup,
           std::vector<QuarklineQ2Indices> const &Q2_lookup,
           size_t const dilE,
           size_t const dilD){

  cmplx result = cmplx(.0,.0);

  const auto &ric0 = 
    ric_lookup[Q2_lookup[lookup[0]].id_ric_lookup].rnd_vec_ids;
  const auto &ric1 = 
    ric_lookup[Q1_lookup[lookup[1]].id_ric_lookup].rnd_vec_ids;
  const auto &ric2 = 
    ric_lookup[rvdaggervr_lookup[lookup[2]].id_ric_lookup].rnd_vec_ids;

  size_t M2_rnd_counter = 0;
  for (const auto &rnd1 : ric1) {

    size_t M1_rnd_counter = 0;
    for (const auto &rnd2 : ric2) {
      for (const auto &rnd0 : ric0) {
        if (rnd0.first == rnd2.second && rnd0.second != rnd2.first) {

          if (rnd1.first != rnd2.second && rnd1.second == rnd2.first &&
              rnd1.first == rnd0.second && rnd1.second != rnd0.first) {

          result += (M1[M1_rnd_counter] *
                                M2[M2_rnd_counter]).trace();
          }
          ++M1_rnd_counter;
        }
      }
    }
    ++M2_rnd_counter;
  }


  return result;
}

cmplx trace_3n(std::vector<Eigen::MatrixXcd> const &L1, 
               std::vector<Eigen::MatrixXcd> const &L2, 
               std::vector<RandomIndexCombinationsQ2> const &ric_lookup,
               std::vector<size_t> const &ric_ids){

   cmplx result = cmplx(0.,0.);

   const auto& ric0 = ric_lookup[ric_ids[0]].rnd_vec_ids;
   const auto& ric1 = ric_lookup[ric_ids[1]].rnd_vec_ids;
   const auto& ric2 = ric_lookup[ric_ids[2]].rnd_vec_ids;

   size_t L1_rnd_counter = 0;
   for (const auto &rnd0 : ric0) {
     for (const auto &rnd1 : ric1) {
       if (rnd0.second == rnd1.first && rnd0.first != rnd1.second) {

        size_t L2_rnd_counter = 0;
        for (const auto &rnd2 : ric2) {
          if (rnd1.second == rnd2.first && rnd2.second == rnd0.first) {
            result += (L1[L1_rnd_counter] *L2[L2_rnd_counter]).trace();
          }
          ++L2_rnd_counter;
        }

        ++L1_rnd_counter;
      }
    }
  }

  return result;
}

/*!
 *  - C40D
 *  - C40V
 *  - C40D
 *  - C40V
 */
compcomp_t trtr(std::vector<cmplx> const &factor1,
          std::vector<cmplx> const &factor2,
          std::vector<RandomIndexCombinationsQ2> const &ric_lookup,
          std::vector<size_t> const &ric_ids){

  compcomp_t result {compcomp_t(.0, .0, .0, .0)};

  const auto &ric0 = ric_lookup[ric_ids[0]].rnd_vec_ids;
  const auto &ric1 = ric_lookup[ric_ids[1]].rnd_vec_ids;

  for (const auto &rnd0 : ric0) {
    for (const auto &rnd1 : ric1) {
      /*! No index in the first trace may appear in the second trace */
      if ((rnd0.first != rnd1.first) && (rnd0.first != rnd1.second) &&
          (rnd0.second != rnd1.first) && (rnd0.second != rnd1.second)) {

        auto const idr0 = &rnd0 - &ric0[0];
        auto const idr1 = &rnd1 - &ric1[0];
  
        result.rere += factor1[idr0].real() * factor2[idr1].real();
        result.reim += factor1[idr0].real() * factor2[idr1].imag();
        result.imre += factor1[idr0].imag() * factor2[idr1].real();
        result.imim += factor1[idr0].imag() * factor2[idr1].imag();
      }
    }
  }

  return result;
}


}  // end of LapH namespace
