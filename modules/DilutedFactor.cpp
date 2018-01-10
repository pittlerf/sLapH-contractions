#include "DilutedFactor.h"

#include <cassert>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>

namespace {

/*! M1xM2 for 4pt functions */
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

/*! M1xM2 for 3pt functions */
void M1xM2(Eigen::MatrixXcd &result, 
           Eigen::MatrixXcd const &M1, 
           std::vector<Eigen::MatrixXcd> const &M2, 
           std::vector<RandomIndexCombinationsQ2> const &ric_lookup,
           std::vector<size_t> const &ric_ids,
           std::pair<size_t, size_t> const & rnd1,
           size_t const dilE,
           size_t const dilD){

  const auto &ric0 = ric_lookup[ric_ids[0]].rnd_vec_ids;
  const auto &ric2 = ric_lookup[ric_ids[2]].rnd_vec_ids;

  size_t M2_rnd_counter = 0;
  for (const auto &rnd2 : ric2) {
    for (const auto &rnd0 : ric0) {
      if (rnd0.first == rnd2.second && rnd0.second != rnd2.first) {

        if (rnd1.first != rnd2.second && rnd1.second == rnd2.first &&
            rnd1.first == rnd0.second && rnd1.second != rnd0.first) {

        result += M2[M2_rnd_counter];
        }
      ++M2_rnd_counter;
      }
    }
  }

  result = (M1 * result);
}

} // end of anonymous namespace 

void Q1(std::vector<Eigen::MatrixXcd> &result, 
                    std::vector<Eigen::MatrixXcd> const &quarklines,
                    std::vector<RandomIndexCombinationsQ2> const &ric_lookup,
                    std::vector<size_t> const &ric_ids,
                    size_t const dilE,
                    size_t const dilD){

  const auto &ric0 = ric_lookup[ric_ids[0]].rnd_vec_ids;

  size_t M2_rnd_counter = 0;
  for (const auto &rnd0 : ric0){
    const auto idr0 = &rnd0 - &ric0[0];
    /*! @Note Allocation should be refactored */
    result.emplace_back(Eigen::MatrixXcd::Zero(dilE * dilD, dilE * dilD));

    result[M2_rnd_counter] = quarklines[idr0];
    ++M2_rnd_counter;
  }
}

cmplx trace_3pt(std::vector<Eigen::MatrixXcd> const &M2,
                std::vector<Eigen::MatrixXcd> const &M1,
                std::vector<RandomIndexCombinationsQ2> const &ric_lookup,
                std::vector<size_t> const &ric_ids,
                size_t const dilE,
                size_t const dilD) {
  Eigen::MatrixXcd M3 = Eigen::MatrixXcd::Zero(dilE * dilD, dilE * dilD);
  cmplx result = cmplx(.0,.0);

  const auto &ric1 = ric_lookup[ric_ids[1]].rnd_vec_ids;

  size_t M2_rnd_counter = 0;
  for (const auto &rnd1 : ric1) {

    // setting matrix values to zero
    M3.setZero(dilE * 4, dilE * 4);

    M1xM2(M3, M2[M2_rnd_counter], M1, ric_lookup, ric_ids, rnd1, dilE, dilD);

    result += M3.trace();

    ++M2_rnd_counter;
  }

  return result;
}

cmplx trace(std::vector<Eigen::MatrixXcd> const &M1,
            std::vector<Eigen::MatrixXcd> const &M2,
            std::vector<RandomIndexCombinationsQ2> const &ric_lookup,
            std::vector<size_t> const &ric_ids,
            size_t const dilE,
            size_t const dilD) {
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


/*! corr0, corrC */
std::vector<cmplx> trace(std::vector<Eigen::MatrixXcd> const &quarkline1,
                         std::vector<Eigen::MatrixXcd> const &quarkline2,
                         std::vector<RandomIndexCombinationsQ2> const &ric_lookup,
                         std::vector<size_t> const &ric_ids) {
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
    result[idr0] += (quarkline1[idr0] * quarkline2[idr1]).trace();
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
                std::vector<size_t> const &ric_ids) {
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
