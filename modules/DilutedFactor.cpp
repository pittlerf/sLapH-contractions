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

bool has_intersection(SmallVectorRndId const &left, SmallVectorRndId const &right) {
  SmallVectorRndId intersection;

  std::set_intersection(std::begin(left),
                        std::end(left),
                        std::begin(right),
                        std::end(right),
                        std::back_inserter(intersection));
  return intersection.size() > 0;
}

void merge_append(SmallVectorRndId &data, SmallVectorRndId const &addition) {
  auto const old_end = std::end(data);
  std::copy(std::begin(addition), std::end(addition), std::back_inserter(data));
  std::inplace_merge(std::begin(data), old_end, std::end(data));
}

void merge_push_back(SmallVectorRndId &data, RndId const &addition) {
  auto const old_end = std::end(data);
  data.push_back(addition);
  std::inplace_merge(std::begin(data), old_end, std::end(data));
}

} // end of anonymous namespace 

std::vector<DilutedFactor> operator*(std::vector<DilutedFactor> const &left_vec,
                                     std::vector<DilutedFactor> const &right_vec) {
  //! @TODO Pull out this magic number.
  auto constexpr rnd_vec_count = 5;
  auto constexpr dilD = 4;

  std::vector<DilutedFactor> result_vec;

  for (auto const &left : left_vec) {
    auto const inner_rnd_id = left.ric.second;

    for (auto const &right : right_vec) {
      // We want to make the inner and outer indices differ. The inner indices need to
      // match because the product would not make sense otherwise.
      bool const is_allowed =
          inner_rnd_id == right.ric.first && left.ric.first != right.ric.second;
      if (!is_allowed) {
          continue;
      }

      // We also need to be careful to not combine factors which have common used random
      // vector indices.
      if (has_intersection(left.used_rnd_ids, right.used_rnd_ids)) {
        continue;
      }

      // We want to keep track of the indices that have been contracted away. These are
      // all the ones from the left factor, all the ones from the right factor and the one
      // that we are contracting over right now.
      SmallVectorRndId used;
      used.push_back(inner_rnd_id);
      merge_append(used, left.used_rnd_ids);
      merge_append(used, right.used_rnd_ids);

      result_vec.push_back({Eigen::MatrixXcd{left.data * right.data},
                            std::make_pair(left.ric.first, right.ric.second),
                            used});
    }
  }

  return result_vec;
}

template<>
void check_random_combinations<QuarkLineType::Q2V>(std::string const &diagram,
                               std::vector<size_t> const &lookup,
                               std::vector<RandomIndexCombinationsQ2> const &ric_lookup,
                               std::vector<VdaggerVRandomLookup> const &rvdaggervr_lookup,
                               std::vector<QuarklineQ2Indices> const &Q2_lookup){
  const auto &ric0 = ric_lookup[Q2_lookup[lookup[0]].id_ric_lookup].rnd_vec_ids;
  const auto &ric1 = ric_lookup[rvdaggervr_lookup[lookup[1]].id_ric_lookup].rnd_vec_ids;
  const auto &ric2 = ric_lookup[Q2_lookup[lookup[2]].id_ric_lookup].rnd_vec_ids;
  const auto &ric3 = ric_lookup[rvdaggervr_lookup[lookup[3]].id_ric_lookup].rnd_vec_ids;

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
  const auto &ric0 = ric_lookup[Q2_lookup[lookup[0]].id_ric_lookup].rnd_vec_ids;
  const auto &ric1 = ric_lookup[rvdaggervr_lookup[lookup[3]].id_ric_lookup].rnd_vec_ids;
  const auto &ric2 = ric_lookup[Q2_lookup[lookup[2]].id_ric_lookup].rnd_vec_ids;
  const auto &ric3 = ric_lookup[rvdaggervr_lookup[lookup[1]].id_ric_lookup].rnd_vec_ids;

  if (ric0.size() != ric1.size() || ric0.size() != ric2.size() ||
      ric0.size() != ric3.size()) {
    std::string error_message =
        std::string("rnd combinations are not the same in ") + diagram;
    throw std::length_error(error_message);
   }
}

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
          //! @TODO: gamma hardcoded
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

std::vector<DilutedTrace> factor_to_trace(std::vector<DilutedFactor> const &vec) {
  std::vector<DilutedTrace> result_vec;

  for (auto const &elem : vec) {
    // We only want to use diagonal elements.
    if (elem.ric.first != elem.ric.second) {
      continue;
    }

    DilutedTrace result = {elem.data.trace(), elem.used_rnd_ids};

    auto const outer_rnd_id = elem.ric.first;
    merge_push_back(result.used_rnd_ids, outer_rnd_id);

    result_vec.push_back(result);
  }

  return result_vec;

}

std::vector<DilutedTrace> factor_to_trace(std::vector<DilutedFactor> const &left_vec,
                                          std::vector<DilutedFactor> const &right_vec) {
  //! @TODO Pull out this magic number.
  auto constexpr rnd_vec_count = 5;

  std::vector<DilutedTrace> result_vec;

  for (auto const &left : left_vec) {
    auto const outer_rnd_id = left.ric.first;
    auto const inner_rnd_id = left.ric.second;

    Eigen::MatrixXcd right_sum(
        Eigen::MatrixXcd::Zero(left.data.rows(), left.data.cols()));

    for (auto const &right : right_vec) {
      // We want to make the inner and outer indices match. The inner indices need to
      // match because the product would not make sense otherwise. The outer indices must
      // match since we want to be able to take the trace over the result. The second
      // condition is where this differs from the other multiplication operator.
      bool const is_allowed =
          inner_rnd_id == right.ric.first && outer_rnd_id == right.ric.second;
      if (!is_allowed) {
          continue;
      }

      // We also need to be careful to not combine factors which have common used random
      // vector indices.
      if (has_intersection(left.used_rnd_ids, right.used_rnd_ids)) {
        continue;
      }

      // We want to keep track of the indices that have been contracted away. These are
      // all the ones from the left factor, all the ones from the right factor and the one
      // that we are contracting over right now.
      SmallVectorRndId used;
      merge_push_back(used, inner_rnd_id);
      merge_push_back(used, outer_rnd_id);
      merge_append(used, left.used_rnd_ids);
      merge_append(used, right.used_rnd_ids);

      // The right sides that we encounter at this point have the same left and right
      // random vector indices. They may differ in the set of used random vector indices.
      // But since we do not plan to contract the result with more DilutedFactor
      // instances, we do not care to preserve the information about the used random
      // vector indices. Therefore we can sum all these elements up to have less
      // multiplications to do.
      result_vec.push_back({DilutedTrace::Data{(left.data * right.data).trace()}, used});
    }
  }

  return result_vec;

}

cmplx trace(std::vector<DilutedFactor> const &left_vec,
            std::vector<DilutedFactor> const &right_vec) {
  //! @TODO Pull out this magic number.
  auto constexpr rnd_vec_count = 5;

  cmplx result(0.0, 0.0);

  for (auto const &left : left_vec) {
    auto const outer_rnd_id = left.ric.first;
    auto const inner_rnd_id = left.ric.second;

    Eigen::MatrixXcd right_sum(
        Eigen::MatrixXcd::Zero(left.data.rows(), left.data.cols()));

    for (auto const &right : right_vec) {
      // We want to make the inner and outer indices match. The inner indices need to
      // match because the product would not make sense otherwise. The outer indices must
      // match since we want to be able to take the trace over the result. The second
      // condition is where this differs from the other multiplication operator.
      bool const is_allowed =
          inner_rnd_id == right.ric.first && outer_rnd_id == right.ric.second;
      if (!is_allowed) {
          continue;
      }

      // We also need to be careful to not combine factors which have common used random
      // vector indices.
      if (has_intersection(left.used_rnd_ids, right.used_rnd_ids)) {
        continue;
      }

      // The right sides that we encounter at this point have the same left and right
      // random vector indices. They may differ in the set of used random vector indices.
      // But since we do not plan to contract the result with more DilutedFactor
      // instances, we do not care to preserve the information about the used random
      // vector indices. Therefore we can sum all these elements up to have less
      // multiplications to do.
      right_sum += right.data;
    }

    auto const &product = left.data * right_sum;
    result += product.trace();
  }

  return result;
}



compcomp_t inner_product(std::vector<DilutedTrace> const &left_vec,
                         std::vector<DilutedTrace> const &right_vec) {
  //! @TODO Pull out this magic number.
  auto constexpr rnd_vec_count = 5;

  compcomp_t result(0.0, 0.0, 0.0, 0.0);

  for (auto const &left : left_vec) {
    cmplx right_sum(0.0, 0.0);

    for (auto const &right : right_vec) {
      // We also need to be careful to not combine factors which have common used random
      // vector indices.
      if (has_intersection(left.used_rnd_ids, right.used_rnd_ids)) {
        continue;
      }

      // The right sides that we encounter at this point have the same left and right
      // random vector indices. They may differ in the set of used random vector indices.
      // But since we do not plan to contract the result with more DilutedFactor
      // instances, we do not care to preserve the information about the used random
      // vector indices. Therefore we can sum all these elements up to have less
      // multiplications to do.
      right_sum += right.data;
    }

    result.rere += left.data.real() * right_sum.real();
    result.reim += left.data.real() * right_sum.imag();
    result.imre += left.data.imag() * right_sum.real();
    result.imim += left.data.imag() * right_sum.imag();
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

/*! corrC */
template <>
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
    size_t const dilD) {
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
      throw std::runtime_error("something wrong with random vectors in build_corrC");
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
