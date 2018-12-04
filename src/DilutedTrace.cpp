#include "DilutedTrace.hpp"

ComplexProduct inner_product(DilutedTraces const &left_vec,
                             DilutedTraces const &right_vec) {
  assert(left_vec.size() > 0);
  assert(right_vec.size() > 0);

  int num_summands = 0;

  ComplexProduct result = complex_product_zero;

  for (auto const &left : left_vec.traces) {
    Complex right_sum(0.0, 0.0);

    for (auto const &right : right_vec.traces) {
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
      ++num_summands;
    }

    result += make_complex_product(left.data, left_vec.ignore_imag) *
              make_complex_product(right_sum, right_vec.ignore_imag);
  }

  if (num_summands == 0) {
    throw std::runtime_error(
        "compcomp_t inner_product(vector<DilutedTrace>, "
        "vector<DilutedTrace>) has an empty result.");
  }

  return result / static_cast<double>(num_summands);
}

std::vector<DilutedTrace> factor_to_trace(std::vector<DilutedFactor> const &left_vec,
                                          std::vector<DilutedFactor> const &right_vec) {
  assert(left_vec.size() > 0);
  assert(right_vec.size() > 0);

  std::vector<DilutedTrace> result_vec;

  LT_ULTRA_FINE_DECLARE;

  for (auto const &left : left_vec) {
    auto const outer_rnd_id = left.ric.first;
    auto const inner_rnd_id = left.ric.second;

    LT_ULTRA_FINE_START;
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
      result_vec.push_back(
          {typename DilutedTrace::Data{(left.data * right.data).trace()}, used});
    }
    LT_ULTRA_FINE_STOP;
    LT_ULTRA_FINE_PRINT("[DilutedFactor::factor_to_trace] multiply_trace");
  }  // for(left_vec)

  if (result_vec.size() == 0) {
    throw std::runtime_error(
        "vector<DilutedTrace> factor_to_trace(vector<DilutedFactor>, "
        "vector<DilutedFactor>) has an empty result.");
  }

  return result_vec;
}

std::vector<DilutedTrace> factor_to_trace(std::vector<DilutedFactor> const &vec) {
  assert(vec.size() > 0);

  std::vector<DilutedTrace> result_vec;

  for (auto const &elem : vec) {
    // We only want to use diagonal elements.
    if (elem.ric.first != elem.ric.second) {
      continue;
    }

    SmallVectorRndId used;
    std::copy(std::begin(elem.used_rnd_ids),
              std::end(elem.used_rnd_ids),
              std::back_inserter(used));

    auto const outer_rnd_id = elem.ric.first;
    merge_push_back(used, outer_rnd_id);

    DilutedTrace result = {elem.data.trace(), used};

    result_vec.push_back(result);
  }

  if (result_vec.size() == 0) {
    throw std::runtime_error(
        "vector<DilutedTrace> factor_to_trace(vector<DilutedFactor>) has an empty "
        "result.");
  }

  return result_vec;
}
