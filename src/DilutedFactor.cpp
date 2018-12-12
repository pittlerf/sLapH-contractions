#include "DilutedFactor.hpp"

std::vector<DilutedFactor> operator*(std::vector<DilutedFactor> const &left_vec,
                                     std::vector<DilutedFactor> const &right_vec) {
  assert(left_vec.size() > 0);
  assert(right_vec.size() > 0);

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
      if ((left.used_rnd_ids & right.used_rnd_ids) != 0u) {
        continue;
      }

      // We want to keep track of the indices that have been contracted away. These are
      // all the ones from the left factor, all the ones from the right factor and the one
      // that we are contracting over right now.
      UsedRnd used = left.used_rnd_ids | right.used_rnd_ids | (1u << inner_rnd_id);

      result_vec.push_back({left.data * right.data,
                            std::make_pair(left.ric.first, right.ric.second),
                            used});
    }
  }

  if (result_vec.size() == 0) {
    throw std::runtime_error(
        "vector<DilutedFactor> operator*(vector<DilutedFactor>, vector<DilutedFactor>) "
        "has an empty result.");
  }

  return result_vec;
}

ComplexProduct trace(std::vector<DilutedFactor> const &left_vec,
                     std::vector<DilutedFactor> const &right_vec) {
  assert(left_vec.size() > 0);
  assert(right_vec.size() > 0);

  Complex result(0.0, 0.0);

  int num_summands = 0;

  LT_ULTRA_FINE_DECLARE;

  for (auto const &left : left_vec) {
    auto const outer_rnd_id = left.ric.first;
    auto const inner_rnd_id = left.ric.second;

    LT_ULTRA_FINE_START;
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
      if ((left.used_rnd_ids & right.used_rnd_ids) == 0u) {
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
    LT_ULTRA_FINE_STOP;
    LT_ULTRA_FINE_PRINT("[DilutedFactor::trace] right_sum");

    LT_ULTRA_FINE_START;
    auto const &product = left.data * right_sum;
    result += product.trace();
    LT_ULTRA_FINE_STOP;
    LT_ULTRA_FINE_PRINT("[DilutedFactor::trace] product_trace");

  }  // for(left_vec)

  if (num_summands == 0) {
    throw std::runtime_error(
        "cmplx trace(vector<DilutedFactor>, vector<DilutedFactor>) has used zero "
        "summands.");
  }

  return make_complex_product(result, false) / static_cast<double>(num_summands);
}
