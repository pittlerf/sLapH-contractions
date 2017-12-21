#include "QuarkLineBlock2.h"

#include <utility>

namespace {
std::complex<double> const I(0.0, 1.0);
}

template <QuarkLineType qlt>
QuarkLineBlock2<qlt>::QuarkLineBlock2(
    RandomVector const &random_vector,
    Perambulator const &perambulator,
    OperatorsForMesons const &_meson_operator,
    size_t const dilT,
    size_t const dilE,
    size_t const nev,
    typename QuarkLineIndices<qlt>::type const &quarkline_indices,
    std::vector<RandomIndexCombinationsQ2> const &ric_lookup)
    : peram(perambulator), rnd_vec(random_vector), 
      meson_operator(_meson_operator), dilT(dilT), dilE(dilE), nev(nev) {}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

template <>
void QuarkLineBlock2<QuarkLineType::Q1>::build_Q1_one_t(
    int const t1,
    int const t2_block,
    typename QuarkLineIndices<QuarkLineType::Q1>::type const &quarkline_indices,
    std::vector<RandomIndexCombinationsQ2> const &ric_lookup) {
  int const eigenspace_dirac_size = dilD * dilE;
  auto const time_key = std::make_pair(t1, t2_block);

  // We have already built this, therefore there is no need to do it again.
  if (Ql.count(time_key) > 0) {
    return;
  }

  for (auto const &op : quarkline_indices) {
    auto const offset = op.offset.first;
    for (auto const &rnd_id : op.rnd_vec_ids) {
      auto const rid1 = rnd_id.first - offset;
      auto const rid2 = rnd_id.second - offset;

      auto const gamma_id = op.gamma[0];
      Eigen::MatrixXcd matrix =
          Eigen::MatrixXcd::Zero(eigenspace_dirac_size, eigenspace_dirac_size);
      
      Eigen::MatrixXcd vdv;
      if (op.need_vdaggerv_daggering == false)
        vdv = meson_operator.return_vdaggerv(op.id_vdaggerv, t1);
      else
        vdv = meson_operator.return_vdaggerv(op.id_vdaggerv, t1).adjoint();

      Eigen::MatrixXcd rvdaggerv = Eigen::MatrixXcd::Zero(4*dilE, nev);
  
      for(size_t block = 0; block < 4; block++){
      for(size_t vec_i = 0; vec_i < nev; ++vec_i) {
        size_t blk =  block + vec_i * 4 + 4 * nev * t1;
        rvdaggerv.block(vec_i%dilE + dilE*block, 0, 1, nev) += 
             vdv.row(vec_i) * std::conj(rnd_vec(rnd_id.first, blk));
      }}

      for (int row = 0; row < 4; row++) {
        for (int col = 0; col < 4; col++) {
          matrix.block(row * dilE, col * dilE, dilE, dilE) =
              gamma_vec[gamma_id].value[row] *
              rvdaggerv.block(row * dilE, 0, dilE, nev) *
              peram[rnd_id.second].block((t1 * 4 + gamma_vec[gamma_id].row[row]) * nev,
                                         (t2_block * 4 + col) * dilE,
                                         nev,
                                         dilE);
        }
      }
      /*! @bug by using rid1 and rid2 instead of rnd_id.first and rnd_id.second,
       *       different quark flavors are not seperated
       */
      Ql[time_key][{op.id}].push_back({matrix, std::make_pair(rid1, rid2), {}});
    }
  }
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
/*! @todo Think about better names for time indices */
template <>
void QuarkLineBlock2<QuarkLineType::Q1>::build_block_pair(
    DilutionIterator const &block_pair,
    typename QuarkLineIndices<QuarkLineType::Q1>::type const &quarkline_indices,
    std::vector<RandomIndexCombinationsQ2> const &ric_lookup) {
  for (auto const slice_pair_one_sink : block_pair.one_sink_slice()) {
    std::cout << slice_pair_one_sink << std::endl;
    build_Q1_one_t(slice_pair_one_sink.source(),
                   slice_pair_one_sink.sink_block(),
                   quarkline_indices,
                   ric_lookup);

    // A `DilutionIterator` contains a source and a sink block, possibly different ones.
    // For the block diagram we need to have the `Q1` objects that start and end in the
    // same block as well. The above `build_Q1_one_t` call just builds them between
    // different blocks, therefore we also need this call. The `build_Q1_one_t` will not
    // build anything that is already there, therefore this does not do any damage.

    // TODO However, this also means that we are doing more work then needed. We should
    // think about keeping those elements with the same source and sink block around
    // longer, because we are going to need them a bunch of times later on.
    build_Q1_one_t(slice_pair_one_sink.source(),
                   slice_pair_one_sink.source_block(),
                   quarkline_indices,
                   ric_lookup);
  }
}

#if 0

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// rvdaggervr is calculated by multiplying vdaggerv with the same quantum
// numbers with random vectors from right and left.
template <>
void QuarkLineBlock2<QuarkLineType::Q0>::build_block_pair(
    RandomVector const &rnd_vec,
    OperatorsForMesons const &meson_operator,
    DilutionIterator const &block_pair,
    typename QuarkLineIndices<QuarkLineType::Q0>::type const &quarkline_indices,
    std::vector<RandomIndexCombinationsQ2> const &ric_lookup) {
  for (auto const slice_pair : block_pair.one_sink_slice()) {
    auto const t1 = slice_pair.source();

    Ql_id.push_front(std::make_pair(t1, -1));

    // Effectively this is a right rotation.
    std::rotate(Ql.rbegin(), Ql.rbegin() + 1, Ql.rend());

    for (const auto &op : quarkline_indices) {
      Eigen::MatrixXcd vdv;
      if (op.need_vdaggerv_daggering == false)
        vdv = meson_operator.return_vdaggerv(op.id_vdaggerv, t1);
      else
        vdv = meson_operator.return_vdaggerv(op.id_vdaggerv, t1).adjoint();

      size_t rnd_counter = 0;
      int check = -1;
      Eigen::MatrixXcd M;  // Intermediate memory

      for (const auto &rnd_id : ric_lookup[op.id_ric_lookup].rnd_vec_ids) {
        if (check != rnd_id.second) {  // this avoids recomputation
          M = Eigen::MatrixXcd::Zero(nev, 4 * dilE);
          for (size_t block = 0; block < 4; block++) {
            for (size_t vec_i = 0; vec_i < nev; vec_i++) {
              size_t blk = block + (vec_i + nev * t1) * 4;
              M.block(0, vec_i % dilE + dilE * block, nev, 1) +=
                  vdv.col(vec_i) * rnd_vec(rnd_id.second, blk);
            }
          }
        }
        for (size_t block_x = 0; block_x < 4; block_x++) {
          for (size_t block_y = 0; block_y < 4; block_y++) {
            for (size_t vec_y = 0; vec_y < nev; ++vec_y) {
              size_t blk = block_y + (vec_y + nev * t1) * 4;
              Ql[0][op.id][rnd_counter].block(
                  dilE * block_y + vec_y % dilE, dilE * block_x, 1, dilE) +=
                  M.block(vec_y, dilE * block_x, 1, dilE) *
                  std::conj(rnd_vec(rnd_id.first, blk));
            }
          }
        }
        check = rnd_id.second;
        rnd_counter++;
      }
    }
  }
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
template <>
void QuarkLineBlock2<QuarkLineType::Q2V>::build_block_pair(
    Perambulator const &peram,
    OperatorsForMesons const &meson_operator,
    DilutionIterator const &block_pair,
    typename QuarkLineIndices<QuarkLineType::Q2V>::type const &quarkline_indices,
    std::vector<RandomIndexCombinationsQ2> const &ric_lookup) {
  for (auto const slice_pair : block_pair.one_sink_slice()) {
    auto const t1 = slice_pair.source();
    auto const b2 = slice_pair.sink_block();

    Ql_id.push_front(std::make_pair(t1, b2));

    // Effectively this is a right rotation.
    std::rotate(Ql.rbegin(), Ql.rbegin() + 1, Ql.rend());

    for (const auto &qll : quarkline_indices) {
      size_t rnd_counter = 0;
      int check = -1;
      Eigen::MatrixXcd M = Eigen::MatrixXcd::Zero(dilD * dilE, 4 * nev);
      for (const auto &rnd_id : ric_lookup[qll.id_ric_lookup].rnd_vec_ids) {
        if (check != rnd_id.first) {  // this avoids recomputation
          for (int row = 0; row < dilD; row++) {
            for (int col = 0; col < 4; col++) {
              if (!qll.need_vdaggerv_dag)
                M.block(col * dilE, row * nev, dilE, nev) =
                    peram[rnd_id.first]
                        .block((t1 * 4 + row) * nev, (b2 * dilD + col) * dilE, nev, dilE)
                        .adjoint() *
                    meson_operator.return_vdaggerv(qll.id_vdaggerv, t1);
              else
                M.block(col * dilE, row * nev, dilE, nev) =
                    peram[rnd_id.first]
                        .block((t1 * 4 + row) * nev, (b2 * dilD + col) * dilE, nev, dilE)
                        .adjoint() *
                    meson_operator.return_vdaggerv(qll.id_vdaggerv, t1).adjoint();
              // gamma_5 trick
              if (((row + col) == 3) || (abs(row - col) > 1))
                M.block(col * dilE, row * nev, dilE, nev) *= -1.;
            }
          }
        }
        Ql[0][qll.id][rnd_counter].setZero(dilD * dilE, dilD * dilE);

        const size_t gamma_id = qll.gamma[0];

        for (size_t block_dil = 0; block_dil < dilD; block_dil++) {
          const cmplx value = gamma_vec[gamma_id].value[block_dil];
          const size_t gamma_index = gamma_vec[gamma_id].row[block_dil];
          for (int row = 0; row < dilD; row++) {
            for (int col = 0; col < dilD; col++) {
              Ql[0][qll.id][rnd_counter].block(row * dilE, col * dilE, dilE, dilE) +=
                  value * M.block(row * dilE, block_dil * nev, dilE, nev) *
                  peram[rnd_id.second].block(
                      (t1 * 4 + gamma_index) * nev, (b2 * dilD + col) * dilE, nev, dilE);
            }
          }
        }
        check = rnd_id.first;
        rnd_counter++;
      }
    }
  }
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
template <>
void QuarkLineBlock2<QuarkLineType::Q2L>::build_block_pair(
    Perambulator const &peram,
    OperatorsForMesons const &meson_operator,
    DilutionIterator const &block_pair,
    typename QuarkLineIndices<QuarkLineType::Q2L>::type const &quarkline_indices,
    std::vector<RandomIndexCombinationsQ2> const &ric_lookup) {
  for (auto const slice_pair : block_pair.one_sink_slice()) {
    auto const t1 = slice_pair.source();
    auto const b1 = slice_pair.source_block();
    auto const b2 = slice_pair.sink_block();

    Ql_id.push_front(std::make_pair(t1, b2));

    // Effectively this is a right rotation.
    std::rotate(Ql.rbegin(), Ql.rbegin() + 1, Ql.rend());

    for (const auto &qll : quarkline_indices) {
      size_t rnd_counter = 0;
      int check = -1;
      Eigen::MatrixXcd M = Eigen::MatrixXcd::Zero(dilD * dilE, 4 * nev);
      for (const auto &rnd_id : ric_lookup[qll.id_ric_lookup].rnd_vec_ids) {
        if (check != rnd_id.first) {  // this avoids recomputation
          for (int row = 0; row < dilD; row++) {
            for (int col = 0; col < 4; col++) {
              if (!qll.need_vdaggerv_dag)
                M.block(col * dilE, row * nev, dilE, nev) =
                    peram[rnd_id.first]
                        .block((t1 * 4 + row) * nev, (b1 * dilD + col) * dilE, nev, dilE)
                        .adjoint() *
                    meson_operator.return_vdaggerv(qll.id_vdaggerv, t1);
              else
                M.block(col * dilE, row * nev, dilE, nev) =
                    peram[rnd_id.first]
                        .block((t1 * 4 + row) * nev, (b1 * dilD + col) * dilE, nev, dilE)
                        .adjoint() *
                    meson_operator.return_vdaggerv(qll.id_vdaggerv, t1).adjoint();
              // gamma_5 trick
              if (((row + col) == 3) || (abs(row - col) > 1))
                M.block(col * dilE, row * nev, dilE, nev) *= -1.;
            }
          }
        }
        Ql[0][qll.id][rnd_counter].setZero(dilD * dilE, dilD * dilE);
        const size_t gamma_id = qll.gamma[0];
        for (size_t block_dil = 0; block_dil < dilD; block_dil++) {
          const cmplx value = gamma_vec[gamma_id].value[block_dil];
          const size_t gamma_index = gamma_vec[gamma_id].row[block_dil];
          for (int row = 0; row < dilD; row++) {
            for (int col = 0; col < dilD; col++) {
              Ql[0][qll.id][rnd_counter].block(row * dilE, col * dilE, dilE, dilE) +=
                  value * M.block(row * dilE, block_dil * nev, dilE, nev) *
                  peram[rnd_id.second].block(
                      (t1 * 4 + gamma_index) * nev, (b2 * dilD + col) * dilE, nev, dilE);
            }
          }
        }
        check = rnd_id.first;
        rnd_counter++;
      }
    }
  }
}

#endif

//template class QuarkLineBlock2<QuarkLineType::Q0>;
template class QuarkLineBlock2<QuarkLineType::Q1>;
//template class QuarkLineBlock2<QuarkLineType::Q2L>;
//template class QuarkLineBlock2<QuarkLineType::Q2V>;
