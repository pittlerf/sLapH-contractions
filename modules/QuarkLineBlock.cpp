#include "QuarkLineBlock.h"

namespace {  // some internal namespace

static const std::complex<double> I(0.0, 1.0);
// Look-up table for gamma matrices. For every Gamma structure (currently 0-15)
// the four non-zero values are specified.

/*! @todo Refactor that into physical quantum number class along with momentum */
static void create_gamma(std::vector<LapH::gamma_lookup> &gamma, const int i) {
  try {
    switch (i) {
      case 0:  // gamma_0
        gamma[0].row[0] = 2;
        gamma[0].value[0] = 1;
        gamma[0].row[1] = 3;
        gamma[0].value[1] = 1;
        gamma[0].row[2] = 0;
        gamma[0].value[2] = 1;
        gamma[0].row[3] = 1;
        gamma[0].value[3] = 1;
        break;

      case 1:  // gamma_1
        gamma[1].row[0] = 3;
        gamma[1].value[0] = -I;
        gamma[1].row[1] = 2;
        gamma[1].value[1] = -I;
        gamma[1].row[2] = 1;
        gamma[1].value[2] = I;
        gamma[1].row[3] = 0;
        gamma[1].value[3] = I;
        break;

      case 2:  // gamma_2
        gamma[2].row[0] = 3;
        gamma[2].value[0] = -1;
        gamma[2].row[1] = 2;
        gamma[2].value[1] = 1;
        gamma[2].row[2] = 1;
        gamma[2].value[2] = 1;
        gamma[2].row[3] = 0;
        gamma[2].value[3] = -1;
        break;

      case 3:  // gamma_3
        gamma[3].row[0] = 2;
        gamma[3].value[0] = -I;
        gamma[3].row[1] = 3;
        gamma[3].value[1] = I;
        gamma[3].row[2] = 0;
        gamma[3].value[2] = I;
        gamma[3].row[3] = 1;
        gamma[3].value[3] = -I;
        break;

      case 4:  // unity
        gamma[4].row[0] = 0;
        gamma[4].value[0] = 1;
        gamma[4].row[1] = 1;
        gamma[4].value[1] = 1;
        gamma[4].row[2] = 2;
        gamma[4].value[2] = 1;
        gamma[4].row[3] = 3;
        gamma[4].value[3] = 1;
        break;

      case 5:  // gamma_5
        gamma[5].row[0] = 0;
        gamma[5].value[0] = 1;
        gamma[5].row[1] = 1;
        gamma[5].value[1] = 1;
        gamma[5].row[2] = 2;
        gamma[5].value[2] = -1;
        gamma[5].row[3] = 3;
        gamma[5].value[3] = -1;
        break;

      case 6:  // gamma_0 * gamma_5
        gamma[6].row[0] = 2;
        gamma[6].value[0] = -1;
        gamma[6].row[1] = 3;
        gamma[6].value[1] = -1;
        gamma[6].row[2] = 0;
        gamma[6].value[2] = 1;
        gamma[6].row[3] = 1;
        gamma[6].value[3] = 1;
        break;

      case 7:  // gamma_1 * gamma_5
        gamma[7].row[0] = 3;
        gamma[7].value[0] = I;
        gamma[7].row[1] = 2;
        gamma[7].value[1] = I;
        gamma[7].row[2] = 1;
        gamma[7].value[2] = I;
        gamma[7].row[3] = 0;
        gamma[7].value[3] = I;
        break;

      case 8:  // gamma_2 * gamma_5
        gamma[8].row[0] = 3;
        gamma[8].value[0] = 1;
        gamma[8].row[1] = 2;
        gamma[8].value[1] = -1;
        gamma[8].row[2] = 1;
        gamma[8].value[2] = 1;
        gamma[8].row[3] = 0;
        gamma[8].value[3] = -1;
        break;

      case 9:  // gamma_3 * gamma_5
        gamma[9].row[0] = 2;
        gamma[9].value[0] = I;
        gamma[9].row[1] = 3;
        gamma[9].value[1] = -I;
        gamma[9].row[2] = 0;
        gamma[9].value[2] = I;
        gamma[9].row[3] = 1;
        gamma[9].value[3] = -I;
        break;

      case 10:  // gamma_0 * gamma_1
        gamma[10].row[0] = 1;
        gamma[10].value[0] = I;
        gamma[10].row[1] = 0;
        gamma[10].value[1] = I;
        gamma[10].row[2] = 3;
        gamma[10].value[2] = -I;
        gamma[10].row[3] = 2;
        gamma[10].value[3] = -I;
        break;

      case 11:  // gamma_0 * gamma_2
        gamma[11].row[0] = 1;
        gamma[11].value[0] = 1;
        gamma[11].row[1] = 0;
        gamma[11].value[1] = -1;
        gamma[11].row[2] = 3;
        gamma[11].value[2] = -1;
        gamma[11].row[3] = 2;
        gamma[11].value[3] = 1;
        break;

      case 12:  // gamma_0 * gamma_3
        gamma[12].row[0] = 0;
        gamma[12].value[0] = I;
        gamma[12].row[1] = 1;
        gamma[12].value[1] = -I;
        gamma[12].row[2] = 2;
        gamma[12].value[2] = -I;
        gamma[12].row[3] = 3;
        gamma[12].value[3] = I;
        break;

      case 13:  // gamma_0 * gamma_1 * gamma_5
        gamma[13].row[0] = 1;
        gamma[13].value[0] = I;
        gamma[13].row[1] = 0;
        gamma[13].value[1] = I;
        gamma[13].row[2] = 3;
        gamma[13].value[2] = I;
        gamma[13].row[3] = 2;
        gamma[13].value[3] = I;
        break;

      case 14:  // gamma_0 * gamma_2 * gamma_5
        gamma[14].row[0] = 1;
        gamma[14].value[0] = 1;
        gamma[14].row[1] = 0;
        gamma[14].value[1] = -1;
        gamma[14].row[2] = 3;
        gamma[14].value[2] = 1;
        gamma[14].row[3] = 2;
        gamma[14].value[3] = -1;
        break;

      case 15:  // gamma_0 * gamma_3 * gamma_5
        gamma[15].row[0] = 0;
        gamma[15].value[0] = I;
        gamma[15].row[1] = 1;
        gamma[15].value[1] = -I;
        gamma[15].row[2] = 2;
        gamma[15].value[2] = I;
        gamma[15].row[3] = 3;
        gamma[15].value[3] = -I;
        break;

      default:
        printf("Dirac component %d not found in BasicOperator::create_gamma\n", i);
        exit(0);
    }
    return;
  } catch (std::exception &e) {
    std::cout << e.what() << "in: BasicOperator::create_gamma\n";
    exit(0);
  }
}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
}  // internal namespace ends here

namespace LapH {

template <QuarkLineType qlt>
QuarkLineBlock<qlt>::QuarkLineBlock(
    const size_t dilT,
    const size_t dilE,
    const size_t nev,
    const typename QuarkLineIndices<qlt>::type &quarkline_indices,
    const std::vector<RandomIndexCombinationsQ2> &ric_lookup)
    : dilT(dilT), dilE(dilE), nev(nev) {
  int const eigenspace_dirac_size = dilD * dilE;
  int const from_source_or_sink_block = 2;
  int const to_source_or_sink_block = 2;
  int const quarklines_per_block_combination =
      from_source_or_sink_block * to_source_or_sink_block * dilT;

  Ql.resize(quarklines_per_block_combination);
  for (int qline_id = 0; qline_id < quarklines_per_block_combination; ++qline_id) {
    Ql[qline_id].resize(quarkline_indices.size());
    for (int op_id = 0; op_id < quarkline_indices.size(); ++op_id) {
      int nb_rnd = ric_lookup[quarkline_indices[op_id].id_ric_lookup].rnd_vec_ids.size();
      Ql[qline_id][op_id].resize(nb_rnd);
      for (int rnd_id = 0; rnd_id < nb_rnd; ++rnd_id) {
        Ql[qline_id][op_id][rnd_id] =
            Eigen::MatrixXcd::Zero(eigenspace_dirac_size, eigenspace_dirac_size);
      }
    }
  }

  Ql_id.set_capacity(quarklines_per_block_combination);

  // creating gamma matrices
  gamma.resize(16);
  for (int i = 0; i < gamma.size(); ++i) create_gamma(gamma, i);

  std::cout << "\tQuarklines initialised" << std::endl;
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

template <>
void QuarkLineBlock<QuarkLineType::Q1>::build_Q1_one_t(
    const Perambulator &peram,
    const OperatorsForMesons &meson_operator,
    const int t1,
    const int t2_block,
    const typename QuarkLineIndices<QuarkLineType::Q1>::type &quarkline_indices,
    const std::vector<RandomIndexCombinationsQ2> &ric_lookup) {
  Ql_id.push_front(std::make_pair(t1, t2_block));

  // Effectively this is a right rotation.
  std::rotate(Ql.rbegin(), Ql.rbegin() + 1, Ql.rend());

  for (const auto &op : quarkline_indices) {
    const size_t offset = ric_lookup[op.id_ric_lookup].offset.first;
    size_t rnd_counter = 0;
    for (const auto &rnd_id : ric_lookup[op.id_ric_lookup].rnd_vec_ids) {
      const size_t rid1 = rnd_id.first - offset;
      //! @todo: hard coded! VERY BAD!!!
      const size_t gamma_id = op.gamma[0];
      for (int row = 0; row < 4; row++) {
        for (int col = 0; col < 4; col++) {
          Ql[0][op.id][rnd_counter].block(row * dilE, col * dilE, dilE, dilE) =
              gamma[gamma_id].value[row] *
              meson_operator.return_rvdaggerv(op.id_rvdaggerv, t1, rid1)
                  .block(row * dilE, 0, dilE, nev) *
              peram[rnd_id.second].block((t1 * 4 + gamma[gamma_id].row[row]) * nev,
                                         (t2_block * 4 + col) * dilE,
                                         nev,
                                         dilE);
        }
      }
      rnd_counter++;
    }
  }

}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
/*! @todo Think about better names for time indices */
template <>
void QuarkLineBlock<QuarkLineType::Q1>::build_block_pair(
    Perambulator const &peram,
    OperatorsForMesons const &meson_operator,
    DilutionIterator const &block_pair,
    typename QuarkLineIndices<QuarkLineType::Q1>::type const &quarkline_indices,
    std::vector<RandomIndexCombinationsQ2> const &ric_lookup) {
  for (auto const slice_pair_one_sink : block_pair.one_sink_slice()) {
    build_Q1_one_t(peram,
                   meson_operator,
                   slice_pair_one_sink.source(),
                   block_pair.sink(),
                   quarkline_indices,
                   ric_lookup);
  }
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
template <>
void QuarkLineBlock<QuarkLineType::Q2V>::build_block_pair(
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
          const cmplx value = gamma[gamma_id].value[block_dil];
          const size_t gamma_index = gamma[gamma_id].row[block_dil];
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
void QuarkLineBlock<QuarkLineType::Q2L>::build_block_pair(
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
                        .block((t1 * 4 + row) * nev,
                               (b1 * dilD + col) * dilE,
                               nev,
                               dilE)
                        .adjoint() *
                    meson_operator.return_vdaggerv(qll.id_vdaggerv, t1);
              else
                M.block(col * dilE, row * nev, dilE, nev) =
                    peram[rnd_id.first]
                        .block((t1 * 4 + row) * nev,
                               (b1 * dilD + col) * dilE,
                               nev,
                               dilE)
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
          const cmplx value = gamma[gamma_id].value[block_dil];
          const size_t gamma_index = gamma[gamma_id].row[block_dil];
          for (int row = 0; row < dilD; row++) {
            for (int col = 0; col < dilD; col++) {
              Ql[0][qll.id][rnd_counter].block(row * dilE, col * dilE, dilE, dilE) +=
                  value * M.block(row * dilE, block_dil * nev, dilE, nev) *
                  peram[rnd_id.second].block((t1 * 4 + gamma_index) * nev,
                                             (b2 * dilD + col) * dilE,
                                             nev,
                                             dilE);
            }
          }
        }
        check = rnd_id.first;
        rnd_counter++;
      }
    }
  }
}

template class QuarkLineBlock<QuarkLineType::Q1>;
template class QuarkLineBlock<QuarkLineType::Q2L>;
template class QuarkLineBlock<QuarkLineType::Q2V>;

}  // end of LapH namespace

