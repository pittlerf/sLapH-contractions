#include "Diagram.h"

void C4cB::contract_slice_pair(std::vector<cmplx> &c,
                               BlockIterator const &slice_pair,
                               QuarkLineBlock2<QuarkLineType::Q0> &q0,
                               QuarkLineBlock2<QuarkLineType::Q1> &q1,
                               QuarkLineBlock2<QuarkLineType::Q2> &q2) {
  OperatorToFactorMap<2, 1> L1;
  OperatorToFactorMap<2, 1> L2;
  for (const auto &ids : quantum_num_ids) {
    multiply<1, 1, 0, 0>(
        L1,
        ids[0],
        q0[{slice_pair.source()}],
        q2[{slice_pair.source_block(), slice_pair.source(), slice_pair.sink_block()}]);

    multiply<1, 1, 0, 0>(
        L2,
        ids[1],
        q0[{slice_pair.sink()}],
        q2[{slice_pair.sink_block(), slice_pair.sink(), slice_pair.source_block()}]);
  }

  for (int i = 0; i != quantum_num_ids.size(); ++i) {
    auto const &ids = quantum_num_ids[i];
    c[i] += trace(L1[ids[0]], L2[ids[1]]);
  }
}
