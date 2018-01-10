#include "Diagram.h"

/*****************************************************************************/
/*                                    C3c                                    */
/*****************************************************************************/

C3c::C3c(std::vector<CorrInfo> const &corr_lookup) : Diagram(corr_lookup) {
  quantum_num_ids_.reserve(corr_lookup.size());

  for (const auto &c_look : corr_lookup) {
    quantum_num_ids_.push_back(
        make_tuple(std::array<size_t, 2>{c_look.lookup[2], c_look.lookup[0]},
                   std::array<size_t, 1>{c_look.lookup[1]}));
  }
}

void C3c::contract_impl(std::vector<cmplx> &c,
                        BlockIterator const &slice_pair,
                        QuarkLineBlockCollection &q) {
  OperatorToFactorMap<2, 1> L1;
  for (const auto &ids : quantum_num_ids_) {
    multiply<1, 1, 0, 0>(
        L1,
        std::get<0>(ids),
        q.q0[{slice_pair.source()}],
        q.q2l[{slice_pair.source_block(), slice_pair.source(), slice_pair.sink_block()}]);
  }

  for (int i = 0; i != quantum_num_ids_.size(); ++i) {
    auto const &ids = quantum_num_ids_[i];
    c[i] +=
        trace(L1[std::get<0>(ids)],
              q.q1[{slice_pair.sink(), slice_pair.source_block()}].at(std::get<1>(ids)));
  }
}

/*****************************************************************************/
/*                                    C30                                    */
/*****************************************************************************/

C30::C30(std::vector<CorrInfo> const &corr_lookup) : Diagram(corr_lookup) {
  quantum_num_ids_.reserve(corr_lookup.size());

  for (const auto &c_look : corr_lookup) {
    quantum_num_ids_.push_back(
        make_tuple(std::array<size_t, 2>{c_look.lookup[2], c_look.lookup[0]},
                   std::array<size_t, 1>{c_look.lookup[1]}));
  }
}

void C30::contract_impl(std::vector<cmplx> &c,
                        BlockIterator const &slice_pair,
                        QuarkLineBlockCollection &q) {
  OperatorToFactorMap<2, 1> L1;
  for (const auto &ids : quantum_num_ids_) {
    multiply<1, 1, 0, 0>(L1,
                         std::get<0>(ids),
                         q.q1[{slice_pair.source(), slice_pair.source_block()}],
                         q.q1[{slice_pair.source(), slice_pair.sink_block()}]);
  }

  for (int i = 0; i != quantum_num_ids_.size(); ++i) {
    auto const &ids = quantum_num_ids_[i];
    c[i] +=
        trace(L1[std::get<0>(ids)],
              q.q1[{slice_pair.sink(), slice_pair.source_block()}].at(std::get<1>(ids)));
  }
}

/*****************************************************************************/
/*                                   C4cB                                    */
/*****************************************************************************/

C4cB::C4cB(std::vector<CorrInfo> const &corr_lookup) : Diagram(corr_lookup) {
  quantum_num_ids_.reserve(corr_lookup.size());

  for (const auto &c_look : corr_lookup) {
    quantum_num_ids_.push_back(
        {std::array<size_t, 2>{c_look.lookup[3], c_look.lookup[0]},
         std::array<size_t, 2>{c_look.lookup[1], c_look.lookup[2]}});
  }
}

void C4cB::contract_impl(std::vector<cmplx> &c,
                         BlockIterator const &slice_pair,
                         QuarkLineBlockCollection &q) {
  OperatorToFactorMap<2, 1> L1;
  OperatorToFactorMap<2, 1> L2;
  for (const auto &ids : quantum_num_ids_) {
    multiply<1, 1, 0, 0>(
        L1,
        ids[0],
        q.q0[{slice_pair.source()}],
        q.q2l[{slice_pair.source_block(), slice_pair.source(), slice_pair.sink_block()}]);

    multiply<1, 1, 0, 0>(
        L2,
        ids[1],
        q.q0[{slice_pair.sink()}],
        q.q2l[{slice_pair.sink_block(), slice_pair.sink(), slice_pair.source_block()}]);
  }

  for (int i = 0; i != quantum_num_ids_.size(); ++i) {
    auto const &ids = quantum_num_ids_[i];
    c[i] += trace(L1[ids[0]], L2[ids[1]]);
  }
}

/*****************************************************************************/
/*                                   C40B                                    */
/*****************************************************************************/

C40B::C40B(std::vector<CorrInfo> const &corr_lookup) : Diagram(corr_lookup) {
  quantum_num_ids_.reserve(corr_lookup.size());

  for (const auto &c_look : corr_lookup) {
    quantum_num_ids_.push_back(
        {std::array<size_t, 2>{c_look.lookup[3], c_look.lookup[0]},
         std::array<size_t, 2>{c_look.lookup[1], c_look.lookup[2]}});
  }
}

void C40B::contract_impl(std::vector<cmplx> &c,
                         BlockIterator const &slice_pair,
                         QuarkLineBlockCollection &q) {
  OperatorToFactorMap<2, 1> L1;
  OperatorToFactorMap<2, 1> L2;
  for (const auto &ids : quantum_num_ids_) {
    multiply<1, 1, 0, 0>(L1,
                         ids[0],
                         q.q1[{slice_pair.source(), slice_pair.source_block()}],
                         q.q1[{slice_pair.source(), slice_pair.sink_block()}]);

    multiply<1, 1, 0, 0>(L2,
                         ids[1],
                         q.q1[{slice_pair.sink(), slice_pair.sink_block()}],
                         q.q1[{slice_pair.sink(), slice_pair.source_block()}]);
  }

  for (int i = 0; i != quantum_num_ids_.size(); ++i) {
    auto const &ids = quantum_num_ids_[i];
    c[i] += trace(L1[ids[0]], L2[ids[1]]);
  }
}

/*****************************************************************************/
/*                                   C4cC                                    */
/*****************************************************************************/

C4cC::C4cC(std::vector<CorrInfo> const &corr_lookup) : Diagram(corr_lookup) {
  quantum_num_ids_.reserve(corr_lookup.size());

  for (const auto &c_look : corr_lookup) {
    quantum_num_ids_.push_back(
        {std::array<size_t, 2>{c_look.lookup[3], c_look.lookup[0]},
         std::array<size_t, 2>{c_look.lookup[1], c_look.lookup[2]}});
  }
}

void C4cC::contract_impl(std::vector<cmplx> &c,
                         BlockIterator const &slice_pair,
                         QuarkLineBlockCollection &q) {
  OperatorToFactorMap<2, 1> L1;
  OperatorToFactorMap<2, 1> L2;
  for (const auto &ids : quantum_num_ids_) {
    multiply<1, 1, 0, 0>(
        L1,
        ids[0],
        q.q0[{slice_pair.sink()}],
        q.q2v[{slice_pair.sink_block(), slice_pair.source(), slice_pair.sink_block()}]);
    multiply<1, 1, 0, 0>(
        L2,
        ids[1],
        q.q0[{slice_pair.sink()}],
        q.q2v[{slice_pair.sink_block(), slice_pair.source(), slice_pair.sink_block()}]);
  }

  for (int i = 0; i != quantum_num_ids_.size(); ++i) {
    auto const &ids = quantum_num_ids_[i];
    c[i] += trace(L1[ids[0]], L2[ids[1]]);
  }
}

/*****************************************************************************/
/*                                   C40C                                    */
/*****************************************************************************/

C40C::C40C(std::vector<CorrInfo> const &corr_lookup) : Diagram(corr_lookup) {
  quantum_num_ids_.reserve(corr_lookup.size());

  for (const auto &c_look : corr_lookup) {
    quantum_num_ids_.push_back(
        {std::array<size_t, 2>{c_look.lookup[3], c_look.lookup[0]},
         std::array<size_t, 2>{c_look.lookup[1], c_look.lookup[2]}});
  }
}

void C40C::contract_impl(std::vector<cmplx> &c,
                         BlockIterator const &slice_pair,
                         QuarkLineBlockCollection &q) {
  OperatorToFactorMap<2, 1> L1;
  OperatorToFactorMap<2, 1> L2;
  for (const auto &ids : quantum_num_ids_) {
    multiply<1, 1, 0, 0>(L1,
                         ids[0],
                         q.q1[{slice_pair.sink(), slice_pair.source_block()}],
                         q.q1[{slice_pair.source(), slice_pair.sink_block()}]);

    multiply<1, 1, 0, 0>(L2,
                         ids[1],
                         q.q1[{slice_pair.sink(), slice_pair.source_block()}],
                         q.q1[{slice_pair.source(), slice_pair.sink_block()}]);
  }

  for (int i = 0; i != quantum_num_ids_.size(); ++i) {
    auto const &ids = quantum_num_ids_[i];
    c[i] += trace(L1[ids[0]], L2[ids[1]]);
  }
}
