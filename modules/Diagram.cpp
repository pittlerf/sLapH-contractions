#include "Diagram.h"

/*****************************************************************************/
/*                                    C2c                                    */
/*****************************************************************************/

void C2c::assemble_impl(std::vector<Complex> &c,
                        BlockIterator const &slice_pair,
                        DiagramParts &q) {
  for (int i = 0; i != ssize(corr_lookup()); ++i) {
    auto const &c_look = corr_lookup()[i];

    auto const &x = q.trQ0Q2[c_look.lookup[0]][slice_pair.source()][slice_pair.sink()];
    c[i] += std::accumulate(std::begin(x), std::end(x), Complex(0.0, 0.0)) /
            static_cast<double>(x.size());
  }
}

/*****************************************************************************/
/*                                    C20                                    */
/*****************************************************************************/

void C20::assemble_impl(std::vector<Complex> &c,
                        BlockIterator const &slice_pair,
                        DiagramParts &q) {
  for (int i = 0; i != ssize(corr_lookup()); ++i) {
    auto const &c_look = corr_lookup()[i];

    auto const &x = q.trQ1Q1[c_look.lookup[0]][slice_pair.source()][slice_pair.sink()];
    c[i] += std::accumulate(std::begin(x), std::end(x), Complex(0.0, 0.0)) /
            static_cast<double>(x.size());
  }
}

/*****************************************************************************/
/*                                    C20V                                   */
/*****************************************************************************/

void C20V::assemble_impl(std::vector<ComplexProduct> &c,
                         BlockIterator const &slice_pair,
                         DiagramParts &q) {
  for (int i = 0; i != ssize(corr_lookup()); ++i) {
    auto const &c_look = corr_lookup()[i];

    c[i] += inner_product(q.trQ1[c_look.lookup[0]][slice_pair.source()],
                          q.trQ1[c_look.lookup[1]][slice_pair.sink()]);
  }
}

/*****************************************************************************/
/*                                    C3c                                    */
/*****************************************************************************/

C3c::C3c(std::vector<DiagramIndex> const &corr_lookup,
         std::string const &output_path,
         std::string const &output_filename,
         int const Lt)
    : DiagramNumeric<Complex>(corr_lookup, output_path, output_filename, Lt) {
  quantum_num_ids_.reserve(corr_lookup.size());

  for (const auto &c_look : corr_lookup) {
    quantum_num_ids_.push_back(
        make_tuple(std::array<ssize_t, 2>{c_look.lookup[2], c_look.lookup[0]},
                   std::array<ssize_t, 1>{c_look.lookup[1]}));
  }
}

void C3c::assemble_impl(std::vector<Complex> &c,
                        BlockIterator const &slice_pair,
                        DiagramParts &q) {
  DilutedFactors<2, 1> L1;
  for (const auto &ids : quantum_num_ids_) {
    multiply<1, 1, 0, 0>(
        L1,
        std::get<0>(ids),
        q.q0[{slice_pair.source()}],
        q.q2l[{slice_pair.source_block(), slice_pair.source(), slice_pair.sink_block()}]);
  }

  for (int i = 0; i != ssize(quantum_num_ids_); ++i) {
    auto const &ids = quantum_num_ids_[i];
    c[i] +=
        trace(L1[std::get<0>(ids)],
              q.q1[{slice_pair.sink(), slice_pair.source_block()}].at(std::get<1>(ids)));
  }
}

/*****************************************************************************/
/*                                    C30                                    */
/*****************************************************************************/

C30::C30(std::vector<DiagramIndex> const &corr_lookup,
         std::string const &output_path,
         std::string const &output_filename,
         int const Lt)
    : DiagramNumeric<Complex>(corr_lookup, output_path, output_filename, Lt) {
  quantum_num_ids_.reserve(corr_lookup.size());

  for (const auto &c_look : corr_lookup) {
    quantum_num_ids_.push_back(
        make_tuple(std::array<ssize_t, 2>{c_look.lookup[2], c_look.lookup[0]},
                   std::array<ssize_t, 1>{c_look.lookup[1]}));
  }
}

void C30::assemble_impl(std::vector<Complex> &c,
                        BlockIterator const &slice_pair,
                        DiagramParts &q) {
  DilutedFactors<2, 1> L1;
  for (const auto &ids : quantum_num_ids_) {
    multiply<1, 1, 0, 0>(L1,
                         std::get<0>(ids),
                         q.q1[{slice_pair.source(), slice_pair.source_block()}],
                         q.q1[{slice_pair.source(), slice_pair.sink_block()}]);
  }

  for (int i = 0; i != ssize(quantum_num_ids_); ++i) {
    auto const &ids = quantum_num_ids_[i];
    c[i] +=
        trace(L1[std::get<0>(ids)],
              q.q1[{slice_pair.sink(), slice_pair.source_block()}].at(std::get<1>(ids)));
  }
}

/*****************************************************************************/
/*                                   C30V                                    */
/*****************************************************************************/

void C30V::assemble_impl(std::vector<ComplexProduct> &c,
                         BlockIterator const &slice_pair,
                         DiagramParts &q) {
  for (int i = 0; i != ssize(corr_lookup()); ++i) {
    auto const &c_look = corr_lookup()[i];

    assert(c_look.lookup[0] < q.trQ1Q1.shape()[0]);
    assert(slice_pair.source() < q.trQ1Q1.shape()[1]);
    assert(slice_pair.source() < q.trQ1Q1.shape()[2]);

    assert(c_look.lookup[1] < q.trQ1.shape()[0]);
    assert(slice_pair.sink() < q.trQ1.shape()[1]);

    c[i] += inner_product(
        q.trQ1Q1[c_look.lookup[0]][slice_pair.source()][slice_pair.source()],
        q.trQ1[c_look.lookup[1]][slice_pair.sink()]);
  }
}

/*****************************************************************************/
/*                                   C4cD                                    */
/*****************************************************************************/

void C4cD::assemble_impl(std::vector<ComplexProduct> &c,
                         BlockIterator const &slice_pair,
                         DiagramParts &q) {
  for (int i = 0; i != ssize(corr_lookup()); ++i) {
    auto const &c_look = corr_lookup()[i];

    c[i] +=
        inner_product(q.trQ0Q2[c_look.lookup[0]][slice_pair.source()][slice_pair.sink()],
                      q.trQ0Q2[c_look.lookup[1]][slice_pair.source()][slice_pair.sink()]);
  }
}

/*****************************************************************************/
/*                                   C40D                                    */
/*****************************************************************************/

void C40D::assemble_impl(std::vector<ComplexProduct> &c,
                         BlockIterator const &slice_pair,
                         DiagramParts &q) {
  for (int i = 0; i != ssize(corr_lookup()); ++i) {
    auto const &c_look = corr_lookup()[i];

    c[i] +=
        inner_product(q.trQ1Q1[c_look.lookup[0]][slice_pair.source()][slice_pair.sink()],
                      q.trQ1Q1[c_look.lookup[1]][slice_pair.source()][slice_pair.sink()]);
  }
}

/*****************************************************************************/
/*                                   C4cV                                    */
/*****************************************************************************/

void C4cV::assemble_impl(std::vector<ComplexProduct> &c,
                         BlockIterator const &slice_pair,
                         DiagramParts &q) {
  for (int i = 0; i != ssize(corr_lookup()); ++i) {
    auto const &c_look = corr_lookup()[i];

    c[i] += inner_product(
        q.trQ0Q2[c_look.lookup[0]][slice_pair.source()][slice_pair.source()],
        q.trQ0Q2[c_look.lookup[1]][slice_pair.sink()][slice_pair.sink()]);
  }
}

/*****************************************************************************/
/*                                   C40V                                    */
/*****************************************************************************/

void C40V::assemble_impl(std::vector<ComplexProduct> &c,
                         BlockIterator const &slice_pair,
                         DiagramParts &q) {
  for (int i = 0; i != ssize(corr_lookup()); ++i) {
    auto const &c_look = corr_lookup()[i];

    c[i] += inner_product(
        q.trQ1Q1[c_look.lookup[0]][slice_pair.source()][slice_pair.source()],
        q.trQ1Q1[c_look.lookup[1]][slice_pair.sink()][slice_pair.sink()]);
  }
}

/*****************************************************************************/
/*                                   C4cB                                    */
/*****************************************************************************/

C4cB::C4cB(std::vector<DiagramIndex> const &corr_lookup,
           std::string const &output_path,
           std::string const &output_filename,
           int const Lt)
    : DiagramNumeric<Complex>(corr_lookup, output_path, output_filename, Lt) {
  quantum_num_ids_.reserve(corr_lookup.size());

  for (const auto &c_look : corr_lookup) {
    quantum_num_ids_.push_back(
        {std::array<ssize_t, 2>{c_look.lookup[3], c_look.lookup[0]},
         std::array<ssize_t, 2>{c_look.lookup[1], c_look.lookup[2]}});
  }
}

void C4cB::assemble_impl(std::vector<Complex> &c,
                         BlockIterator const &slice_pair,
                         DiagramParts &q) {
  DilutedFactors<2, 1> L1;
  DilutedFactors<2, 1> L2;
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

  for (int i = 0; i != ssize(quantum_num_ids_); ++i) {
    auto const &ids = quantum_num_ids_[i];
    c[i] += trace(L1[ids[0]], L2[ids[1]]);
  }
}

/*****************************************************************************/
/*                                   C40B                                    */
/*****************************************************************************/

C40B::C40B(std::vector<DiagramIndex> const &corr_lookup,
           std::string const &output_path,
           std::string const &output_filename,
           int const Lt)
    : DiagramNumeric<Complex>(corr_lookup, output_path, output_filename, Lt) {
  quantum_num_ids_.reserve(corr_lookup.size());

  for (const auto &c_look : corr_lookup) {
    quantum_num_ids_.push_back(
        {std::array<ssize_t, 2>{c_look.lookup[3], c_look.lookup[0]},
         std::array<ssize_t, 2>{c_look.lookup[1], c_look.lookup[2]}});
  }
}

void C40B::assemble_impl(std::vector<Complex> &c,
                         BlockIterator const &slice_pair,
                         DiagramParts &q) {
  DilutedFactors<2, 1> L1;
  DilutedFactors<2, 1> L2;
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

  for (int i = 0; i != ssize(quantum_num_ids_); ++i) {
    auto const &ids = quantum_num_ids_[i];
    c[i] += trace(L1[ids[0]], L2[ids[1]]);
  }
}

/*****************************************************************************/
/*                                   C4cC                                    */
/*****************************************************************************/

C4cC::C4cC(std::vector<DiagramIndex> const &corr_lookup,
           std::string const &output_path,
           std::string const &output_filename,
           int const Lt)
    : DiagramNumeric<Complex>(corr_lookup, output_path, output_filename, Lt) {
  quantum_num_ids_.reserve(corr_lookup.size());

  for (const auto &c_look : corr_lookup) {
    quantum_num_ids_.push_back(
        {std::array<ssize_t, 2>{c_look.lookup[3], c_look.lookup[0]},
         std::array<ssize_t, 2>{c_look.lookup[1], c_look.lookup[2]}});
  }
}

void C4cC::assemble_impl(std::vector<Complex> &c,
                         BlockIterator const &slice_pair,
                         DiagramParts &q) {
  DilutedFactors<2, 1> L1;
  DilutedFactors<2, 1> L2;
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

  for (int i = 0; i != ssize(quantum_num_ids_); ++i) {
    auto const &ids = quantum_num_ids_[i];
    c[i] += trace(L1[ids[0]], L2[ids[1]]);
  }
}

/*****************************************************************************/
/*                                   C40C                                    */
/*****************************************************************************/

C40C::C40C(std::vector<DiagramIndex> const &corr_lookup,
           std::string const &output_path,
           std::string const &output_filename,
           int const Lt)
    : DiagramNumeric<Complex>(corr_lookup, output_path, output_filename, Lt) {
  quantum_num_ids_.reserve(corr_lookup.size());

  for (const auto &c_look : corr_lookup) {
    quantum_num_ids_.push_back(
        {std::array<ssize_t, 2>{c_look.lookup[3], c_look.lookup[0]},
         std::array<ssize_t, 2>{c_look.lookup[1], c_look.lookup[2]}});
  }
}

void C40C::assemble_impl(std::vector<Complex> &c,
                         BlockIterator const &slice_pair,
                         DiagramParts &q) {
  DilutedFactors<2, 1> L1;
  DilutedFactors<2, 1> L2;
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

  for (int i = 0; i != ssize(quantum_num_ids_); ++i) {
    auto const &ids = quantum_num_ids_[i];
    c[i] += trace(L1[ids[0]], L2[ids[1]]);
  }
}
