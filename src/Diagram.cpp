#include "Diagram.hpp"

#include "local_timer.hpp"

#include <omp.h>
#include <boost/range/adaptor/indexed.hpp>

ComplexProduct resolve_request(std::vector<TraceRequest> const &trace_requests,
                               BlockIterator const &slice_pair,
                               DiagramParts &q) {
  assert(ssize(trace_requests) == 1);
  auto const &trace_request0 = trace_requests.at(0);
  auto const &locations0 = trace_request0.locations;

  auto const &x0 = q.trace_factories[trace_request0.tr_name]
                       ->get(slice_pair, locations0)
                       .at(trace_request0.tr_id);

  return make_complex_product(
      std::accumulate(std::begin(x0), std::end(x0), Complex{0.0, 0.0}) /
      static_cast<double>(x0.size()), false);
}

ComplexProduct resolve_request2(std::vector<TraceRequest> const &trace_requests,
                                BlockIterator const &slice_pair,
                                DiagramParts &q) {
  assert(ssize(trace_requests) == 2);
  auto const &trace_request0 = trace_requests.at(0);
  auto const &locations0 = trace_request0.locations;
  auto const &x0 = q.trace_factories[trace_request0.tr_name]
                       ->get(slice_pair, locations0)
                       .at(trace_request0.tr_id);
  DilutedTraces t0{x0, false};

  auto const &trace_request1 = trace_requests.at(1);
  auto const &locations1 = trace_request1.locations;
  auto const &x1 = q.trace_factories[trace_request1.tr_name]
                       ->get(slice_pair, locations1)
                       .at(trace_request1.tr_id);
  DilutedTraces t1{x1, false};

  return inner_product(t0, t1);
}

/*****************************************************************************/
/*                                    C2c                                    */
/*****************************************************************************/

void C2c::assemble_impl(std::vector<ComplexProduct> &c,
                        BlockIterator const &slice_pair,
                        DiagramParts &q) {
  assert(correlator_requests().size() == corr_lookup().size());
  for (auto const &request : correlator_requests() | boost::adaptors::indexed(0)) {
    c.at(request.index()) +=
        resolve_request(request.value().trace_requests, slice_pair, q);
  }
}

/*****************************************************************************/
/*                                    C20                                    */
/*****************************************************************************/

void C20::assemble_impl(std::vector<ComplexProduct> &c,
                        BlockIterator const &slice_pair,
                        DiagramParts &q) {
  assert(correlator_requests().size() == corr_lookup().size());
  for (auto const &request : correlator_requests() | boost::adaptors::indexed(0)) {
    c.at(request.index()) +=
        resolve_request(request.value().trace_requests, slice_pair, q);
  }
}

/*****************************************************************************/
/*                                    C20V                                   */
/*****************************************************************************/

void C20V::assemble_impl(std::vector<ComplexProduct> &c,
                         BlockIterator const &slice_pair,
                         DiagramParts &q) {
  assert(correlator_requests().size() == corr_lookup().size());
  for (auto const &request : correlator_requests() | boost::adaptors::indexed(0)) {
    c.at(request.index()) +=
        resolve_request2(request.value().trace_requests, slice_pair, q);
  }
}

/*****************************************************************************/
/*                                    C3c                                    */
/*****************************************************************************/

C3c::C3c(std::vector<DiagramIndex> const &corr_lookup,
         std::vector<CorrelatorRequest> const &corr_requests,
         std::string const &output_path,
         std::string const &output_filename,
         int const Lt)
    : Diagram(
          corr_lookup, corr_requests, output_path, output_filename, Lt) {
  quantum_num_ids_.reserve(corr_lookup.size());

  for (const auto &c_look : corr_lookup) {
    quantum_num_ids_.push_back(
        make_tuple(std::array<ssize_t, 2>{c_look.lookup[2], c_look.lookup[0]},
                   std::array<ssize_t, 1>{c_look.lookup[1]}));
  }
}

void C3c::assemble_impl(std::vector<ComplexProduct> &c,
                        BlockIterator const &slice_pair,
                        DiagramParts &q) {
  DilutedFactorsMap<2> L1;

  LT_DIAGRAMS_DECLARE;
  LT_DIAGRAMS_START;
  for (const auto &ids : quantum_num_ids_) {
    multiply<1, 1>(
        L1,
        std::get<0>(ids),
        q.q0[{slice_pair.source()}],
        q.q2[{slice_pair.source_block(), slice_pair.source(), slice_pair.sink_block()}]);
  }
  LT_DIAGRAMS_STOP;
  LT_DIAGRAMS_PRINT("[C3c::assemble_impl] multiply");

  LT_DIAGRAMS_START;
  for (int i = 0; i != ssize(quantum_num_ids_); ++i) {
    auto const &ids = quantum_num_ids_[i];
    c[i] +=
        trace(L1[std::get<0>(ids)],
              q.q1[{slice_pair.sink(), slice_pair.source_block()}].at(std::get<1>(ids)));
  }
  LT_DIAGRAMS_STOP;
  LT_DIAGRAMS_PRINT("[C3c::assemble_impl] trace");
}

/*****************************************************************************/
/*                                    C30                                    */
/*****************************************************************************/

void C30::assemble_impl(std::vector<ComplexProduct> &c,
                        BlockIterator const &slice_pair,
                        DiagramParts &q) {
  assert(correlator_requests().size() == corr_lookup().size());
  for (auto const &request : correlator_requests() | boost::adaptors::indexed(0)) {
    c.at(request.index()) +=
        resolve_request(request.value().trace_requests, slice_pair, q);
  }
}

/*****************************************************************************/
/*                                   C30V                                    */
/*****************************************************************************/

void C30V::assemble_impl(std::vector<ComplexProduct> &c,
                         BlockIterator const &slice_pair,
                         DiagramParts &q) {
  LT_DIAGRAMS_DECLARE;
  LT_DIAGRAMS_START;
  assert(correlator_requests().size() == corr_lookup().size());

  for (auto const &request : correlator_requests() | boost::adaptors::indexed(0)) {
    c.at(request.index()) +=
        resolve_request2(request.value().trace_requests, slice_pair, q);
  }
  LT_DIAGRAMS_STOP;
  LT_DIAGRAMS_PRINT("[C30::assemble_impl] inner_product");
}

/*****************************************************************************/
/*                                   C4cD                                    */
/*****************************************************************************/

void C4cD::assemble_impl(std::vector<ComplexProduct> &c,
                         BlockIterator const &slice_pair,
                         DiagramParts &q) {
  LT_DIAGRAMS_DECLARE;
  LT_DIAGRAMS_START;
  assert(correlator_requests().size() == corr_lookup().size());
  for (auto const &request : correlator_requests() | boost::adaptors::indexed(0)) {
    c.at(request.index()) +=
        resolve_request2(request.value().trace_requests, slice_pair, q);
  }
  LT_DIAGRAMS_STOP;
  LT_DIAGRAMS_PRINT("[C4cD::assemble_impl] inner_product");
}

/*****************************************************************************/
/*                                   C40D                                    */
/*****************************************************************************/

void C40D::assemble_impl(std::vector<ComplexProduct> &c,
                         BlockIterator const &slice_pair,
                         DiagramParts &q) {
  LT_DIAGRAMS_DECLARE;
  LT_DIAGRAMS_START;
  assert(correlator_requests().size() == corr_lookup().size());
  for (auto const &request : correlator_requests() | boost::adaptors::indexed(0)) {
    c.at(request.index()) +=
        resolve_request2(request.value().trace_requests, slice_pair, q);
  }
  LT_DIAGRAMS_STOP;
  LT_DIAGRAMS_PRINT("[C40D::assemble_impl] inner_product");
}

/*****************************************************************************/
/*                                   C4cV                                    */
/*****************************************************************************/

void C4cV::assemble_impl(std::vector<ComplexProduct> &c,
                         BlockIterator const &slice_pair,
                         DiagramParts &q) {
  LT_DIAGRAMS_DECLARE;
  LT_DIAGRAMS_START;
  assert(correlator_requests().size() == corr_lookup().size());
  for (auto const &request : correlator_requests() | boost::adaptors::indexed(0)) {
    c.at(request.index()) +=
        resolve_request2(request.value().trace_requests, slice_pair, q);
  }
  LT_DIAGRAMS_STOP;
  LT_DIAGRAMS_PRINT("[C4cV::assemble_impl] inner_product");
}

/*****************************************************************************/
/*                                   C40V                                    */
/*****************************************************************************/

void C40V::assemble_impl(std::vector<ComplexProduct> &c,
                         BlockIterator const &slice_pair,
                         DiagramParts &q) {
  LT_DIAGRAMS_DECLARE;
  LT_DIAGRAMS_START;
  assert(correlator_requests().size() == corr_lookup().size());
  for (auto const &request : correlator_requests() | boost::adaptors::indexed(0)) {
    c.at(request.index()) +=
        resolve_request2(request.value().trace_requests, slice_pair, q);
  }
  LT_DIAGRAMS_STOP;
  LT_DIAGRAMS_PRINT("[C40V::assemble_impl] inner_product");
}

/*****************************************************************************/
/*                                   C4cB                                    */
/*****************************************************************************/

C4cB::C4cB(std::vector<DiagramIndex> const &corr_lookup,
           std::vector<CorrelatorRequest> const &corr_requests,
           std::string const &output_path,
           std::string const &output_filename,
           int const Lt)
    : Diagram(
          corr_lookup, corr_requests, output_path, output_filename, Lt) {
  quantum_num_ids_.reserve(corr_lookup.size());

  for (const auto &c_look : corr_lookup) {
    quantum_num_ids_.push_back(
        {std::array<ssize_t, 2>{c_look.lookup[3], c_look.lookup[0]},
         std::array<ssize_t, 2>{c_look.lookup[1], c_look.lookup[2]}});
  }
}

void C4cB::assemble_impl(std::vector<ComplexProduct> &c,
                         BlockIterator const &slice_pair,
                         DiagramParts &q) {
  LT_DIAGRAMS_DECLARE;
  LT_DIAGRAMS_START;
  DilutedFactorsMap<2> L1;
  DilutedFactorsMap<2> L2;
  for (const auto &ids : quantum_num_ids_) {
    multiply<1, 1>(
        L1,
        ids[0],
        q.q0[{slice_pair.source()}],
        q.q2[{slice_pair.source_block(), slice_pair.source(), slice_pair.sink_block()}]);

    multiply<1, 1>(
        L2,
        ids[1],
        q.q0[{slice_pair.sink()}],
        q.q2[{slice_pair.sink_block(), slice_pair.sink(), slice_pair.source_block()}]);
  }
  LT_DIAGRAMS_STOP;
  LT_DIAGRAMS_PRINT("[C4cB::assemble_impl] multiply");

  LT_DIAGRAMS_START;
  for (int i = 0; i != ssize(quantum_num_ids_); ++i) {
    auto const &ids = quantum_num_ids_[i];
    c[i] += trace(L1[ids[0]], L2[ids[1]]);
  }
  LT_DIAGRAMS_STOP;
  LT_DIAGRAMS_PRINT("[C4cB::assemble_impl] trace");
}

/*****************************************************************************/
/*                                   C40B                                    */
/*****************************************************************************/

C40B::C40B(std::vector<DiagramIndex> const &corr_lookup,
           std::vector<CorrelatorRequest> const &corr_requests,
           std::string const &output_path,
           std::string const &output_filename,
           int const Lt)
    : Diagram(
          corr_lookup, corr_requests, output_path, output_filename, Lt) {
  quantum_num_ids_.reserve(corr_lookup.size());

  for (const auto &c_look : corr_lookup) {
    quantum_num_ids_.push_back(
        {std::array<ssize_t, 2>{c_look.lookup[3], c_look.lookup[0]},
         std::array<ssize_t, 2>{c_look.lookup[1], c_look.lookup[2]}});
  }
}

void C40B::assemble_impl(std::vector<ComplexProduct> &c,
                         BlockIterator const &slice_pair,
                         DiagramParts &q) {
  LT_DIAGRAMS_DECLARE;
  LT_DIAGRAMS_START;
  DilutedFactorsMap<2> L1;
  DilutedFactorsMap<2> L2;
  for (const auto &ids : quantum_num_ids_) {
    multiply<1, 1>(L1,
                   ids[0],
                   q.q1[{slice_pair.source(), slice_pair.source_block()}],
                   q.q1[{slice_pair.source(), slice_pair.sink_block()}]);

    multiply<1, 1>(L2,
                   ids[1],
                   q.q1[{slice_pair.sink(), slice_pair.sink_block()}],
                   q.q1[{slice_pair.sink(), slice_pair.source_block()}]);
  }
  LT_DIAGRAMS_STOP;
  LT_DIAGRAMS_PRINT("[C40B::assemble_impl] multiply");

  LT_DIAGRAMS_START;
  for (int i = 0; i != ssize(quantum_num_ids_); ++i) {
    auto const &ids = quantum_num_ids_[i];
    c[i] += trace(L1[ids[0]], L2[ids[1]]);
  }
  LT_DIAGRAMS_STOP;
  LT_DIAGRAMS_PRINT("[C40B::assemble_impl] trace");
}

/*****************************************************************************/
/*                                   C4cC                                    */
/*****************************************************************************/

C4cC::C4cC(std::vector<DiagramIndex> const &corr_lookup,
           std::vector<CorrelatorRequest> const &corr_requests,
           std::string const &output_path,
           std::string const &output_filename,
           int const Lt)
    : Diagram(
          corr_lookup, corr_requests, output_path, output_filename, Lt) {
  quantum_num_ids_.reserve(corr_lookup.size());

  for (const auto &c_look : corr_lookup) {
    quantum_num_ids_.push_back(
        {std::array<ssize_t, 2>{c_look.lookup[3], c_look.lookup[0]},
         std::array<ssize_t, 2>{c_look.lookup[1], c_look.lookup[2]}});
  }
}

void C4cC::assemble_impl(std::vector<ComplexProduct> &c,
                         BlockIterator const &slice_pair,
                         DiagramParts &q) {
  LT_DIAGRAMS_DECLARE;
  LT_DIAGRAMS_START;
  DilutedFactorsMap<2> L1;
  DilutedFactorsMap<2> L2;
  for (const auto &ids : quantum_num_ids_) {
    multiply<1, 1>(
        L1,
        ids[0],
        q.q0[{slice_pair.sink()}],
        q.q2[{slice_pair.sink_block(), slice_pair.source(), slice_pair.sink_block()}]);
    multiply<1, 1>(
        L2,
        ids[1],
        q.q0[{slice_pair.sink()}],
        q.q2[{slice_pair.sink_block(), slice_pair.source(), slice_pair.sink_block()}]);
  }
  LT_DIAGRAMS_STOP;
  LT_DIAGRAMS_PRINT("[C4cC::assemble_impl] multiply");

  LT_DIAGRAMS_START;
  for (int i = 0; i != ssize(quantum_num_ids_); ++i) {
    auto const &ids = quantum_num_ids_[i];
    c[i] += trace(L1[ids[0]], L2[ids[1]]);
  }
  LT_DIAGRAMS_STOP;
  LT_DIAGRAMS_PRINT("[C4cC::assemble_impl] trace");
}

/*****************************************************************************/
/*                                   C40C                                    */
/*****************************************************************************/

C40C::C40C(std::vector<DiagramIndex> const &corr_lookup,
           std::vector<CorrelatorRequest> const &corr_requests,
           std::string const &output_path,
           std::string const &output_filename,
           int const Lt)
    : Diagram(
          corr_lookup, corr_requests, output_path, output_filename, Lt) {
  quantum_num_ids_.reserve(corr_lookup.size());

  for (const auto &c_look : corr_lookup) {
    quantum_num_ids_.push_back(
        {std::array<ssize_t, 2>{c_look.lookup[3], c_look.lookup[0]},
         std::array<ssize_t, 2>{c_look.lookup[1], c_look.lookup[2]}});
  }
}

void C40C::assemble_impl(std::vector<ComplexProduct> &c,
                         BlockIterator const &slice_pair,
                         DiagramParts &q) {
  LT_DIAGRAMS_DECLARE;
  LT_DIAGRAMS_START;
  DilutedFactorsMap<2> L1;
  DilutedFactorsMap<2> L2;
  for (const auto &ids : quantum_num_ids_) {
    multiply<1, 1>(L1,
                   ids[0],
                   q.q1[{slice_pair.sink(), slice_pair.source_block()}],
                   q.q1[{slice_pair.source(), slice_pair.sink_block()}]);

    multiply<1, 1>(L2,
                   ids[1],
                   q.q1[{slice_pair.sink(), slice_pair.source_block()}],
                   q.q1[{slice_pair.source(), slice_pair.sink_block()}]);
  }
  LT_DIAGRAMS_STOP;
  LT_DIAGRAMS_PRINT("[C40C::assemble_impl] multiply");

  LT_DIAGRAMS_START;
  for (int i = 0; i != ssize(quantum_num_ids_); ++i) {
    auto const &ids = quantum_num_ids_[i];
    c[i] += trace(L1[ids[0]], L2[ids[1]]);
  }
  LT_DIAGRAMS_STOP;
  LT_DIAGRAMS_PRINT("[C40C::assemble_impl] trace");
}
