#include "Diagram.hpp"

#include "local_timer.hpp"

#include <omp.h>
#include <boost/range/adaptor/indexed.hpp>

ComplexProduct resolve_request1(std::vector<TraceRequest> const &trace_requests,
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
          static_cast<double>(x0.size()),
      false);
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

ComplexProduct resolve_request(std::vector<TraceRequest> const &trace_requests,
                               BlockIterator const &slice_pair,
                               DiagramParts &q) {
  if (ssize(trace_requests) == 1) {
    return resolve_request1(trace_requests, slice_pair, q);
  } else if (ssize(trace_requests) == 2) {
    return resolve_request2(trace_requests, slice_pair, q);
  } else {
    throw std::runtime_error("This many traces are not implemented yet.");
  }
}

void Diagram::assemble_impl(std::vector<ComplexProduct> &c,
                            BlockIterator const &slice_pair,
                            DiagramParts &q) {
  assert(correlator_requests().size() == corr_lookup().size());
  for (auto const &request : correlator_requests() | boost::adaptors::indexed(0)) {
    c.at(request.index()) +=
        resolve_request(request.value().trace_requests, slice_pair, q);
  }
}
