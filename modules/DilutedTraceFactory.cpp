#include "DilutedTraceFactory.h"

template <>
void DilutedTraceCollection<DilutedFactorType::Q0, DilutedFactorType::Q2, 2>::build(
    DilutedFactorFactory<DilutedFactorType::Q0> &q0,
    DilutedFactorFactory<DilutedFactorType::Q2> &q2v,
                 DiagramIndex const &c_look,
                 int const t1,
                 int const t2,
                 int const b1,
                 int const b2) {
  tr[{t1,t2}][c_look.id] = factor_to_trace(
      q0[{t2}].at({c_look.lookup[1]}), q2v[{b2, t1, b2}].at({c_look.lookup[0]}));
}

// dummy must be passed to keep the function call consistent, but it is a reference to 
// the same object.
template <>
void DilutedTraceCollection<DilutedFactorType::Q1, DilutedFactorType::Q1, 2>::build(
    DilutedFactorFactory<DilutedFactorType::Q1> &q1,
    DilutedFactorFactory<DilutedFactorType::Q1> &dummy,
                 DiagramIndex const &c_look,
                 int const t1,
                 int const t2,
                 int const b1,
                 int const b2) {

  tr[{t1,t2}][c_look.id] = factor_to_trace(q1[{t1, b2}].at({c_look.lookup[0]}), 
      q1[{t2, b1}].at({c_look.lookup[1]}));
}

template class DilutedTraceCollection<DilutedFactorType::Q0, DilutedFactorType::Q2, 2>;
template class DilutedTraceCollection<DilutedFactorType::Q1, DilutedFactorType::Q1, 2>;
