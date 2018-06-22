#include "DilutedTraceFactory.h"

template <>
void DilutedTraceCollection<DilutedFactorType::Q0, DilutedFactorType::Q2, 2>::build(
                 int const t1,
                 int const t2,
                 int const b1,
                 int const b2) {

  for (const auto &c_look : diagram_index_collection) {
  tr[{t1,t2}][c_look.id] = factor_to_trace(
      df1[{t2}].at({c_look.lookup[1]}), df2[{b2, t1, b2}].at({c_look.lookup[0]}));
  }
}

template <>
void DilutedTraceCollection<DilutedFactorType::Q1, DilutedFactorType::Q1, 2>::build(
                 int const t1,
                 int const t2,
                 int const b1,
                 int const b2) {

  for (const auto &c_look : diagram_index_collection) {
    tr[{t1,t2}][c_look.id] = factor_to_trace(df1[{t1, b2}].at({c_look.lookup[0]}), 
      df2[{t2, b1}].at({c_look.lookup[1]}));
  }
}

template class DilutedTraceCollection<DilutedFactorType::Q0, DilutedFactorType::Q2, 2>;
template class DilutedTraceCollection<DilutedFactorType::Q1, DilutedFactorType::Q1, 2>;
