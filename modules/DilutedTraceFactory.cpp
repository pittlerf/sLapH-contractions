#include "DilutedTraceFactory.h"

template <>
void DilutedTraceCollection<DilutedFactorType::Q0, DilutedFactorType::Q2, 2>::build(
                 DilutionIterator const &block_pair) {

  for (auto const slice_pair : block_pair){
    auto t1 = slice_pair.source();
    auto t2 = slice_pair.sink();
    auto b1 = dilution_scheme.time_to_block(t1);
    auto b2 = dilution_scheme.time_to_block(t2);

    for (const auto &c_look : diagram_index_collection) {
      tr[{t1,t2}][c_look.id] = factor_to_trace(
        df1[{t2}].at({c_look.lookup[1]}), df2[{b2, t1, b2}].at({c_look.lookup[0]}));
    }

    t1 = slice_pair.source();
    t2 = slice_pair.source();
    b1 = dilution_scheme.time_to_block(t1);
    b2 = dilution_scheme.time_to_block(t2);

    for (const auto &c_look : diagram_index_collection) {
      tr[{t1,t2}][c_look.id] = factor_to_trace(
        df1[{t2}].at({c_look.lookup[1]}), df2[{b2, t1, b2}].at({c_look.lookup[0]}));
    }

  }
}

template <>
void DilutedTraceCollection<DilutedFactorType::Q1, DilutedFactorType::Q1, 2>::build(
                 DilutionIterator const &block_pair) {

  for (auto const slice_pair : block_pair){
    auto t1 = slice_pair.source();
    auto t2 = slice_pair.sink();
    auto b1 = dilution_scheme.time_to_block(t1);
    auto b2 = dilution_scheme.time_to_block(t2);

    for (const auto &c_look : diagram_index_collection) {
      tr[{t1,t2}][c_look.id] = factor_to_trace(df1[{t1, b2}].at({c_look.lookup[0]}), 
        df2[{t2, b1}].at({c_look.lookup[1]}));
    }

    t1 = slice_pair.source();
    t2 = slice_pair.source();
    b1 = dilution_scheme.time_to_block(t1);
    b2 = dilution_scheme.time_to_block(t2);

    for (const auto &c_look : diagram_index_collection) {
      tr[{t1,t2}][c_look.id] = factor_to_trace(df1[{t1, b2}].at({c_look.lookup[0]}), 
        df2[{t2, b1}].at({c_look.lookup[1]}));
    }
  }
}

template class DilutedTraceCollection<DilutedFactorType::Q0, DilutedFactorType::Q2, 2>;
template class DilutedTraceCollection<DilutedFactorType::Q1, DilutedFactorType::Q1, 2>;
