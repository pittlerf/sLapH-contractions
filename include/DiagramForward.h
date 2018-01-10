#pragma once

#include "typedefs.h"

class DiagramComp;
class DiagramCompComp;

template <typename Numeric>
struct DiagramTraits {};

template <>
struct DiagramTraits<cmplx> {
  using Diagram = DiagramComp;
};

template <>
struct DiagramTraits<compcomp_t> {
  using Diagram = DiagramCompComp;
};
