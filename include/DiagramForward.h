#pragma once

#include "typedefs.h"

template <typename Numeric>
class DiagramNumeric;

template <typename Numeric>
struct DiagramTraits {
  using Diagram = DiagramNumeric<Numeric>;
};

class DiagramParts;
