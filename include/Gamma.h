#pragma once

#include "typedefs.h"

#include <vector>

namespace LapH {

struct gamma_lookup {
  std::array<int, 4> row;
  std::array<cmplx, 4> value;
};

void create_gamma(std::vector<gamma_lookup> &gamma, int const ignored);
}
