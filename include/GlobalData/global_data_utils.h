#ifndef GLOBALDATA_UTILS_H_
#define GLOBALDATA_UTILS_H_

#include <array>
#include <iostream>
#include <string>

#include "global_data_typedefs.h"
#include "typedefs.h"

namespace global_data_utils {

  // functions for input handling
  quark make_quark (const std::string& quark_string);
  void quark_check(quark quarks);
  Operator_list make_operator_list(const std::string& operator_string);

} // end of namespace global_data_utils

#endif // GLOBALDATA_UTILS_H_
