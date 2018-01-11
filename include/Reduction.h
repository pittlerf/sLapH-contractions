#pragma once

#include "DiagramForward.h"
#include "h5-wrapper.h"

template <typename Numeric>
class Reduction {
 public:
  using MyDiagram = typename DiagramTraits<Numeric>::Diagram;

  Reduction(MyDiagram &diagram,
            std::string const &output_path,
            std::string const &output_filename,
            int const Lt)
      : diagram_(diagram),
        filehandle_(
            output_path, diagram_.name(), output_filename, comp_type_factory<Numeric>()),
        correlator_(diagram_.corr_lookup().size(), std::vector<Numeric>(Lt, Numeric{})),
        Lt_(Lt) {}

  std::vector<std::vector<Numeric>> &correlator() { return correlator_; }

 private:
  MyDiagram &diagram_;
  WriteHDF5Correlator filehandle_;
  std::vector<std::vector<Numeric>> correlator_;
  int const Lt_;
};
