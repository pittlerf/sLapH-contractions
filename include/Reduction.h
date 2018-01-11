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

  void reduce(std::vector<std::vector<Numeric>> const &c) {
#pragma omp critical
    {
      for (int i = 0; i != correlator_.size(); ++i) {
        for (size_t t = 0; t < Lt_; t++) {
          correlator_[i][t] += c[t][i];
        }
      }
    }
  }

  void write() {
    for (int i = 0; i != correlator_.size(); ++i) {
      for (auto &corr : correlator_[i]) {
        corr /= Lt_;
      }
      // write data to file
      filehandle_.write(correlator_[i], diagram_.corr_lookup()[i]);
    }
  }

  std::vector<std::vector<Numeric>> &correlator() { return correlator_; }



 private:
  MyDiagram &diagram_;
  WriteHDF5Correlator filehandle_;
  std::vector<std::vector<Numeric>> correlator_;
  int const Lt_;
};
