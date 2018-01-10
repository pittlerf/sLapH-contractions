#pragma once

#include "QuarkLineBlock2.h"
#include "typedefs.h"

class Diagram {
  public:
   virtual void contract_slice_pair(std::vector<cmplx> &c,
                                    BlockIterator const &slice_pair,
                                    QuarkLineBlock2<QuarkLineType::Q0> &q0,
                                    QuarkLineBlock2<QuarkLineType::Q1> &q1,
                                    QuarkLineBlock2<QuarkLineType::Q2> &q2) = 0;
};

class C4cB : public Diagram {
 public:
  C4cB(std::vector<std::array<std::array<size_t, 2>, 2>> const &quantum_num_ids)
      : quantum_num_ids(quantum_num_ids) {}

  void contract_slice_pair(std::vector<cmplx> &c,
                           BlockIterator const &slice_pair,
                           QuarkLineBlock2<QuarkLineType::Q0> &q0,
                           QuarkLineBlock2<QuarkLineType::Q1> &q1,
                           QuarkLineBlock2<QuarkLineType::Q2> &q2) override;

 private:
  std::vector<std::array<std::array<size_t, 2>, 2>> const &quantum_num_ids;
};
