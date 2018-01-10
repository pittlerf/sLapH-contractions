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

   virtual ~Diagram() {}
};

class C4cB : public Diagram {
 public:
  C4cB(std::vector<CorrInfo> const &corr_lookup) {
    quantum_num_ids_.reserve(corr_lookup.size());

    for (const auto &c_look : corr_lookup) {
      quantum_num_ids_.push_back(
          {std::array<size_t, 2>{c_look.lookup[3], c_look.lookup[0]},
           std::array<size_t, 2>{c_look.lookup[1], c_look.lookup[2]}});
    }
  }

  void contract_slice_pair(std::vector<cmplx> &c,
                           BlockIterator const &slice_pair,
                           QuarkLineBlock2<QuarkLineType::Q0> &q0,
                           QuarkLineBlock2<QuarkLineType::Q1> &q1,
                           QuarkLineBlock2<QuarkLineType::Q2> &q2) override;

 private:
  std::vector<std::array<std::array<size_t, 2>, 2>> quantum_num_ids_;
};
