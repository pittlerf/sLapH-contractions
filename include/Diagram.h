#pragma once

#include "QuarkLineBlock2.h"
#include "typedefs.h"

struct QuarkLineBlockCollection {
  QuarkLineBlock2<QuarkLineType::Q0> &q0;
  QuarkLineBlock2<QuarkLineType::Q1> &q1;
  QuarkLineBlock2<QuarkLineType::Q2> &q2l;
  QuarkLineBlock2<QuarkLineType::Q2> &q2v;
};

class Diagram {
 public:
  void contract(std::vector<cmplx> &c,
                BlockIterator const &slice_pair,
                QuarkLineBlockCollection &q) {
    contract_impl(c, slice_pair, q);
  }

  virtual char const *name() const = 0;

  virtual ~Diagram() {}

 private:
  virtual void contract_impl(std::vector<cmplx> &c,
                             BlockIterator const &slice_pair,
                             QuarkLineBlockCollection &q) = 0;
};

class C4cB : public Diagram {
 public:
  C4cB(std::vector<CorrInfo> const &corr_lookup);

  char const *name() const override { return "C4+B"; }

 private:
  void contract_impl(std::vector<cmplx> &c,
                     BlockIterator const &slice_pair,
                     QuarkLineBlockCollection &q) override;

  std::vector<std::array<std::array<size_t, 2>, 2>> quantum_num_ids_;
};
