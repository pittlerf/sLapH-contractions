#pragma once

#include "Correlators.h"
#include "QuarkLineBlock2.h"
#include "Reduction.h"
#include "typedefs.h"

struct QuarkLineBlockCollection {
  QuarkLineBlockCollection(RandomVector const &random_vector,
                           Perambulator const &perambulator,
                           OperatorsForMesons const &meson_operator,
                           size_t const dilT,
                           size_t const dilE,
                           size_t const nev,
                           DilutedFactorLookup const &dil_fac_lookup,
                           DilutedTraceCollection<2> &corr0,
                           DilutedTraceCollection<2> &corrC,
                           DilutedTraceCollection2<1> &corr_part_trQ1)
      : q0(random_vector,
           perambulator,
           meson_operator,
           dilT,
           dilE,
           nev,
           dil_fac_lookup.Q0),
        q1(random_vector,
           perambulator,
           meson_operator,
           dilT,
           dilE,
           nev,
           dil_fac_lookup.Q1),
        q2l(random_vector,
            perambulator,
            meson_operator,
            dilT,
            dilE,
            nev,
            dil_fac_lookup.Q2L),
        q2v(random_vector,
            perambulator,
            meson_operator,
            dilT,
            dilE,
            nev,
            dil_fac_lookup.Q2V),
        corr0(corr0),
        corrC(corrC),
        corr_part_trQ1(corr_part_trQ1) {}

  void clear() {
    q0.clear();
    q1.clear();
    q2l.clear();
    q2v.clear();
  }

  QuarkLineBlock2<QuarkLineType::Q0> q0;
  QuarkLineBlock2<QuarkLineType::Q1> q1;
  QuarkLineBlock2<QuarkLineType::Q2> q2l;
  QuarkLineBlock2<QuarkLineType::Q2> q2v;

  DilutedTraceCollection<2> &corr0;
  DilutedTraceCollection<2> &corrC;
  DilutedTraceCollection2<1> &corr_part_trQ1;
};

class Diagram {
 public:
  Diagram(std::vector<CorrInfo> const &corr_lookup) : corr_lookup_(corr_lookup) {}

  virtual ~Diagram() {}

  virtual char const *name() const = 0;

  std::vector<CorrInfo> const &corr_lookup() const { return corr_lookup_; }

  virtual void contract(int const t,
                        BlockIterator const &slice_pair,
                        QuarkLineBlockCollection &q) = 0;

  virtual void reduce() = 0;

  virtual void write() = 0;

 private:
  std::vector<CorrInfo> const &corr_lookup_;
};

template <typename Numeric_>
class DiagramNumeric : public Diagram {
 public:
  using Numeric = Numeric_;

  DiagramNumeric(std::vector<CorrInfo> const &_corr_lookup,
                 std::string const &output_path,
                 std::string const &output_filename,
                 int const Lt)
      : Diagram(_corr_lookup),
        output_path_(output_path),
        output_filename_(output_filename),
        Lt_(Lt),
        correlator_(corr_lookup().size(), std::vector<Numeric>(Lt, Numeric{})),
        c_(omp_get_max_threads(),
           std::vector<std::vector<Numeric>>(
               Lt, std::vector<Numeric>(corr_lookup().size(), Numeric{}))) {}

  void contract(int const t,
                BlockIterator const &slice_pair,
                QuarkLineBlockCollection &q) override {
    int const tid = omp_get_thread_num();
    contract_impl(c_.at(tid).at(t), slice_pair, q);
  }

  void reduce() override {
    int const tid = omp_get_thread_num();

#pragma omp critical
    {
      for (int i = 0; i != correlator_.size(); ++i) {
        for (size_t t = 0; t < Lt_; t++) {
          correlator_[i][t] += c_[tid][t][i];
        }
      }
    }
  }

  void write() override {
    assert(output_path_ != "");
    assert(output_filename_ != "");

    WriteHDF5Correlator filehandle(
        output_path_, name(), output_filename_, comp_type_factory<Numeric>());

    for (int i = 0; i != correlator_.size(); ++i) {
      for (auto &corr : correlator_[i]) {
        corr /= Lt_;
      }
      // write data to file
      filehandle.write(correlator_[i], corr_lookup()[i]);
    }
  }

 private:
  virtual void contract_impl(std::vector<Numeric> &c,
                             BlockIterator const &slice_pair,
                             QuarkLineBlockCollection &q) = 0;

  std::string const &output_path_;
  std::string const &output_filename_;

  int const Lt_;

  std::vector<std::vector<Numeric>> correlator_;
  std::vector<std::vector<std::vector<Numeric>>> c_;
};

/*****************************************************************************/
/*                                    C2                                     */
/*****************************************************************************/

/*! Build charged 2pt correlation function
 *  @f{align}{
 *    C = \langle \gamma_5 D_\mathtt{Q0}^{-1}(t|t')^\dagger \gamma_5  \Gamma_\mathtt{Op0}
 *                D_\mathtt{Q1}^{-1}(t|t') \Gamma_\mathtt{Op1} \rangle
 *  @f}
 */
class C2c : public DiagramNumeric<cmplx> {
 public:
  C2c(std::vector<CorrInfo> const &corr_lookup,
      std::string const &output_path,
      std::string const &output_filename,
      int const Lt);

  char const *name() const override { return "C2+"; }

 private:
  void contract_impl(std::vector<cmplx> &c,
                     BlockIterator const &slice_pair,
                     QuarkLineBlockCollection &q) override;
};

/*! Build neutral 2pt correlation function
 *  @f{align}{
 *    C = \langle D_\mathtt{Q0}^{-1}(t'|t) \Gamma_\mathtt{Op0}
 *                D_\mathtt{Q1}^{-1}(t|t') \Gamma_\mathtt{Op1} \rangle
 *  @f}
 */
class C20 : public DiagramNumeric<cmplx> {
 public:
  C20(std::vector<CorrInfo> const &corr_lookup,
      std::string const &output_path,
      std::string const &output_filename,
      int const Lt);

  char const *name() const override { return "C20"; }

 private:
  void contract_impl(std::vector<cmplx> &c,
                     BlockIterator const &slice_pair,
                     QuarkLineBlockCollection &q) override;
};

/*! Build neutral 2pt correlation function
 *  @f{align}{
 *    C = \langle D_\mathtt{Q0}^{-1}(t|t) \Gamma_\mathtt{Op0} \rangle \cdot
 *        \langle D_\mathtt{Q1}^{-1}(t'|t') \Gamma_\mathtt{Op1} \rangle
 *  @f}
 */
class C20V : public DiagramNumeric<compcomp_t> {
 public:
  C20V(std::vector<CorrInfo> const &corr_lookup,
       std::string const &output_path,
       std::string const &output_filename,
       int const Lt);

  char const *name() const override { return "C20V"; }

 private:
  void contract_impl(std::vector<compcomp_t> &c,
                     BlockIterator const &slice_pair,
                     QuarkLineBlockCollection &q) override;
};

/*****************************************************************************/
/*                                    C3                                     */
/*****************************************************************************/

/*! Build neutral 3pt correlation function
 *  @f{align}{
 *    C = \langle \gamma_5 D_\mathtt{Q0}^{-1}(t|t)^\dagger \gamma_5 \Gamma_\mathtt{Op0}
 *                D_\mathtt{Q1}^{-1}(t|t') \Gamma_\mathtt{Op1}
 *                D_\mathtt{Q2}^{-1}(t'|t) \Gamma_\mathtt{Op2} \rangle
 *  @f}
 */
class C3c : public DiagramNumeric<cmplx> {
 public:
  C3c(std::vector<CorrInfo> const &corr_lookup,
      std::string const &output_path,
      std::string const &output_filename,
      int const Lt);

  char const *name() const override { return "C3+"; }

 private:
  void contract_impl(std::vector<cmplx> &c,
                     BlockIterator const &slice_pair,
                     QuarkLineBlockCollection &q) override;

  std::vector<std::tuple<std::array<size_t, 2>, std::array<size_t, 1>>> quantum_num_ids_;
};

/*! Build neutral 3pt correlation function
 *  @f{align}{
 *    C = \langle D_\mathtt{Q0}^{-1}(t|t) \Gamma_\mathtt{Op0}
 *                D_\mathtt{Q1}^{-1}(t|t') \Gamma_\mathtt{Op1}
 *                D_\mathtt{Q2}^{-1}(t'|t) \Gamma_\mathtt{Op2} \rangle
 *  @f}
 */
class C30 : public DiagramNumeric<cmplx> {
 public:
  C30(std::vector<CorrInfo> const &corr_lookup,
      std::string const &output_path,
      std::string const &output_filename,
      int const Lt);

  char const *name() const override { return "C30"; }

 private:
  void contract_impl(std::vector<cmplx> &c,
                     BlockIterator const &slice_pair,
                     QuarkLineBlockCollection &q) override;

  std::vector<std::tuple<std::array<size_t, 2>, std::array<size_t, 1>>> quantum_num_ids_;
};

class C30V : public DiagramNumeric<compcomp_t> {
 public:
  C30V(std::vector<CorrInfo> const &corr_lookup,
       std::string const &output_path,
       std::string const &output_filename,
       int const Lt);

  char const *name() const override { return "C30V"; }

 private:
  void contract_impl(std::vector<compcomp_t> &c,
                     BlockIterator const &slice_pair,
                     QuarkLineBlockCollection &q) override;
};

/*****************************************************************************/
/*                                   C4                                     */
/*****************************************************************************/

/*! Build charged 4pt correlation function: Direct diagram
 *  @f{align}{
 *    C = \langle \gamma_5 D_\mathtt{Q0}^{-1}(t|t')^\dagger \gamma_5 \Gamma_\mathtt{Op0}
 *                D_\mathtt{Q1}^{-1}(t|t') \Gamma_\mathtt{Op1} \rangle \cdot
 *        \langle \gamma_5 D_\mathtt{Q2}^{-1}(t|t')^\dagger \gamma_5 \Gamma_\mathtt{Op2}
 *                D_\mathtt{Q3}^{-1}(t|t') \Gamma_\mathtt{Op3} \rangle
 *  @f}
 */
class C4cD : public DiagramNumeric<compcomp_t> {
 public:
  C4cD(std::vector<CorrInfo> const &corr_lookup,
       std::string const &output_path,
       std::string const &output_filename,
       int const Lt);

  char const *name() const override { return "C4+D"; }

 private:
  void contract_impl(std::vector<compcomp_t> &c,
                     BlockIterator const &slice_pair,
                     QuarkLineBlockCollection &q) override;
};

/*! Build neutral 4pt correlation function: Direct diagram
 *  @f{align}{
 *    C = \langle D_\mathtt{Q0}^{-1}(t'|t) \Gamma_\mathtt{Op0}
 *                D_\mathtt{Q1}^{-1}(t|t') \Gamma_\mathtt{Op1} \rangle \cdot
 *        \langle D_\mathtt{Q2}^{-1}(t'|t) \Gamma_\mathtt{Op2}
 *                D_\mathtt{Q3}^{-1}(t|t') \Gamma_\mathtt{Op3} \rangle
 *  @f}
 */
class C40D : public DiagramNumeric<compcomp_t> {
 public:
  C40D(std::vector<CorrInfo> const &corr_lookup,
       std::string const &output_path,
       std::string const &output_filename,
       int const Lt);

  char const *name() const override { return "C40D"; }

 private:
  void contract_impl(std::vector<compcomp_t> &c,
                     BlockIterator const &slice_pair,
                     QuarkLineBlockCollection &q) override;
};

/*! Build charged 4pt correlation function: Vacuum diagram
 *  @f{align}{
 *    C = \langle \gamma_5 D_\mathtt{Q0}^{-1}(t|t)^\dagger \gamma_5 \Gamma_\mathtt{Op0}
 *                D_\mathtt{Q1}^{-1}(t|t) \Gamma_\mathtt{Op1} \rangle \cdot
 *        \langle \gamma_5 D_\mathtt{Q2}^{-1}(t'|t')^\dagger \gamma_5 \Gamma_\mathtt{Op2}
 *                D_\mathtt{Q3}^{-1}(t'|t') \Gamma_\mathtt{Op3} \rangle
 *  @f}
 */
class C4cV : public DiagramNumeric<compcomp_t> {
 public:
  C4cV(std::vector<CorrInfo> const &corr_lookup,
       std::string const &output_path,
       std::string const &output_filename,
       int const Lt);

  char const *name() const override { return "C4+V"; }

 private:
  void contract_impl(std::vector<compcomp_t> &c,
                     BlockIterator const &slice_pair,
                     QuarkLineBlockCollection &q) override;
};

/*! Build neutral 4pt correlation function: Vacuum diagram
 *  @f{align}{
 *    C = \langle D_\mathtt{Q0}^{-1}(t|t) \Gamma_\mathtt{Op0}
 *                D_\mathtt{Q1}^{-1}(t|t) \Gamma_\mathtt{Op1} \rangle \cdot
 *        \langle D_\mathtt{Q2}^{-1}(t'|t') \Gamma_\mathtt{Op2}
 *                D_\mathtt{Q3}^{-1}(t'|t') \Gamma_\mathtt{Op3} \rangle
 *  @f}
 */
class C40V : public DiagramNumeric<compcomp_t> {
 public:
  C40V(std::vector<CorrInfo> const &corr_lookup,
       std::string const &output_path,
       std::string const &output_filename,
       int const Lt);

  char const *name() const override { return "C40V"; }

 private:
  void contract_impl(std::vector<compcomp_t> &c,
                     BlockIterator const &slice_pair,
                     QuarkLineBlockCollection &q) override;
};

/*! Build charged 4pt correlation function: Box diagram
 *  @f{align}{
 *    C = \langle \gamma_5 D_\mathtt{Q0}^{-1}(t|t)^\dagger \gamma_5 \Gamma_\mathtt{Op0}
 *                D_\mathtt{Q1}^{-1}(t|t') \Gamma_\mathtt{Op1}
 *                \gamma_5 D_\mathtt{Q2}^{-1}(t'|t')^\dagger \gamma_5 \Gamma_\mathtt{Op2}
 *                D_\mathtt{Q3}^{-1}(t'|t) \Gamma_\mathtt{Op3} \rangle
 *  @f}
 */
class C4cB : public DiagramNumeric<cmplx> {
 public:
  C4cB(std::vector<CorrInfo> const &corr_lookup,
       std::string const &output_path,
       std::string const &output_filename,
       int const Lt);

  char const *name() const override { return "C4+B"; }

 private:
  void contract_impl(std::vector<cmplx> &c,
                     BlockIterator const &slice_pair,
                     QuarkLineBlockCollection &q) override;

  std::vector<std::array<std::array<size_t, 2>, 2>> quantum_num_ids_;
};

/*! Build neutral 4pt correlation function: Box diagram
 *  @f{align}{
 *    C = \langle D_\mathtt{Q0}^{-1}(t|t) \Gamma_\mathtt{Op0}
 *                D_\mathtt{Q1}^{-1}(t|t') \Gamma_\mathtt{Op1}
 *                D_\mathtt{Q2}^{-1}(t'|t') \Gamma_\mathtt{Op2}
 *                D_\mathtt{Q3}^{-1}(t'|t) \Gamma_\mathtt{Op3} \rangle
 *  @f}
 */
class C40B : public DiagramNumeric<cmplx> {
 public:
  C40B(std::vector<CorrInfo> const &corr_lookup,
       std::string const &output_path,
       std::string const &output_filename,
       int const Lt);

  char const *name() const override { return "C40B"; }

 private:
  void contract_impl(std::vector<cmplx> &c,
                     BlockIterator const &slice_pair,
                     QuarkLineBlockCollection &q) override;

  std::vector<std::array<std::array<size_t, 2>, 2>> quantum_num_ids_;
};

/*! Build charged 4pt correlation function: Cross diagram
 *  @f{align}{
 *    C = \langle \gamma_5 D_\mathtt{Q0}^{-1}(t|t')^\dagger \gamma_5 \Gamma_\mathtt{Op0}
 *                D_\mathtt{Q1}^{-1}(t|t') \Gamma_\mathtt{Op1}
 *                \gamma_5 D_\mathtt{Q2}^{-1}(t|t')^\dagger \gamma_5 \Gamma_\mathtt{Op2}
 *                D_\mathtt{Q3}^{-1}(t|t') \Gamma_\mathtt{Op3} \rangle
 *  @f}
 */
class C4cC : public DiagramNumeric<cmplx> {
 public:
  C4cC(std::vector<CorrInfo> const &corr_lookup,
       std::string const &output_path,
       std::string const &output_filename,
       int const Lt);

  char const *name() const override { return "C4+C"; }

 private:
  void contract_impl(std::vector<cmplx> &c,
                     BlockIterator const &slice_pair,
                     QuarkLineBlockCollection &q) override;

  std::vector<std::array<std::array<size_t, 2>, 2>> quantum_num_ids_;
};

/*! Build neutral 4pt correlation function: Cross diagram
 *  @f{align}{
 *    C = \langle D_\mathtt{Q0}^{-1}(t'|t) \Gamma_\mathtt{Op0}
 *                D_\mathtt{Q1}^{-1}(t|t') \Gamma_\mathtt{Op1}
 *                D_\mathtt{Q2}^{-1}(t'|t) \Gamma_\mathtt{Op2}
 *                D_\mathtt{Q3}^{-1}(t|t') \Gamma_\mathtt{Op3} \rangle
 *  @f}
 */
class C40C : public DiagramNumeric<cmplx> {
 public:
  C40C(std::vector<CorrInfo> const &corr_lookup,
       std::string const &output_path,
       std::string const &output_filename,
       int const Lt);

  char const *name() const override { return "C40C"; }

 private:
  void contract_impl(std::vector<cmplx> &c,
                     BlockIterator const &slice_pair,
                     QuarkLineBlockCollection &q) override;

  std::vector<std::array<std::array<size_t, 2>, 2>> quantum_num_ids_;
};
