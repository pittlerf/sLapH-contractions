#pragma once

#include "Correlators.h"
#include "QuarkLineBlock2.h"
#include "Reduction.h"
#include "typedefs.h"

struct QuarkLineBlockCollection {
  QuarkLineBlockCollection(RandomVector const &random_vector,
                           Perambulator const &perambulator,
                           OperatorFactory const &meson_operator,
                           size_t const dilT,
                           size_t const dilE,
                           size_t const nev,
                           size_t const Lt,
                           DilutedFactorLookup const &dil_fac_lookup,
                           DiagramIndicesCollection const &corr_lookup)
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
            dil_fac_lookup.Q2V) {
    corrC.resize(boost::extents[corr_lookup.corrC.size()][Lt][Lt]);
    corr0.resize(boost::extents[corr_lookup.corr0.size()][Lt][Lt]);
    corr_part_trQ1.resize(boost::extents[corr_lookup.trQ1.size()][Lt]);
  }

  void clear() {
    q0.clear();
    q1.clear();
    q2l.clear();
    q2v.clear();
  }

  DilutedFactorFactory<DilutedFactorType::Q0> q0;
  DilutedFactorFactory<DilutedFactorType::Q1> q1;
  DilutedFactorFactory<DilutedFactorType::Q2> q2l;
  DilutedFactorFactory<DilutedFactorType::Q2> q2v;

  //< Temporal memory for tr(rVdaggerV*Q1*rVdaggerV*Q1)
  DilutedTracesTwoTimes<2> corr0;

  //< Temporal memory for tr(Q2V*rVdaggerVr)
  DilutedTracesTwoTimes<2> corrC;

  //< Temporal memory for tr(Q1)
  DilutedTraceOneTime<1> corr_part_trQ1;
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

#pragma omp critical(Diagram_reduce)
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
class C2c : public DiagramNumeric<Complex> {
 public:
  using DiagramNumeric<Complex>::DiagramNumeric;

  char const *name() const override { return "C2+"; }

 private:
  void contract_impl(std::vector<Complex> &c,
                     BlockIterator const &slice_pair,
                     QuarkLineBlockCollection &q) override;
};

/*! Build neutral 2pt correlation function
 *  @f{align}{
 *    C = \langle D_\mathtt{Q0}^{-1}(t'|t) \Gamma_\mathtt{Op0}
 *                D_\mathtt{Q1}^{-1}(t|t') \Gamma_\mathtt{Op1} \rangle
 *  @f}
 */
class C20 : public DiagramNumeric<Complex> {
 public:
  using DiagramNumeric<Complex>::DiagramNumeric;

  char const *name() const override { return "C20"; }

 private:
  void contract_impl(std::vector<Complex> &c,
                     BlockIterator const &slice_pair,
                     QuarkLineBlockCollection &q) override;
};

/*! Build neutral 2pt correlation function
 *  @f{align}{
 *    C = \langle D_\mathtt{Q0}^{-1}(t|t) \Gamma_\mathtt{Op0} \rangle \cdot
 *        \langle D_\mathtt{Q1}^{-1}(t'|t') \Gamma_\mathtt{Op1} \rangle
 *  @f}
 */
class C20V : public DiagramNumeric<ComplexProduct> {
 public:
  using DiagramNumeric<ComplexProduct>::DiagramNumeric;

  char const *name() const override { return "C20V"; }

 private:
  void contract_impl(std::vector<ComplexProduct> &c,
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
class C3c : public DiagramNumeric<Complex> {
 public:
  C3c(std::vector<CorrInfo> const &corr_lookup,
      std::string const &output_path,
      std::string const &output_filename,
      int const Lt);

  char const *name() const override { return "C3+"; }

 private:
  void contract_impl(std::vector<Complex> &c,
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
class C30 : public DiagramNumeric<Complex> {
 public:
  C30(std::vector<CorrInfo> const &corr_lookup,
      std::string const &output_path,
      std::string const &output_filename,
      int const Lt);

  char const *name() const override { return "C30"; }

 private:
  void contract_impl(std::vector<Complex> &c,
                     BlockIterator const &slice_pair,
                     QuarkLineBlockCollection &q) override;

  std::vector<std::tuple<std::array<size_t, 2>, std::array<size_t, 1>>> quantum_num_ids_;
};

class C30V : public DiagramNumeric<ComplexProduct> {
 public:
  using DiagramNumeric<ComplexProduct>::DiagramNumeric;

  char const *name() const override { return "C30V"; }

 private:
  void contract_impl(std::vector<ComplexProduct> &c,
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
class C4cD : public DiagramNumeric<ComplexProduct> {
 public:
  using DiagramNumeric<ComplexProduct>::DiagramNumeric;

  char const *name() const override { return "C4+D"; }

 private:
  void contract_impl(std::vector<ComplexProduct> &c,
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
class C40D : public DiagramNumeric<ComplexProduct> {
 public:
  using DiagramNumeric<ComplexProduct>::DiagramNumeric;

  char const *name() const override { return "C40D"; }

 private:
  void contract_impl(std::vector<ComplexProduct> &c,
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
class C4cV : public DiagramNumeric<ComplexProduct> {
 public:
  using DiagramNumeric<ComplexProduct>::DiagramNumeric;

  char const *name() const override { return "C4+V"; }

 private:
  void contract_impl(std::vector<ComplexProduct> &c,
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
class C40V : public DiagramNumeric<ComplexProduct> {
 public:
  using DiagramNumeric<ComplexProduct>::DiagramNumeric;

  char const *name() const override { return "C40V"; }

 private:
  void contract_impl(std::vector<ComplexProduct> &c,
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
class C4cB : public DiagramNumeric<Complex> {
 public:
  C4cB(std::vector<CorrInfo> const &corr_lookup,
       std::string const &output_path,
       std::string const &output_filename,
       int const Lt);

  char const *name() const override { return "C4+B"; }

 private:
  void contract_impl(std::vector<Complex> &c,
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
class C40B : public DiagramNumeric<Complex> {
 public:
  C40B(std::vector<CorrInfo> const &corr_lookup,
       std::string const &output_path,
       std::string const &output_filename,
       int const Lt);

  char const *name() const override { return "C40B"; }

 private:
  void contract_impl(std::vector<Complex> &c,
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
class C4cC : public DiagramNumeric<Complex> {
 public:
  C4cC(std::vector<CorrInfo> const &corr_lookup,
       std::string const &output_path,
       std::string const &output_filename,
       int const Lt);

  char const *name() const override { return "C4+C"; }

 private:
  void contract_impl(std::vector<Complex> &c,
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
class C40C : public DiagramNumeric<Complex> {
 public:
  C40C(std::vector<CorrInfo> const &corr_lookup,
       std::string const &output_path,
       std::string const &output_filename,
       int const Lt);

  char const *name() const override { return "C40C"; }

 private:
  void contract_impl(std::vector<Complex> &c,
                     BlockIterator const &slice_pair,
                     QuarkLineBlockCollection &q) override;

  std::vector<std::array<std::array<size_t, 2>, 2>> quantum_num_ids_;
};
