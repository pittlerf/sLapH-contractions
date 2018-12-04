#pragma once

#include "Correlators.hpp"
#include "DilutedFactorY.hpp"
#include "DilutedTraceFactory.hpp"
#include "h5-wrapper.hpp"

#include <omp.h>

#include <mutex>

struct DiagramParts {
  DiagramParts(RandomVector const &random_vector,
               Perambulator const &perambulator,
               OperatorFactory const &meson_operator,
               DilutionScheme const &dilution_scheme,
               ssize_t const dilT,
               ssize_t const dilE,
               ssize_t const nev,
               ssize_t const Lt,
               DilutedFactorIndicesCollection const &dil_fac_lookup,
               DiagramIndicesCollection const &corr_lookup,
               TraceIndicesCollection const &trace_indices_map)
      : q0(random_vector,
           perambulator,
           meson_operator,
           dilT,
           dilE,
           nev,
           dil_fac_lookup.at("Q0")),
        q1(random_vector,
           perambulator,
           meson_operator,
           dilT,
           dilE,
           nev,
           dil_fac_lookup.at("Q1")),
        q2(random_vector,
           perambulator,
           meson_operator,
           dilT,
           dilE,
           nev,
           dil_fac_lookup.at("Q2")),
        trQ1(q1, trace_indices_map.at("trQ1"), dilution_scheme),
        trQ1Q1(q1, q1, trace_indices_map.at("trQ1Q1"), dilution_scheme),
        trQ0Q2(q0, q2, trace_indices_map.at("trQ0Q2"), dilution_scheme),
        trQ1Q1Q1(q1, q1, q1, trace_indices_map.at("trQ1Q1Q1"), dilution_scheme),
        trQ1Q0Q2(q1, q0, q2, trace_indices_map.at("trQ1Q0Q2"), dilution_scheme),
        trQ1Q1Q1Q1(q1, q1, q1, q1, trace_indices_map.at("trQ1Q1Q1Q1"), dilution_scheme) {
    trace_factories["trQ1"] = &trQ1;
    trace_factories["trQ1Q1"] = &trQ1Q1;
    trace_factories["trQ0Q2"] = &trQ0Q2;
    trace_factories["trQ1Q1Q1"] = &trQ1Q1Q1;
    trace_factories["trQ1Q0Q2"] = &trQ1Q0Q2;
    trace_factories["trQ1Q1Q1Q1"] = &trQ1Q1Q1Q1;
  }

  void clear() {
    q0.clear();
    q1.clear();
    q2.clear();

    for (auto const &elem : trace_factories) {
      elem.second->clear();
    }
  }

  DilutedFactorFactory<DilutedFactorType::Q0> q0;
  DilutedFactorFactory<DilutedFactorType::Q1> q1;
  DilutedFactorFactory<DilutedFactorType::Q2> q2;

  //< Temporal memory for tr(Q1)
  DilutedTrace1Factory<DilutedFactorType::Q1> trQ1;

  //< Temporal memory for tr(rVdaggerV*Q1*rVdaggerV*Q1)
  DilutedTrace2Factory<DilutedFactorType::Q1, DilutedFactorType::Q1> trQ1Q1;

  //< Temporal memory for tr(Q2*rVdaggerVr)
  DilutedTrace2Factory<DilutedFactorType::Q0, DilutedFactorType::Q2> trQ0Q2;

  DilutedTrace3Factory<DilutedFactorType::Q1,
                       DilutedFactorType::Q1,
                       DilutedFactorType::Q1>
      trQ1Q1Q1;

  DilutedTrace3Factory<DilutedFactorType::Q1,
                       DilutedFactorType::Q0,
                       DilutedFactorType::Q2>
      trQ1Q0Q2;

  DilutedTrace4Factory<DilutedFactorType::Q1,
                       DilutedFactorType::Q1,
                       DilutedFactorType::Q1,
                       DilutedFactorType::Q1>
      trQ1Q1Q1Q1;

  std::map<std::string, AbstractDilutedTraceFactory *> trace_factories;
};

class Diagram {
 public:
  Diagram(std::vector<DiagramIndex> const &corr_lookup,
          std::vector<CorrelatorRequest> const &corr_requests,
          std::string const &output_path,
          std::string const &output_filename,
          int const Lt)
      : corr_lookup_(corr_lookup),
        corr_requests_(corr_requests),
        output_path_(output_path),
        output_filename_(output_filename),
        Lt_(Lt),
        correlator_(Lt,
                    std::vector<ComplexProduct>(corr_lookup.size(), ComplexProduct{})),
        c_(omp_get_max_threads(),
           std::vector<ComplexProduct>(corr_lookup.size(), ComplexProduct{})),
        mutexes_(Lt) {}

  virtual ~Diagram() {}

  virtual char const *name() const = 0;

  std::vector<DiagramIndex> const &corr_lookup() const { return corr_lookup_; }

  std::vector<CorrelatorRequest> const &correlator_requests() const {
    return corr_requests_;
  }

  void assemble(int const t, BlockIterator const &slice_pair, DiagramParts &q) {
    int const tid = omp_get_thread_num();

    for (int i = 0; i != ssize(corr_lookup()); ++i) {
      c_[tid][i] = ComplexProduct{};
    }

    assemble_impl(c_.at(tid), slice_pair, q);

    {
      std::lock_guard<std::mutex> lock(mutexes_[t]);

      for (int i = 0; i != ssize(corr_lookup()); ++i) {
        correlator_[t][i] += c_[tid][i];
      }
    }
  }

  void write() {
    assert(output_path_ != "");
    assert(output_filename_ != "");

    WriteHDF5Correlator filehandle(
        output_path_, name(), output_filename_, make_comp_type<ComplexProduct>());

    std::vector<ComplexProduct> one_corr(Lt_);

    for (int i = 0; i != ssize(corr_lookup()); ++i) {
      for (int t = 0; t < Lt_; ++t) {
        one_corr[t] = correlator_[t][i] / static_cast<double>(Lt_);
      }
      // Write data to file.
      filehandle.write(one_corr, corr_lookup()[i]);
    }
  }

 private:
  virtual void assemble_impl(std::vector<ComplexProduct> &c,
                             BlockIterator const &slice_pair,
                             DiagramParts &q) = 0;

  std::vector<DiagramIndex> const &corr_lookup_;
  std::vector<CorrelatorRequest> const &corr_requests_;

  std::string const &output_path_;
  std::string const &output_filename_;

  int const Lt_;

  /** OpenMP-shared correlators, indices are (1) time and (2) correlator id. */
  std::vector<std::vector<ComplexProduct>> correlator_;

  /** OpenMP-shared correlators, indices are (1) thread id and (2) correlator id. */
  std::vector<std::vector<ComplexProduct>> c_;

  std::vector<std::mutex> mutexes_;
};

class GeneralDiagram : public Diagram {
 public:
  GeneralDiagram(std::vector<DiagramIndex> const &corr_lookup,
                 std::vector<CorrelatorRequest> const &corr_requests,
                 std::string const &output_path,
                 std::string const &output_filename,
                 int const Lt,
                 char const *name)
      : Diagram(corr_lookup, corr_requests, output_path, output_filename, Lt),
        name_(name) {}

  char const *name() const override { return name_; }

 private:
  void assemble_impl(std::vector<ComplexProduct> &c,
                     BlockIterator const &slice_pair,
                     DiagramParts &q) override;

  char const *name_;
};

/*****************************************************************************/
/*                                   C4                                     */
/*****************************************************************************/

class C4cB : public Diagram {
 public:
  C4cB(std::vector<DiagramIndex> const &corr_lookup,
       std::vector<CorrelatorRequest> const &corr_requests,
       std::string const &output_path,
       std::string const &output_filename,
       int const Lt);

  char const *name() const override { return "C4cB"; }

 private:
  void assemble_impl(std::vector<ComplexProduct> &c,
                     BlockIterator const &slice_pair,
                     DiagramParts &q) override;

  std::vector<std::array<std::array<ssize_t, 2>, 2>> quantum_num_ids_;
};

class C4cC : public Diagram {
 public:
  C4cC(std::vector<DiagramIndex> const &corr_lookup,
       std::vector<CorrelatorRequest> const &corr_requests,
       std::string const &output_path,
       std::string const &output_filename,
       int const Lt);

  char const *name() const override { return "C4cC"; }

 private:
  void assemble_impl(std::vector<ComplexProduct> &c,
                     BlockIterator const &slice_pair,
                     DiagramParts &q) override;

  std::vector<std::array<std::array<ssize_t, 2>, 2>> quantum_num_ids_;
};
