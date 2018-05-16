#include "Correlators.h"

//#define DILUTION_ITERATOR_PRINT

#include "Diagram.h"
#include "DilutedFactor.h"
#include "QuarkLineBlock2.h"
#include "StopWatch.h"
#include "dilution-iterator.h"
#include "typedefs.h"

#include "local_timer.h"

#include <iomanip>

#include <omp.h>

int get_time_delta(BlockIterator const &slice_pair, int const Lt) {
  return abs((slice_pair.sink() - slice_pair.source() - Lt) % Lt);
}

void build_trQ0Q2(DiagramParts &q,
                 DiagramIndex const &c_look,
                 int const t1,
                 int const t2,
                 int const b1,
                 int const b2) {
  q.trQ0Q2[c_look.id][t1][t2] = factor_to_trace(
      q.q0[{t2}].at({c_look.lookup[1]}), q.q2v[{b2, t1, b2}].at({c_look.lookup[0]}));
}

void build_trQ1Q1(DiagramParts &q,
                 DiagramIndex const &c_look,
                 int const t1,
                 int const t2,
                 int const b1,
                 int const b2) {
  q.trQ1Q1[c_look.id][t1][t2] = factor_to_trace(q.q1[{t1, b2}].at({c_look.lookup[0]}),
                                                q.q1[{t2, b1}].at({c_look.lookup[1]}));
}

void build_trQ1(DiagramParts &q,
                DiagramIndex const &c_look,
                int const t,
                int const b) {
  q.trQ1[c_look.id][t] = factor_to_trace(q.q1[{t, b}].at({c_look.lookup[0]}));
}



/******************************************************************************/
/*!
 *  @param quarklines       Instance of Quarklines. Contains prebuilt
 *                          combinations of operators and perambulators
 *  @param meson_operator   Instance of OperatorsForMesons. Contains
 *                          operators (@f$ V^\dagger V $f$) with momenta
 *                          and with/without dilution.
 *  @param perambulators    Instance of Perambulator class. Contains
 *                          Perambulator data
 *  @param operator_lookup
 *  @param corr_lookup
 *  @param quark_lookup
 *
 *  If a diagram is not specified in the infile, corr_lookup contains an empty
 *  vector for this diagram and the build function immediately returns
 */
void contract(const ssize_t Lt,
              const ssize_t dilT,
              const ssize_t dilE,
              const ssize_t nev,
              OperatorFactory const &meson_operator,
              RandomVector const &randomvectors,
              Perambulator const &perambulators,
              OperatorLookup const &operator_lookup,
              DiagramIndicesCollection const &corr_lookup,
              DilutedFactorIndicesCollection const &quark_lookup,
              std::string const output_path,
              std::string const output_filename) {
  DilutedFactorLookup const dil_fac_lookup(
      {quark_lookup.Q0, quark_lookup.Q1, quark_lookup.Q2V, quark_lookup.Q2L});

  // XXX If we had C++14, we could do `make_unique`.
  std::vector<std::unique_ptr<Diagram>> diagrams;

  if (!corr_lookup.C2c.empty())
    diagrams.emplace_back(new C2c(corr_lookup.C2c, output_path, output_filename, Lt));
  if (!corr_lookup.C20.empty())
    diagrams.emplace_back(new C20(corr_lookup.C20, output_path, output_filename, Lt));
  if (!corr_lookup.C20V.empty())
    diagrams.emplace_back(new C20V(corr_lookup.C20V, output_path, output_filename, Lt));
  
  if (!corr_lookup.C3c.empty())
    diagrams.emplace_back(new C3c(corr_lookup.C3c, output_path, output_filename, Lt));
  if (!corr_lookup.C30.empty())
    diagrams.emplace_back(new C30(corr_lookup.C30, output_path, output_filename, Lt));
  if (!corr_lookup.C30V.empty())
    diagrams.emplace_back(new C30V(corr_lookup.C30V, output_path, output_filename, Lt));
  
  if (!corr_lookup.C4cB.empty())
    diagrams.emplace_back(new C4cB(corr_lookup.C4cB, output_path, output_filename, Lt));
  if (!corr_lookup.C4cC.empty())
    diagrams.emplace_back(new C4cC(corr_lookup.C4cC, output_path, output_filename, Lt));
  if (!corr_lookup.C40B.empty())
    diagrams.emplace_back(new C40B(corr_lookup.C40B, output_path, output_filename, Lt));
  if (!corr_lookup.C40C.empty())
    diagrams.emplace_back(new C40C(corr_lookup.C40C, output_path, output_filename, Lt));
  
  if (!corr_lookup.C4cD.empty())
    diagrams.emplace_back(new C4cD(corr_lookup.C4cD, output_path, output_filename, Lt));
  if (!corr_lookup.C4cV.empty())
    diagrams.emplace_back(new C4cV(corr_lookup.C4cV, output_path, output_filename, Lt));
  if (!corr_lookup.C40D.empty())
    diagrams.emplace_back(new C40D(corr_lookup.C40D, output_path, output_filename, Lt));
  if (!corr_lookup.C40V.empty())
    diagrams.emplace_back(new C40V(corr_lookup.C40V, output_path, output_filename, Lt));

  DilutionScheme const dilution_scheme(Lt, dilT, DilutionType::block);

  StopWatch swatch("All contractions");

#pragma omp parallel
  {
    swatch.start();

    LT_CORRELATOR_DECLARE;

    DiagramParts q(randomvectors,
                   perambulators,
                   meson_operator,
                   dilT,
                   dilE,
                   nev,
                   Lt,
                   dil_fac_lookup,
                   corr_lookup);

#pragma omp for schedule(dynamic)
    for (int b = 0; b < dilution_scheme.size(); ++b) {
#pragma omp critical(cout)
      {
        std::cout << "Thread " << std::setw(3) << omp_get_thread_num() << " of "
                  << std::setw(3) << omp_get_num_threads() << " starts with block pair "
                  << std::setw(5) << b << " of " << std::setw(5) << dilution_scheme.size()
                  << "." << std::endl;
      }

      auto const block_pair = dilution_scheme[b];

      LT_CORRELATOR_START;
      // Build trQ0Q2.
      for (auto const slice_pair : block_pair) {
        for (const auto &c_look : corr_lookup.trQ0Q2) {
          build_trQ0Q2(q,
                      c_look,
                      slice_pair.source(),
                      slice_pair.sink(),
                      slice_pair.source_block(),
                      slice_pair.sink_block());

          build_trQ0Q2(q,
                      c_look,
                      slice_pair.source(),
                      slice_pair.source(),
                      slice_pair.source_block(),
                      slice_pair.source_block());
        }
      }
      LT_CORRELATOR_STOP;
      LT_CORRELATOR_PRINT("[contract] build_trQ0Q2");

      LT_CORRELATOR_START;
      // Build trQ1Q1.
      for (auto const slice_pair : block_pair) {
        for (const auto &c_look : corr_lookup.trQ1Q1) {
          build_trQ1Q1(q,
                      c_look,
                      slice_pair.source(),
                      slice_pair.sink(),
                      slice_pair.source_block(),
                      slice_pair.sink_block());

          build_trQ1Q1(q,
                      c_look,
                      slice_pair.source(),
                      slice_pair.source(),
                      slice_pair.source_block(),
                      slice_pair.source_block());
        }
      }
      LT_CORRELATOR_STOP;
      LT_CORRELATOR_PRINT("[contract] build_trQ1Q1");

      LT_CORRELATOR_START;
      // Build tr(Q1).
      for (auto const slice_pair : block_pair.one_sink_slice()) {
        for (const auto &c_look : corr_lookup.trQ1) {
          build_trQ1(q,
                      c_look,
                      slice_pair.source(),
                      slice_pair.source_block());
        }
      }
      LT_CORRELATOR_STOP;
      LT_CORRELATOR_PRINT("[contract] build_trQ1");

      // Build the diagrams.
      for (auto &diagram : diagrams) {
        if (diagram->corr_lookup().empty()) {
          continue;
        }
        LT_CORRELATOR_START;
        for (auto const slice_pair : block_pair) {
          int const t = get_time_delta(slice_pair, Lt);

          diagram->assemble(t, slice_pair, q);
        }  // End of slice pair loop.
        LT_CORRELATOR_STOP;
        LT_CORRELATOR_PRINT( std::string("[contract] ") + diagram->name() );
      }    // End of diagram loop.

      q.clear();
    }  // End of block pair loop.

    swatch.stop();
  }  // End of parallel section.

  swatch.print();

  for (auto &diagram : diagrams) {
    diagram->write();
  }
}
