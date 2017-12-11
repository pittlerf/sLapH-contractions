#include "Correlators.h"

//#define DILUTION_ITERATOR_PRINT

#include "DilutedFactor.h"
#include "QuarkLineBlock.h"
#include "QuarkLineBlock2.h"
#include "StopWatch.h"
#include "dilution-iterator.h"
#include "typedefs.h"

namespace {

/*! Creates compound datatype to write complex numbers from LapH::complex_t
 *  vectors to HDF5 file
 *
 *  @Returns cmplx_w   HDF5 compound datatype for complex numbers
 */
H5::CompType comp_type_factory_tr() {
  H5::CompType cmplx_w(2 * sizeof(double));
  auto type = H5::PredType::NATIVE_DOUBLE;
  cmplx_w.insertMember("re", HOFFSET(LapH::complex_t, re), type);
  cmplx_w.insertMember("im", HOFFSET(LapH::complex_t, im), type);

  return cmplx_w;
}

/*! Creates compound datatype to write complex numbers from LapH::compcomp_t
 *  vectors to HDF5 file
 *
 *  @Returns cmplx_w   HDF5 compound datatype for structs of four doubles
 */
H5::CompType comp_type_factory_trtr() {
  H5::CompType cmplxcmplx_w(4 * sizeof(double));
  auto type = H5::PredType::NATIVE_DOUBLE;
  cmplxcmplx_w.insertMember("rere", HOFFSET(LapH::compcomp_t, rere), type);
  cmplxcmplx_w.insertMember("reim", HOFFSET(LapH::compcomp_t, reim), type);
  cmplxcmplx_w.insertMember("imre", HOFFSET(LapH::compcomp_t, imre), type);
  cmplxcmplx_w.insertMember("imim", HOFFSET(LapH::compcomp_t, imim), type);

  return cmplxcmplx_w;
}

/*! Class to write correlations function to files in hdf5 format
 *
 *  @Warning  Dependency inversion principle is violated: Class depends on the
 *            concrete implementation of hdf5
 */
class WriteHDF5Correlator {

public:
  WriteHDF5Correlator(const std::string output_path, const std::string diagram,
                      const std::string output_filename,
                      const H5::CompType &_comp_type)
      : comp_type(_comp_type) {

    create_folder_for_hdf5_file(output_path.c_str());

    const H5std_string file_name(
        (output_path + "/" + diagram + output_filename).c_str());
    open_or_create_hdf5_file(file_name);
  }

  /*! Writes data to file
   *
   *  @Param corr       The data to write
   *  @Param corr_info  Contains the hdf5_dataset_name
   *
   *  @todo   It is sufficient to pass the hdf5_dataset_name
   *  @todo   corr_datatype would be better as class template typename
   *
   *  @remark The type corr_datatype is always either LapH::complex_t or
   *          LapH::compcomp_t. The function body is identical for both types
   *          as everything is specified by corr_info. Thus the template
   *          overload
   */
  template <typename corr_datatype>
  void write(const std::vector<corr_datatype> &corr,
             const CorrInfo &corr_info) {

    // Turn off the autoprinting when an exception occurs because fuck you
    // thats why
    H5::Exception::dontPrint();
    // That's right bitches!

    // create the dataset to write data ----------------------------------------
    H5::Group group;
    H5::DataSet dset;
    H5std_string dataset_name((corr_info.hdf5_dataset_name).c_str());
    hsize_t dim(corr.size());
    H5::DataSpace dspace(1, &dim);

    // actual write
    try {
      dset = file.createDataSet(dataset_name, comp_type, dspace);
      dset.write(&corr[0], comp_type);

      // closing of dset is delegated to destructor ~DataSet

    } catch (H5::Exception &e) {
      e.printError();
    }
  }

private:
  /*! Checks whether output path exists and if not creates it
   *
   *  @param[in] path Path where hdf5 file shall be written
   */
  void create_folder_for_hdf5_file(const char *path) {
    if (access(path, 0) != 0) {
      std::cout << "\tdirectory " << path
                << " does not exist and will be created";
      boost::filesystem::path dir(path);
      if (boost::filesystem::create_directories(dir))
        std::cout << "\tSuccess" << std::endl;
      else
        std::cout << "\tFailure" << std::endl;
    }
  }

  /*! Opens output file in Truncation mode by default
   *
   *  @param[in]  name String containing path+filename of the desired file
   *  @param[out] file File pointer to hdf5 file
   *
   *  @todo If overwrite=no flag is set, rather than H5F_ACC_TRUNC should be
   *        replaced by H5F_ACC_EXCL
   *
   *  @throws H5::exception if something goes wrong
   */
  void open_or_create_hdf5_file(const H5std_string &name) {

    file = H5::H5File(name, H5F_ACC_TRUNC);
  }

  /*! The hdf5 file pointer */
  H5::H5File file;
  /*! The hdf5 compound datatype.
   *
   *  @see  H5::CompType comp_type_factory_tr()
   *  @see  H5::CompType comp_type_factory_trtr()
   */
  H5::CompType comp_type;

}; // end of class WriteHDF5Correlator

} // end of anonymous namespace

int get_time_delta(BlockIterator const &slice_pair, int const Lt) {
  return abs((slice_pair.sink() - slice_pair.source() - Lt) % Lt);
}

/******************************************************************************/
/******************************************************************************/

LapH::Correlators::Correlators(const size_t Lt, const size_t dilT,
                               const size_t dilE, const size_t nev,
                               const CorrelatorLookup &corr_lookup,
                               OperatorLookup const &operator_lookup,
                               QuarklineLookup const &quark_lookup)
    : Lt(Lt), dilT(dilT), dilE(dilE), nev(nev),
      dil_fac_lookup(DilutedFactorLookup(operator_lookup.rvdaggervr_lookuptable,
                                         quark_lookup.Q1, quark_lookup.Q2V,
                                         quark_lookup.Q2L)),
      ric_lookup(operator_lookup.ricQ2_lookup) {}

/*!
 *  @deprecated
 */
void LapH::Correlators::build_C1(OperatorsForMesons const &meson_operator,
                                 Perambulator const &perambulators,
                                 std::vector<CorrInfo> const &corr_lookup,
                                 std::string const output_path,
                                 std::string const output_filename) {

  if (corr_lookup.size() == 0)
    return;

  // every element of corr_lookup contains the same filename. Wlog choose the
  // first element
  WriteHDF5Correlator filehandle(output_path, "C1", output_filename,
                                 comp_type_factory_tr());

  StopWatch swatch("C1", 1);

  std::vector<vec> correlator(corr_lookup.size(), vec(Lt, cmplx(.0, .0)));

  DilutionScheme const dilution_scheme(Lt, dilT, DilutionType::block);

  //#pragma omp parallel
  {
    swatch.start();

    QuarkLineBlock<QuarkLineType::Q1> quarklines(dilT, dilE, nev,
                                                 dil_fac_lookup.Q1, ric_lookup);

    std::vector<vec> C(corr_lookup.size(), vec(Lt, cmplx(.0, .0)));

    /*! @todo can DilutionIterator also give just a single block? */
    //#pragma omp for schedule(dynamic)
    for (int t = 0; t < Lt; t++) {
      int const b = t / dilT;
      quarklines.build_Q1_one_t(perambulators, meson_operator, t, b,
                                dil_fac_lookup.Q1, ric_lookup);

      for (const auto &c_look : corr_lookup) {
        const auto &ric =
            ric_lookup[dil_fac_lookup.Q1[c_look.lookup[0]].id_ric_lookup]
                .rnd_vec_ids;

        std::vector<cmplx> correlator(Lt, cmplx(.0, .0));
        for (const auto &id : ric) {
          auto const idr0 = &id - &ric[0];
          C[c_look.id][t] += quarklines(t, b, c_look.lookup[0], idr0).trace();
        }
      } // loop over operators ends here
    }   // loop over time ends here

    //#pragma omp critical
    {
      for (const auto &c_look : corr_lookup)
        for (size_t t = 0; t < Lt; t++)
          correlator[c_look.id][t] += C[c_look.id][t];
    }
    swatch.stop();
  } // parallel part ends here

  // normalisation
  for (const auto &c_look : corr_lookup) {
    for (auto &corr : correlator[c_look.id]) {
      // TODO: Hard Coded atm - Be carefull
      corr /= 5 * Lt;
    }
    // write data to file
    filehandle.write(correlator[c_look.id], c_look);
  }

  swatch.print();
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void LapH::Correlators::build_corr0(OperatorsForMesons const &meson_operator,
                                    Perambulator const &perambulators,
                                    std::vector<CorrInfo> const &corr_lookup) {
  if (corr_lookup.size() == 0)
    return;

  StopWatch swatch("corr0");

  corr0.resize(boost::extents[corr_lookup.size()][Lt][Lt]);

  DilutionScheme const dilution_scheme(Lt, dilT, DilutionType::block);

#pragma omp parallel
  {
    swatch.start();

    QuarkLineBlock<QuarkLineType::Q1> quarklines_local(
        dilT, dilE, nev, dil_fac_lookup.Q1, ric_lookup);

#pragma omp for schedule(dynamic)
    for (int b = 0; b < dilution_scheme.size(); ++b) {
      // Notation is that `t1` is the source and `t2` the sink. Both will be
      // done eventually, so this is symmetric.

      auto const block_pair = dilution_scheme[b];

      quarklines_local.build_block_pair(perambulators, meson_operator,
                                        block_pair, dil_fac_lookup.Q1,
                                        ric_lookup);

      for (auto const slice_pair : block_pair) {
        for (const auto &c_look : corr_lookup) {

          const auto &ric0 =
              ric_lookup[dil_fac_lookup.Q1[c_look.lookup[0]].id_ric_lookup]
                  .rnd_vec_ids;
          const auto &ric1 =
              ric_lookup[dil_fac_lookup.Q1[c_look.lookup[1]].id_ric_lookup]
                  .rnd_vec_ids;
          if (ric0.size() != ric1.size()) {
            std::cout << "rnd combinations are not the same in build_corr0"
                      << std::endl;
            exit(1);
          }

          std::vector<size_t> random_index_combination_ids =
              std::vector<size_t>(
                  {dil_fac_lookup.Q1[c_look.lookup[0]].id_ric_lookup,
                   dil_fac_lookup.Q1[c_look.lookup[1]].id_ric_lookup});

          corr0[c_look.id][slice_pair.source()][slice_pair.sink()] = trace(
              quarklines_local(slice_pair.source(), slice_pair.sink_block(),
                               c_look.lookup[0]),
              quarklines_local(slice_pair.sink(), slice_pair.source_block(),
                               c_look.lookup[1]),
              ric_lookup, random_index_combination_ids);
        }
      }
    }
    swatch.stop();
  }

  swatch.print();
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void LapH::Correlators::build_C20(std::vector<CorrInfo> const &corr_lookup,
                                  std::string const output_path,
                                  std::string const output_filename) {
  if (corr_lookup.empty())
    return;

  StopWatch swatch("C20");
  swatch.start();

  // every element of corr_lookup contains the same filename. Wlog choose the
  // first element
  WriteHDF5Correlator filehandle(output_path, "C20", output_filename,
                                 comp_type_factory_tr());

  for (const auto &c_look : corr_lookup) {
    std::vector<cmplx> correlator(Lt, cmplx(.0, .0));

    DilutionScheme const dilution_scheme(Lt, dilT, DilutionType::block);
    for (auto const block_pair : dilution_scheme) {
      for (auto const slice_pair : block_pair) {
        int const t = get_time_delta(slice_pair, Lt);
        /*! @todo hidden because range based but this is a loop over random */
        correlator[t] += std::accumulate(
            corr0[c_look.lookup[0]][slice_pair.source()][slice_pair.sink()]
                .begin(),
            corr0[c_look.lookup[0]][slice_pair.source()][slice_pair.sink()]
                .end(),
            cmplx(0.0, 0.0));
      }
    }
    // normalisation
    for (auto &corr : correlator) {
      corr /= Lt * corr0[c_look.lookup[0]][0][0].size();
    }

    // write data to file
    filehandle.write(correlator, c_look);
  }

  swatch.stop();
  swatch.print();
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void LapH::Correlators::build_C40D(CorrelatorLookup const &corr_lookup,
                                   std::string const output_path,
                                   std::string const output_filename) {
  if (corr_lookup.C40D.empty())
    return;

  StopWatch swatch("C40D", 1);
  swatch.start();

  // every element of corr_lookup contains the same filename. Wlog choose the
  // first element
  WriteHDF5Correlator filehandle(output_path, "C40D", output_filename,
                                 comp_type_factory_trtr());

  DilutionScheme const dilution_scheme(Lt, dilT, DilutionType::block);

  for (const auto &c_look : corr_lookup.C40D) {

    std::vector<LapH::compcomp_t> correlator(Lt,
                                             LapH::compcomp_t(.0, .0, .0, .0));

    const size_t id0 = corr_lookup.corr0[c_look.lookup[0]].lookup[0];
    const size_t id1 = corr_lookup.corr0[c_look.lookup[1]].lookup[0];

    for (auto const block_pair : dilution_scheme) {
      for (auto const slice_pair : block_pair) {
        int const t = get_time_delta(slice_pair, Lt);

        std::vector<size_t> random_index_combination_ids =
            std::vector<size_t>({dil_fac_lookup.Q1[id0].id_ric_lookup,
                                 dil_fac_lookup.Q1[id1].id_ric_lookup});

        /*! @todo Write move assignment for compcomp_t and give trtr return
         * parameter */
        correlator[t] += trtr(
            corr0[c_look.lookup[0]][slice_pair.source()][slice_pair.sink()],
            corr0[c_look.lookup[1]][slice_pair.source()][slice_pair.sink()],
            ric_lookup, random_index_combination_ids);
      }
    }

    // normalisation
    for (auto &corr : correlator) {
      corr.rere /= (5U * 4U * 3U * 2U) * Lt;
      corr.reim /= (5U * 4U * 3U * 2U) * Lt;
      corr.imre /= (5U * 4U * 3U * 2U) * Lt;
      corr.imim /= (5U * 4U * 3U * 2U) * Lt;
    }

    // write data to file
    filehandle.write(correlator, c_look);
  }

  swatch.stop();
  swatch.print();
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void LapH::Correlators::build_C40V(CorrelatorLookup const &corr_lookup,
                                   std::string const output_path,
                                   std::string const output_filename) {

  if (corr_lookup.C40V.empty())
    return;

  StopWatch swatch("C40V");
  swatch.start();

  // every element of corr_lookup contains the same filename. Wlog choose the
  // first element
  WriteHDF5Correlator filehandle(output_path, "C40V", output_filename,
                                 comp_type_factory_trtr());

  DilutionScheme const dilution_scheme(Lt, dilT, DilutionType::block);

  for (const auto &c_look : corr_lookup.C40V) {

    std::vector<LapH::compcomp_t> correlator(Lt,
                                             LapH::compcomp_t(.0, .0, .0, .0));

    const size_t id0 = corr_lookup.corr0[c_look.lookup[0]].lookup[0];
    const size_t id1 = corr_lookup.corr0[c_look.lookup[1]].lookup[0];

    for (auto const block_pair : dilution_scheme) {
      for (auto const slice_pair : block_pair) {
        int const t = get_time_delta(slice_pair, Lt);

        std::vector<size_t> random_index_combination_ids =
            std::vector<size_t>({dil_fac_lookup.Q1[id0].id_ric_lookup,
                                 dil_fac_lookup.Q1[id1].id_ric_lookup});

        /*! @todo Write move assignment for compcomp_t and give trtr return
         * parameter */
        correlator[t] += trtr(
            corr0[c_look.lookup[0]][slice_pair.source()][slice_pair.source()],
            corr0[c_look.lookup[1]][slice_pair.sink()][slice_pair.sink()],
            ric_lookup, random_index_combination_ids);
      }
    }

    // normalisation
    for (auto &corr : correlator) {
      corr.rere /= (5U * 4U * 3U * 2U) * Lt;
      corr.reim /= (5U * 4U * 3U * 2U) * Lt;
      corr.imre /= (5U * 4U * 3U * 2U) * Lt;
      corr.imim /= (5U * 4U * 3U * 2U) * Lt;
    }

    // write data to file
    filehandle.write(correlator, c_look);
  }

  swatch.stop();
  swatch.print();
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void LapH::Correlators::build_corrC(RandomVector const &randomvectors,
                                    Perambulator const &perambulators,
                                    OperatorsForMesons const &meson_operator,
                                    std::vector<CorrInfo> const &corr_lookup) {
  if (corr_lookup.size() == 0)
    return;

  StopWatch swatch("corrC");

  corrC.resize(boost::extents[corr_lookup.size()][Lt][Lt]);

  DilutionScheme const dilution_scheme(Lt, dilT, DilutionType::block);

#pragma omp parallel
  {
    swatch.start();

    QuarkLineBlock<QuarkLineType::Q2V> quarklines_Q2V(
        dilT, dilE, nev, dil_fac_lookup.Q2V, ric_lookup);
    QuarkLineBlock<QuarkLineType::Q0> quarklines_Q0(
        dilT, dilE, nev, dil_fac_lookup.Q0, ric_lookup);

#pragma omp for schedule(dynamic)
    for (int b = 0; b < dilution_scheme.size(); ++b) {
      auto const block_pair = dilution_scheme[b];

      quarklines_Q2V.build_block_pair(perambulators, meson_operator, block_pair,
                                      dil_fac_lookup.Q2V, ric_lookup);
      quarklines_Q0.build_block_pair(randomvectors, meson_operator, block_pair,
                                     dil_fac_lookup.Q0, ric_lookup);

      for (auto const slice_pair : block_pair) {
        auto const t1 = slice_pair.source();
        auto const t2 = slice_pair.sink();

        // building correlator
        for (const auto &c_look : corr_lookup) {
          const auto &ric0 =
              ric_lookup[dil_fac_lookup.Q2V[c_look.lookup[0]].id_ric_lookup]
                  .rnd_vec_ids;
          const auto &ric1 =
              ric_lookup[ // just for checking
                  dil_fac_lookup.Q0[c_look.lookup[1]].id_ric_lookup]
                  .rnd_vec_ids;
          if (ric0.size() != ric1.size()) {
            std::cout << "rnd combinations are not the same in build_corrC"
                      << std::endl;
            exit(0);
          }

          std::vector<size_t> random_index_combination_ids{
              dil_fac_lookup.Q2V[c_look.lookup[0]].id_ric_lookup,
              dil_fac_lookup.Q0[c_look.lookup[1]].id_ric_lookup};

          corrC[c_look.id][t1][t2] =
              trace<QuarkLineType::Q2V, QuarkLineType::Q0>(quarklines_Q2V,
                                                           quarklines_Q0,
                                                           slice_pair.source(),
                                                           slice_pair.sink_block(),
                                                           slice_pair.sink(),
                                                           c_look.lookup,
                                                           ric_lookup,
                                                           random_index_combination_ids,
                                                           c_look.gamma[0],
                                                           dilE,
                                                           4);
        }
      }
    }
    swatch.stop();
  }

  swatch.print();
}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void LapH::Correlators::build_C2c(std::vector<CorrInfo> const &corr_lookup,
                                  std::string const output_path,
                                  std::string const output_filename) {
  if (corr_lookup.empty())
    return;

  StopWatch swatch("C2c");
  swatch.start();

  // every element of corr_lookup contains the same filename. Wlog choose the
  // first element
  WriteHDF5Correlator filehandle(output_path, "C2+", output_filename,
                                 comp_type_factory_tr());

  for (auto const &c_look : corr_lookup) {
    std::vector<cmplx> correlator(Lt, cmplx(.0, .0));
    //    if(c_look.outfile.find("Check") == 0){
    //      for(int t1 = 0; t1 < Lt; t1++){
    //        for(const auto& corr : corrC[c_look.lookup[0]][t1][t1]){
    //          correlator[t1] += corr;
    //        }
    //      }
    //      // normalisation
    //      for(auto& corr : correlator){
    //        corr /= corrC[c_look.lookup[0]][0][0].size();
    //      }
    //      // write data to file
    //      filehandle.write(correlator, c_look);
    //    }
    //    else{
    DilutionScheme const dilution_scheme(Lt, dilT, DilutionType::block);
    for (auto const block_pair : dilution_scheme) {
      for (auto const slice_pair : block_pair) {
        int const t = get_time_delta(slice_pair, Lt);
        for (const auto &corr :
             corrC[c_look.lookup[0]][slice_pair.source()][slice_pair.sink()]) {
          correlator[t] += corr;
        }
      }
    }
    // normalisation
    for (auto &corr : correlator) {
      corr /= Lt * corrC[c_look.lookup[0]][0][0].size();
    }
    // write data to file
    filehandle.write(correlator, c_look);
    //    }
  }

  swatch.stop();
  swatch.print();
}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void LapH::Correlators::build_C4cD(CorrelatorLookup const &corr_lookup,
                                   std::string const output_path,
                                   std::string const output_filename) {
  if (corr_lookup.C4cD.empty())
    return;

  StopWatch swatch("C4cD");
  swatch.start();

  // every element of corr_lookup contains the same filename. Wlog choose the
  // first element
  WriteHDF5Correlator filehandle(output_path, "C4+D", output_filename,
                                 comp_type_factory_trtr());

  DilutionScheme const dilution_scheme(Lt, dilT, DilutionType::block);

  for (const auto &c_look : corr_lookup.C4cD) {

    std::vector<LapH::compcomp_t> correlator(Lt,
                                             LapH::compcomp_t(.0, .0, .0, .0));

    const size_t id0 = corr_lookup.corrC[c_look.lookup[0]].lookup[0];
    const size_t id1 = corr_lookup.corrC[c_look.lookup[1]].lookup[0];

    for (auto const block_pair : dilution_scheme) {
      for (auto const slice_pair : block_pair) {
        int const t = get_time_delta(slice_pair, Lt);

        std::vector<size_t> random_index_combination_ids =
            std::vector<size_t>({dil_fac_lookup.Q2V[id0].id_ric_lookup,
                                 dil_fac_lookup.Q2V[id1].id_ric_lookup});

        /*! @todo Write move assignment for compcomp_t and give trtr return
         * parameter */
        correlator[t] += trtr(
            corrC[c_look.lookup[0]][slice_pair.source()][slice_pair.sink()],
            corrC[c_look.lookup[1]][slice_pair.source()][slice_pair.sink()],
            ric_lookup, random_index_combination_ids);
      }
    }

    // normalisation
    for (auto &corr : correlator) {
      corr.rere /= (5U * 4U * 3U * 2U) * Lt;
      corr.reim /= (5U * 4U * 3U * 2U) * Lt;
      corr.imre /= (5U * 4U * 3U * 2U) * Lt;
      corr.imim /= (5U * 4U * 3U * 2U) * Lt;
    }

    // write data to file
    filehandle.write(correlator, c_look);
  }

  swatch.stop();
  swatch.print();
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void LapH::Correlators::build_C4cV(CorrelatorLookup const &corr_lookup,
                                   std::string const output_path,
                                   std::string const output_filename) {
  if (corr_lookup.C4cV.empty())
    return;

  StopWatch swatch("C4cV");
  swatch.start();

  // every element of corr_lookup contains the same filename. Wlog choose the
  // first element
  WriteHDF5Correlator filehandle(output_path, "C4+V", output_filename,
                                 comp_type_factory_trtr());

  DilutionScheme const dilution_scheme(Lt, dilT, DilutionType::block);

  for (const auto &c_look : corr_lookup.C4cV) {

    std::vector<LapH::compcomp_t> correlator(Lt,
                                             LapH::compcomp_t(.0, .0, .0, .0));

    const size_t id0 = corr_lookup.corrC[c_look.lookup[0]].lookup[0];
    const size_t id1 = corr_lookup.corrC[c_look.lookup[1]].lookup[0];

    for (auto const block_pair : dilution_scheme) {
      for (auto const slice_pair : block_pair) {
        int const t = get_time_delta(slice_pair, Lt);

        std::vector<size_t> random_index_combination_ids =
            std::vector<size_t>({dil_fac_lookup.Q2V[id0].id_ric_lookup,
                                 dil_fac_lookup.Q2V[id1].id_ric_lookup});

        /*! @todo Write move assignment for compcomp_t and give trtr return
         * parameter */
        correlator[t] += trtr(
            corrC[c_look.lookup[0]][slice_pair.source()][slice_pair.source()],
            corrC[c_look.lookup[1]][slice_pair.sink()][slice_pair.sink()],
            ric_lookup, random_index_combination_ids);
      }
    }

    // normalisation
    for (auto &corr : correlator) {
      corr.rere /= (5U * 4U * 3U * 2U) * Lt;
      corr.reim /= (5U * 4U * 3U * 2U) * Lt;
      corr.imre /= (5U * 4U * 3U * 2U) * Lt;
      corr.imim /= (5U * 4U * 3U * 2U) * Lt;
    }

    // write data to file
    filehandle.write(correlator, c_look);
  }

  swatch.stop();
  swatch.print();
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void LapH::Correlators::build_C4cC(RandomVector const &randomvectors,
                                   OperatorsForMesons const &meson_operator,
                                   Perambulator const &perambulators,
                                   std::vector<CorrInfo> const &corr_lookup,
                                   std::string const output_path,
                                   std::string const output_filename) {
  if (corr_lookup.empty())
    return;

  StopWatch swatch("C4cC");

  // every element of corr_lookup contains the same filename. Wlog choose the
  // first element
  WriteHDF5Correlator filehandle(output_path, "C4+C", output_filename,
                                 comp_type_factory_tr());

  std::vector<vec> correlator(corr_lookup.size(), vec(Lt, cmplx(.0, .0)));

  DilutionScheme const dilution_scheme(Lt, dilT, DilutionType::block);

// This is necessary to ensure the correct summation of the correlation function
#pragma omp parallel
  {
    swatch.start();
    std::vector<vec> C(corr_lookup.size(), vec(Lt, cmplx(.0, .0)));
    // building the quark line directly frees up a lot of memory
    QuarkLineBlock<QuarkLineType::Q0> quarkline_Q0(
        dilT, dilE, nev, dil_fac_lookup.Q0, ric_lookup);
    QuarkLineBlock<QuarkLineType::Q2V> quarkline_Q2V(
        dilT, dilE, nev, dil_fac_lookup.Q2V, ric_lookup);

    // creating memory arrays M1, M2 for intermediate storage of Quarklines
    // ------
    std::vector<std::vector<Eigen::MatrixXcd>> M1, M2;
    std::vector<std::array<size_t, 3>> M1_look;
    std::vector<std::array<size_t, 3>> M2_look;
    size_t M1_counter = 0;
    size_t M2_counter = 0;
    for (const auto &c_look : corr_lookup) {

      try {
        check_random_combinations<QuarkLineType::Q2V>(
            std::string("C4cC"), c_look.lookup, ric_lookup, dil_fac_lookup.Q0,
            dil_fac_lookup.Q2V);
      } catch (const std::length_error &le) {
        std::cerr << "Length error: " << le.what() << '\n';
      }

      size_t norm = 0;

      // creating lookup_table for M1 -----------------------------------------
      const size_t id0 = c_look.lookup[3];
      const size_t id1 = c_look.lookup[0];
      auto it1 = std::find_if(M1_look.begin(), M1_look.end(),
                              [&id0, &id1](std::array<size_t, 3> check) {
                                return (id0 == check[1] && id1 == check[2]);
                              });
      if (!(it1 != M1_look.end())) {
        M1.emplace_back(std::vector<Eigen::MatrixXcd>());

        M1_look.emplace_back(std::array<size_t, 3>({{M1_counter, id0, id1}}));

        M1_counter++;
      }

      // creating lookup_table for M2 -----------------------------------------
      const size_t id2 = c_look.lookup[1];
      const size_t id3 = c_look.lookup[2];
      auto it2 = std::find_if(M2_look.begin(), M2_look.end(),
                              [&id2, &id3](std::array<size_t, 3> check) {
                                return (id2 == check[1] && id3 == check[2]);
                              });
      if (!(it2 != M2_look.end())) {
        M2.emplace_back(std::vector<Eigen::MatrixXcd>());
        M2_look.emplace_back(std::array<size_t, 3>({{M2_counter, id2, id3}}));

        M2_counter++;
      }
    } // first run over lookuptable ends here - memory and new lookuptable
// are generated ------------------------------------------------------------

#pragma omp for schedule(dynamic)
    for (int b = 0; b < dilution_scheme.size(); ++b) {
      auto const block_pair = dilution_scheme[b];
      // creating quarklines
      quarkline_Q0.build_block_pair(randomvectors, meson_operator, block_pair,
                                  dil_fac_lookup.Q0, ric_lookup);
      quarkline_Q2V.build_block_pair(perambulators, meson_operator, block_pair,
                                  dil_fac_lookup.Q2V, ric_lookup);

      for (auto const slice_pair : block_pair) {
        int const t = get_time_delta(slice_pair, Lt);

        /*! Optimization by saving all products Q2V \cdot rVdaggerVr. Reuse not
         *  optimal as M1 and M2 don't share objects but time indices are
         *  identical in cross diagram
         */

        // build M1
        // ----------------------------------------------------------------
        for (const auto &look : M1_look) {

          std::vector<size_t> random_index_combination_ids =
              std::vector<size_t>(
                  {dil_fac_lookup.Q0 [look[1]].id_ric_lookup,
                   dil_fac_lookup.Q2L[look[2]].id_ric_lookup});

          rVdaggerVrxQ2(M1[look[0]], 
              quarkline_Q0(slice_pair.sink(), -1, look[1]),
              quarkline_Q2V(slice_pair.source(), slice_pair.sink_block(), look[2]),
              ric_lookup, random_index_combination_ids, dilE, 4);
        }
        // build M2
        // ----------------------------------------------------------------
        for (const auto &look : M2_look) {

          std::vector<size_t> random_index_combination_ids =
              std::vector<size_t>(
                  {dil_fac_lookup.Q0 [look[1]].id_ric_lookup,
                   dil_fac_lookup.Q2L[look[2]].id_ric_lookup});

          rVdaggerVrxQ2(M2[look[0]], 
              quarkline_Q0(slice_pair.sink(), -1, look[1]),
              quarkline_Q2V(slice_pair.source(), slice_pair.sink_block(), look[2]),
              ric_lookup, random_index_combination_ids, dilE, 4);
        }

        /*! Optimization by summing M2 over all random vectors before
         *  multiplying with M1
         */
        // Final summation for correlator
        // ------------------------------------------

        for (const auto &c_look : corr_lookup) {
          const size_t id0 = c_look.lookup[3];
          const size_t id1 = c_look.lookup[0];
          // TODO: This should be an access operator
          auto it1 = std::find_if(M1_look.begin(), M1_look.end(),
                                  [&id0, &id1](std::array<size_t, 3> check) {
                                    return (id0 == check[1] && id1 == check[2]);
                                  });
          const size_t id2 = c_look.lookup[1];
          const size_t id3 = c_look.lookup[2];
          auto it2 = std::find_if(M2_look.begin(), M2_look.end(),
                                  [&id2, &id3](std::array<size_t, 3> check) {
                                    return (id2 == check[1] && id3 == check[2]);
                                  });

          std::vector<size_t> random_index_combination_ids =
              std::vector<size_t>(
                {dil_fac_lookup.Q0 [c_look.lookup[3]].id_ric_lookup,
                 dil_fac_lookup.Q2L[c_look.lookup[0]].id_ric_lookup,
                 dil_fac_lookup.Q0 [c_look.lookup[1]].id_ric_lookup,
                 dil_fac_lookup.Q2L[c_look.lookup[2]].id_ric_lookup});

          // M1 and M2 implicitly contain time indices. Thus += over time is
          // necessary
          C[c_look.id][t] += trace(
              M1[(*it1)[0]], M2[(*it2)[0]], ric_lookup,
              random_index_combination_ids, dilE, 4);

        } // loop over operators ends here
      }
    } // loops over time end here
#pragma omp critical
    {
      for (const auto &c_look : corr_lookup)
        for (size_t t = 0; t < Lt; t++)
          correlator[c_look.id][t] += C[c_look.id][t];
    }

    swatch.stop();
  } // parallel part ends here

  // normalisation
  for (const auto &c_look : corr_lookup) {
    for (auto &corr : correlator[c_look.id]) {
      corr /= (6 * 5 * 4 * 3) * Lt; // TODO: Hard Coded atm - Be carefull
    }
    // write data to file
    filehandle.write(correlator[c_look.id], c_look);
  }

  swatch.print();
}

/*****************************************************************************/
/*                                 build_C3c                                 */
/*****************************************************************************/

void LapH::Correlators::build_C3c(RandomVector const &randomvectors,
                                  OperatorsForMesons const &meson_operator,
                                  Perambulator const &perambulators,
                                  std::vector<CorrInfo> const &corr_lookup,
                                  std::string const output_path,
                                  std::string const output_filename) {
  if (corr_lookup.empty())
    return;

  // every element of corr_lookup contains the same filename. Wlog choose the
  // first element
  WriteHDF5Correlator filehandle(output_path, "C3+", output_filename,
                                 comp_type_factory_tr());

  StopWatch swatch("C3c");

  DilutionScheme const dilution_scheme(Lt, dilT, DilutionType::block);

  std::vector<vec> correlator(corr_lookup.size(), vec(Lt, cmplx(.0, .0)));

// This is necessary to ensure the correct summation of the correlation function
#pragma omp parallel
  {
    swatch.start();
    std::vector<vec> C(corr_lookup.size(), vec(Lt, cmplx(.0, .0)));
    // building the quark line directly frees up a lot of memory
    QuarkLineBlock<QuarkLineType::Q0> quarkline_Q0(
        dilT, dilE, nev, dil_fac_lookup.Q0, ric_lookup);
    QuarkLineBlock<QuarkLineType::Q2L> quarkline_Q2L(
        dilT, dilE, nev, dil_fac_lookup.Q2L, ric_lookup);
    QuarkLineBlock<QuarkLineType::Q1> quarkline_Q1(
        dilT, dilE, nev, dil_fac_lookup.Q1, ric_lookup);

    // creating memory arrays M1, M2 for intermediate storage of Quarklines
    // ------
    std::vector<std::vector<Eigen::MatrixXcd>> M1, M2;

    std::vector<std::array<size_t, 3>> M1_look;
    std::vector<std::array<size_t, 2>> M2_look;

    size_t M1_counter = 0;
    size_t M2_counter = 0;

    for (const auto &c_look : corr_lookup) {
      const auto &ric0 =
          ric_lookup[dil_fac_lookup.Q2L[c_look.lookup[0]].id_ric_lookup]
              .rnd_vec_ids;
      const auto &ric1 =
          ric_lookup[dil_fac_lookup.Q1[c_look.lookup[1]].id_ric_lookup]
              .rnd_vec_ids;
      const auto &ric2 =
          ric_lookup[dil_fac_lookup.Q0[c_look.lookup[2]].id_ric_lookup]
              .rnd_vec_ids;
      if (ric0.size() != ric1.size() || ric0.size() != ric2.size()) {
        std::cout << "rnd combinations are not the same in build_C3+"
                  << std::endl;
        exit(0);
      }

      // creating memeory for M1
      // -------------------------------------------------
      /*! @warning For some reason indices 2 and 0 are interchanged here */
      const size_t id0 = c_look.lookup[2];
      const size_t id2 = c_look.lookup[0];
      auto it1 = std::find_if(M1_look.begin(), M1_look.end(),
                              [&id0, &id2](std::array<size_t, 3> check) {
                                return (id0 == check[1] && id2 == check[2]);
                              });

      if (!(it1 != M1_look.end())) {
        M1.emplace_back(std::vector<Eigen::MatrixXcd>());
        M1_look.emplace_back(std::array<size_t, 3>({{M1_counter, id0, id2}}));
        ++M1_counter;
      }

      // creating memory for M2
      // -------------------------------------------------
      const size_t id1 = c_look.lookup[1];
      auto it2 = std::find_if(
          M2_look.begin(), M2_look.end(),
          [&id1](std::array<size_t, 2> check) { return (id1 == check[1]); });

      /*! No reuse between different gamma structures as for every entry of
       * corr_lookup the whole thing is rebuilt
       */
      if (!(it2 != M2_look.end())) {
        M2.emplace_back(std::vector<Eigen::MatrixXcd>());
        M2_look.emplace_back(std::array<size_t, 2>({{M2_counter, id1}}));
        ++M2_counter;
      }

    } // first run over lookup table ends here - memory and new lookup table
// are generated ------------------------------------------------------------

#pragma omp for schedule(dynamic)
    for (int b = 0; b < dilution_scheme.size(); ++b) {
      auto const block_pair = dilution_scheme[b];
      // creating quarklines
      quarkline_Q0.build_block_pair(randomvectors, meson_operator, block_pair,
                                      dil_fac_lookup.Q0, ric_lookup);
      quarkline_Q2L.build_block_pair(perambulators, meson_operator, block_pair,
                                      dil_fac_lookup.Q2L, ric_lookup);
      quarkline_Q1.build_block_pair(perambulators, meson_operator, block_pair,
                                     dil_fac_lookup.Q1, ric_lookup);

      for (auto const slice_pair : block_pair) {
        int const t = get_time_delta(slice_pair, Lt);
        // build M1
        // ----------------------------------------------------------------
        for (const auto &look : M1_look) {
          
          std::vector<size_t> random_index_combination_ids =
              std::vector<size_t>(
                  {dil_fac_lookup.Q0 [look[1]].id_ric_lookup,
                   dil_fac_lookup.Q2L[look[2]].id_ric_lookup});

          rVdaggerVrxQ2(M1[look[0]], 
              quarkline_Q0(slice_pair.source(), -1, look[1]),
              quarkline_Q2L(slice_pair.source(), slice_pair.sink_block(), look[2]),
              ric_lookup, random_index_combination_ids, dilE, 4);
        }

        // build M2
        // ----------------------------------------------------------------
        for (const auto &look : M2_look) {

          std::vector<size_t> random_index_combination_ids =
              std::vector<size_t>(
                  {dil_fac_lookup.Q1[look[1]].id_ric_lookup});

          Q1(M2[look[0]], 
             quarkline_Q1(slice_pair.sink(), slice_pair.source_block(), look[1]),
             ric_lookup, 
             random_index_combination_ids,
             dilE, 4);
        }

        // Final summation for correlator
        // ------------------------------------------
        for (const auto &c_look : corr_lookup) {
          const size_t id0 = c_look.lookup[2];
          const size_t id1 = c_look.lookup[1];
          const size_t id2 = c_look.lookup[0];
          auto it1 = std::find_if(M1_look.begin(), M1_look.end(),
                                  [&id0, &id2](std::array<size_t, 3> check) {
                                    return (id0 == check[1] && id2 == check[2]);
                                  });
          auto it2 = std::find_if(M2_look.begin(), M2_look.end(),
                                  [&id1](std::array<size_t, 2> check) {
                                    return (id1 == check[1]);
                                  });

          std::vector<size_t> random_index_combination_ids =
              std::vector<size_t>(
                {dil_fac_lookup.Q2L[id0].id_ric_lookup,
                 dil_fac_lookup.Q1 [id1].id_ric_lookup,
                 dil_fac_lookup.Q0 [id2].id_ric_lookup});

          C[c_look.id][t] += trace_3pt(
              M2[(*it2)[0]], M1[(*it1)[0]], ric_lookup,
              random_index_combination_ids, dilE, 4);

        }
      }
    } // loops over time end here
#pragma omp critical
    {
      for (const auto &c_look : corr_lookup)
        for (size_t t = 0; t < Lt; t++)
          correlator[c_look.id][t] += C[c_look.id][t];
    }
    swatch.stop();
  } // parallel part ends here

  // normalisation
  for (const auto &c_look : corr_lookup) {
    for (auto &corr : correlator[c_look.id]) {
      /*! @todo Hard Coded atm - Be careful */
      corr /= (6 * 5 * 4) * Lt;
    }
    // write data to file
    filehandle.write(correlator[c_look.id], c_look);
  }
  swatch.print();
}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void LapH::Correlators::build_C4cB(RandomVector const &randomvectors,
                                   OperatorsForMesons const &meson_operator,
                                   Perambulator const &perambulators,
                                   std::vector<CorrInfo> const &corr_lookup,
                                   std::string const output_path,
                                   std::string const output_filename) {
  if (corr_lookup.empty())
    return;

  StopWatch swatch("C4cB");

  WriteHDF5Correlator filehandle(
      output_path, "C4+B", output_filename, comp_type_factory_tr());

  std::vector<vec> correlator(corr_lookup.size(), vec(Lt, cmplx(.0, .0)));

  DilutionScheme const dilution_scheme(Lt, dilT, DilutionType::block);

// This is necessary to ensure the correct summation of the correlation function
#pragma omp parallel
  {
    swatch.start();
    std::vector<vec> C(corr_lookup.size(), vec(Lt, cmplx(.0, .0)));
    // building the quark line directly frees up a lot of memory
    QuarkLineBlock<QuarkLineType::Q0> quarkline_Q0(
        dilT, dilE, nev, dil_fac_lookup.Q0, ric_lookup);

    QuarkLineBlock<QuarkLineType::Q2L> quarkline_Q2L(
        dilT, dilE, nev, dil_fac_lookup.Q2L, ric_lookup);

    // Create lookup tables for M1, M2 to save intermediate results that can be
    // reused (rVdaggerVr \cdot Q2V).
    std::vector<std::vector<Eigen::MatrixXcd>> M1, M2;
    std::vector<std::array<size_t, 3>> M1_look;
    std::vector<std::array<size_t, 3>> M2_look;

    //! @TODO (Martin Ueding): Remove these and use `M1.size()` instead.
    size_t M1_counter = 0;
    size_t M2_counter = 0;

    for (const auto &c_look : corr_lookup) {
      try {
        check_random_combinations<QuarkLineType::Q2L>(std::string("C4cB"),
                                                      c_look.lookup,
                                                      ric_lookup,
                                                      dil_fac_lookup.Q0,
                                                      dil_fac_lookup.Q2L);
      } catch (const std::length_error &le) {
        std::cerr << "Length error: " << le.what() << '\n';
      }

      // Create lookup table for M1.
      auto const id3 = c_look.lookup[3];
      auto const id0 = c_look.lookup[0];
      auto it1 = std::find_if(
          M1_look.begin(), M1_look.end(), [&id3, &id0](std::array<size_t, 3> check) {
            return (id3 == check[1] && id0 == check[2]);
          });
      if (!(it1 != M1_look.end())) {
        M1.emplace_back(std::vector<Eigen::MatrixXcd>());

        M1_look.emplace_back(std::array<size_t, 3>({{M1_counter, id3, id0}}));
        M1_counter++;
      }

      // Create lookup table for M2.
      auto const id1 = c_look.lookup[1];
      auto const id2 = c_look.lookup[2];
      auto it2 = std::find_if(
          M2_look.begin(), M2_look.end(), [&id1, &id2](std::array<size_t, 3> check) {
            return (id1 == check[1] && id2 == check[2]);
          });
      if (!(it2 != M2_look.end())) {
        M2.emplace_back(std::vector<Eigen::MatrixXcd>());
        M2_look.emplace_back(std::array<size_t, 3>({{M2_counter, id1, id2}}));
        M2_counter++;
      }
    }

#pragma omp for schedule(dynamic)
    // Perform contraction here
    for (int b = 0; b < dilution_scheme.size(); ++b) {
      auto const block_pair = dilution_scheme[b];
      // Create quarklines for all time combinations in block_pair
      quarkline_Q0.build_block_pair(
          randomvectors, meson_operator, block_pair, dil_fac_lookup.Q0, ric_lookup);
      quarkline_Q2L.build_block_pair(
          perambulators, meson_operator, block_pair, dil_fac_lookup.Q2L, ric_lookup);

      for (auto const slice_pair : block_pair) {
        int const t = get_time_delta(slice_pair, Lt);

        // Calculate M1
        for (const auto &look : M1_look) {
          std::vector<size_t> random_index_combination_ids{
              dil_fac_lookup.Q0[look[1]].id_ric_lookup,
              dil_fac_lookup.Q2L[look[2]].id_ric_lookup};

          rVdaggerVrxQ2(
              M1[look[0]],
              quarkline_Q0(slice_pair.source(), -1, look[1]),
              quarkline_Q2L(slice_pair.source(), slice_pair.sink_block(), look[2]),
              ric_lookup,
              random_index_combination_ids,
              dilE,
              4);
        }

        // Calculate M2
        for (const auto &look : M2_look) {
          std::vector<size_t> random_index_combination_ids{
              dil_fac_lookup.Q0[look[1]].id_ric_lookup,
              dil_fac_lookup.Q2L[look[2]].id_ric_lookup};

          rVdaggerVrxQ2(
              M2[look[0]],
              quarkline_Q0(slice_pair.sink(), -1, look[1]),
              quarkline_Q2L(slice_pair.sink(), slice_pair.source_block(), look[2]),
              ric_lookup,
              random_index_combination_ids,
              dilE,
              4);
        }

        // Calculate final sum now only depending on M1 and M2
        for (const auto &c_look : corr_lookup) {
          /*! @TODO These are essentially getters for M1_look and M2_look. 
           *        Simplify that!
           */
          auto const id3 = c_look.lookup[3];
          auto const id0 = c_look.lookup[0];
          auto it1 = std::find_if(
              M1_look.begin(), M1_look.end(), [&id3, &id0](std::array<size_t, 3> check) {
                return (id3 == check[1] && id0 == check[2]);
              });
          auto const id1 = c_look.lookup[1];
          auto const id2 = c_look.lookup[2];
          auto it2 = std::find_if(
              M2_look.begin(), M2_look.end(), [&id1, &id2](std::array<size_t, 3> check) {
                return (id1 == check[1] && id2 == check[2]);
              });

          std::vector<size_t> random_index_combination_ids{
              dil_fac_lookup.Q0[c_look.lookup[3]].id_ric_lookup,
              dil_fac_lookup.Q2L[c_look.lookup[0]].id_ric_lookup,
              dil_fac_lookup.Q0[c_look.lookup[1]].id_ric_lookup,
              dil_fac_lookup.Q2L[c_look.lookup[2]].id_ric_lookup};

          // Perform trace in random index-, eigen- and dirac space
          // += because multiple block pairs contribute to the same t
          // Magic number 4 is number_of_blocks in Dilution of Dirac space
          C[c_look.id][t] += trace(M1[(*it1)[0]],
                                   M2[(*it2)[0]],
                                   ric_lookup,
                                   random_index_combination_ids,
                                   dilE,
                                   4);

        } // loop over operators ends here
      }
    } // loops over time end here
#pragma omp critical
    {
      for (const auto &c_look : corr_lookup)
        for (size_t t = 0; t < Lt; t++)
          correlator[c_look.id][t] += C[c_look.id][t];
    }
    swatch.stop();
  } // parallel part ends here

  // normalisation
  for (const auto &c_look : corr_lookup) {
    for (auto &corr : correlator[c_look.id]) {
      // @todo Hard Coded atm - Be careful
      corr /= (6 * 5 * 4 * 3) * Lt; 
    }
    // write data to file
    filehandle.write(correlator[c_look.id], c_look);
  }
  swatch.print();
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void LapH::Correlators::build_C30(OperatorsForMesons const &meson_operator,
                                  Perambulator const &perambulators,
                                  std::vector<CorrInfo> const &corr_lookup,
                                  std::string const output_path,
                                  std::string const output_filename) {
  if (corr_lookup.empty())
    return;

  StopWatch swatch("C30");
  swatch.start();

  // every element of corr_lookup contains the same filename. Wlog choose the
  // first element
  WriteHDF5Correlator filehandle(output_path, "C30", output_filename,
                                 comp_type_factory_tr());

  std::vector<vec> correlator(corr_lookup.size(), vec(Lt, cmplx(.0, .0)));

  DilutionScheme const dilution_scheme(Lt, dilT, DilutionType::block);

#pragma omp parallel
  {
    swatch.start();
    // std::vector<vec> C(corr_lookup.size(), vec(Lt, cmplx(.0, .0)));
    std::vector<vec> C(corr_lookup.size(), vec(Lt, cmplx(.0, .0)));
    // building the quark line directly frees up a lot of memory
    QuarkLineBlock<QuarkLineType::Q1> quarklines(dilT, dilE, nev,
                                                 dil_fac_lookup.Q1, ric_lookup);

    // creating memory arrays M1, M2 for intermediate storage of Quarklines
    // ------
    std::vector<std::vector<Eigen::MatrixXcd>> L1, L2;
    std::vector<std::array<size_t, 3>> L1_look;
    std::vector<std::array<size_t, 2>> L2_look;
    size_t L1_counter = 0;
    size_t L2_counter = 0;

    for (const auto &c_look : corr_lookup) {

      const auto &ric0 =
          ric_lookup[dil_fac_lookup.Q1[c_look.lookup[0]].id_ric_lookup]
              .rnd_vec_ids;
      const auto &ric1 =
          ric_lookup[dil_fac_lookup.Q1[c_look.lookup[1]].id_ric_lookup]
              .rnd_vec_ids;
      const auto &ric2 =
          ric_lookup[dil_fac_lookup.Q1[c_look.lookup[2]].id_ric_lookup]
              .rnd_vec_ids;
      if (ric0.size() != ric1.size() || ric0.size() != ric2.size()) {
        std::cout << "rnd combinations are not the same in build_C30"
                  << std::endl;
        exit(0);
      }

      // creating memeory for L1
      // -------------------------------------------------
      const size_t id0 = c_look.lookup[0];
      const size_t id1 = c_look.lookup[1];
      auto it1 = std::find_if(L1_look.begin(), L1_look.end(),
                              [&id0, &id1](std::array<size_t, 3> check) {
                                return (id0 == check[1] && id1 == check[2]);
                              });
      if (!(it1 != L1_look.end())) {
        L1.emplace_back(std::vector<Eigen::MatrixXcd>());

        L1_look.emplace_back(std::array<size_t, 3>({{L1_counter, id0, id1}}));
        L1_counter++;
      }
      // creating memeory for L2
      // -------------------------------------------------
      const size_t id2 = c_look.lookup[2];
      auto it2 = std::find_if(
          L2_look.begin(), L2_look.end(), [&id2](std::array<size_t, 2> check) {
            return (id2 == check[1]);
          });
      if (!(it2 != L2_look.end())) {
        L2.emplace_back(std::vector<Eigen::MatrixXcd>());
        L2_look.emplace_back(std::array<size_t, 2>({{L2_counter, id2}}));
        L2_counter++;
      }
    }  // first run over lookuptable ends here - memory and new lookuptable
// are generated ------------------------------------------------------------

#pragma omp for schedule(dynamic)
    for (int b = 0; b < dilution_scheme.size(); ++b) {
      auto const block_pair = dilution_scheme[b];

      quarklines.build_block_pair(perambulators, meson_operator, block_pair,
                                  dil_fac_lookup.Q1, ric_lookup);

      for (auto const slice_pair : block_pair) {
        int const t = get_time_delta(slice_pair, Lt);

        // build L1
        // ----------------------------------------------------------------
        for (const auto &look : L1_look) {

          std::vector<size_t> random_index_combination_ids =
              std::vector<size_t>(
                  {dil_fac_lookup.Q1[look[1]].id_ric_lookup,
                   dil_fac_lookup.Q1[look[2]].id_ric_lookup});

          Q1xQ1(L1[look[0]], 
                quarklines(slice_pair.source(), slice_pair.sink_block(), look[1]),
                quarklines(slice_pair.sink(), slice_pair.source_block(), look[2]),
                ric_lookup, random_index_combination_ids, dilE, 4);
        }

        // build L2
        // ----------------------------------------------------------------
        for (const auto &look : L2_look) {
          std::vector<size_t> random_index_combination_ids =
              std::vector<size_t>(
                  {dil_fac_lookup.Q1[look[1]].id_ric_lookup});

          Q1(L2[look[0]], 
             quarklines(slice_pair.source(), slice_pair.source_block(), look[1]),
             ric_lookup, 
             random_index_combination_ids,
             dilE, 4);
        }

        for (const auto &c_look : corr_lookup) {

          const size_t id0 = c_look.lookup[0];
          const size_t id1 = c_look.lookup[1];
          const size_t id2 = c_look.lookup[2];

          auto it1 = std::find_if(L1_look.begin(), L1_look.end(),
                                  [&id0, &id1](std::array<size_t, 3> check) {
                                    return (id0 == check[1] && id1 == check[2]);
                                  });
          auto it2 = std::find_if(L2_look.begin(), L2_look.end(),
                                  [&id2](std::array<size_t, 2> check) {
                                    return (id2 == check[1]);
                                  });

          std::vector<size_t> random_index_combination_ids =
              std::vector<size_t>(
                  {dil_fac_lookup.Q1[c_look.lookup[0]].id_ric_lookup,
                   dil_fac_lookup.Q1[c_look.lookup[1]].id_ric_lookup,
                   dil_fac_lookup.Q1[c_look.lookup[2]].id_ric_lookup});

          C[c_look.id][t] += trace_3pt(
              L2[(*it2)[0]], L1[(*it1)[0]], ric_lookup,
              random_index_combination_ids, dilE, 4);

        }
      }
    } // loop over time ends here

#pragma omp critical
    {
      for (const auto &c_look : corr_lookup)
        for (size_t t = 0; t < Lt; t++)
          correlator[c_look.id][t] += C[c_look.id][t];
    }
    swatch.stop();
  } // parallel part ends here

  // normalisation
  for (const auto &c_look : corr_lookup) {
    for (auto &corr : correlator[c_look.id]) {
      corr /= (5 * 4 * 3) * Lt; // TODO: Hard Coded atm - Be carefull
    }
    // write data to file
    filehandle.write(correlator[c_look.id], c_look);
  }
  swatch.print();
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
/*!
 *  @note Not optimised.
 */
void LapH::Correlators::build_C40C(OperatorsForMesons const &meson_operator,
                                   Perambulator const &perambulators,
                                   std::vector<CorrInfo> const &corr_lookup,
                                   std::string const output_path,
                                   std::string const output_filename) {
  if (corr_lookup.empty())
    return;

  StopWatch swatch("C40C");
  swatch.start();

  // every element of corr_lookup contains the same filename. Wlog choose the
  // first element
  WriteHDF5Correlator filehandle(output_path, "C40C", output_filename,
                                 comp_type_factory_tr());

  std::vector<vec> correlator(corr_lookup.size(), vec(Lt, cmplx(.0, .0)));

  DilutionScheme const dilution_scheme(Lt, dilT, DilutionType::block);

#pragma omp parallel
  {
    swatch.start();
    // std::vector<vec> C(corr_lookup.size(), vec(Lt, cmplx(.0, .0)));
    std::vector<vec> C(corr_lookup.size(), vec(Lt, cmplx(.0, .0)));
    // building the quark line directly frees up a lot of memory
    QuarkLineBlock<QuarkLineType::Q1> quarklines(dilT, dilE, nev,
                                                 dil_fac_lookup.Q1, ric_lookup);

    // creating memory arrays M1, M2 for intermediate storage of Quarklines
    // ------
    std::vector<std::vector<Eigen::MatrixXcd>> L1, L2;
    std::vector<std::array<size_t, 3>> L1_look;
    std::vector<std::array<size_t, 3>> L2_look;
    size_t L1_counter = 0;
    size_t L2_counter = 0;

    for (const auto &c_look : corr_lookup) {

      //const auto &ric0 =
      //    ric_lookup[dil_fac_lookup.Q1[c_look.lookup[0]].id_ric_lookup]
      //        .rnd_vec_ids;
      //const auto &ric1 =
      //    ric_lookup[dil_fac_lookup.Q1[c_look.lookup[1]].id_ric_lookup]
      //        .rnd_vec_ids;
      //const auto &ric2 =
      //    ric_lookup[dil_fac_lookup.Q1[c_look.lookup[2]].id_ric_lookup]
      //        .rnd_vec_ids;
      //const auto &ric3 =
      //    ric_lookup[dil_fac_lookup.Q1[c_look.lookup[3]].id_ric_lookup]
      //        .rnd_vec_ids;
      //    if(ric0.size() != ric1.size() || ric0.size() != ric2.size() ||
      //       ric0.size() != ric3.size()){
      //      std::cout << "rnd combinations are not the same in build_corr0"
      //                << std::endl;
      //      exit(0);
      //    }

      //      try {
      //        check_random_combinations<QuarkLineType::Q2L>(
      //            std::string("C4cB"), c_look.lookup, ric_lookup,
      //            dil_fac_lookup.Q0, dil_fac_lookup.Q2L);
      //      }
      //      catch (const std::length_error& le) {
      //        std::cerr << "Length error: " << le.what() << '\n';
      //      }

      // creating memeory for L1
      // -------------------------------------------------
      const size_t id0 = c_look.lookup[3];
      const size_t id1 = c_look.lookup[0];
      auto it1 = std::find_if(L1_look.begin(), L1_look.end(),
                              [&id0, &id1](std::array<size_t, 3> check) {
                                return (id0 == check[1] && id1 == check[2]);
                              });
      if (!(it1 != L1_look.end())) {
        L1.emplace_back(std::vector<Eigen::MatrixXcd>());

        L1_look.emplace_back(std::array<size_t, 3>({{L1_counter, id0, id1}}));
        L1_counter++;
      }
      // creating memeory for L2
      // -------------------------------------------------
      const size_t id2 = c_look.lookup[1];
      const size_t id3 = c_look.lookup[2];
      auto it2 = std::find_if(L2_look.begin(), L2_look.end(),
                              [&id2, &id3](std::array<size_t, 3> check) {
                                return (id2 == check[1] && id3 == check[2]);
                              });
      if (!(it2 != L2_look.end())) {
        L2.emplace_back(std::vector<Eigen::MatrixXcd>());
        L2_look.emplace_back(std::array<size_t, 3>({{L2_counter, id2, id3}}));
        L2_counter++;
      }
    } // first run over lookuptable ends here - memory and new lookuptable
// are generated ------------------------------------------------------------

#pragma omp for schedule(dynamic)
    for (int b = 0; b < dilution_scheme.size(); ++b) {
      auto const block_pair = dilution_scheme[b];

      quarklines.build_block_pair(perambulators, meson_operator, block_pair,
                                  dil_fac_lookup.Q1, ric_lookup);

      for (auto const slice_pair : block_pair) {
        int const t = get_time_delta(slice_pair, Lt);

        // build L1
        // ----------------------------------------------------------------
        for (const auto &look : L1_look) {
          std::vector<size_t> random_index_combination_ids =
              std::vector<size_t>(
                  {dil_fac_lookup.Q1[look[1]].id_ric_lookup,
                   dil_fac_lookup.Q1[look[2]].id_ric_lookup});

          Q1xQ1(L1[look[0]], 
                quarklines(slice_pair.sink(), slice_pair.source_block(), look[1]),
                quarklines(slice_pair.source(), slice_pair.sink_block(), look[2]),
                ric_lookup, random_index_combination_ids, dilE, 4);
        }

        // build L2
        // ----------------------------------------------------------------
        for (const auto &look : L2_look) {
          std::vector<size_t> random_index_combination_ids =
              std::vector<size_t>(
                  {dil_fac_lookup.Q1[look[1]].id_ric_lookup,
                   dil_fac_lookup.Q1[look[2]].id_ric_lookup});

          Q1xQ1(L2[look[0]], 
                quarklines(slice_pair.sink(), slice_pair.source_block(), look[1]),
                quarklines(slice_pair.source(), slice_pair.sink_block(), look[2]),
                ric_lookup, random_index_combination_ids, dilE, 4);
        }

        for (const auto &c_look : corr_lookup) {

          const size_t id0 = c_look.lookup[3];
          const size_t id1 = c_look.lookup[0];
          const size_t id2 = c_look.lookup[1];
          const size_t id3 = c_look.lookup[2];

          auto it1 = std::find_if(L1_look.begin(), L1_look.end(),
                                  [&id0, &id1](std::array<size_t, 3> check) {
                                    return (id0 == check[1] && id1 == check[2]);
                                  });
          auto it2 = std::find_if(L2_look.begin(), L2_look.end(),
                                  [&id2, &id3](std::array<size_t, 3> check) {
                                    return (id2 == check[1] && id3 == check[2]);
                                  });

          std::vector<size_t> random_index_combination_ids =
              std::vector<size_t>(
                {dil_fac_lookup.Q1[c_look.lookup[3]].id_ric_lookup,
                 dil_fac_lookup.Q1[c_look.lookup[0]].id_ric_lookup,
                 dil_fac_lookup.Q1[c_look.lookup[1]].id_ric_lookup,
                 dil_fac_lookup.Q1[c_look.lookup[2]].id_ric_lookup});

          C[c_look.id][t] += trace(
              L1[(*it1)[0]], L2[(*it2)[0]], ric_lookup,
              random_index_combination_ids, dilE, 4);

        }
      }
    }
#pragma omp critical
    {
      for (const auto &c_look : corr_lookup)
        for (size_t t = 0; t < Lt; t++)
          correlator[c_look.id][t] += C[c_look.id][t];
    }
    swatch.stop();
  } // parallel part ends here

  // normalisation
  for (const auto &c_look : corr_lookup) {
    for (auto &corr : correlator[c_look.id]) {
      corr /= (5 * 4 * 3 * 2) * Lt; // TODO: Hard Coded atm - Be carefull
    }
    // write data to file
    filehandle.write(correlator[c_look.id], c_look);
  }
  swatch.print();
}

void make_Q1_Q1_map(LapH::OperatorToFactorMap<2> &L,
                    size_t const op_id0,
                    size_t const op_id1,
                    std::vector<QuarklineQ1Indices> const &Q1_lookup,
                    LapH::OperatorToFactorMap<1> const &factor0,
                    LapH::OperatorToFactorMap<1> const &factor1,
                    std::vector<RandomIndexCombinationsQ2> const &ric_lookup,
                    int const dilE,
                    int const dilD) {
  typename LapH::OperatorToFactorMap<2>::key_type const key = {op_id0, op_id1};

  if (L.count(key) == 0) {
    L[key] = factor0.at({op_id0}) * factor1.at({op_id1});
  }
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
/*!
 *  @note Not optimised.
 */
void LapH::Correlators::build_C40B(OperatorsForMesons const &meson_operator,
                                   Perambulator const &perambulators,
                                   std::vector<CorrInfo> const &corr_lookup,
                                   std::string const output_path,
                                   std::string const output_filename) {

  if (corr_lookup.empty())
    return;

  StopWatch swatch("C40B");

  // every element of corr_lookup contains the same filename. Wlog choose the
  // first element
  WriteHDF5Correlator filehandle(output_path, "C40B", output_filename,
                                 comp_type_factory_tr());

  std::vector<vec> correlator(corr_lookup.size(), vec(Lt, cmplx(.0, .0)));

  DilutionScheme const dilution_scheme(Lt, dilT, DilutionType::block);

#pragma omp parallel
  {
    swatch.start();
    // std::vector<vec> C(corr_lookup.size(), vec(Lt, cmplx(.0, .0)));
    std::vector<vec> C(corr_lookup.size(), vec(Lt, cmplx(.0, .0)));
    // building the quark line directly frees up a lot of memory
    QuarkLineBlock2<QuarkLineType::Q1> quarklines(
        dilT, dilE, nev, dil_fac_lookup.Q1, ric_lookup);

      //const auto &ric0 =
      //    ric_lookup[dil_fac_lookup.Q1[c_look.lookup[0]].id_ric_lookup]
      //        .rnd_vec_ids;
      //const auto &ric1 =
      //    ric_lookup[dil_fac_lookup.Q1[c_look.lookup[1]].id_ric_lookup]
      //        .rnd_vec_ids;
      //const auto &ric2 =
      //    ric_lookup[dil_fac_lookup.Q1[c_look.lookup[2]].id_ric_lookup]
      //        .rnd_vec_ids;
      //const auto &ric3 =
      //    ric_lookup[dil_fac_lookup.Q1[c_look.lookup[3]].id_ric_lookup]
      //        .rnd_vec_ids;
      //    if(ric0.size() != ric1.size() || ric0.size() != ric2.size() ||
      //       ric0.size() != ric3.size()){
      //      std::cout << "rnd combinations are not the same in build_corr0"
      //                << std::endl;
      //      exit(0);
      //    }

      //      try {
      //        check_random_combinations<QuarkLineType::Q2L>(
      //            std::string("C4cB"), c_look.lookup, ric_lookup,
      //            dil_fac_lookup.Q0, dil_fac_lookup.Q2L);
      //      }
      //      catch (const std::length_error& le) {
      //        std::cerr << "Length error: " << le.what() << '\n';
      //      }


#pragma omp for schedule(dynamic)
    for (int b = 0; b < dilution_scheme.size(); ++b) {
      auto const block_pair = dilution_scheme[b];

      std::cout << block_pair << std::endl;

      quarklines.build_block_pair(perambulators, meson_operator, block_pair,
                                  dil_fac_lookup.Q1, ric_lookup);

      for (auto const slice_pair : block_pair) {
        int const t = get_time_delta(slice_pair, Lt);

        OperatorToFactorMap<2> L1;
        for (const auto &c_look : corr_lookup) {
          auto const id0 = static_cast<size_t>(c_look.lookup[3]);
          auto const id1 = static_cast<size_t>(c_look.lookup[0]);
          make_Q1_Q1_map(L1,
                         id0,
                         id1,
                         dil_fac_lookup.Q1,
                         quarklines(slice_pair.source(), slice_pair.source_block()),
                         quarklines(slice_pair.source(), slice_pair.sink_block()),
                         ric_lookup,
                         dilE,
                         4);
        }

        OperatorToFactorMap<2> L2;
        for (const auto &c_look : corr_lookup) {
          auto const id0 = static_cast<size_t>(c_look.lookup[1]);
          auto const id1 = static_cast<size_t>(c_look.lookup[2]);
          make_Q1_Q1_map(L2,
                         id0,
                         id1,
                         dil_fac_lookup.Q1,
                         quarklines(slice_pair.sink(), slice_pair.sink_block()),
                         quarklines(slice_pair.sink(), slice_pair.source_block()),
                         ric_lookup,
                         dilE,
                         4);
        }

        for (const auto &c_look : corr_lookup) {
          const size_t id0 = c_look.lookup[3];
          const size_t id1 = c_look.lookup[0];
          const size_t id2 = c_look.lookup[1];
          const size_t id3 = c_look.lookup[2];

          typename OperatorToFactorMap<2>::key_type const key1 = {id0, id1};
          typename OperatorToFactorMap<2>::key_type const key2 = {id2, id3};

          C[c_look.id][t] += trace(L1[key1], L2[key2]);
        }
      }
    }
#pragma omp critical
    {
      for (const auto &c_look : corr_lookup)
        for (size_t t = 0; t < Lt; t++)
          correlator[c_look.id][t] += C[c_look.id][t];
    }
    swatch.stop();
  } // parallel part ends here

  // normalisation
  for (const auto &c_look : corr_lookup) {
    for (auto &corr : correlator[c_look.id]) {
      corr /= (5 * 4 * 3 * 2) * Lt; // TODO: Hard Coded atm - Be carefull
    }
    // write data to file
    filehandle.write(correlator[c_look.id], c_look);
  }
  swatch.print();
}

/******************************************************************************/
/*!
 *  @param quarklines       Instance of Quarklines. Contains prebuilt
 *                          combinations of operators and perambulators
 *  @param meson_operator   Instance of LapH::OperatorsForMesons. Contains
 *                          operators (@f$ V^\dagger V $f$) with momenta
 *                          and with/without dilution.
 *  @param perambulators    Instance of LapH::Perambulator class. Contains
 *                          Perambulator data
 *  @param operator_lookup
 *  @param corr_lookup
 *  @param quark_lookup
 *
 *  If a diagram is not specified in the infile, corr_lookup contains an empty
 *  vector for this diagram and the build function immediately returns
 */
// void LapH::Correlators::contract (Quarklines& quarklines,
void LapH::Correlators::contract(OperatorsForMesons const &meson_operator,
                                 RandomVector const &randomvectors,
                                 Perambulator const &perambulators,
                                 OperatorLookup const &operator_lookup,
                                 CorrelatorLookup const &corr_lookup,
                                 QuarklineLookup const &quark_lookup,
                                 std::string const output_path,
                                 std::string const output_filename) {

  // 1. Build all functions which need corrC and free it afterwards.
  build_corrC(randomvectors, perambulators, meson_operator, corr_lookup.corrC);
  build_C2c(corr_lookup.C2c, output_path, output_filename);
  build_C4cD(corr_lookup, output_path, output_filename);
  build_C4cV(corr_lookup, output_path, output_filename);

  // 2. Build all functions which need corr0 and free it afterwards.
  build_corr0(meson_operator, perambulators, corr_lookup.corr0);
  // in C3c, also corr0 is build, since this is much faster
  build_C3c(randomvectors, meson_operator, perambulators, corr_lookup.C3c, output_path,
            output_filename);
  build_C20(corr_lookup.C20, output_path, output_filename);
  build_C40D(corr_lookup, output_path, output_filename);
  build_C40V(corr_lookup, output_path, output_filename);

  // 3. Build all other correlation functions.
  //  build_C1(meson_operator, perambulators, operator_lookup, corr_lookup.C1,
  //                                           quark_lookup, output_path,
  //                                           output_filename);
  build_C4cC(randomvectors, meson_operator, perambulators, corr_lookup.C4cC, output_path,
             output_filename);
  build_C4cB(randomvectors, meson_operator, perambulators, corr_lookup.C4cB, output_path,
             output_filename);
  build_C30(meson_operator, perambulators, corr_lookup.C30, output_path,
            output_filename);
  build_C40C(meson_operator, perambulators, corr_lookup.C40C, output_path,
             output_filename);
  build_C40B(meson_operator, perambulators, corr_lookup.C40B, output_path,
             output_filename);
}
