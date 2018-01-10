#include "Correlators.h"

//#define DILUTION_ITERATOR_PRINT

#include "DilutedFactor.h"
#include "QuarkLineBlock2.h"
#include "StopWatch.h"
#include "dilution-iterator.h"
#include "h5-wrapper.h"
#include "typedefs.h"

template <int n1, int n2, size_t rvecs1, size_t rvecs2>
void multiply(OperatorToFactorMap<n1 + n2, rvecs1 + rvecs2 + 1> &L,
              std::array<size_t, n1 + n2> const &key,
              OperatorToFactorMap<n1, rvecs1> const &factor0,
              OperatorToFactorMap<n2, rvecs2> const &factor1) {
  if (L.count(key) == 0) {
    std::array<size_t, n1> key1;
    std::array<size_t, n2> key2;

    std::copy_n(std::begin(key) + 0, n1, std::begin(key1));
    std::copy_n(std::begin(key) + n1, n2, std::begin(key2));

#pragma omp critical(cout)
    {
//      std::cout << "multiply\t";
//      print(key);

//      MU_DEBUG(factor0.size());
//      for (auto const &elem : factor0) {
//        print(elem.first);
//      }
//      MU_DEBUG(factor1.size());
//      for (auto const &elem : factor1) {
//        print(elem.first);
//      }
    }

    auto const &f0 = factor0.at(key1);
    auto const &f1 = factor1.at(key2);

    L[key] = f0 * f1;
  }
}

namespace {

/*! Creates compound datatype to write complex numbers from complex_t
 *  vectors to HDF5 file
 *
 *  @Returns cmplx_w   HDF5 compound datatype for complex numbers
 */
H5::CompType comp_type_factory_tr() {
  H5::CompType cmplx_w(2 * sizeof(double));
  auto type = H5::PredType::NATIVE_DOUBLE;
  cmplx_w.insertMember("re", HOFFSET(complex_t, re), type);
  cmplx_w.insertMember("im", HOFFSET(complex_t, im), type);

  return cmplx_w;
}

/*! Creates compound datatype to write complex numbers from compcomp_t
 *  vectors to HDF5 file
 *
 *  @Returns cmplx_w   HDF5 compound datatype for structs of four doubles
 */
H5::CompType comp_type_factory_trtr() {
  H5::CompType cmplxcmplx_w(4 * sizeof(double));
  auto type = H5::PredType::NATIVE_DOUBLE;
  cmplxcmplx_w.insertMember("rere", HOFFSET(compcomp_t, rere), type);
  cmplxcmplx_w.insertMember("reim", HOFFSET(compcomp_t, reim), type);
  cmplxcmplx_w.insertMember("imre", HOFFSET(compcomp_t, imre), type);
  cmplxcmplx_w.insertMember("imim", HOFFSET(compcomp_t, imim), type);

  return cmplxcmplx_w;
}

/*! Class to write correlations function to files in hdf5 format
 *
 *  @Warning  Dependency inversion principle is violated: Class depends on the
 *            concrete implementation of hdf5
 */
class WriteHDF5Correlator {
 public:
  WriteHDF5Correlator(const std::string output_path,
                      const std::string diagram,
                      const std::string output_filename,
                      const H5::CompType &_comp_type)
      : comp_type(_comp_type) {
    create_folder_for_hdf5_file(output_path.c_str());

    const H5std_string file_name((output_path + "/" + diagram + output_filename).c_str());
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
   *  @remark The type corr_datatype is always either complex_t or
   *          compcomp_t. The function body is identical for both types
   *          as everything is specified by corr_info. Thus the template
   *          overload
   */
  template <typename corr_datatype>
  void write(const std::vector<corr_datatype> &corr, const CorrInfo &corr_info) {
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
      std::cout << "\tdirectory " << path << " does not exist and will be created";
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

};  // end of class WriteHDF5Correlator

}  // end of anonymous namespace

int get_time_delta(BlockIterator const &slice_pair, int const Lt) {
  return abs((slice_pair.sink() - slice_pair.source() - Lt) % Lt);
}

/******************************************************************************/
/******************************************************************************/

Correlators::Correlators(const size_t Lt,
                         const size_t dilT,
                         const size_t dilE,
                         const size_t nev,
                         const CorrelatorLookup &corr_lookup,
                         OperatorLookup const &operator_lookup,
                         QuarklineLookup const &quark_lookup)
    : Lt(Lt),
      dilT(dilT),
      dilE(dilE),
      nev(nev),
      dil_fac_lookup(
          {quark_lookup.Q0, quark_lookup.Q1, quark_lookup.Q2V, quark_lookup.Q2L}),
      ric_lookup(operator_lookup.ricQ2_lookup) {}

/*!
 *  @deprecated
 */
void Correlators::build_part_trQ1(RandomVector const &randomvectors,
                                  OperatorsForMesons const &meson_operator,
                                  Perambulator const &perambulators,
                                  std::vector<CorrInfo> const &corr_lookup,
                                  std::string const output_path,
                                  std::string const output_filename) {
  if (corr_lookup.empty())
    return;

  StopWatch swatch("C1");

  corr_part_trQ1.resize(boost::extents[corr_lookup.size()][Lt]);

  DilutionScheme const dilution_scheme(Lt, dilT, DilutionType::block);

#pragma omp parallel
  {
    swatch.start();

    QuarkLineBlock2<QuarkLineType::Q1> quarklines(
        randomvectors, perambulators, meson_operator, dilT, dilE, nev, dil_fac_lookup.Q1);

#pragma omp for schedule(dynamic)
    for (int t = 0; t < Lt; ++t) {
      auto const b = dilution_scheme.time_to_block(t);

      for (const auto &c_look : corr_lookup) {
        corr_part_trQ1[c_look.id][t] =
            factor_to_trace(quarklines[{t, b}].at({c_look.lookup[0]}));
      }
    }

    swatch.stop();
  }  // parallel part ends here

  HDF5Handle handle(output_path, "C1", output_filename);

  // normalisation
  for (const auto &c_look : corr_lookup) {
    for (auto &corr_t : corr_part_trQ1[c_look.id]) {
      for (auto &diluted_trace : corr_t) {
        // TODO: Hard Coded atm - Be carefull
        diluted_trace.data /= Lt;
      }
    }

    auto group = handle.create_group(c_look.hdf5_dataset_name);
    std::cout << "Going to write" << std::endl;
    write_heterogenious(group, corr_part_trQ1[c_look.id]);
  }
  swatch.print();
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void Correlators::build_C20V(std::vector<CorrInfo> const &corr_lookup,
                             std::string const output_path,
                             std::string const output_filename) {
  if (corr_lookup.empty())
    return;

  StopWatch swatch("C20V");
  swatch.start();

  // every element of corr_lookup contains the same filename. Wlog choose the
  // first element
  WriteHDF5Correlator filehandle(
      output_path, "C20V", output_filename, comp_type_factory_trtr());

  DilutionScheme const dilution_scheme(Lt, dilT, DilutionType::block);

  for (const auto &c_look : corr_lookup) {
    std::vector<compcomp_t> correlator(Lt, compcomp_t(.0, .0, .0, .0));

    for (auto const block_pair : dilution_scheme) {
      for (auto const slice_pair : block_pair) {
        int const t = get_time_delta(slice_pair, Lt);
        correlator[t] +=
            inner_product(corr_part_trQ1[c_look.lookup[0]][slice_pair.source()],
                          corr_part_trQ1[c_look.lookup[1]][slice_pair.sink()]);
      }
    }

    // normalisation
    for (auto &corr : correlator) {
      corr /= Lt;
    }

    // write data to file
    filehandle.write(correlator, c_look);
  }

  swatch.stop();
  swatch.print();
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void Correlators::build_corr0(RandomVector const &randomvectors,
                              OperatorsForMesons const &meson_operator,
                              Perambulator const &perambulators,
                              std::vector<CorrInfo> const &corr_lookup) {
  if (corr_lookup.empty())
    return;

  StopWatch swatch("corr0");

  corr0.resize(boost::extents[corr_lookup.size()][Lt][Lt]);

  DilutionScheme const dilution_scheme(Lt, dilT, DilutionType::block);

#pragma omp parallel
  {
    swatch.start();

    QuarkLineBlock2<QuarkLineType::Q1> quarklines(
        randomvectors, perambulators, meson_operator, dilT, dilE, nev, dil_fac_lookup.Q1);

#pragma omp for schedule(dynamic)
    for (int b = 0; b < dilution_scheme.size(); ++b) {
      // Notation is that `t1` is the source and `t2` the sink. Both will be
      // done eventually, so this is symmetric.

      auto const block_pair = dilution_scheme[b];

      for (auto const slice_pair : block_pair) {
        for (const auto &c_look : corr_lookup) {
          corr0[c_look.id][slice_pair.source()][slice_pair.sink()] = factor_to_trace(
              quarklines[{slice_pair.source(), slice_pair.sink_block()}].at(
                  {c_look.lookup[0]}),
              quarklines[{slice_pair.sink(), slice_pair.source_block()}].at(
                  {c_look.lookup[1]}));
        }
      }
    }
    swatch.stop();
  }

  swatch.print();
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void Correlators::build_C20(std::vector<CorrInfo> const &corr_lookup,
                            std::string const output_path,
                            std::string const output_filename) {
  if (corr_lookup.empty())
    return;

  StopWatch swatch("C20");
  swatch.start();

  // every element of corr_lookup contains the same filename. Wlog choose the
  // first element
  WriteHDF5Correlator filehandle(
      output_path, "C20", output_filename, comp_type_factory_tr());

  for (const auto &c_look : corr_lookup) {
    std::vector<cmplx> correlator(Lt, cmplx(.0, .0));

    DilutionScheme const dilution_scheme(Lt, dilT, DilutionType::block);
    for (auto const block_pair : dilution_scheme) {
      for (auto const slice_pair : block_pair) {
        int const t = get_time_delta(slice_pair, Lt);
        /*! @todo hidden because range based but this is a loop over random */

        auto const &c = corr0[c_look.lookup[0]][slice_pair.source()][slice_pair.sink()];
        correlator[t] += std::accumulate(std::begin(c), std::end(c), cmplx(0.0, 0.0));
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
void Correlators::build_C40D(std::vector<CorrInfo> const &corr_lookup,
                             std::string const output_path,
                             std::string const output_filename) {
  if (corr_lookup.empty())
    return;

  StopWatch swatch("C40D", 1);
  swatch.start();

  // every element of corr_lookup contains the same filename. Wlog choose the
  // first element
  WriteHDF5Correlator filehandle(
      output_path, "C40D", output_filename, comp_type_factory_trtr());

  DilutionScheme const dilution_scheme(Lt, dilT, DilutionType::block);

  for (const auto &c_look : corr_lookup) {
    std::vector<compcomp_t> correlator(Lt, compcomp_t(.0, .0, .0, .0));

    for (auto const block_pair : dilution_scheme) {
      for (auto const slice_pair : block_pair) {
        int const t = get_time_delta(slice_pair, Lt);

        /*! @todo Write move assignment for compcomp_t and give trtr return
         * parameter */
        correlator[t] += inner_product(
            corr0[c_look.lookup[0]][slice_pair.source()][slice_pair.sink()],
            corr0[c_look.lookup[1]][slice_pair.source()][slice_pair.sink()]);
      }
    }

    // normalisation
    for (auto &corr : correlator) {
      corr /= Lt;
    }

    // write data to file
    filehandle.write(correlator, c_look);
  }

  swatch.stop();
  swatch.print();
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void Correlators::build_C40V(std::vector<CorrInfo> const &corr_lookup,
                             std::string const output_path,
                             std::string const output_filename) {
  if (corr_lookup.empty())
    return;

  StopWatch swatch("C40V");
  swatch.start();

  // every element of corr_lookup contains the same filename. Wlog choose the
  // first element
  WriteHDF5Correlator filehandle(
      output_path, "C40V", output_filename, comp_type_factory_trtr());

  DilutionScheme const dilution_scheme(Lt, dilT, DilutionType::block);

  for (const auto &c_look : corr_lookup) {
    std::vector<compcomp_t> correlator(Lt, compcomp_t(.0, .0, .0, .0));

    for (auto const block_pair : dilution_scheme) {
      for (auto const slice_pair : block_pair) {
        int const t = get_time_delta(slice_pair, Lt);
        correlator[t] += inner_product(
            corr0[c_look.lookup[0]][slice_pair.source()][slice_pair.source()],
            corr0[c_look.lookup[1]][slice_pair.sink()][slice_pair.sink()]);
      }
    }

    // normalisation
    for (auto &corr : correlator) {
      corr /= Lt;
    }

    // write data to file
    filehandle.write(correlator, c_look);
  }

  swatch.stop();
  swatch.print();
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void Correlators::build_corrC(RandomVector const &randomvectors,
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

    QuarkLineBlock2<QuarkLineType::Q0> quarklines_Q0(
        randomvectors, perambulators, meson_operator, dilT, dilE, nev, dil_fac_lookup.Q0);
    QuarkLineBlock2<QuarkLineType::Q2> quarklines_Q2V(randomvectors,
                                                      perambulators,
                                                      meson_operator,
                                                      dilT,
                                                      dilE,
                                                      nev,
                                                      dil_fac_lookup.Q2V);

#pragma omp for schedule(dynamic)
    for (int b = 0; b < dilution_scheme.size(); ++b) {
      auto const block_pair = dilution_scheme[b];

      for (auto const slice_pair : block_pair) {
        for (const auto &c_look : corr_lookup) {
          corrC[c_look.id][slice_pair.source()][slice_pair.sink()] =
              factor_to_trace(quarklines_Q0[{slice_pair.sink()}].at({c_look.lookup[1]}),
                              quarklines_Q2V[{slice_pair.sink_block(),
                                              slice_pair.source(),
                                              slice_pair.sink_block()}]
                                  .at({c_look.lookup[0]}));
        }
      }
    }
    swatch.stop();
  }

  swatch.print();
}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void Correlators::build_C2c(std::vector<CorrInfo> const &corr_lookup,
                            std::string const output_path,
                            std::string const output_filename) {
  if (corr_lookup.empty())
    return;

  StopWatch swatch("C2c");
  swatch.start();

  // every element of corr_lookup contains the same filename. Wlog choose the
  // first element
  WriteHDF5Correlator filehandle(
      output_path, "C2+", output_filename, comp_type_factory_tr());

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
        auto const &c = corrC[c_look.lookup[0]][slice_pair.source()][slice_pair.sink()];
        correlator[t] += std::accumulate(std::begin(c), std::end(c), cmplx(0.0, 0.0));
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
void Correlators::build_C4cD(CorrelatorLookup const &corr_lookup,
                             std::string const output_path,
                             std::string const output_filename) {
  if (corr_lookup.C4cD.empty())
    return;

  StopWatch swatch("C4cD");
  swatch.start();

  // every element of corr_lookup contains the same filename. Wlog choose the
  // first element
  WriteHDF5Correlator filehandle(
      output_path, "C4+D", output_filename, comp_type_factory_trtr());

  DilutionScheme const dilution_scheme(Lt, dilT, DilutionType::block);

  for (const auto &c_look : corr_lookup.C4cD) {
    std::vector<compcomp_t> correlator(Lt, compcomp_t(.0, .0, .0, .0));

    for (auto const block_pair : dilution_scheme) {
      for (auto const slice_pair : block_pair) {
        int const t = get_time_delta(slice_pair, Lt);

        /*! @todo Write move assignment for compcomp_t and give trtr return
         * parameter */
        correlator[t] += inner_product(
            corrC[c_look.lookup[0]][slice_pair.source()][slice_pair.sink()],
            corrC[c_look.lookup[1]][slice_pair.source()][slice_pair.sink()]);
      }
    }

    // normalisation
    for (auto &corr : correlator) {
      corr /= Lt;
    }

    // write data to file
    filehandle.write(correlator, c_look);
  }

  swatch.stop();
  swatch.print();
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void Correlators::build_C4cV(CorrelatorLookup const &corr_lookup,
                             std::string const output_path,
                             std::string const output_filename) {
  if (corr_lookup.C4cV.empty())
    return;

  StopWatch swatch("C4cV");
  swatch.start();

  // every element of corr_lookup contains the same filename. Wlog choose the
  // first element
  WriteHDF5Correlator filehandle(
      output_path, "C4+V", output_filename, comp_type_factory_trtr());

  DilutionScheme const dilution_scheme(Lt, dilT, DilutionType::block);

  for (const auto &c_look : corr_lookup.C4cV) {
    std::vector<compcomp_t> correlator(Lt, compcomp_t(.0, .0, .0, .0));

    for (auto const block_pair : dilution_scheme) {
      for (auto const slice_pair : block_pair) {
        int const t = get_time_delta(slice_pair, Lt);

        /*! @todo Write move assignment for compcomp_t and give trtr return
         * parameter */
        correlator[t] += inner_product(
            corrC[c_look.lookup[0]][slice_pair.source()][slice_pair.source()],
            corrC[c_look.lookup[1]][slice_pair.sink()][slice_pair.sink()]);
      }
    }

    // normalisation
    for (auto &corr : correlator) {
      corr /= Lt;
    }

    // write data to file
    filehandle.write(correlator, c_look);
  }

  swatch.stop();
  swatch.print();
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void Correlators::build_C4cC(RandomVector const &randomvectors,
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
  WriteHDF5Correlator filehandle(
      output_path, "C4+C", output_filename, comp_type_factory_tr());

  std::vector<std::vector<cmplx>> correlator(corr_lookup.size(),
                                             std::vector<cmplx>(Lt, cmplx(.0, .0)));

  DilutionScheme const dilution_scheme(Lt, dilT, DilutionType::block);

  std::vector<std::array<std::array<size_t, 2>, 2>> quantum_num_ids;
  quantum_num_ids.reserve(corr_lookup.size());
  for (const auto &c_look : corr_lookup) {
    quantum_num_ids.push_back(
        {std::array<size_t, 2>{c_look.lookup[3], c_look.lookup[0]},
         std::array<size_t, 2>{c_look.lookup[1], c_look.lookup[2]}});
  }

#pragma omp parallel
  {
    swatch.start();
    std::vector<std::vector<cmplx>> C(corr_lookup.size(),
                                      std::vector<cmplx>(Lt, cmplx(.0, .0)));

    // building the quark line directly frees up a lot of memory
    QuarkLineBlock2<QuarkLineType::Q0> quarkline_Q0(
        randomvectors, perambulators, meson_operator, dilT, dilE, nev, dil_fac_lookup.Q0);

    QuarkLineBlock2<QuarkLineType::Q2> quarkline_Q2V(randomvectors,
                                                     perambulators,
                                                     meson_operator,
                                                     dilT,
                                                     dilE,
                                                     nev,
                                                     dil_fac_lookup.Q2V);

#pragma omp for schedule(dynamic)
    // Perform contraction here
    for (int b = 0; b < dilution_scheme.size(); ++b) {
      auto const block_pair = dilution_scheme[b];
      // Create quarklines for all time combinations in block_pair

      for (auto const slice_pair : block_pair) {
        int const t = get_time_delta(slice_pair, Lt);

        OperatorToFactorMap<2, 1> L1;
        OperatorToFactorMap<2, 1> L2;
        for (const auto &ids : quantum_num_ids) {
          multiply<1, 1, 0, 0>(L1,
                               ids[0],
                               quarkline_Q0[{slice_pair.sink()}],
                               quarkline_Q2V[{slice_pair.sink_block(),
                                              slice_pair.source(),
                                              slice_pair.sink_block()}]);

          multiply<1, 1, 0, 0>(L2,
                               ids[1],
                               quarkline_Q0[{slice_pair.sink()}],
                               quarkline_Q2V[{slice_pair.sink_block(),
                                              slice_pair.source(),
                                              slice_pair.sink_block()}]);
        }

        for (int i = 0; i != quantum_num_ids.size(); ++i) {
          auto const &ids = quantum_num_ids[i];
          C[i][t] += trace(L1[ids[0]], L2[ids[1]]);
        }
      }

      quarkline_Q0.clear();
      quarkline_Q2V.clear();
    }  // loops over time end here
#pragma omp critical
    {
      for (int i = 0; i != quantum_num_ids.size(); ++i) {
        for (size_t t = 0; t < Lt; t++)
          correlator[i][t] += C[i][t];
      }
    }
    swatch.stop();
  }  // parallel part ends here

  // normalisation
  for (int i = 0; i != quantum_num_ids.size(); ++i) {
    for (auto &corr : correlator[i]) {
      corr /= 3 * Lt;
    }
    // write data to file
    filehandle.write(correlator[i], corr_lookup[i]);
  }
  swatch.print();
}

/*****************************************************************************/
/*                                 build_C3c                                 */
/*****************************************************************************/
void Correlators::build_C3c(RandomVector const &randomvectors,
                            OperatorsForMesons const &meson_operator,
                            Perambulator const &perambulators,
                            std::vector<CorrInfo> const &corr_lookup,
                            std::string const output_path,
                            std::string const output_filename) {
  if (corr_lookup.empty())
    return;

  StopWatch swatch("C3c");

  // every element of corr_lookup contains the same filename. Wlog choose the
  // first element
  WriteHDF5Correlator filehandle(
      output_path, "C3+", output_filename, comp_type_factory_tr());

  std::vector<std::vector<cmplx>> correlator(corr_lookup.size(),
                                             std::vector<cmplx>(Lt, cmplx(.0, .0)));

  DilutionScheme const dilution_scheme(Lt, dilT, DilutionType::block);

  std::vector<std::tuple<std::array<size_t, 2>, std::array<size_t, 1>>> quantum_num_ids;
  quantum_num_ids.reserve(corr_lookup.size());
  for (const auto &c_look : corr_lookup) {
    quantum_num_ids.push_back(
        make_tuple(std::array<size_t, 2>{c_look.lookup[2], c_look.lookup[0]},
                   std::array<size_t, 1>{c_look.lookup[1]}));
  }

#pragma omp parallel
  {
    swatch.start();
    // std::vector<std::vector<cmplx>> C(corr_lookup.size(), std::vector<cmplx>(Lt,
    // cmplx(.0, .0)));
    std::vector<std::vector<cmplx>> C(corr_lookup.size(),
                                      std::vector<cmplx>(Lt, cmplx(.0, .0)));

    QuarkLineBlock2<QuarkLineType::Q0> quarklines_Q0(
        randomvectors, perambulators, meson_operator, dilT, dilE, nev, dil_fac_lookup.Q0);

    QuarkLineBlock2<QuarkLineType::Q1> quarklines_Q1(
        randomvectors, perambulators, meson_operator, dilT, dilE, nev, dil_fac_lookup.Q1);

    QuarkLineBlock2<QuarkLineType::Q2> quarklines_Q2(randomvectors,
                                                     perambulators,
                                                     meson_operator,
                                                     dilT,
                                                     dilE,
                                                     nev,
                                                     dil_fac_lookup.Q2L);

#pragma omp for schedule(dynamic)
    for (int b = 0; b < dilution_scheme.size(); ++b) {
      auto const block_pair = dilution_scheme[b];

      std::cout << block_pair << std::endl;

      for (auto const slice_pair : block_pair) {
        int const t = get_time_delta(slice_pair, Lt);

        OperatorToFactorMap<2, 1> L1;
        for (const auto &ids : quantum_num_ids) {
          multiply<1, 1, 0, 0>(L1,
                               std::get<0>(ids),
                               quarklines_Q0[{slice_pair.source()}],
                               quarklines_Q2[{slice_pair.source_block(),
                                              slice_pair.source(),
                                              slice_pair.sink_block()}]);
        }

        for (int i = 0; i != quantum_num_ids.size(); ++i) {
          auto const &ids = quantum_num_ids[i];
          C[i][t] +=
              trace(L1[std::get<0>(ids)],
                    quarklines_Q1[{slice_pair.sink(), slice_pair.source_block()}].at(
                        std::get<1>(ids)));
        }
      }

      quarklines_Q0.clear();
      quarklines_Q1.clear();
      quarklines_Q2.clear();
    }
#pragma omp critical
    {
      for (int i = 0; i != quantum_num_ids.size(); ++i) {
        for (size_t t = 0; t < Lt; t++)
          correlator[i][t] += C[i][t];
      }
    }
    swatch.stop();
  }  // parallel part ends here

  // normalisation
  for (int i = 0; i != quantum_num_ids.size(); ++i) {
    for (auto &corr : correlator[i]) {
      corr /= 2* Lt;
    }
    // write data to file
    filehandle.write(correlator[i], corr_lookup[i]);
  }
  swatch.print();
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void Correlators::build_C4cB(RandomVector const &randomvectors,
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

  std::vector<std::vector<cmplx>> correlator(corr_lookup.size(),
                                             std::vector<cmplx>(Lt, cmplx(.0, .0)));

  DilutionScheme const dilution_scheme(Lt, dilT, DilutionType::block);

  std::vector<std::array<std::array<size_t, 2>, 2>> quantum_num_ids;
  quantum_num_ids.reserve(corr_lookup.size());
  for (const auto &c_look : corr_lookup) {
    quantum_num_ids.push_back(
        {std::array<size_t, 2>{c_look.lookup[3], c_look.lookup[0]},
         std::array<size_t, 2>{c_look.lookup[1], c_look.lookup[2]}});
  }

// This is necessary to ensure the correct summation of the correlation function
#pragma omp parallel
  {
    swatch.start();
    std::vector<std::vector<cmplx>> C(corr_lookup.size(),
                                      std::vector<cmplx>(Lt, cmplx(.0, .0)));

    // building the quark line directly frees up a lot of memory
    QuarkLineBlock2<QuarkLineType::Q0> quarkline_Q0(
        randomvectors, perambulators, meson_operator, dilT, dilE, nev, dil_fac_lookup.Q0);

    QuarkLineBlock2<QuarkLineType::Q2> quarkline_Q2L(randomvectors,
                                                     perambulators,
                                                     meson_operator,
                                                     dilT,
                                                     dilE,
                                                     nev,
                                                     dil_fac_lookup.Q2L);

#pragma omp for schedule(dynamic)
    // Perform contraction here
    for (int b = 0; b < dilution_scheme.size(); ++b) {
      auto const block_pair = dilution_scheme[b];
      // Create quarklines for all time combinations in block_pair

      for (auto const slice_pair : block_pair) {
        int const t = get_time_delta(slice_pair, Lt);

        OperatorToFactorMap<2, 1> L1;
        OperatorToFactorMap<2, 1> L2;
        for (const auto &ids : quantum_num_ids) {
          multiply<1, 1, 0, 0>(L1,
                               ids[0],
                               quarkline_Q0[{slice_pair.source()}],
                               quarkline_Q2L[{slice_pair.source_block(),
                                              slice_pair.source(),
                                              slice_pair.sink_block()}]);

          multiply<1, 1, 0, 0>(L2,
                               ids[1],
                               quarkline_Q0[{slice_pair.sink()}],
                               quarkline_Q2L[{slice_pair.sink_block(),
                                              slice_pair.sink(),
                                              slice_pair.source_block()}]);
        }

        for (int i = 0; i != quantum_num_ids.size(); ++i) {
          auto const &ids = quantum_num_ids[i];
          C[i][t] += trace(L1[ids[0]], L2[ids[1]]);
        }
      }

      quarkline_Q0.clear();
      quarkline_Q2L.clear();
    }  // loops over time end here
#pragma omp critical
    {
      for (int i = 0; i != quantum_num_ids.size(); ++i) {
        for (size_t t = 0; t < Lt; t++)
          correlator[i][t] += C[i][t];
      }
    }
    swatch.stop();
  }  // parallel part ends here

  // normalisation
  for (int i = 0; i != quantum_num_ids.size(); ++i) {
    for (auto &corr : correlator[i]) {
      corr /= 3 * Lt;
    }
    // write data to file
    filehandle.write(correlator[i], corr_lookup[i]);
  }
  swatch.print();
}

void Correlators::build_C30(RandomVector const &randomvectors,
                            OperatorsForMesons const &meson_operator,
                            Perambulator const &perambulators,
                            std::vector<CorrInfo> const &corr_lookup,
                            std::string const output_path,
                            std::string const output_filename) {
  if (corr_lookup.empty())
    return;

  StopWatch swatch("C30");

  // every element of corr_lookup contains the same filename. Wlog choose the
  // first element
  WriteHDF5Correlator filehandle(
      output_path, "C30", output_filename, comp_type_factory_tr());

  std::vector<std::vector<cmplx>> correlator(corr_lookup.size(),
                                             std::vector<cmplx>(Lt, cmplx(.0, .0)));

  DilutionScheme const dilution_scheme(Lt, dilT, DilutionType::block);

  std::vector<std::tuple<std::array<size_t, 2>, std::array<size_t, 1>>> quantum_num_ids;
  quantum_num_ids.reserve(corr_lookup.size());
  for (const auto &c_look : corr_lookup) {
    quantum_num_ids.push_back(
        make_tuple(std::array<size_t, 2>{c_look.lookup[2], c_look.lookup[0]},
                   std::array<size_t, 1>{c_look.lookup[1]}));
  }

#pragma omp parallel
  {
    swatch.start();
    std::vector<std::vector<cmplx>> C(corr_lookup.size(),
                                      std::vector<cmplx>(Lt, cmplx(.0, .0)));
    // building the quark line directly frees up a lot of memory
    QuarkLineBlock2<QuarkLineType::Q1> quarklines(
        randomvectors, perambulators, meson_operator, dilT, dilE, nev, dil_fac_lookup.Q1);

#pragma omp for schedule(dynamic)
    for (int b = 0; b < dilution_scheme.size(); ++b) {
      auto const block_pair = dilution_scheme[b];

      std::cout << block_pair << std::endl;

      for (auto const slice_pair : block_pair) {
        int const t = get_time_delta(slice_pair, Lt);

        OperatorToFactorMap<2, 1> L1;
        for (const auto &ids : quantum_num_ids) {
          multiply<1, 1, 0, 0>(
              L1,
              std::get<0>(ids),
              quarklines[{slice_pair.source(), slice_pair.source_block()}],
              quarklines[{slice_pair.source(), slice_pair.sink_block()}]);
        }

        for (int i = 0; i != quantum_num_ids.size(); ++i) {
          auto const &ids = quantum_num_ids[i];
          C[i][t] += trace(L1[std::get<0>(ids)],
                           quarklines[{slice_pair.sink(), slice_pair.source_block()}].at(
                               std::get<1>(ids)));
        }
      }

      quarklines.clear();
    }
#pragma omp critical
    {
      for (int i = 0; i != quantum_num_ids.size(); ++i) {
        for (size_t t = 0; t < Lt; t++)
          correlator[i][t] += C[i][t];
      }
    }
    swatch.stop();
  }  // parallel part ends here

  // normalisation
  for (int i = 0; i != quantum_num_ids.size(); ++i) {
    for (auto &corr : correlator[i]) {
      corr /= Lt;
    }
    // write data to file
    filehandle.write(correlator[i], corr_lookup[i]);
  }
  swatch.print();
}

void Correlators::build_C40C(RandomVector const &randomvectors,
                             OperatorsForMesons const &meson_operator,
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
  WriteHDF5Correlator filehandle(
      output_path, "C40C", output_filename, comp_type_factory_tr());

  std::vector<std::vector<cmplx>> correlator(corr_lookup.size(),
                                             std::vector<cmplx>(Lt, cmplx(.0, .0)));

  DilutionScheme const dilution_scheme(Lt, dilT, DilutionType::block);

  std::vector<std::array<std::array<size_t, 2>, 2>> quantum_num_ids;
  quantum_num_ids.reserve(corr_lookup.size());
  for (const auto &c_look : corr_lookup) {
    quantum_num_ids.push_back(
        {std::array<size_t, 2>{c_look.lookup[3], c_look.lookup[0]},
         std::array<size_t, 2>{c_look.lookup[1], c_look.lookup[2]}});
  }

#pragma omp parallel
  {
    swatch.start();
    // std::vector<std::vector<cmplx>> C(corr_lookup.size(), std::vector<cmplx>(Lt,
    // cmplx(.0, .0)));
    std::vector<std::vector<cmplx>> C(corr_lookup.size(),
                                      std::vector<cmplx>(Lt, cmplx(.0, .0)));
    // building the quark line directly frees up a lot of memory
    QuarkLineBlock2<QuarkLineType::Q1> quarklines(
        randomvectors, perambulators, meson_operator, dilT, dilE, nev, dil_fac_lookup.Q1);

#pragma omp for schedule(dynamic)
    for (int b = 0; b < dilution_scheme.size(); ++b) {
      auto const block_pair = dilution_scheme[b];

      std::cout << block_pair << std::endl;

      for (auto const slice_pair : block_pair) {
        int const t = get_time_delta(slice_pair, Lt);

        OperatorToFactorMap<2, 1> L1;
        OperatorToFactorMap<2, 1> L2;
        for (const auto &ids : quantum_num_ids) {
          multiply<1, 1, 0, 0>(
              L1,
              ids[0],
              quarklines[{slice_pair.sink(), slice_pair.source_block()}],
              quarklines[{slice_pair.source(), slice_pair.sink_block()}]);

          multiply<1, 1, 0, 0>(
              L2,
              ids[1],
              quarklines[{slice_pair.sink(), slice_pair.source_block()}],
              quarklines[{slice_pair.source(), slice_pair.sink_block()}]);
        }

        for (int i = 0; i != quantum_num_ids.size(); ++i) {
          auto const &ids = quantum_num_ids[i];
          C[i][t] += trace(L1[ids[0]], L2[ids[1]]);
        }
      }

      quarklines.clear();
    }
#pragma omp critical
    {
      for (int i = 0; i != quantum_num_ids.size(); ++i) {
        for (size_t t = 0; t < Lt; t++)
          correlator[i][t] += C[i][t];
      }
    }
    swatch.stop();
  }  // parallel part ends here

  // normalisation
  for (int i = 0; i != quantum_num_ids.size(); ++i) {
    for (auto &corr : correlator[i]) {
      corr /= Lt;
    }
    // write data to file
    filehandle.write(correlator[i], corr_lookup[i]);
  }
  swatch.print();
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
/*!
 *  @note Not optimised.
 */
void Correlators::build_C40B(RandomVector const &randomvectors,
                             OperatorsForMesons const &meson_operator,
                             Perambulator const &perambulators,
                             std::vector<CorrInfo> const &corr_lookup,
                             std::string const output_path,
                             std::string const output_filename) {
  if (corr_lookup.empty())
    return;

  StopWatch swatch("C40B");

  // every element of corr_lookup contains the same filename. Wlog choose the
  // first element
  WriteHDF5Correlator filehandle(
      output_path, "C40B", output_filename, comp_type_factory_tr());

  std::vector<std::vector<cmplx>> correlator(corr_lookup.size(),
                                             std::vector<cmplx>(Lt, cmplx(.0, .0)));

  DilutionScheme const dilution_scheme(Lt, dilT, DilutionType::block);

  std::vector<std::array<std::array<size_t, 2>, 2>> quantum_num_ids;
  quantum_num_ids.reserve(corr_lookup.size());
  for (const auto &c_look : corr_lookup) {
    quantum_num_ids.push_back(
        {std::array<size_t, 2>{c_look.lookup[3], c_look.lookup[0]},
         std::array<size_t, 2>{c_look.lookup[1], c_look.lookup[2]}});
  }

#pragma omp parallel
  {
    swatch.start();
    // std::vector<std::vector<cmplx>> C(corr_lookup.size(), std::vector<cmplx>(Lt,
    // cmplx(.0, .0)));
    std::vector<std::vector<cmplx>> C(corr_lookup.size(),
                                      std::vector<cmplx>(Lt, cmplx(.0, .0)));
    // building the quark line directly frees up a lot of memory
    QuarkLineBlock2<QuarkLineType::Q1> quarklines(
        randomvectors, perambulators, meson_operator, dilT, dilE, nev, dil_fac_lookup.Q1);

#pragma omp for schedule(dynamic)
    for (int b = 0; b < dilution_scheme.size(); ++b) {
      auto const block_pair = dilution_scheme[b];

      std::cout << block_pair << std::endl;

      for (auto const slice_pair : block_pair) {
        int const t = get_time_delta(slice_pair, Lt);

        OperatorToFactorMap<2, 1> L1;
        OperatorToFactorMap<2, 1> L2;
        for (const auto &ids : quantum_num_ids) {
          multiply<1, 1, 0, 0>(
              L1,
              ids[0],
              quarklines[{slice_pair.source(), slice_pair.source_block()}],
              quarklines[{slice_pair.source(), slice_pair.sink_block()}]);

          multiply<1, 1, 0, 0>(
              L2,
              ids[1],
              quarklines[{slice_pair.sink(), slice_pair.sink_block()}],
              quarklines[{slice_pair.sink(), slice_pair.source_block()}]);
        }

        for (int i = 0; i != quantum_num_ids.size(); ++i) {
          auto const &ids = quantum_num_ids[i];
          C[i][t] += trace(L1[ids[0]], L2[ids[1]]);
        }
      }

      quarklines.clear();
    }
#pragma omp critical
    {
      for (int i = 0; i != quantum_num_ids.size(); ++i) {
        for (size_t t = 0; t < Lt; t++)
          correlator[i][t] += C[i][t];
      }
    }
    swatch.stop();
  }  // parallel part ends here

  // normalisation
  for (int i = 0; i != quantum_num_ids.size(); ++i) {
    for (auto &corr : correlator[i]) {
      corr /= Lt;
    }
    // write data to file
    filehandle.write(correlator[i], corr_lookup[i]);
  }
  swatch.print();
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
// void Correlators::contract (Quarklines& quarklines,
void Correlators::contract(OperatorsForMesons const &meson_operator,
                           RandomVector const &randomvectors,
                           Perambulator const &perambulators,
                           OperatorLookup const &operator_lookup,
                           CorrelatorLookup const &corr_lookup,
                           QuarklineLookup const &quark_lookup,
                           std::string const output_path,
                           std::string const output_filename) {
  build_part_trQ1(randomvectors,
                  meson_operator,
                  perambulators,
                  corr_lookup.C1,
                  output_path,
                  output_filename);
  build_C20V(corr_lookup.C20V, output_path, output_filename);

  // 1. Build all functions which need corrC and free it afterwards.
  build_corrC(randomvectors, perambulators, meson_operator, corr_lookup.corrC);
  build_C2c(corr_lookup.C2c, output_path, output_filename);
  build_C4cD(corr_lookup, output_path, output_filename);
  build_C4cV(corr_lookup, output_path, output_filename);

  // 2. Build all functions which need corr0 and free it afterwards.
  build_corr0(randomvectors, meson_operator, perambulators, corr_lookup.corr0);
  // in C3c, also corr0 is build, since this is much faster
  build_C3c(randomvectors,
            meson_operator,
            perambulators,
            corr_lookup.C3c,
            output_path,
            output_filename);
  build_C20(corr_lookup.C20, output_path, output_filename);
  build_C40D(corr_lookup.C40D, output_path, output_filename);
  build_C40V(corr_lookup.C40V, output_path, output_filename);

  // 3. Build all other correlation functions.
  build_C4cC(randomvectors,
             meson_operator,
             perambulators,
             corr_lookup.C4cC,
             output_path,
             output_filename);
  build_C4cB(randomvectors,
             meson_operator,
             perambulators,
             corr_lookup.C4cB,
             output_path,
             output_filename);
  build_C30(randomvectors,
            meson_operator,
            perambulators,
            corr_lookup.C30,
            output_path,
            output_filename);
  build_C40C(randomvectors,
             meson_operator,
             perambulators,
             corr_lookup.C40C,
             output_path,
             output_filename);
  build_C40B(randomvectors,
             meson_operator,
             perambulators,
             corr_lookup.C40B,
             output_path,
             output_filename);
}
