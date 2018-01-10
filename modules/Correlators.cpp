#include "Correlators.h"

//#define DILUTION_ITERATOR_PRINT

#include "DilutedFactor.h"
#include "QuarkLineBlock2.h"
#include "StopWatch.h"
#include "dilution-iterator.h"
#include "h5-wrapper.h"
#include "typedefs.h"
#include "Diagram.h"

int get_time_delta(BlockIterator const &slice_pair, int const Lt) {
  return abs((slice_pair.sink() - slice_pair.source() - Lt) % Lt);
}

namespace {

template <typename Numeric>
H5::CompType comp_type_factory();

/*! Creates compound datatype to write complex numbers from complex_t
 *  vectors to HDF5 file
 *
 *  @Returns cmplx_w   HDF5 compound datatype for complex numbers
 */
template <>
H5::CompType comp_type_factory<cmplx>() {
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
template <>
H5::CompType comp_type_factory<compcomp_t>() {
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
   *  @see  H5::CompType comp_type_factory<compcomp_t>()
   */
  H5::CompType comp_type;

};  // end of class WriteHDF5Correlator

}  // end of anonymous namespace

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
          {quark_lookup.Q0, quark_lookup.Q1, quark_lookup.Q2V, quark_lookup.Q2L}) {}

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
template <typename Numeric>
void build_diagram(typename DiagramTraits<Numeric>::Diagram &diagram,
                   RandomVector const &randomvectors,
                   OperatorsForMesons const &meson_operator,
                   Perambulator const &perambulators,
                   std::string const output_path,
                   std::string const output_filename,
                   const size_t Lt,
                   const size_t dilT,
                   const size_t dilE,
                   const size_t nev,
                   DilutedFactorLookup const &dil_fac_lookup,
                   DilutedTraceCollection<2> &corr0,
                   DilutedTraceCollection<2> &corrC,
                   DilutedTraceCollection2<1> &corr_part_trQ1) {
  if (diagram.corr_lookup().empty())
    return;

  StopWatch swatch(diagram.name());

  WriteHDF5Correlator filehandle(
      output_path, diagram.name(), output_filename, comp_type_factory<Numeric>());

  std::vector<std::vector<Numeric>> correlator(diagram.corr_lookup().size(),
                                             std::vector<Numeric>(Lt, Numeric{}));

  DilutionScheme const dilution_scheme(Lt, dilT, DilutionType::block);

// This is necessary to ensure the correct summation of the correlation function
#pragma omp parallel
  {
    swatch.start();
    std::vector<std::vector<Numeric>> C(
        Lt, std::vector<Numeric>(diagram.corr_lookup().size(), Numeric{}));

    // building the quark line directly frees up a lot of memory
    QuarkLineBlock2<QuarkLineType::Q0> quarkline_Q0(
        randomvectors, perambulators, meson_operator, dilT, dilE, nev, dil_fac_lookup.Q0);

    QuarkLineBlock2<QuarkLineType::Q1> quarkline_Q1(
        randomvectors, perambulators, meson_operator, dilT, dilE, nev, dil_fac_lookup.Q1);

    QuarkLineBlock2<QuarkLineType::Q2> quarkline_Q2L(randomvectors,
                                                     perambulators,
                                                     meson_operator,
                                                     dilT,
                                                     dilE,
                                                     nev,
                                                     dil_fac_lookup.Q2L);

    QuarkLineBlock2<QuarkLineType::Q2> quarkline_Q2V(randomvectors,
                                                     perambulators,
                                                     meson_operator,
                                                     dilT,
                                                     dilE,
                                                     nev,
                                                     dil_fac_lookup.Q2V);

    QuarkLineBlockCollection part_collection = {quarkline_Q0,
                                                quarkline_Q1,
                                                quarkline_Q2L,
                                                quarkline_Q2V,
                                                corr0,
                                                corrC,
                                                corr_part_trQ1};

#pragma omp for schedule(dynamic)
    // Perform contraction here
    for (int b = 0; b < dilution_scheme.size(); ++b) {
      auto const block_pair = dilution_scheme[b];
      for (auto const slice_pair : block_pair) {
        int const t = get_time_delta(slice_pair, Lt);

        diagram.contract(C[t], slice_pair, part_collection);
      }

      quarkline_Q0.clear();
      quarkline_Q1.clear();
      quarkline_Q2L.clear();
      quarkline_Q2V.clear();

    }  // loops over time end here
#pragma omp critical
    {
      for (int i = 0; i != correlator.size(); ++i) {
        for (size_t t = 0; t < Lt; t++) {
          correlator[i][t] += C[t][i];
        }
      }
    }
    swatch.stop();
  }  // parallel part ends here

  // normalisation
  for (int i = 0; i != correlator.size(); ++i) {
    for (auto &corr : correlator[i]) {
      corr /= Lt;
    }
    // write data to file
    filehandle.write(correlator[i], diagram.corr_lookup()[i]);
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

  // 1. Build all functions which need corrC and free it afterwards.
  build_corrC(randomvectors, perambulators, meson_operator, corr_lookup.corrC);

  // 2. Build all functions which need corr0 and free it afterwards.
  build_corr0(randomvectors, meson_operator, perambulators, corr_lookup.corr0);

  // 3. Build all other correlation functions.

  // XXX If we had C++14, we could do `make_unique`.
  std::vector<std::unique_ptr<Diagram>> diagrams;

  diagrams.emplace_back(new C2c(corr_lookup.C2c));
  diagrams.emplace_back(new C20(corr_lookup.C20));
  diagrams.emplace_back(new C20V(corr_lookup.C20V));

  diagrams.emplace_back(new C3c(corr_lookup.C3c));
  diagrams.emplace_back(new C30(corr_lookup.C30));

  diagrams.emplace_back(new C4cB(corr_lookup.C4cB));
  diagrams.emplace_back(new C40B(corr_lookup.C40B));
  diagrams.emplace_back(new C4cC(corr_lookup.C4cC));
  diagrams.emplace_back(new C40C(corr_lookup.C40C));

  diagrams.emplace_back(new C4cD(corr_lookup.C4cD));
  diagrams.emplace_back(new C40D(corr_lookup.C40D));
  diagrams.emplace_back(new C4cV(corr_lookup.C4cV));
  diagrams.emplace_back(new C40V(corr_lookup.C40V));

  for (auto &diagram : diagrams) {
    diagram->build(randomvectors,
                   meson_operator,
                   perambulators,
                   output_path,
                   output_filename,
                   Lt,
                   dilT,
                   dilE,
                   nev,
                   dil_fac_lookup,
                   corr0,
                   corrC,
                   corr_part_trQ1);
  }
}

template void build_diagram<cmplx>(typename DiagramTraits<cmplx>::Diagram &diagram,
                                   RandomVector const &randomvectors,
                                   OperatorsForMesons const &meson_operator,
                                   Perambulator const &perambulators,
                                   std::string const output_path,
                                   std::string const output_filename,
                                   const size_t Lt,
                                   const size_t dilT,
                                   const size_t dilE,
                                   const size_t nev,
                                   DilutedFactorLookup const &dil_fac_lookup,
                                   DilutedTraceCollection<2> &corr0,
                                   DilutedTraceCollection<2> &corrC,
                                   DilutedTraceCollection2<1> &corr_part_trQ1);

template void build_diagram<compcomp_t>(
    typename DiagramTraits<compcomp_t>::Diagram &diagram,
    RandomVector const &randomvectors,
    OperatorsForMesons const &meson_operator,
    Perambulator const &perambulators,
    std::string const output_path,
    std::string const output_filename,
    const size_t Lt,
    const size_t dilT,
    const size_t dilE,
    const size_t nev,
    DilutedFactorLookup const &dil_fac_lookup,
    DilutedTraceCollection<2> &corr0,
    DilutedTraceCollection<2> &corrC,
    DilutedTraceCollection2<1> &corr_part_trQ1);
