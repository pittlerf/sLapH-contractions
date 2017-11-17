#include "Correlators.h"

//#define DILUTION_ITERATOR_PRINT

#include "QuarkLineBlock.h"
#include "dilution-iterator.h"
#include "StopWatch.h"

namespace
{

/*! Creates compound datatype to write complex numbers from LapH::complex_t 
 *  vectors to HDF5 file
 *
 *  @Returns cmplx_w   HDF5 compound datatype for complex numbers
 */
H5::CompType comp_type_factory_tr(){
  H5::CompType cmplx_w(2*sizeof(double));
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
H5::CompType comp_type_factory_trtr(){
  H5::CompType cmplxcmplx_w(4*sizeof(double));
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
class WriteHDF5Correlator{

public:
  WriteHDF5Correlator(const std::string output_path, 
                      const std::string diagram,
                      const std::string output_filename, 
                      const H5::CompType& _comp_type) : comp_type(_comp_type){

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
  template < typename corr_datatype >
  void write(const std::vector<corr_datatype>& corr, const CorrInfo& corr_info){

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
    try
    {
      dset = file.createDataSet(dataset_name, comp_type, dspace);
      dset.write(&corr[0], comp_type);
  
      // closing of dset is delegated to destructor ~DataSet
  
    } 
    catch(H5::Exception& e){
      e.printError();
    }
  }

private:

  /*! Checks whether output path exists and if not creates it 
   * 
   *  @param[in] path Path where hdf5 file shall be written
   */
  void create_folder_for_hdf5_file(const char* path){
    if(access(path, 0 ) != 0) {
      std::cout << "\tdirectory " << path 
                << " does not exist and will be created";
      boost::filesystem::path dir(path);
      if(boost::filesystem::create_directories(dir))
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
  void open_or_create_hdf5_file(const H5std_string& name){

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

/******************************************************************************/
/******************************************************************************/

/*!
 *  @deprecated
 */
void LapH::Correlators::build_C1(const Quarklines& quarklines,
                    const std::vector<CorrInfo>& corr_lookup,
                    const QuarklineLookup& quark_lookup,
                    const std::vector<RandomIndexCombinationsQ2>& ric_lookup) {

  if(corr_lookup.empty())
    return;

  StopWatch swatch("C1", 1);
  swatch.start();

  // every element of corr_lookup contains the same filename. Wlog choose the 
  // first element
//  WriteHDF5Correlator filehandle(corr_lookup[0], comp_type_factory_tr() );

  for(const auto& c_look : corr_lookup){
    const auto& ric = ric_lookup[quark_lookup.Q1[c_look.lookup[0]].
                                                     id_ric_lookup].rnd_vec_ids;
    for(const auto& id : ric){

      std::vector<cmplx> correlator(Lt, cmplx(.0,.0));
      for(size_t t = 0; t < Lt; t++){
        correlator[t] +=
         quarklines.return_Q1(t, t/dilT, c_look.lookup[0], &id-&ric[0]).trace();
      }

      // write data to file
//      filehandle.write(correlator, c_look);
    }

  }
  swatch.stop();
  swatch.print();
}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------




// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void LapH::Correlators::build_corr0(const OperatorsForMesons& meson_operator,
                                    const Perambulator& perambulators,
                                    const std::vector<CorrInfo>& corr_lookup,
                                    const QuarklineLookup& quark_lookup,
                                    const OperatorLookup& operator_lookup) {
  if (corr_lookup.size() == 0) return;

  StopWatch swatch("corr0");

  corr0.resize(boost::extents[corr_lookup.size()][Lt][Lt]);

  DilutionScheme dilution_scheme(Lt, dilT, DilutionType::block);

#pragma omp parallel
  {
    swatch.start();

    QuarkLineBlock<QuarkLineType::Q1> quarklines_local(
        dilT, dilE, nev, quark_lookup.Q1, operator_lookup.ricQ2_lookup);

#pragma omp for schedule(dynamic)
    for (int b = 0; b < dilution_scheme.size(); ++b) {
      // Notation is that `t1` is the source and `t2` the sink. Both will be
      // done eventually, so this is symmetric.

      auto const block_pair = dilution_scheme[b];

      for (auto const slice_pair_one_sink : block_pair.one_sink_slice()) {
        quarklines_local.build_Q1_one_t(perambulators,
                                        meson_operator,
                                        slice_pair_one_sink.source(),
                                        block_pair.sink(),
                                        quark_lookup.Q1,
                                        operator_lookup.ricQ2_lookup);
      }

      for (auto const slice_pair : block_pair) {
        for (const auto& c_look : corr_lookup) {
          const auto& ric0 =
              operator_lookup
                  .ricQ2_lookup[quark_lookup.Q1[c_look.lookup[0]].id_ric_lookup]
                  .rnd_vec_ids;
          const auto& ric1 =
              operator_lookup
                  .ricQ2_lookup[quark_lookup.Q1[c_look.lookup[1]].id_ric_lookup]
                  .rnd_vec_ids;
          if (ric0.size() != ric1.size()) {
            std::cout << "rnd combinations are not the same in build_corr0" << std::endl;
            exit(1);
          }
          corr0[c_look.id][slice_pair.source()][slice_pair.sink()].resize(ric0.size());
          for (auto& corr : corr0[c_look.id][slice_pair.source()][slice_pair.sink()])
            corr = cmplx(0.0, 0.0);
          for (const auto& rnd : ric0) {
            const auto id = &rnd - &ric0[0];
            const auto it1 = std::find_if(
                ric1.begin(), ric1.end(), [&](std::pair<size_t, size_t> pair) {
                  return (pair == std::make_pair(rnd.second, rnd.first));
                });
            if (it1 == ric1.end()) {
              std::cout << "something wrong with random vectors in build_corr0"
                        << std::endl;
              exit(1);
            }

            /*! @todo How do I properly get the block indices for sink? */
            corr0[c_look.id][slice_pair.source()][slice_pair.sink()][id] +=
                (quarklines_local(
                     slice_pair.source(), slice_pair.sink_block(), c_look.lookup[0], id) *
                 quarklines_local(slice_pair.sink(),
                                  slice_pair.source_block(),
                                  c_look.lookup[1],
                                  it1 - ric1.begin()))
                    .trace();
          }
        }
      }
    }
    swatch.stop();
  }

  swatch.print();
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void LapH::Correlators::build_C20(const std::vector<CorrInfo>& corr_lookup, 
                                  const std::string output_path,
                                  const std::string output_filename) {

  if(corr_lookup.empty())
    return;

  StopWatch swatch("C20");
  swatch.start();


  // every element of corr_lookup contains the same filename. Wlog choose the 
  // first element
  WriteHDF5Correlator filehandle(output_path, "C20", output_filename, comp_type_factory_tr() );

  for(const auto& c_look : corr_lookup){
    std::vector<cmplx> correlator(Lt, cmplx(.0,.0));
    for(int t1 = 0; t1 < Lt; t1++){
    for(int t2 = 0; t2 < Lt; t2++){
      int t = abs((t2 - t1 - (int)Lt) % (int)Lt);
      for(const auto& corr : corr0[c_look.lookup[0]][t1][t2])
        correlator[t] += corr;
    }}
    // normalisation
    for(auto& corr : correlator){
      corr /= Lt*corr0[c_look.lookup[0]][0][0].size();
    }

    // write data to file
    filehandle.write(correlator, c_look);
  }

  swatch.stop();
  swatch.print();
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void LapH::Correlators::build_C40D(const OperatorLookup &operator_lookup,
                                   const CorrelatorLookup &corr_lookup,
                                   const QuarklineLookup &quark_lookup,
                                   const std::string output_path,
                                   const std::string output_filename) {
  if (corr_lookup.C40D.empty()) return;

  StopWatch swatch("C40D", 1);
  swatch.start();

  // every element of corr_lookup contains the same filename. Wlog choose the
  // first element
  WriteHDF5Correlator filehandle(
      output_path, "C40D", output_filename, comp_type_factory_trtr());

  for (const auto &c_look : corr_lookup.C40D) {
    std::vector<LapH::compcomp_t> correlator(Lt, LapH::compcomp_t(.0, .0, .0, .0));
    const size_t id0 = corr_lookup.corr0[c_look.lookup[0]].lookup[0];
    const size_t id1 = corr_lookup.corr0[c_look.lookup[1]].lookup[0];
    const auto &ric0 =
        operator_lookup.ricQ2_lookup[quark_lookup.Q1[id0].id_ric_lookup].rnd_vec_ids;
    const auto &ric1 =
        operator_lookup.ricQ2_lookup[quark_lookup.Q1[id1].id_ric_lookup].rnd_vec_ids;
    size_t norm = 0;

    DilutionScheme dilution_scheme(Lt, dilT, DilutionType::block);
    for (auto const block_pair : dilution_scheme) {
      for (auto const slice_pair : block_pair) {
        int const t =
            abs((slice_pair.sink() - slice_pair.source() - static_cast<int>(Lt)) %
                static_cast<int>(Lt));

        for (const auto &rnd0 : ric0) {
          for (const auto &rnd1 : ric1) {
            if ((rnd0.first != rnd1.first) && (rnd0.first != rnd1.second) &&
                (rnd0.second != rnd1.first) && (rnd0.second != rnd1.second)) {
              auto const &factor1 =
                  corr0[c_look.lookup[0]][slice_pair.source()][slice_pair.sink()].at(
                      &rnd0 - &ric0[0]);
              auto const &factor2 =
                  corr0[c_look.lookup[1]][slice_pair.source()][slice_pair.sink()].at(
                      &rnd1 - &ric1[0]);

              correlator[t].rere += factor1.real() * factor2.real();
              correlator[t].reim += factor1.real() * factor2.imag();
              correlator[t].imre += factor1.imag() * factor2.real();
              correlator[t].imim += factor1.imag() * factor2.imag();
              norm++;
            }
          }
        }
      }
    }

    // normalisation
    for (auto &corr1 : correlator) {
      corr1.rere /= norm / Lt;
      corr1.reim /= norm / Lt;
      corr1.imre /= norm / Lt;
      corr1.imim /= norm / Lt;
    }

    // write data to file
    filehandle.write(correlator, c_look);
  }

  swatch.stop();
  swatch.print();
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void LapH::Correlators::build_C40V(const OperatorLookup& operator_lookup, 
                                   const CorrelatorLookup& corr_lookup,
                                   const QuarklineLookup& quark_lookup, 
                                   const std::string output_path,
                                   const std::string output_filename) {

  if(corr_lookup.C40V.empty())
    return;

  StopWatch swatch("C40V");

  // every element of corr_lookup contains the same filename. Wlog choose the 
  // first element
  WriteHDF5Correlator filehandle(output_path, "C40V", output_filename, comp_type_factory_trtr() );

  swatch.start();

  for(const auto& c_look : corr_lookup.C40V){
    std::vector<LapH::compcomp_t> correlator(Lt, LapH::compcomp_t(.0,.0,.0,.0));
    const size_t id0 = corr_lookup.corr0[c_look.lookup[0]].lookup[0];
    const size_t id1 = corr_lookup.corr0[c_look.lookup[1]].lookup[0];
    const auto& ric0 = operator_lookup.ricQ2_lookup[quark_lookup.Q1[id0].
                                                     id_ric_lookup].rnd_vec_ids;
    const auto& ric1 = operator_lookup.ricQ2_lookup[quark_lookup.Q1[id1].
                                                     id_ric_lookup].rnd_vec_ids;
    size_t norm = 0;
    DilutionScheme dilution_scheme(Lt, dilT, DilutionType::block);
    for (auto const block_pair : dilution_scheme) {
      for (auto const slice_pair : block_pair) {
        int const t =
            abs((slice_pair.sink() - slice_pair.source() - static_cast<int>(Lt)) %
                static_cast<int>(Lt));
        for (const auto &rnd0 : ric0) {
          for (const auto &rnd1 : ric1) {
            if ((rnd0.first != rnd1.first) && (rnd0.first != rnd1.second) &&
                (rnd0.second != rnd1.first) && (rnd0.second != rnd1.second)) {
              auto const &factor1 =
                  corr0[c_look.lookup[0]][slice_pair.source()][slice_pair.source()].at(
                      &rnd0 - &ric0[0]);
              auto const &factor2 =
                  corr0[c_look.lookup[1]][slice_pair.sink()][slice_pair.sink()].at(
                      &rnd1 - &ric1[0]);

              correlator[t].rere += factor1.real() * factor2.real();
              correlator[t].reim += factor1.real() * factor2.imag();
              correlator[t].imre += factor1.imag() * factor2.real();
              correlator[t].imim += factor1.imag() * factor2.imag();
              norm++;
            }
          }
        }
      }
    }
    // normalisation
    for (auto &corr1 : correlator) {
      corr1.rere /= norm / Lt;
      corr1.reim /= norm / Lt;
      corr1.imre /= norm / Lt;
      corr1.imim /= norm / Lt;
    }
    // write data to file
    filehandle.write(correlator, c_look);
  }

  swatch.stop();
  swatch.print();
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void LapH::Correlators::build_corrC(const Perambulator &perambulators,
                                    const OperatorsForMesons &meson_operator,
                                    const OperatorLookup &operator_lookup,
                                    const std::vector<CorrInfo> &corr_lookup,
                                    const QuarklineLookup &quark_lookup) {
  if (corr_lookup.size() == 0) return;

  StopWatch swatch("corrC");

  corrC.resize(boost::extents[corr_lookup.size()][Lt][Lt]);

#pragma omp parallel
  {
    swatch.start();

    // building the quark line directly frees up a lot of memory
    QuarkLine_one_t<QuarkLineType::Q2V> quarklines(
        dilT, dilE, nev, quark_lookup.Q2V, operator_lookup.ricQ2_lookup);
#pragma omp for schedule(dynamic)
    for (int t1_i = 0; t1_i < Lt / dilT; t1_i++) {
      for (int t2_i = t1_i; t2_i < Lt / dilT; t2_i++) {
        quarklines.build(perambulators,
                         meson_operator,
                         t1_i,
                         t2_i,
                         quark_lookup.Q2V,
                         operator_lookup.ricQ2_lookup);
        for (int dir = 0; dir < 2; dir++) {
          if ((t1_i == t2_i) && (dir == 1)) continue;

          int t1_min, t2_min, t1_max, t2_max;
          if (dir == 0) {
            t1_min = dilT * t1_i;
            t1_max = dilT * (t1_i + 1);
            t2_min = dilT * t2_i;
            t2_max = dilT * (t2_i + 1);
          } else {
            t1_min = dilT * t2_i;
            t1_max = dilT * (t2_i + 1);
            t2_min = dilT * t1_i;
            t2_max = dilT * (t1_i + 1);
          }

          for (int t1 = t1_min; t1 < t1_max; t1++) {
            for (int t2 = t2_min; t2 < t2_max; t2++) {
              // quarkline indices
              int id_Q2L_1;

              if (t1_i == t2_i) {
                if (t1_min != 0)
                  id_Q2L_1 = t1 % t1_min;
                else {
                  id_Q2L_1 = t1;
                }
              } else {
                if (t1_min != 0)
                  id_Q2L_1 = (dir)*dilT + t1 % t1_min;
                else
                  id_Q2L_1 = ((dir)*dilT + t1);
              }

              // building correlator -------------------------------------------------
              for (const auto &c_look : corr_lookup) {
                const auto &ric0 =
                    operator_lookup
                        .ricQ2_lookup[quark_lookup.Q2V[c_look.lookup[0]].id_ric_lookup]
                        .rnd_vec_ids;
                const auto &ric1 =
                    operator_lookup
                        .ricQ2_lookup[  // just for checking
                            operator_lookup.rvdaggervr_lookuptable[c_look.lookup[1]]
                                .id_ricQ_lookup]
                        .rnd_vec_ids;
                if (ric0.size() != ric1.size()) {
                  std::cout << "rnd combinations are not the same in build_corrC"
                            << std::endl;
                  exit(0);
                }
                corrC[c_look.id][t1][t2].resize(ric0.size());
                for (auto &corr : corrC[c_look.id][t1][t2]) corr = cmplx(0.0, 0.0);
                for (const auto &rnd : ric0) {
                  const auto id = &rnd - &ric0[0];
                  for (size_t block = 0; block < 4; block++) {
                    const auto gamma_index =
                        quarklines.return_gamma_row(c_look.gamma[0], block);
                    corrC[c_look.id][t1][t2][id] +=
                        quarklines.return_gamma_val(c_look.gamma[0], block) *
                        (quarklines.return_Ql(id_Q2L_1, c_look.lookup[0], id)
                             .block(block * dilE, gamma_index * dilE, dilE, dilE) *
                         meson_operator.return_rvdaggervr(c_look.lookup[1], t2, id)
                             .block(gamma_index * dilE, block * dilE, dilE, dilE))
                            .trace();
                  }
                }
              }
            }
          }  // t1, t2 end here
        }    // dir (directions) end here
      }
    }  // block times end here
    swatch.stop();
  }  // omp parall ends here

  swatch.print();
}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void LapH::Correlators::build_C2c(const std::vector<CorrInfo>& corr_lookup, 
                                  const std::string output_path,
                                  const std::string output_filename) {
  if(corr_lookup.empty())
    return;

  StopWatch swatch("C2c");

  // every element of corr_lookup contains the same filename. Wlog choose the 
  // first element
  WriteHDF5Correlator filehandle(output_path, "C2+", output_filename, comp_type_factory_tr() );

  swatch.start();

  for(const auto& c_look : corr_lookup){
    std::vector<cmplx> correlator(Lt, cmplx(.0,.0));
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
      for(int t1 = 0; t1 < Lt; t1++){
      for(int t2 = 0; t2 < Lt; t2++){
        int t = abs((t2 - t1 - (int)Lt) % (int)Lt);
        for(const auto& corr : corrC[c_look.lookup[0]][t1][t2]){
          correlator[t] += corr;
        }
      }}
      // normalisation
      for(auto& corr : correlator){
        corr /= Lt*corrC[c_look.lookup[0]][0][0].size();
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
void LapH::Correlators::build_C4cD(const OperatorLookup& operator_lookup, 
                                   const CorrelatorLookup& corr_lookup,
                                   const QuarklineLookup& quark_lookup, 
                                   const std::string output_path,
                                   const std::string output_filename) {

  if(corr_lookup.C4cD.empty())
    return;

  StopWatch swatch("C4cD");
  swatch.start();

  // every element of corr_lookup contains the same filename. Wlog choose the 
  // first element
  WriteHDF5Correlator filehandle(output_path, "C4+D", output_filename, comp_type_factory_trtr() );

  for(const auto& c_look : corr_lookup.C4cD){
    std::vector<LapH::compcomp_t> correlator(Lt, LapH::compcomp_t(.0,.0,.0,.0));
    const size_t id0 = corr_lookup.corrC[c_look.lookup[0]].lookup[0];
    const size_t id1 = corr_lookup.corrC[c_look.lookup[1]].lookup[0];
    const auto& ric0 = operator_lookup.ricQ2_lookup[quark_lookup.Q2V[id0].
                                                     id_ric_lookup].rnd_vec_ids;
    const auto& ric1 = operator_lookup.ricQ2_lookup[quark_lookup.Q2V[id1].
                                                     id_ric_lookup].rnd_vec_ids;

    size_t norm = 0;
    for(int t1 = 0; t1 < Lt; t1++){
    for(int t2 = 0; t2 < Lt; t2++){
      int t = abs((t2 - t1 - (int)Lt) % (int)Lt);

      for(const auto& rnd0 : ric0)
      for(const auto& rnd1 : ric1)

      if((rnd0.first != rnd1.first) && (rnd0.first != rnd1.second) &&
         (rnd0.second != rnd1.first) && (rnd0.second != rnd1.second)){

        correlator[t].rere += 
                   corrC[c_look.lookup[0]][t1][t2].at(&rnd0 - &ric0[0]).real() *
                   corrC[c_look.lookup[1]][t1][t2].at(&rnd1 - &ric1[0]).real();
        correlator[t].reim += 
                   corrC[c_look.lookup[0]][t1][t2].at(&rnd0 - &ric0[0]).real() *
                   corrC[c_look.lookup[1]][t1][t2].at(&rnd1 - &ric1[0]).imag();
        correlator[t].imre += 
                   corrC[c_look.lookup[0]][t1][t2].at(&rnd0 - &ric0[0]).imag() *
                   corrC[c_look.lookup[1]][t1][t2].at(&rnd1 - &ric1[0]).real();
        correlator[t].imim += 
                   corrC[c_look.lookup[0]][t1][t2].at(&rnd0 - &ric0[0]).imag() *
                   corrC[c_look.lookup[1]][t1][t2].at(&rnd1 - &ric1[0]).imag();
        norm++;
      }
    }}
    // normalisation
    for(auto& corr1 : correlator){
      corr1.rere /= norm/Lt;
      corr1.reim /= norm/Lt;
      corr1.imre /= norm/Lt;
      corr1.imim /= norm/Lt;
    }
    // write data to file
    filehandle.write(correlator, c_look);
  }

  swatch.stop();
  swatch.print();
}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void LapH::Correlators::build_C4cV(const OperatorLookup& operator_lookup, 
                                   const CorrelatorLookup& corr_lookup,
                                   const QuarklineLookup& quark_lookup, 
                                   const std::string output_path,
                                   const std::string output_filename) {
  if(corr_lookup.C4cV.empty())
    return;

  StopWatch swatch("C4cV");
  swatch.start();

  // every element of corr_lookup contains the same filename. Wlog choose the 
  // first element
  WriteHDF5Correlator filehandle(output_path, "C4+V", output_filename, comp_type_factory_trtr() );

  for(const auto& c_look : corr_lookup.C4cV){
    std::vector<LapH::compcomp_t> correlator(Lt, LapH::compcomp_t(.0,.0,.0,.0));
    const size_t id0 = corr_lookup.corrC[c_look.lookup[0]].lookup[0];
    const size_t id1 = corr_lookup.corrC[c_look.lookup[1]].lookup[0];
    const auto& ric0 = operator_lookup.ricQ2_lookup[quark_lookup.Q2V[id0].
                                                     id_ric_lookup].rnd_vec_ids;
    const auto& ric1 = operator_lookup.ricQ2_lookup[quark_lookup.Q2V[id1].
                                                     id_ric_lookup].rnd_vec_ids;

    size_t norm = 0;
    for(int t1 = 0; t1 < Lt; t1++){
    for(int t2 = 0; t2 < Lt; t2++){
      int t = abs((t2 - t1 - (int)Lt) % (int)Lt);

      for(const auto& rnd0 : ric0)
      for(const auto& rnd1 : ric1)

      if((rnd0.first != rnd1.first) && (rnd0.first != rnd1.second) &&
         (rnd0.second != rnd1.first) && (rnd0.second != rnd1.second)){

        correlator[t].rere += 
                   corrC[c_look.lookup[0]][t1][t1].at(&rnd0 - &ric0[0]).real() *
                   corrC[c_look.lookup[1]][t2][t2].at(&rnd1 - &ric1[0]).real();
        correlator[t].reim += 
                   corrC[c_look.lookup[0]][t1][t1].at(&rnd0 - &ric0[0]).real() *
                   corrC[c_look.lookup[1]][t2][t2].at(&rnd1 - &ric1[0]).imag();
        correlator[t].imre += 
                   corrC[c_look.lookup[0]][t1][t1].at(&rnd0 - &ric0[0]).imag() *
                   corrC[c_look.lookup[1]][t2][t2].at(&rnd1 - &ric1[0]).real();
        correlator[t].imim += 
                   corrC[c_look.lookup[0]][t1][t1].at(&rnd0 - &ric0[0]).imag() *
                   corrC[c_look.lookup[1]][t2][t2].at(&rnd1 - &ric1[0]).imag();
        norm++;
      }
    }}
    // normalisation
    for(auto& corr1 : correlator){
      corr1.rere /= norm/Lt;
      corr1.reim /= norm/Lt;
      corr1.imre /= norm/Lt;
      corr1.imim /= norm/Lt;
    }
    // write data to file
    filehandle.write(correlator, c_look);
  }

  swatch.stop();
  swatch.print();
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void LapH::Correlators::build_C4cC(const OperatorsForMesons& meson_operator,
                                   const Perambulator& perambulators,
                                   const OperatorLookup& operator_lookup,
                                   const std::vector<CorrInfo>& corr_lookup,
                                   const QuarklineLookup& quark_lookup, 
                                   const std::string output_path,
                                   const std::string output_filename) {
  if(corr_lookup.empty())
    return;

  StopWatch swatch("C4cC");

  // every element of corr_lookup contains the same filename. Wlog choose the
  // first element
  WriteHDF5Correlator filehandle(output_path, "C4+C", output_filename, comp_type_factory_tr() );

  std::vector<vec> correlator(corr_lookup.size(), vec(Lt, cmplx(.0,.0)));

// This is necessary to ensure the correct summation of the correlation function
#pragma omp parallel
{
  swatch.start();
  std::vector<vec> C(corr_lookup.size(), vec(Lt, cmplx(.0,.0)));
  // building the quark line directly frees up a lot of memory
  QuarkLine_one_t<QuarkLineType::Q2V> quarklines(dilT, dilE, nev, quark_lookup.Q2V, 
                        operator_lookup.ricQ2_lookup);
  // creating memory arrays M1, M2 for intermediate storage of Quarklines ------
  std::vector<std::vector<Eigen::MatrixXcd> > M1, M2;
  std::vector<std::array<size_t, 3> > M1_look;
  std::vector<std::array<size_t, 3> > M2_look;
  size_t M1_counter = 0;
  size_t M2_counter = 0;
  for(const auto& c_look : corr_lookup){
    const auto& ric0 = operator_lookup.ricQ2_lookup[
                  quark_lookup.Q2V[c_look.lookup[0]].id_ric_lookup].rnd_vec_ids;
    const auto& ric1 = operator_lookup.ricQ2_lookup[//needed only for checking
               operator_lookup.rvdaggervr_lookuptable[c_look.lookup[1]].
                                                    id_ricQ_lookup].rnd_vec_ids;
    const auto& ric2 = operator_lookup.ricQ2_lookup[
                  quark_lookup.Q2V[c_look.lookup[2]].id_ric_lookup].rnd_vec_ids;
    const auto& ric3 = operator_lookup.ricQ2_lookup[//needed only for checking
               operator_lookup.rvdaggervr_lookuptable[c_look.lookup[3]].
                                                    id_ricQ_lookup].rnd_vec_ids;
    if(ric0.size() != ric1.size() || ric0.size() != ric2.size() || 
       ric0.size() != ric3.size()){
      std::cout << "rnd combinations are not the same in C4+C" << std::endl;
    }

    size_t norm = 0;

    // creating memory for M1 -------------------------------------------------
    const size_t id0 = c_look.lookup[0];
    const size_t id1 = c_look.lookup[1];
    auto it1 = std::find_if(M1_look.begin(), M1_look.end(), 
                           [&id0, &id1](std::array<size_t, 3> check){
                             return (id0 == check[1] && id1 == check[2]);
                           });
    if(!(it1 != M1_look.end())){
      M1.emplace_back(std::vector<Eigen::MatrixXcd>());
      for(const auto& rnd0 : ric0)
      for(const auto& rnd1 : ric1)
      if(rnd0.first != rnd1.first && rnd0.second == rnd1.second)
        M1[M1_counter].emplace_back(Eigen::MatrixXcd::Zero(4*dilE, 4*dilE));
      M1_look.emplace_back(std::array<size_t, 3>({{M1_counter, id0, id1}}));
      M1_counter++;
    }
    // creating memeory for M2 -------------------------------------------------
    const size_t id2 = c_look.lookup[2];
    const size_t id3 = c_look.lookup[3];
    auto it2 = std::find_if(M2_look.begin(), M2_look.end(), 
                            [&id2, &id3](std::array<size_t, 3> check){
                              return (id2 == check[1] && id3 == check[2]);
                            });
    if(!(it2 != M2_look.end())){
      M2.emplace_back(std::vector<Eigen::MatrixXcd>());
      for(const auto& rnd2 : ric2)
      for(const auto& rnd3 : ric3)
      if(rnd2.first != rnd3.first && rnd2.second == rnd3.second)
        M2[M2_counter].emplace_back(Eigen::MatrixXcd::Zero(4*dilE, 4*dilE));
      M2_look.emplace_back(std::array<size_t, 3>({{M2_counter, id2, id3}}));
      M2_counter++;
    }
  }// first run over lookuptable ends here - memory and new lookuptable 
   // are generated ------------------------------------------------------------

  #pragma omp for schedule(dynamic)
  for(int t1_i = 0; t1_i < Lt/dilT; t1_i++){
  for(int t2_i = t1_i; t2_i < Lt/dilT; t2_i++){
    // creating quarklines
    quarklines.build(perambulators, meson_operator, t1_i, t2_i,
                               quark_lookup.Q2V, operator_lookup.ricQ2_lookup);

  for(int bla = 0; bla < 2; bla++){

  if((t1_i == t2_i) && (bla == 1))
    continue;

  int t1_min, t2_min, t1_max, t2_max;
  if(bla == 0){
    t1_min = dilT*t1_i;
    t1_max = dilT*(t1_i+1);
    t2_min = dilT*t2_i;
    t2_max = dilT*(t2_i+1);
  }
  else{
    t1_min = dilT*t2_i;
    t1_max = dilT*(t2_i+1);
    t2_min = dilT*t1_i;
    t2_max = dilT*(t1_i+1);
  }

  for(int t1 = t1_min; t1 < t1_max; t1++){
  for(int t2 = t2_min; t2 < t2_max; t2++){
    int t = abs((t2 - t1 - (int)Lt) % (int)Lt);

    // quarkline indices
    int id_Q2V_1, id_Q2V_2;

//    if(t1_i == t2_i){                                                                  
//        if(t1_min != 0)                                                                  
//          id_Q2V_1 = t1%t1_min;                                                       
//        else{                                                                         
//          id_Q2V_1 = t1;                                                                 
//        }                                                                                
//        if(t1_min != 0)                                                                  
//          id_Q2V_2 = t1%t1_min;                                                       
//        else{                                                                            
//          id_Q2V_2 = t1;                                                                 
//        }                                                                         
//    }                                                                           
//    else{                                                                              
    if(t1_min != 0)                                                               
      id_Q2V_1 = (bla)*dilT + t1%t1_min;                                         
    else                                                                  
      id_Q2V_1 = (bla)*dilT + t1;
    if(t1_min != 0)                                                               
      id_Q2V_2 = (bla)*dilT + t1%t1_min;                           
    else                                                                     
      id_Q2V_2 = (bla)*dilT + t1;               
//    }                                                  

    // build M1 ----------------------------------------------------------------
    for(const auto& look : M1_look){
      const auto& ric0 = operator_lookup.ricQ2_lookup[
                  quark_lookup.Q2V[look[1]].id_ric_lookup].rnd_vec_ids;
      const auto& ric1 = operator_lookup.ricQ2_lookup[//needed only for checking
                 operator_lookup.rvdaggervr_lookuptable[look[2]].
                                                    id_ricQ_lookup].rnd_vec_ids;
      size_t M1_rnd_counter = 0;
      for(const auto& rnd0 : ric0){
      for(const auto& rnd1 : ric1){
      if(rnd0.first != rnd1.first && rnd0.second == rnd1.second){
        const size_t idr0 = &rnd0 - &ric0[0];
        const size_t idr1 = &rnd1 - &ric1[0];
        for(size_t col = 0; col < 4; col++){
          const cmplx value = quarklines.return_gamma_val(5, col); // TODO: gamma hardcoded
          const size_t gamma_index = quarklines.return_gamma_row(
                                                          5, col); // TODO: gamma hardcoded
//          const cmplx value = quarklines.return_gamma_val(c_look.gamma[0], col);
//          const size_t gamma_index = quarklines.return_gamma_row(
//                                                          c_look.gamma[0], col);
          M1[look[0]][M1_rnd_counter].block(0, col*dilE, 4*dilE, dilE) = value *
            quarklines.return_Ql(id_Q2V_1, look[1], idr0).
                               block(0, gamma_index*dilE, 4*dilE, dilE) *
            meson_operator.return_rvdaggervr(look[2], t2, idr1).
                                block(gamma_index*dilE, col*dilE, dilE, dilE);
        }
        M1_rnd_counter++;
      }}}
    }
    // build M2 ----------------------------------------------------------------
    for(const auto& look : M2_look){
      const auto& ric2 = operator_lookup.ricQ2_lookup[
                  quark_lookup.Q2V[look[1]].id_ric_lookup].rnd_vec_ids;
      const auto& ric3 = operator_lookup.ricQ2_lookup[//needed only for checking
                 operator_lookup.rvdaggervr_lookuptable[look[2]].
                                                    id_ricQ_lookup].rnd_vec_ids;
      size_t M2_rnd_counter = 0;
      for(const auto& rnd2 : ric2){
      for(const auto& rnd3 : ric3){
        if(rnd2.first != rnd3.first && rnd2.second == rnd3.second){
        const size_t idr2 = &rnd2 - &ric2[0];
        const size_t idr3 = &rnd3 - &ric3[0];
        for(size_t col = 0; col < 4; col++){
          const cmplx value = quarklines.return_gamma_val(5, col); // TODO: gamma hardcoded
          const size_t gamma_index = quarklines.return_gamma_row(
                                                          5, col); // TODO: gamma hardcoded
//          const cmplx value = quarklines.return_gamma_val(c_look.gamma[1], col);
//          const size_t gamma_index = quarklines.return_gamma_row(
//                                                            c_look.gamma[1], col);

          M2[look[0]][M2_rnd_counter].block(0, col*dilE, 4*dilE, dilE) = value * 
            quarklines.return_Ql(id_Q2V_2, look[1], idr2).
                               block(0, gamma_index*dilE, 4*dilE, dilE) *
            meson_operator.return_rvdaggervr(look[2], t2, idr3).
                                block(gamma_index*dilE, col*dilE, dilE, dilE);

        }
        M2_rnd_counter++;
      }}}
    }
    // Final summation for correlator ------------------------------------------
    Eigen::MatrixXcd M3 = Eigen::MatrixXcd::Zero(4*dilE, 4*dilE);
    for(const auto& c_look : corr_lookup){
      const auto& ric0 = operator_lookup.ricQ2_lookup[
                  quark_lookup.Q2V[c_look.lookup[0]].id_ric_lookup].rnd_vec_ids;
      const auto& ric1 = operator_lookup.ricQ2_lookup[//needed only for checking
                 operator_lookup.rvdaggervr_lookuptable[c_look.lookup[1]].
                                                    id_ricQ_lookup].rnd_vec_ids;
      const auto& ric2 = operator_lookup.ricQ2_lookup[
                  quark_lookup.Q2V[c_look.lookup[2]].id_ric_lookup].rnd_vec_ids;
      const auto& ric3 = operator_lookup.ricQ2_lookup[//needed only for checking
                 operator_lookup.rvdaggervr_lookuptable[c_look.lookup[3]].
                                                    id_ricQ_lookup].rnd_vec_ids;
      const size_t id0 = c_look.lookup[0];
      const size_t id1 = c_look.lookup[1];
      auto it1 = std::find_if(M1_look.begin(), M1_look.end(), 
                             [&id0, &id1](std::array<size_t, 3> check){
                               return (id0 == check[1] && id1 == check[2]);
                             });
      const size_t id2 = c_look.lookup[2];
      const size_t id3 = c_look.lookup[3];
      auto it2 = std::find_if(M2_look.begin(), M2_look.end(), 
                              [&id2, &id3](std::array<size_t, 3> check){
                                return (id2 == check[1] && id3 == check[2]);
                              });
      size_t M1_rnd_counter = 0;
      for(const auto& rnd0 : ric0){
      for(const auto& rnd1 : ric1){

//      if(rnd0.first == rnd1.first && rnd0.second != rnd1.second){
      if(rnd0.first != rnd1.first && rnd0.second == rnd1.second){
        M3.setZero(4 * dilE, 4 * dilE); // setting matrix values to zero
        size_t M2_rnd_counter = 0;
        for(const auto& rnd2 : ric2){
        for(const auto& rnd3 : ric3){
        if(rnd2.first != rnd3.first && rnd2.second == rnd3.second){
          if(rnd1.first == rnd2.first && rnd0.first == rnd3.first &&
           rnd0.second != rnd3.second){
//        if(rnd2.first == rnd3.first && rnd2.second != rnd3.second){
//          if(rnd0.second == rnd3.second && rnd1.second == rnd2.second &&
//             rnd0.first != rnd2.first){
            M3 += M2[(*it2)[0]][M2_rnd_counter];
          }
          M2_rnd_counter++;
        }}}
        C[c_look.id][t] += (M1[(*it1)[0]][M1_rnd_counter++]*M3).trace();
      }}}
    } // loop over operators ends here
  }}}}} // loops over time end here
  #pragma omp critical
  {
    for(const auto& c_look : corr_lookup)
      for(size_t t = 0; t < Lt; t++)
        correlator[c_look.id][t] += C[c_look.id][t];
  }

  swatch.stop();
}// parallel part ends here

  // normalisation
  for(const auto& c_look : corr_lookup){
    for(auto& corr : correlator[c_look.id]){
      corr /= (6*5*4*3)*Lt; // TODO: Hard Coded atm - Be carefull
    }
    // write data to file
    filehandle.write(correlator[c_look.id], c_look);
  }

  swatch.print();
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
//void LapH::Correlators::build_C4cC(const Quarklines& quarklines, 
//                                   const OperatorsForMesons& meson_operator,
//                                   const OperatorLookup& operator_lookup,
//                                   const std::vector<CorrInfo>& corr_lookup,
//                                   const QuarklineLookup& quark_lookup) {
//
//  std::cout << "\tcomputing C4cC:";
//  clock_t time = clock();
//
//  for(const auto& c_look : corr_lookup){
//    std::vector<cmplx> correlator(Lt, cmplx(.0,.0));
//    const auto& ric0 = operator_lookup.ricQ2_lookup[
//                  quark_lookup.Q2V[c_look.lookup[0]].id_ric_lookup].rnd_vec_ids;
//    const auto& ric1 = operator_lookup.ricQ2_lookup[//needed only for checking
//               operator_lookup.rvdaggervr_lookuptable[c_look.lookup[1]].
//                                                    id_ricQ_lookup].rnd_vec_ids;
//    const auto& ric2 = operator_lookup.ricQ2_lookup[
//                  quark_lookup.Q2V[c_look.lookup[2]].id_ric_lookup].rnd_vec_ids;
//    const auto& ric3 = operator_lookup.ricQ2_lookup[//needed only for checking
//               operator_lookup.rvdaggervr_lookuptable[c_look.lookup[3]].
//                                                    id_ricQ_lookup].rnd_vec_ids;
//    if(ric0.size() != ric1.size() || ric0.size() != ric2.size() || 
//       ric0.size() != ric3.size()){
//      std::cout << "rnd combinations are not the same in C4+C" 
//                << std::endl;
//      }
//
//    size_t norm = 0;
//// This is necessary to ensure the correct summation of the correlation function
//#pragma omp parallel reduction(+:norm)
//{
//    std::vector<cmplx> C(Lt, cmplx(.0,.0));
//    #pragma omp for schedule(dynamic) 
//    for(int t1 = 0; t1 < Lt; t1++){
//    for(int t2 = 0; t2 < Lt; t2++){
//      int t = abs((t2 - t1 - (int)Lt) % (int)Lt);
//      for(const auto& rnd0 : ric0){
//      for(const auto& rnd1 : ric1){
//      if(rnd0.first != rnd1.first && rnd0.second == rnd1.second){
//        const size_t id0 = &rnd0 - &ric0[0];
//        const size_t id1 = &rnd1 - &ric1[0];
//        Eigen::MatrixXcd M1 = Eigen::MatrixXcd::Zero(4 * dilE, 4 * dilE);
//        for(size_t col = 0; col < 4; col++){
//
//          const cmplx value = quarklines.return_gamma_val(c_look.gamma[0], col);
//          const size_t gamma_index = quarklines.return_gamma_row(
//                                                          c_look.gamma[0], col);
//          M1.block(0, col*dilE, 4*dilE, dilE) = value *
//            quarklines.return_Q2V(t1, t2/dilT, c_look.lookup[0], id0).
//                               block(0, gamma_index*dilE, 4*dilE, dilE) *
//            meson_operator.return_rvdaggervr(c_look.lookup[1], t2, id1).
//                                block(gamma_index*dilE, col*dilE, dilE, dilE);
//        }
//        for(const auto& rnd2 : ric2){
//        for(const auto& rnd3 : ric3){
//        if(rnd2.first != rnd3.first && rnd2.second == rnd3.second && 
//           rnd1.first == rnd2.first && rnd0.first == rnd3.first &&
//           rnd0.second != rnd3.second){
//          const size_t id2 = &rnd2 - &ric2[0];
//          const size_t id3 = &rnd3 - &ric3[0];
//          Eigen::MatrixXcd M2 = Eigen::MatrixXcd::Zero(4 * dilE, 4 * dilE);
//          for(size_t col = 0; col < 4; col++){
//            const cmplx value = quarklines.return_gamma_val(c_look.gamma[1], col);
//            const size_t gamma_index = quarklines.return_gamma_row(
//                                                            c_look.gamma[1], col);
//            M2.block(0, col*dilE, 4*dilE, dilE) = value *
//              quarklines.return_Q2V(t1, t2/dilT, c_look.lookup[2], id2).
//                                 block(0, gamma_index*dilE, 4*dilE, dilE) *
//              meson_operator.return_rvdaggervr(c_look.lookup[3], t2, id3).
//                                  block(gamma_index*dilE, col*dilE, dilE, dilE);
//          }
//          C[t] += (M1*M2).trace();
//          norm++;
//        }}}
//      }}}
//    }}
//    #pragma omp critical
//    {
//      for(size_t t = 0; t < Lt; t++)
//        correlator[t] += C[t];
//    }
//}// parallel part ends here
//
//    // normalisation
//    for(auto& corr : correlator)
//      corr /= norm/Lt;
//    // write data to file
//    write_correlators(correlator, c_look);
//  }
//
//  time = clock() - time;
//  std::cout << "\t\t\tSUCCESS - " << ((float) time) / CLOCKS_PER_SEC 
//            << " seconds" << std::endl;
//}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void LapH::Correlators::build_C3c(const OperatorsForMesons& meson_operator,
                                  const Perambulator& perambulators,
                                  const OperatorLookup& operator_lookup,
                                  const std::vector<CorrInfo>& corr_lookup,
                                  const QuarklineLookup& quark_lookup, 
                                  const std::string output_path,
                                  const std::string output_filename) {
  if(corr_lookup.empty())
    return;

  // every element of corr_lookup contains the same filename. Wlog choose the 
  // first element
  WriteHDF5Correlator filehandle(output_path, "C3+", output_filename, comp_type_factory_tr() );

  StopWatch swatch("C3c");


  std::vector<vec> correlator(corr_lookup.size(), vec(Lt, cmplx(.0,.0)));

// This is necessary to ensure the correct summation of the correlation function
#pragma omp parallel
{
  swatch.start();
  std::vector<vec> C(corr_lookup.size(), vec(Lt, cmplx(.0,.0)));
  // building the quark line directly frees up a lot of memory
  QuarkLine_one_t<QuarkLineType::Q2L> quarklines_Q2L(dilT, dilE, nev, quark_lookup.Q2L, 
                        operator_lookup.ricQ2_lookup);
  QuarkLine_one_t<QuarkLineType::Q1> quarklines_Q1(dilT, dilE, nev, quark_lookup.Q1, 
                        operator_lookup.ricQ2_lookup);

  // creating memory arrays M1, M2 for intermediate storage of Quarklines ------
  std::vector<std::vector<Eigen::MatrixXcd> > M1, M2;
  std::vector<std::array<size_t, 3> > M1_look;
  size_t M1_counter = 0;
  for(const auto& c_look : corr_lookup){

    const auto& ric0 = operator_lookup.ricQ2_lookup[
                quark_lookup.Q2L[c_look.lookup[0]].id_ric_lookup].rnd_vec_ids;
    const auto& ric1 = operator_lookup.ricQ2_lookup[
                quark_lookup.Q1[c_look.lookup[1]].id_ric_lookup].rnd_vec_ids;
    const auto& ric2 = operator_lookup.ricQ2_lookup[
                operator_lookup.rvdaggervr_lookuptable[c_look.lookup[2]].
                                                  id_ricQ_lookup].rnd_vec_ids;
    if(ric0.size() != ric1.size() || ric0.size() != ric2.size()){
      std::cout << "rnd combinations are not the same in build_C3+" 
                << std::endl;
      exit(0);
    }

    // creating memeory for M1 -------------------------------------------------
    const size_t id0 = c_look.lookup[0];
    const size_t id2 = c_look.lookup[2];
    auto it1 = std::find_if(M1_look.begin(), M1_look.end(), 
                           [&id0, &id2](std::array<size_t, 3> check){
                             return (id0 == check[1] && id2 == check[2]);
                           });

    if(!(it1 != M1_look.end())){
      M1.emplace_back(std::vector<Eigen::MatrixXcd>());
      for(const auto& rnd0 : ric0)
      for(const auto& rnd2 : ric2)
      if(rnd0.first == rnd2.first && rnd0.second != rnd2.second)
        M1[M1_counter].emplace_back(Eigen::MatrixXcd::Zero(4*dilE, 4*dilE));
      M1_look.emplace_back(std::array<size_t, 3>({{M1_counter, id0, id2}}));
      M1_counter++;
    }
  }// first run over lookuptable ends here - memory and new lookuptable 
   // are generated ------------------------------------------------------------

  #pragma omp for schedule(dynamic)
  for(int t1_i = 0; t1_i < Lt/dilT; t1_i++){
  for(int t2_i = t1_i; t2_i < Lt/dilT; t2_i++){
    // creating quarklines
    quarklines_Q2L.build(perambulators, meson_operator, t1_i, t2_i,
                               quark_lookup.Q2L, operator_lookup.ricQ2_lookup);
    quarklines_Q1.build(perambulators, meson_operator, t1_i, t2_i,
                               quark_lookup.Q1, operator_lookup.ricQ2_lookup);

  for(int dir = 0; dir < 2; dir++){

  if((t1_i == t2_i) && (dir == 1))
    continue;

  int t1_min, t2_min, t1_max, t2_max;
  if(dir == 0){
    t1_min = dilT*t1_i;
    t1_max = dilT*(t1_i+1);
    t2_min = dilT*t2_i;
    t2_max = dilT*(t2_i+1);
  }
  else{
    t1_min = dilT*t2_i;
    t1_max = dilT*(t2_i+1);
    t2_min = dilT*t1_i;
    t2_max = dilT*(t1_i+1);
  }

  for(int t1 = t1_min; t1 < t1_max; t1++){
  for(int t2 = t2_min; t2 < t2_max; t2++){
    int t = abs((t2 - t1 - (int)Lt) % (int)Lt);

    // quarkline indices
    int id_Q2L_1, id_Q2L_2;

    if(t1_i == t2_i){                                                                  
      if(t1_min != 0)                                                                  
        id_Q2L_1 = t1%t1_min;                                                       
      else{                                                                         
        id_Q2L_1 = t1;                                                                 
      }                                                                                
      if(t2_min != 0)                                                                  
        id_Q2L_2 = t2%t2_min;                                                       
      else{                                                                            
        id_Q2L_2 = t2;                                                                 
      }                                                                         
    }                                                                           
    else{                                                                              
      if(t1_min != 0)                                                               
         id_Q2L_1 = (dir)*dilT + t1%t1_min;                                         
       else                                                                  
         id_Q2L_1 = ((dir)*dilT + t1);                                                          
       if(t2_min != 0)                                                               
         id_Q2L_2 = (dir+1)%2*dilT + t2%t2_min;                           
       else                                                                     
         id_Q2L_2 = (dir+1)%2*dilT + t2;               
    }                                                  

    // build M1 ----------------------------------------------------------------
    for(const auto& look : M1_look){
      const auto& ric0 = operator_lookup.ricQ2_lookup[
                  quark_lookup.Q2L[look[1]].id_ric_lookup].rnd_vec_ids;
      const auto& ric2 = operator_lookup.ricQ2_lookup[
                  operator_lookup.rvdaggervr_lookuptable[look[2]].
                                                    id_ricQ_lookup].rnd_vec_ids;
      size_t M1_rnd_counter = 0;
      for(const auto& rnd0 : ric0){
      for(const auto& rnd2 : ric2){
      if(rnd0.first == rnd2.first && rnd0.second != rnd2.second){
        const size_t idr0 = &rnd0 - &ric0[0];
        const size_t idr2 = &rnd2 - &ric2[0];
        for(size_t col = 0; col < 4; col++){
  
          const cmplx value = quarklines_Q2L.return_gamma_val(5, col); // TODO: gamma hardcoded
          const size_t gamma_index = quarklines_Q2L.return_gamma_row(
                                                          5, col); // TODO: gamma hardcoded

          M1[look[0]][M1_rnd_counter].block(col*dilE, 0, dilE, 4*dilE) = value *
              meson_operator.return_rvdaggervr(look[2], t1, idr2).
                                 block(col*dilE, gamma_index*dilE, dilE, dilE) *
              quarklines_Q2L.return_Ql(id_Q2L_1, look[1], idr0).
                                 block(gamma_index*dilE, 0, dilE, 4*dilE);
        }
        M1_rnd_counter++;
      }}}
    }

    // Final summation for correlator ------------------------------------------
    for(const auto& c_look : corr_lookup){

      const auto& ric0 = operator_lookup.ricQ2_lookup[
                  quark_lookup.Q2L[c_look.lookup[0]].id_ric_lookup].rnd_vec_ids;
      const auto& ric1 = operator_lookup.ricQ2_lookup[
                  quark_lookup.Q1[c_look.lookup[1]].id_ric_lookup].rnd_vec_ids;
      const auto& ric2 = operator_lookup.ricQ2_lookup[
                  operator_lookup.rvdaggervr_lookuptable[c_look.lookup[2]].
                                                    id_ricQ_lookup].rnd_vec_ids;

      const size_t id0 = c_look.lookup[0];
      const size_t id1 = c_look.lookup[1];
      const size_t id2 = c_look.lookup[2];
      auto it1 = std::find_if(M1_look.begin(), M1_look.end(), 
                             [&id0, &id2](std::array<size_t, 3> check){
                               return (id0 == check[1] && id2 == check[2]);
                             });
      size_t M1_rnd_counter = 0;
      for(const auto& rnd0 : ric0){
      for(const auto& rnd2 : ric2){
      if(rnd0.first == rnd2.first && rnd0.second != rnd2.second){
        for(const auto& rnd1 : ric1){
        const auto idr1 = &rnd1 - &ric1[0];
        if(rnd1.first != rnd2.first  && rnd1.second == rnd2.second &&
           rnd1.first == rnd0.second && rnd1.second != rnd0.first){
          C[c_look.id][t] += (M1[(*it1)[0]][M1_rnd_counter]*
                    quarklines_Q1.return_Ql(id_Q2L_2, c_look.lookup[1], idr1)).trace();
        }}
        M1_rnd_counter++;
      }}}
    } // loop over operators ends here
//std::cout << "\n\nhier 2\n\n" << std::endl;
  }}}}} // loops over time end here
  #pragma omp critical
  {
    for(const auto& c_look : corr_lookup)
      for(size_t t = 0; t < Lt; t++)
        correlator[c_look.id][t] += C[c_look.id][t];
  }
  swatch.stop();
}// parallel part ends here


  // normalisation
  for(const auto& c_look : corr_lookup){
    for(auto& corr : correlator[c_look.id]){
      corr /= (6*5*4)*Lt; // TODO: Hard Coded atm - Be carefull
    }
    // write data to file
    filehandle.write(correlator[c_look.id], c_look);
  }
  swatch.print();
}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void LapH::Correlators::build_C4cB(const OperatorsForMesons& meson_operator,
                                   const Perambulator& perambulators,
                                   const OperatorLookup& operator_lookup,
                                   const std::vector<CorrInfo>& corr_lookup,
                                   const QuarklineLookup& quark_lookup, 
                                   const std::string output_path,
                                   const std::string output_filename) {
  if(corr_lookup.empty())
    return;

  StopWatch swatch("C4cB");

  // every element of corr_lookup contains the same filename. Wlog choose the 
  // first element
  WriteHDF5Correlator filehandle(output_path, "C4+B", output_filename, comp_type_factory_tr() );

  std::vector<vec> correlator(corr_lookup.size(), vec(Lt, cmplx(.0,.0)));

// This is necessary to ensure the correct summation of the correlation function
#pragma omp parallel
{
  swatch.start();
  std::vector<vec> C(corr_lookup.size(), vec(Lt, cmplx(.0,.0)));
  // building the quark line directly frees up a lot of memory
  QuarkLine_one_t<QuarkLineType::Q2L> quarklines(dilT, dilE, nev, quark_lookup.Q2L, 
                        operator_lookup.ricQ2_lookup);
  // creating memory arrays M1, M2 for intermediate storage of Quarklines ------
  std::vector<std::vector<Eigen::MatrixXcd> > M1, M2;
  std::vector<std::array<size_t, 3> > M1_look;
  std::vector<std::array<size_t, 3> > M2_look;
  size_t M1_counter = 0;
  size_t M2_counter = 0;
  for(const auto& c_look : corr_lookup){
    const auto& ric0 = operator_lookup.ricQ2_lookup[
                quark_lookup.Q2L[c_look.lookup[0]].id_ric_lookup].rnd_vec_ids;
    const auto& ric1 = operator_lookup.ricQ2_lookup[//needed only for checking
               operator_lookup.rvdaggervr_lookuptable[c_look.lookup[3]].
                                                  id_ricQ_lookup].rnd_vec_ids;
    const auto& ric2 = operator_lookup.ricQ2_lookup[
                quark_lookup.Q2L[c_look.lookup[2]].id_ric_lookup].rnd_vec_ids;
    const auto& ric3 = operator_lookup.ricQ2_lookup[//needed only for checking
               operator_lookup.rvdaggervr_lookuptable[c_look.lookup[1]].
                                                    id_ricQ_lookup].rnd_vec_ids;
    if(ric0.size() != ric1.size() || ric0.size() != ric2.size() || 
       ric0.size() != ric3.size()){
      std::cout << "rnd combinations are not the same in build_C4cB" 
                << std::endl;
    }
    // creating memeory for M1 -------------------------------------------------
    const size_t id3 = c_look.lookup[3];
    const size_t id0 = c_look.lookup[0];
    auto it1 = std::find_if(M1_look.begin(), M1_look.end(), 
                           [&id3, &id0](std::array<size_t, 3> check){
                             return (id3 == check[1] && id0 == check[2]);
                           });
    if(!(it1 != M1_look.end())){
      M1.emplace_back(std::vector<Eigen::MatrixXcd>());
      for(const auto& rnd0 : ric0)
      for(const auto& rnd1 : ric1)
      if(rnd0.first == rnd1.first && rnd0.second != rnd1.second)
        M1[M1_counter].emplace_back(Eigen::MatrixXcd::Zero(4*dilE, 4*dilE));
      M1_look.emplace_back(std::array<size_t, 3>({{M1_counter, id3, id0}}));
      M1_counter++;
    }
    // creating memeory for M2 -------------------------------------------------
    const size_t id1 = c_look.lookup[1];
    const size_t id2 = c_look.lookup[2];
    auto it2 = std::find_if(M2_look.begin(), M2_look.end(), 
                            [&id1, &id2](std::array<size_t, 3> check){
                              return (id1 == check[1] && id2 == check[2]);
                            });
    if(!(it2 != M2_look.end())){
      M2.emplace_back(std::vector<Eigen::MatrixXcd>());
      for(const auto& rnd2 : ric2)
      for(const auto& rnd3 : ric3)
      if(rnd2.first == rnd3.first && rnd2.second != rnd3.second)
        M2[M2_counter].emplace_back(Eigen::MatrixXcd::Zero(4*dilE, 4*dilE));
      M2_look.emplace_back(std::array<size_t, 3>({{M2_counter, id1, id2}}));
      M2_counter++;
    }
  }// first run over lookuptable ends here - memory and new lookuptable 
   // are generated ------------------------------------------------------------

  #pragma omp for schedule(dynamic)
  for(int t1_i = 0; t1_i < Lt/dilT; t1_i++){
  for(int t2_i = t1_i; t2_i < Lt/dilT; t2_i++){
    // creating quarklines
    quarklines.build(perambulators, meson_operator, t1_i, t2_i,
                               quark_lookup.Q2L, operator_lookup.ricQ2_lookup);

  for(int bla = 0; bla < 2; bla++){

  if((t1_i == t2_i) && (bla == 1))
    continue;

  int t1_min, t2_min, t1_max, t2_max;
  if(bla == 0){
    t1_min = dilT*t1_i;
    t1_max = dilT*(t1_i+1);
    t2_min = dilT*t2_i;
    t2_max = dilT*(t2_i+1);
  }
  else{
    t1_min = dilT*t2_i;
    t1_max = dilT*(t2_i+1);
    t2_min = dilT*t1_i;
    t2_max = dilT*(t1_i+1);
  }

  for(int t1 = t1_min; t1 < t1_max; t1++){
  for(int t2 = t2_min; t2 < t2_max; t2++){
    int t = abs((t2 - t1 - (int)Lt) % (int)Lt);

    // quarkline indices
    int id_Q2L_1, id_Q2L_2;

    if(t1_i == t2_i){                                                                  
      if(t1_min != 0)                                                                  
        id_Q2L_1 = t1%t1_min;                                                       
      else{                                                                         
        id_Q2L_1 = t1;                                                                 
      }                                                                                
      if(t2_min != 0)                                                                  
        id_Q2L_2 = t2%t2_min;                                                       
      else{                                                                            
        id_Q2L_2 = t2;                                                                 
      }                                                                         
    }                                                                           
    else{                                                                              
      if(t1_min != 0)                                                               
         id_Q2L_1 = (bla)*dilT + t1%t1_min;                                         
       else                                                                  
         id_Q2L_1 = ((bla)*dilT + t1);                                                          
       if(t2_min != 0)                                                               
         id_Q2L_2 = (bla+1)%2*dilT + t2%t2_min;                           
       else                                                                     
         id_Q2L_2 = (bla+1)%2*dilT + t2;               
    }                                                  

    // build M1 ----------------------------------------------------------------
    for(const auto& look : M1_look){
      const auto& ric0 = operator_lookup.ricQ2_lookup[
                  quark_lookup.Q2L[look[2]].id_ric_lookup].rnd_vec_ids;
      const auto& ric1 = operator_lookup.ricQ2_lookup[//needed only for checking
                 operator_lookup.rvdaggervr_lookuptable[look[1]].
                                                    id_ricQ_lookup].rnd_vec_ids;
      size_t M1_rnd_counter = 0;
      for(const auto& rnd0 : ric0){
      for(const auto& rnd1 : ric1){
      if(rnd0.first == rnd1.first && rnd0.second != rnd1.second){
        const size_t idr0 = &rnd0 - &ric0[0];
        const size_t idr1 = &rnd1 - &ric1[0];
        for(size_t col = 0; col < 4; col++){
          const cmplx value = quarklines.return_gamma_val(5, col); // TODO: gamma hardcoded
          const size_t gamma_index = quarklines.return_gamma_row(
                                                          5, col); // TODO: gamma hardcoded
          M1[look[0]][M1_rnd_counter].block(col*dilE, 0, dilE, 4*dilE) = value *
              meson_operator.return_rvdaggervr(look[1], t1, idr1).
                                 block(col*dilE, gamma_index*dilE, dilE, dilE) *
              quarklines.return_Ql(id_Q2L_1, look[2], idr0).
                                 block(gamma_index*dilE, 0, dilE, 4*dilE);
        }
        M1_rnd_counter++;
      }}}
    }
    // build M2 ----------------------------------------------------------------
    for(const auto& look : M2_look){
      const auto& ric2 = operator_lookup.ricQ2_lookup[
                  quark_lookup.Q2L[look[2]].id_ric_lookup].rnd_vec_ids;
      const auto& ric3 = operator_lookup.ricQ2_lookup[//needed only for checking
                 operator_lookup.rvdaggervr_lookuptable[look[1]].
                                                    id_ricQ_lookup].rnd_vec_ids;
      size_t M2_rnd_counter = 0;
      for(const auto& rnd2 : ric2){
      for(const auto& rnd3 : ric3){
      if(rnd2.first == rnd3.first && rnd2.second != rnd3.second){
        const size_t idr2 = &rnd2 - &ric2[0];
        const size_t idr3 = &rnd3 - &ric3[0];
        for(size_t col = 0; col < 4; col++){
          const cmplx value = 
                           quarklines.return_gamma_val(5, col); // TODO: gamma hardcoded
          const size_t gamma_index = quarklines.return_gamma_row(5, col); // TODO: gamma hardcoded
          M2[look[0]][M2_rnd_counter].
            block(col*dilE, 0, dilE, 4*dilE) = value *
            meson_operator.return_rvdaggervr(look[1], t2, idr3).
                                block(col*dilE, gamma_index*dilE, dilE, dilE)*
            quarklines.return_Ql(id_Q2L_2, look[2], idr2).
                               block(gamma_index*dilE, 0, dilE, 4*dilE);
        }
        M2_rnd_counter++;
      }}}
    }
    // Final summation for correlator ------------------------------------------
    Eigen::MatrixXcd M3 = Eigen::MatrixXcd::Zero(4*dilE, 4*dilE);
    for(const auto& c_look : corr_lookup){
      const auto& ric0 = operator_lookup.ricQ2_lookup[
                  quark_lookup.Q2L[c_look.lookup[0]].id_ric_lookup].rnd_vec_ids;
      const auto& ric1 = operator_lookup.ricQ2_lookup[//needed only for checking
                 operator_lookup.rvdaggervr_lookuptable[c_look.lookup[3]].
                                                    id_ricQ_lookup].rnd_vec_ids;
      const auto& ric2 = operator_lookup.ricQ2_lookup[
                  quark_lookup.Q2L[c_look.lookup[2]].id_ric_lookup].rnd_vec_ids;
      const auto& ric3 = operator_lookup.ricQ2_lookup[//needed only for checking
                 operator_lookup.rvdaggervr_lookuptable[c_look.lookup[1]].
                                                    id_ricQ_lookup].rnd_vec_ids;
      const size_t id3 = c_look.lookup[3];
      const size_t id0 = c_look.lookup[0];
      auto it1 = std::find_if(M1_look.begin(), M1_look.end(), 
                             [&id3, &id0](std::array<size_t, 3> check){
                               return (id3 == check[1] && id0 == check[2]);
                             });
      const size_t id1 = c_look.lookup[1];
      const size_t id2 = c_look.lookup[2];
      auto it2 = std::find_if(M2_look.begin(), M2_look.end(), 
                              [&id1, &id2](std::array<size_t, 3> check){
                                return (id1 == check[1] && id2 == check[2]);
                              });
      size_t M1_rnd_counter = 0;
      for(const auto& rnd0 : ric0){
      for(const auto& rnd1 : ric1){
      if(rnd0.first == rnd1.first && rnd0.second != rnd1.second){
        M3.setZero(4 * dilE, 4 * dilE); // setting matrix values to zero
        size_t M2_rnd_counter = 0;
        for(const auto& rnd2 : ric2){
        for(const auto& rnd3 : ric3){
        if(rnd2.first == rnd3.first && rnd2.second != rnd3.second){
          if(rnd0.second == rnd3.second && rnd1.second == rnd2.second &&
             rnd0.first != rnd2.first){
            M3 += M2[(*it2)[0]][M2_rnd_counter];
          }
          M2_rnd_counter++;
        }}}
        C[c_look.id][t] += (M1[(*it1)[0]][M1_rnd_counter++]*M3).trace();
      }}}
    } // loop over operators ends here
  }}}}} // loops over time end here
  #pragma omp critical
  {
    for(const auto& c_look : corr_lookup)
      for(size_t t = 0; t < Lt; t++)
        correlator[c_look.id][t] += C[c_look.id][t];
  }
  swatch.stop();
}// parallel part ends here


  // normalisation
  for(const auto& c_look : corr_lookup){
    for(auto& corr : correlator[c_look.id]){
      corr /= (6*5*4*3)*Lt; // TODO: Hard Coded atm - Be carefull
    }
    // write data to file
    filehandle.write(correlator[c_look.id], c_look);
  }
  swatch.print();
}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void LapH::Correlators::build_C30(const Quarklines& quarklines,
                                  const std::vector<CorrInfo>& corr_lookup,
                                  const QuarklineLookup& quark_lookup,
                                  const std::vector<RandomIndexCombinationsQ2>& ric_lookup, 
                                  const std::string output_path,
                                  const std::string output_filename) {
  if(corr_lookup.empty())
    return;

  StopWatch swatch("C30");
  swatch.start();

  // every element of corr_lookup contains the same filename. Wlog choose the 
  // first element
  WriteHDF5Correlator filehandle(output_path, "C30", output_filename, comp_type_factory_tr() );

  for(const auto& c_look : corr_lookup){
    std::vector<cmplx> correlator(Lt, cmplx(.0,.0));
    const auto& ric0 = ric_lookup[quark_lookup.Q1[c_look.lookup[0]].
                                                     id_ric_lookup].rnd_vec_ids;
    const auto& ric1 = ric_lookup[quark_lookup.Q1[c_look.lookup[1]].
                                                     id_ric_lookup].rnd_vec_ids;
    const auto& ric2 = ric_lookup[quark_lookup.Q1[c_look.lookup[2]].
                                                     id_ric_lookup].rnd_vec_ids;
    if(ric0.size() != ric1.size() || ric0.size() != ric2.size()){
      std::cout << "rnd combinations are not the same in build_C30" 
                << std::endl;
      exit(0);
    }

    size_t norm = 0;
// This is necessary to ensure the correct summation of the correlation function
#pragma omp parallel reduction(+:norm)
{
    std::vector<cmplx> C(Lt, cmplx(.0,.0));
    #pragma omp for schedule(dynamic) 
    for(int t1 = 0; t1 < Lt; t1++){
    for(int t2 = 0; t2 < Lt; t2++){
      int t = abs((t2 - t1 - (int)Lt) % (int)Lt);
      for(const auto& rnd0 : ric0){
      for(const auto& rnd1 : ric1){
      if(rnd0.second == rnd1.first && rnd0.first != rnd1.second){
        const auto L1 =
          quarklines.return_Q1(t1, t2/dilT, c_look.lookup[0], &rnd0-&ric0[0]) *
          quarklines.return_Q1(t2, t1/dilT, c_look.lookup[1], &rnd1-&ric1[0]);
        for(const auto& rnd2 : ric2){
        if(rnd1.second == rnd2.first && rnd2.second == rnd0.first){
          C[t] += (L1*quarklines.return_Q1(t1, t1/dilT, c_look.lookup[2], 
                                                       &rnd2-&ric2[0])).trace();
          norm++;
        }}
      }}}
    }}
    #pragma omp critical
    {
      for(size_t t = 0; t < Lt; t++)
        correlator[t] += C[t];
    }
}// parallel part ends here

    // normalisation
    for(auto& corr : correlator){
      corr /= norm/Lt;
    }
    // write data to file
    filehandle.write(correlator, c_look);
  }
  swatch.stop();
  swatch.print();
}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
/*! 
 *  @note Not optimised.
 */
void LapH::Correlators::build_C40C(const Quarklines& quarklines,
                                   const std::vector<CorrInfo>& corr_lookup,
                                   const QuarklineLookup& quark_lookup,
                                   const std::vector<RandomIndexCombinationsQ2>& ric_lookup, 
                                   const std::string output_path,
                                   const std::string output_filename) {

  if(corr_lookup.empty())
    return;

  StopWatch swatch("C40C");
  swatch.start();

  // every element of corr_lookup contains the same filename. Wlog choose the
  // first element
  WriteHDF5Correlator filehandle(output_path, "C40C", output_filename, comp_type_factory_tr() );

  for(const auto& c_look : corr_lookup){
    std::vector<cmplx> correlator(Lt, cmplx(.0,.0));
    const auto& ric0 = ric_lookup[quark_lookup.Q1[c_look.lookup[0]].
                                                     id_ric_lookup].rnd_vec_ids;
    const auto& ric1 = ric_lookup[quark_lookup.Q1[c_look.lookup[1]].
                                                     id_ric_lookup].rnd_vec_ids;
    const auto& ric2 = ric_lookup[quark_lookup.Q1[c_look.lookup[2]].
                                                     id_ric_lookup].rnd_vec_ids;
    const auto& ric3 = ric_lookup[quark_lookup.Q1[c_look.lookup[3]].
                                                     id_ric_lookup].rnd_vec_ids;
    if(ric0.size() != ric1.size() || ric0.size() != ric2.size() || 
       ric0.size() != ric3.size()){
      std::cout << "rnd combinations are not the same in C40C" 
                << std::endl;
      exit(0);
    }

    size_t norm = 0;
// This is necessary to ensure the correct summation of the correlation function
#pragma omp parallel reduction(+:norm)
{
    std::vector<cmplx> C(Lt, cmplx(.0,.0));
    #pragma omp for schedule(dynamic) 
    for(int t1 = 0; t1 < Lt; t1++){
    for(int t2 = 0; t2 < Lt; t2++){
      int t = abs((t2 - t1 - (int)Lt) % (int)Lt);
      for(const auto& rnd0 : ric0){
      for(const auto& rnd1 : ric1){
      if(rnd0.second == rnd1.first && rnd0.first != rnd1.second){
        const auto L1 =
          quarklines.return_Q1(t1, t2/dilT, c_look.lookup[0], &rnd0-&ric0[0]) *
          quarklines.return_Q1(t2, t1/dilT, c_look.lookup[1], &rnd1-&ric1[0]);
        for(const auto& rnd2 : ric2){
        for(const auto& rnd3 : ric3){
        if(rnd1.second == rnd2.first && rnd2.second == rnd3.first && 
           rnd3.second == rnd0.first && rnd2.first != rnd3.second &&
           rnd0.second != rnd3.first){
          const auto L2 =
            quarklines.return_Q1(t1, t2/dilT, c_look.lookup[2], &rnd2-&ric2[0])*
            quarklines.return_Q1(t2, t1/dilT, c_look.lookup[3], &rnd3-&ric3[0]);
          C[t] += (L1*L2).trace();
          norm++;
        }}}
      }}}
    }}
    #pragma omp critical
    {
      for(size_t t = 0; t < Lt; t++)
        correlator[t] += C[t];
    }
}// parallel part ends here

    // normalisation
    for(auto& corr : correlator){
      corr /= norm/Lt;
    }
    // write data to file
    filehandle.write(correlator, c_look);
  }
  swatch.stop();
  swatch.print();
}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
/*! 
 *  @note Not optimised.
 */
void LapH::Correlators::build_C40B(const Quarklines& quarklines,
                     const std::vector<CorrInfo>& corr_lookup,
                     const QuarklineLookup& quark_lookup,
                     const std::vector<RandomIndexCombinationsQ2>& ric_lookup, 
                     const std::string output_path,
                     const std::string output_filename) {

  if(corr_lookup.empty())
    return;

  StopWatch swatch("C40B");
  swatch.start();

  // every element of corr_lookup contains the same filename. Wlog choose the
  // first element
  WriteHDF5Correlator filehandle(output_path, "C40B", output_filename, comp_type_factory_tr() );

  for(const auto& c_look : corr_lookup){
    std::vector<cmplx> correlator(Lt, cmplx(.0,.0));
    const auto& ric0 = ric_lookup[quark_lookup.Q1[c_look.lookup[0]].
                                                     id_ric_lookup].rnd_vec_ids;
    const auto& ric1 = ric_lookup[quark_lookup.Q1[c_look.lookup[1]].
                                                     id_ric_lookup].rnd_vec_ids;
    const auto& ric2 = ric_lookup[quark_lookup.Q1[c_look.lookup[2]].
                                                     id_ric_lookup].rnd_vec_ids;
    const auto& ric3 = ric_lookup[quark_lookup.Q1[c_look.lookup[3]].
                                                     id_ric_lookup].rnd_vec_ids;
//    if(ric0.size() != ric1.size() || ric0.size() != ric2.size() || 
//       ric0.size() != ric3.size()){
//      std::cout << "rnd combinations are not the same in build_corr0" 
//                << std::endl;
//      exit(0);
//    }

    size_t norm = 0;
// This is necessary to ensure the correct summation of the correlation function
#pragma omp parallel reduction(+:norm)
{
    std::vector<cmplx> C(Lt, cmplx(.0,.0));
    #pragma omp for schedule(dynamic) 
    for(int t1 = 0; t1 < Lt; t1++){
    for(int t2 = 0; t2 < Lt; t2++){
      int t = abs((t2 - t1 - (int)Lt) % (int)Lt);
      for(const auto& rnd0 : ric0){
      for(const auto& rnd1 : ric1){
      if(rnd0.second == rnd1.first && rnd0.first != rnd1.second){
        const auto L1 =
          quarklines.return_Q1(t1, t2/dilT, c_look.lookup[0], &rnd0-&ric0[0]) *
          quarklines.return_Q1(t2, t2/dilT, c_look.lookup[1], &rnd1-&ric1[0]);
        for(const auto& rnd2 : ric2){
        for(const auto& rnd3 : ric3){
        if(rnd1.second == rnd2.first && rnd2.second == rnd3.first && 
           rnd3.second == rnd0.first && rnd2.first != rnd3.second &&
           rnd0.second != rnd3.first){
          const auto L2 =
            quarklines.return_Q1(t2, t1/dilT, c_look.lookup[2], &rnd2-&ric2[0])*
            quarklines.return_Q1(t1, t1/dilT, c_look.lookup[3], &rnd3-&ric3[0]);
          C[t] += (L1*L2).trace();
          norm++;
        }}}
      }}}
    }}
    #pragma omp critical
    {
      for(size_t t = 0; t < Lt; t++)
        correlator[t] += C[t];
    }
}// parallel part ends here

    // normalisation
    for(auto& corr : correlator){
      corr /= norm/Lt;
    }
    // write data to file
    filehandle.write(correlator, c_look);
  }

  swatch.stop();
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
//void LapH::Correlators::contract (Quarklines& quarklines, 
void LapH::Correlators::contract (Quarklines& quarklines, 
                     const OperatorsForMesons& meson_operator,
                     const Perambulator& perambulators,
                     const OperatorLookup& operator_lookup,
                     const CorrelatorLookup& corr_lookup, 
                     const QuarklineLookup& quark_lookup,
                     const std::string output_path,
                     const std::string output_filename) {

  // 1. Build all functions which need corrC and free it afterwards.
  build_corrC(perambulators, meson_operator, operator_lookup, 
              corr_lookup.corrC, quark_lookup);
  build_C2c(corr_lookup.C2c, output_path, output_filename);
  build_C4cD(operator_lookup, corr_lookup, quark_lookup, output_path, output_filename);
  build_C4cV(operator_lookup, corr_lookup, quark_lookup, output_path, output_filename);
  // 2. Build all functions which need corr0 and free it afterwards.
  build_corr0(meson_operator, perambulators, corr_lookup.corr0, 
              quark_lookup, operator_lookup);
  // in C3c, also corr0 is build, since this is much faster
  build_C3c(meson_operator, perambulators, operator_lookup, corr_lookup.C3c, 
                                           quark_lookup, output_path, output_filename);
  build_C20(corr_lookup.C20, output_path, output_filename);
  build_C40D(operator_lookup, corr_lookup, quark_lookup, output_path, output_filename);
  build_C40V(operator_lookup, corr_lookup, quark_lookup, output_path, output_filename);
  // 3. Build all other correlation functions.
//  build_C1(quarklines, corr_lookup.C1, quark_lookup, 
//                                                 operator_lookup.ricQ2_lookup);
  build_C4cC(meson_operator, perambulators, operator_lookup, corr_lookup.C4cC, 
                                           quark_lookup, output_path, output_filename);
//  build_C4cC(quarklines, meson_operator, operator_lookup, corr_lookup.C4cC, 
//                                                                 quark_lookup);
  build_C4cB(meson_operator, perambulators, operator_lookup, corr_lookup.C4cB, 
                                           quark_lookup, output_path, output_filename);
  build_C30(quarklines, corr_lookup.C30, quark_lookup, 
                           operator_lookup.ricQ2_lookup, output_path, output_filename);
  build_C40C(quarklines, corr_lookup.C40C, quark_lookup, 
                           operator_lookup.ricQ2_lookup, output_path, output_filename);
  build_C40B(quarklines, corr_lookup.C40B, quark_lookup, 
                           operator_lookup.ricQ2_lookup, output_path, output_filename);
}

