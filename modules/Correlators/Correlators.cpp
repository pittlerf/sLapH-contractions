#include "Correlators.h"

/*! @TODO Why is the hdf5 stuff not in an unnamed namespace or a seperate 
 *        file? 
 */
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// ugly way to check if group exists and if not to create it -
// Alternative: use c API (don't want to)
inline void open_or_create_hdf5_group(const std::string& GROUP_NAME,
                                      const H5::H5File& file, H5::Group& group){
  try{
    group = file.openGroup(GROUP_NAME.c_str());
  }
  catch(H5::Exception& e){
    group = file.createGroup(GROUP_NAME.c_str());
  }
}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// ugly stays ugly :(
inline void open_or_create_hdf5_file(const H5std_string FILE_NAME, 
                                     H5::H5File&file){

  try{
    file = H5::H5File(FILE_NAME, H5F_ACC_EXCL);
  }
  catch(H5::Exception& e){
    file = H5::H5File(FILE_NAME, H5F_ACC_RDWR);
  }
}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
static void write_correlators(const std::vector<cmplx>& corr, 
                              const CorrInfo& corr_info){
  // check if directory exists
  if(access( corr_info.outpath.c_str(), 0 ) != 0) {
      std::cout << "\tdirectory " << corr_info.outpath.c_str() 
                << " does not exist and will be created";
      boost::filesystem::path dir(corr_info.outpath.c_str());
      if(boost::filesystem::create_directories(dir))
        std::cout << "\tSuccess" << std::endl;
      else
        std::cout << "\tFailure" << std::endl;
  }
  // writing the data ----------------------------------------------------------
  try
  {
    // exceptins will be catched at the end and not printed
    H5::Exception::dontPrint(); 
    // hdf5 data
    H5::H5File file;
    H5::Group group;
    H5::DataSet dset;
    // create new memory data type for writing for COMPLEX numbers -------------
    H5::CompType cmplx_w(sizeof(std::complex<double>));
    auto type = H5::PredType::NATIVE_DOUBLE;
    cmplx_w.insertMember("re", HOFFSET(LapH::complex_t, re), type);
    cmplx_w.insertMember("im", HOFFSET(LapH::complex_t, im), type);
    // open file or create the file if it does not exist -----------------------
    const H5std_string FILE_NAME((corr_info.outpath+corr_info.outfile).c_str());
    open_or_create_hdf5_file(FILE_NAME, file);
    // create the dataset to write data ----------------------------------------
    H5std_string DATASET_NAME((corr_info.hdf5_dataset_name).c_str());
    hsize_t dim(corr.size());
    H5::DataSpace dspace(1, &dim);
    dset = file.createDataSet(DATASET_NAME, cmplx_w, dspace);
    dset.write(&corr[0], cmplx_w);
    dset.close();
    file.close();

  } // end of try block - catching all the bad things --------------------------
  // catch failure caused by the H5File operations
  catch(H5::FileIException error){
     error.printError();
  }
  // catch failure caused by the DataSet operations
  catch(H5::DataSetIException error){
     error.printError();
  }
  // catch failure caused by the DataSpace operations
  catch(H5::DataSpaceIException error){
     error.printError();
  }
  // catch failure caused by the DataSpace operations
  catch(H5::DataTypeIException error){
     error.printError();
  }
  // catch failure caused by the Group operations
  catch(H5::GroupIException error){
     error.printError();
  }
}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// TODO: Bad style: code duplication - see write_correlators
static void write_4pt_correlators(const std::vector<LapH::compcomp_t>& corr, 
                                  const CorrInfo& corr_info){
  // check if directory exists
  if(access( corr_info.outpath.c_str(), 0 ) != 0) {
      std::cout << "\tdirectory " << corr_info.outpath.c_str() 
                << " does not exist and will be created";
      boost::filesystem::path dir(corr_info.outpath.c_str());
      if(!boost::filesystem::create_directories(dir))
        std::cout << "\tSuccess" << std::endl;
      else
        std::cout << "\tFailure" << std::endl;
  }
  // writing the data ----------------------------------------------------------
  try
  {
    // exceptins will be catched at the end and not printed
    H5::Exception::dontPrint(); 
    // hdf5 data
    H5::H5File file;
    H5::Group group;
    H5::DataSet dset;
    // create new memory data type for writing for COMPLEX numbers -------------
    H5::CompType cmplxcmplx_w(4*sizeof(double));
    auto type = H5::PredType::NATIVE_DOUBLE;
    cmplxcmplx_w.insertMember("rere", HOFFSET(LapH::compcomp_t, rere), type);
    cmplxcmplx_w.insertMember("reim", HOFFSET(LapH::compcomp_t, reim), type);
    cmplxcmplx_w.insertMember("imre", HOFFSET(LapH::compcomp_t, imre), type);
    cmplxcmplx_w.insertMember("imim", HOFFSET(LapH::compcomp_t, imim), type);
    // open file or create the file if it does not exist -----------------------
    const H5std_string FILE_NAME((corr_info.outpath+corr_info.outfile).c_str());
    open_or_create_hdf5_file(FILE_NAME, file);
    // create the dataset to write data ----------------------------------------
    H5std_string DATASET_NAME((corr_info.hdf5_dataset_name).c_str());
    hsize_t dim(corr.size());
    H5::DataSpace dspace(1, &dim);
    dset = file.createDataSet(DATASET_NAME, cmplxcmplx_w, dspace);
    dset.write(&corr[0], cmplxcmplx_w);
    dset.close();
    file.close();

  } // end of try block - catching all the bad things --------------------------
  // catch failure caused by the H5File operations
  catch(H5::FileIException error){
     error.printError();
  }
  // catch failure caused by the DataSet operations
  catch(H5::DataSetIException error){
     error.printError();
  }
  // catch failure caused by the DataSpace operations
  catch(H5::DataSpaceIException error){
     error.printError();
  }
  // catch failure caused by the DataSpace operations
  catch(H5::DataTypeIException error){
     error.printError();
  }
  // catch failure caused by the Group operations
  catch(H5::GroupIException error){
     error.printError();
  }
}

/******************************************************************************/
/******************************************************************************/

/*******************************************************************************/
/**
 *
 * @Param quarklines
 * @Param corr_lookup
 * @Param quark_lookup
 * @Param ric_lookup
 */
void LapH::Correlators::build_C1(const Quarklines& quarklines,
                    const std::vector<CorrInfo>& corr_lookup,
                    const QuarklineLookup& quark_lookup,
                    const std::vector<RandomIndexCombinationsQ2>& ric_lookup) {

  std::cout << "\tcomputing C1:";
  clock_t time = clock();

  for(const auto& c_look : corr_lookup){
    const auto& ric = ric_lookup[quark_lookup.Q1[c_look.lookup[0]].
                                                     id_ric_lookup].rnd_vec_ids;
    std::vector<cmplx> correlator(Lt*ric.size(), cmplx(.0,.0));
    for(size_t t = 0; t < Lt; t++){
      for(const auto& id : ric){
        correlator[(&id-&ric[0])*Lt + t] +=
         quarklines.return_Q1(t, t/dilT, c_look.lookup[0], &id-&ric[0]).trace();
      }
    }
    // write data to file
    write_correlators(correlator, c_look);
  }
  time = clock() - time;
  std::cout << "\t\t\tSUCCESS - " << ((float) time) / CLOCKS_PER_SEC 
            << " seconds" << std::endl;
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

  if(corr_lookup.size() == 0)
    return;

  std::cout << "\tcomputing corr0:";
  clock_t time = clock();

  corr0.resize(boost::extents[corr_lookup.size()][Lt][Lt]);
#pragma omp parallel
{
  Quarklines_one_t quarklines_intern(2, dilT, dilE, nev, quark_lookup, 
                               operator_lookup.ricQ2_lookup);

  #pragma omp for schedule(dynamic) 
  for(int t1 = 0; t1 < Lt; t1++){
  for(int t2 = 0; t2 < Lt; t2++){
    quarklines_intern.build_Q1_one_t(perambulators, meson_operator, t1, t2,
                              quark_lookup.Q1, operator_lookup.ricQ2_lookup);

    for(const auto& c_look : corr_lookup){
      const auto& ric0 = operator_lookup.ricQ2_lookup[quark_lookup.Q1[
                                   c_look.lookup[0]].id_ric_lookup].rnd_vec_ids;
      const auto& ric1 = operator_lookup.ricQ2_lookup[quark_lookup.Q1[
                                   c_look.lookup[1]].id_ric_lookup].rnd_vec_ids;
      if(ric0.size() != ric1.size()){
        std::cout << "rnd combinations are not the same in build_corr0" 
                  << std::endl;
        exit(0);
      }
      corr0[c_look.id][t1][t2].resize(ric0.size());
      for(auto& corr : corr0[c_look.id][t1][t2])
        corr = cmplx(0.0,0.0);
      for(const auto& rnd : ric0){
        const auto id = &rnd - &ric0[0];
        const auto it1 = std::find_if(ric1.begin(), ric1.end(),
                                [&](std::pair<size_t, size_t> pair){
                                  return (pair == 
                                       std::make_pair(rnd.second, rnd.first));
                                });
        if(it1 == ric1.end()){
          std::cout << "something wrong with random vectors in build_corr0" 
                    << std::endl;
          exit(0);
        }
        corr0[c_look.id][t1][t2][id] += 
                    (quarklines_intern.return_Q1(0, 0, c_look.lookup[0], id) *
                     quarklines_intern.return_Q1(1, 0, c_look.lookup[1], 
                                                     it1-ric1.begin())).trace();
      }
    }
  }}
}
  time = clock() - time;
  std::cout << "\t\tSUCCESS - " << ((float) time) / CLOCKS_PER_SEC 
            << " seconds" << std::endl;
}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void LapH::Correlators::build_C20(const std::vector<CorrInfo>& corr_lookup) {

  std::cout << "\tcomputing C20:";
  clock_t time = clock();

  for(const auto& c_look : corr_lookup){
    std::vector<cmplx> correlator(Lt, cmplx(.0,.0));
    for(int t1 = 0; t1 < Lt; t1++){
    for(int t2 = 0; t2 < Lt; t2++){
      int t = abs((t2 - t1 - (int)Lt) % (int)Lt);
      for(const auto& corr : corr0[c_look.lookup[0]][t1][t2])
        correlator[t] += corr;
    }}
    // normalisation
    for(auto& corr : correlator)
      corr /= Lt*corr0[c_look.lookup[0]][0][0].size();
    // write data to file
    write_correlators(correlator, c_look);
  }

  time = clock() - time;
  std::cout << "\t\tSUCCESS - " << ((float) time) / CLOCKS_PER_SEC 
            << " seconds" << std::endl;
}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void LapH::Correlators::build_C40D(const OperatorLookup& operator_lookup, 
                                   const CorrelatorLookup& corr_lookup,
                                   const QuarklineLookup& quark_lookup) {

  std::cout << "\tcomputing C40D:";
  clock_t time = clock();

  for(const auto& c_look : corr_lookup.C40D){
    std::vector<LapH::compcomp_t> correlator(Lt, LapH::compcomp_t(.0,.0,.0,.0));
    const size_t id0 = corr_lookup.corr0[c_look.lookup[0]].lookup[0];
    const size_t id1 = corr_lookup.corr0[c_look.lookup[1]].lookup[0];
    const auto& ric0 = operator_lookup.ricQ2_lookup[quark_lookup.Q1[id0].
                                                     id_ric_lookup].rnd_vec_ids;
    const auto& ric1 = operator_lookup.ricQ2_lookup[quark_lookup.Q1[id1].
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
                   corr0[c_look.lookup[0]][t1][t2].at(&rnd0 - &ric0[0]).real() *
                   corr0[c_look.lookup[1]][t1][t2].at(&rnd1 - &ric1[0]).real();
        correlator[t].reim += 
                   corr0[c_look.lookup[0]][t1][t2].at(&rnd0 - &ric0[0]).real() *
                   corr0[c_look.lookup[1]][t1][t2].at(&rnd1 - &ric1[0]).imag();
        correlator[t].imre += 
                   corr0[c_look.lookup[0]][t1][t2].at(&rnd0 - &ric0[0]).imag() *
                   corr0[c_look.lookup[1]][t1][t2].at(&rnd1 - &ric1[0]).real();
        correlator[t].imim += 
                   corr0[c_look.lookup[0]][t1][t2].at(&rnd0 - &ric0[0]).imag() *
                   corr0[c_look.lookup[1]][t1][t2].at(&rnd1 - &ric1[0]).imag();
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
    write_4pt_correlators(correlator, c_look);
  }

  time = clock() - time;
  std::cout << "\t\tSUCCESS - " << ((float) time) / CLOCKS_PER_SEC 
            << " seconds" << std::endl;
}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void LapH::Correlators::build_C40V(const OperatorLookup& operator_lookup, 
                                   const CorrelatorLookup& corr_lookup,
                                   const QuarklineLookup& quark_lookup) {

  std::cout << "\tcomputing C40V:";
  clock_t time = clock();

  for(const auto& c_look : corr_lookup.C40V){
    std::vector<LapH::compcomp_t> correlator(Lt, LapH::compcomp_t(.0,.0,.0,.0));
    const size_t id0 = corr_lookup.corr0[c_look.lookup[0]].lookup[0];
    const size_t id1 = corr_lookup.corr0[c_look.lookup[1]].lookup[0];
    const auto& ric0 = operator_lookup.ricQ2_lookup[quark_lookup.Q1[id0].
                                                     id_ric_lookup].rnd_vec_ids;
    const auto& ric1 = operator_lookup.ricQ2_lookup[quark_lookup.Q1[id1].
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
                   corr0[c_look.lookup[0]][t1][t1].at(&rnd0 - &ric0[0]).real() *
                   corr0[c_look.lookup[1]][t2][t2].at(&rnd1 - &ric1[0]).real();
        correlator[t].reim += 
                   corr0[c_look.lookup[0]][t1][t1].at(&rnd0 - &ric0[0]).real() *
                   corr0[c_look.lookup[1]][t2][t2].at(&rnd1 - &ric1[0]).imag();
        correlator[t].imre += 
                   corr0[c_look.lookup[0]][t1][t1].at(&rnd0 - &ric0[0]).imag() *
                   corr0[c_look.lookup[1]][t2][t2].at(&rnd1 - &ric1[0]).real();
        correlator[t].imim += 
                   corr0[c_look.lookup[0]][t1][t1].at(&rnd0 - &ric0[0]).imag() *
                   corr0[c_look.lookup[1]][t2][t2].at(&rnd1 - &ric1[0]).imag();
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
    write_4pt_correlators(correlator, c_look);
  }

  time = clock() - time;
  std::cout << "\t\tSUCCESS - " << ((float) time) / CLOCKS_PER_SEC 
            << " seconds" << std::endl;
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void LapH::Correlators::build_corrC(const Perambulator& perambulators,
                                    const OperatorsForMesons& meson_operator,
                                    const OperatorLookup& operator_lookup,
                                    const std::vector<CorrInfo>& corr_lookup,
                                    const QuarklineLookup& quark_lookup) {

  if(corr_lookup.size() == 0)
    return;

  std::cout << "\tcomputing corrC:";
  clock_t time = clock();

  corrC.resize(boost::extents[corr_lookup.size()][Lt][Lt]);

#pragma omp parallel
{
  // building the quark line directly frees up a lot of memory
  Quarklines_one_t quarklines(2*dilT, dilT, dilE, nev, quark_lookup, 
                        operator_lookup.ricQ2_lookup);
  #pragma omp for schedule(dynamic)
  for(int t1_i = 0; t1_i < Lt/dilT; t1_i++){
  for(int t2_i = t1_i; t2_i < Lt/dilT; t2_i++){
    quarklines.build_Q2V_one_t(perambulators, meson_operator, t1_i, t2_i,
                              quark_lookup.Q2V, operator_lookup.ricQ2_lookup);
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

        // quarkline indices
        int id_Q2L_1;
    
        if(t1_i == t2_i){   
          if(t1_min != 0)  
            id_Q2L_1 = t1%t1_min;    
          else{        
            id_Q2L_1 = t1;    
          }    
        }      
        else{    
          if(t1_min != 0) 
             id_Q2L_1 = (dir)*dilT + t1%t1_min; 
           else   
             id_Q2L_1 = ((dir)*dilT + t1);     
        }                                                  

        // building correlator -------------------------------------------------
        for(const auto& c_look : corr_lookup){
          const auto& ric0 = operator_lookup.ricQ2_lookup[
                  quark_lookup.Q2V[c_look.lookup[0]].id_ric_lookup].rnd_vec_ids;
          const auto& ric1 = operator_lookup.ricQ2_lookup[//just for checking
                       operator_lookup.rvdaggervr_lookuptable[c_look.lookup[1]].
                                                    id_ricQ_lookup].rnd_vec_ids;
          if(ric0.size() != ric1.size()){
            std::cout << "rnd combinations are not the same in build_corrC" 
                      << std::endl;
            exit(0);
          }
          corrC[c_look.id][t1][t2].resize(ric0.size());
          for(auto& corr : corrC[c_look.id][t1][t2])
            corr = cmplx(0.0,0.0);
          for(const auto& rnd : ric0){
            const auto id = &rnd - &ric0[0];
            for(size_t block = 0; block < 4; block++){
              const auto gamma_index = 
                           quarklines.return_gamma_row(c_look.gamma[0], block);
              corrC[c_look.id][t1][t2][id] += 
                   quarklines.return_gamma_val(c_look.gamma[0], block) *
                   (quarklines.return_Q2V(id_Q2L_1, 0, c_look.lookup[0], id).
                              block(block*dilE, gamma_index*dilE, dilE, dilE) *
                   meson_operator.return_rvdaggervr(c_look.lookup[1], t2, id).
                      block(gamma_index*dilE, block*dilE, dilE, dilE)).trace();
            }
          }
        }
      }}// t1, t2 end here
    }// dir (directions) end here
  }}// block times end here
}// omp parall ends here

  time = clock() - time;
  std::cout << "\t\tSUCCESS - " << ((float) time) / CLOCKS_PER_SEC 
            << " seconds" << std::endl;
}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void LapH::Correlators::build_C2c(const std::vector<CorrInfo>& corr_lookup) {

  std::cout << "\tcomputing C2c:";
  clock_t time = clock();

  for(const auto& c_look : corr_lookup){
    std::vector<cmplx> correlator(Lt, cmplx(.0,.0));
    if(c_look.outfile.find("Check") == 0){
      for(int t1 = 0; t1 < Lt; t1++){
        for(const auto& corr : corrC[c_look.lookup[0]][t1][t1]){
          correlator[t1] += corr;
        }
      }
      // normalisation
      for(auto& corr : correlator)
        corr /= corrC[c_look.lookup[0]][0][0].size();
      // write data to file
      write_correlators(correlator, c_look);
    }
    else{
      for(int t1 = 0; t1 < Lt; t1++){
      for(int t2 = 0; t2 < Lt; t2++){
        int t = abs((t2 - t1 - (int)Lt) % (int)Lt);
        for(const auto& corr : corrC[c_look.lookup[0]][t1][t2]){
          correlator[t] += corr;
        }
      }}
      // normalisation
      for(auto& corr : correlator)
        corr /= Lt*corrC[c_look.lookup[0]][0][0].size();
      // write data to file
      write_correlators(correlator, c_look);
    }
  }

  time = clock() - time;
  std::cout << "\t\tSUCCESS - " << ((float) time) / CLOCKS_PER_SEC 
            << " seconds" << std::endl;
}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void LapH::Correlators::build_C4cD(const OperatorLookup& operator_lookup, 
                                   const CorrelatorLookup& corr_lookup,
                                   const QuarklineLookup& quark_lookup) {

  std::cout << "\tcomputing C4cD:";
  clock_t time = clock();

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
    write_4pt_correlators(correlator, c_look);
  }

  time = clock() - time;
  std::cout << "\t\tSUCCESS - " << ((float) time) / CLOCKS_PER_SEC 
            << " seconds" << std::endl;
}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void LapH::Correlators::build_C4cV(const OperatorLookup& operator_lookup, 
                                   const CorrelatorLookup& corr_lookup,
                                   const QuarklineLookup& quark_lookup) {

  std::cout << "\tcomputing C4cV:";
  clock_t time = clock();

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
    write_4pt_correlators(correlator, c_look);
  }

  time = clock() - time;
  std::cout << "\t\tSUCCESS - " << ((float) time) / CLOCKS_PER_SEC 
            << " seconds" << std::endl;
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void LapH::Correlators::build_C4cC(const OperatorsForMesons& meson_operator,
                                   const Perambulator& perambulators,
                                   const OperatorLookup& operator_lookup,
                                   const std::vector<CorrInfo>& corr_lookup,
                                   const QuarklineLookup& quark_lookup) {

  std::cout << "\tcomputing C4cC:";
  clock_t time = clock();
  
  std::vector<vec> correlator(corr_lookup.size(), vec(Lt, cmplx(.0,.0)));

// This is necessary to ensure the correct summation of the correlation function
#pragma omp parallel
{
  std::vector<vec> C(corr_lookup.size(), vec(Lt, cmplx(.0,.0)));
  // building the quark line directly frees up a lot of memory
  Quarklines_one_t quarklines(2*dilT, dilT, dilE, nev, quark_lookup, 
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
    quarklines.build_Q2V_one_t(perambulators, meson_operator, t1_i, t2_i,
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
            quarklines.return_Q2V(id_Q2V_1, 0, look[1], idr0).
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
            quarklines.return_Q2V(id_Q2V_2, 0, look[1], idr2).
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
}// parallel part ends here

  // normalisation
  for(const auto& c_look : corr_lookup){
    for(auto& corr : correlator[c_look.id])
      corr /= (6*5*4*3)*Lt; // TODO: Hard Coded atm - Be carefull
    // write data to file
    write_correlators(correlator[c_look.id], c_look);
  }

  time = clock() - time;
  std::cout << "\t\t\tSUCCESS - " << ((float) time) / CLOCKS_PER_SEC 
            << " seconds" << std::endl;
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
                                  const QuarklineLookup& quark_lookup) {
  if(corr_lookup.size() == 0)
    return;

  std::cout << "\tcomputing C3c:";
  clock_t time = clock();

  std::vector<vec> correlator(corr_lookup.size(), vec(Lt, cmplx(.0,.0)));

// This is necessary to ensure the correct summation of the correlation function
#pragma omp parallel
{
  std::vector<vec> C(corr_lookup.size(), vec(Lt, cmplx(.0,.0)));
  // building the quark line directly frees up a lot of memory
  Quarklines_one_t quarklines(2*dilT, dilT, dilE, nev, quark_lookup, 
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
    quarklines.build_Q2L_one_t(perambulators, meson_operator, t1_i, t2_i,
                               quark_lookup.Q2L, operator_lookup.ricQ2_lookup);
    quarklines.build_Q1_mult_t(perambulators, meson_operator, t1_i, t2_i,
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
  
          const cmplx value = quarklines.return_gamma_val(5, col); // TODO: gamma hardcoded
          const size_t gamma_index = quarklines.return_gamma_row(
                                                          5, col); // TODO: gamma hardcoded

          M1[look[0]][M1_rnd_counter].block(col*dilE, 0, dilE, 4*dilE) = value *
              meson_operator.return_rvdaggervr(look[2], t1, idr2).
                                 block(col*dilE, gamma_index*dilE, dilE, dilE) *
              quarklines.return_Q2L(id_Q2L_1, 0, look[1], idr0).
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
                    quarklines.return_Q1(id_Q2L_2, 0, c_look.lookup[1], idr1)).trace();
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
}// parallel part ends here


  // normalisation
  for(const auto& c_look : corr_lookup){
    for(auto& corr : correlator[c_look.id])
      //corr /= (3)*Lt; // TODO: Hard Coded atm - Be carefull
      corr /= (6*5*4)*Lt; // TODO: Hard Coded atm - Be carefull
    // write data to file
    write_correlators(correlator[c_look.id], c_look);
  }
  time = clock() - time;
  std::cout << "\t\t\tSUCCESS - " << ((float) time) / CLOCKS_PER_SEC 
            << " seconds" << std::endl;
}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void LapH::Correlators::build_C4cB(const OperatorsForMesons& meson_operator,
                                   const Perambulator& perambulators,
                                   const OperatorLookup& operator_lookup,
                                   const std::vector<CorrInfo>& corr_lookup,
                                   const QuarklineLookup& quark_lookup) {

  if(corr_lookup.size() == 0)
    return;

  std::cout << "\tcomputing C4cB:";
  clock_t time = clock();

  std::vector<vec> correlator(corr_lookup.size(), vec(Lt, cmplx(.0,.0)));

// This is necessary to ensure the correct summation of the correlation function
#pragma omp parallel
{
  std::vector<vec> C(corr_lookup.size(), vec(Lt, cmplx(.0,.0)));
  // building the quark line directly frees up a lot of memory
  Quarklines_one_t quarklines(2*dilT, dilT, dilE, nev, quark_lookup, 
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
    quarklines.build_Q2L_one_t(perambulators, meson_operator, t1_i, t2_i,
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
              quarklines.return_Q2L(id_Q2L_1, 0, look[2], idr0).
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
            quarklines.return_Q2L(id_Q2L_2, 0, look[2], idr2).
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
}// parallel part ends here


  // normalisation
  for(const auto& c_look : corr_lookup){
    for(auto& corr : correlator[c_look.id])
      corr /= (6*5*4*3)*Lt; // TODO: Hard Coded atm - Be carefull
    // write data to file
    write_correlators(correlator[c_look.id], c_look);
  }
  time = clock() - time;
  std::cout << "\t\t\tSUCCESS - " << ((float) time) / CLOCKS_PER_SEC 
            << " seconds" << std::endl;
}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void LapH::Correlators::build_C30(const Quarklines& quarklines,
                     const std::vector<CorrInfo>& corr_lookup,
                     const QuarklineLookup& quark_lookup,
                     const std::vector<RandomIndexCombinationsQ2>& ric_lookup) {

  std::cout << "\tcomputing C30:";
  clock_t time = clock();

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
    for(auto& corr : correlator)
      corr /= norm/Lt;
    // write data to file
    write_correlators(correlator, c_look);
  }
  time = clock() - time;
  std::cout << "\t\tSUCCESS - " << ((float) time) / CLOCKS_PER_SEC 
            << " seconds" << std::endl;
}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void LapH::Correlators::build_C40C(const Quarklines& quarklines,
                     const std::vector<CorrInfo>& corr_lookup,
                     const QuarklineLookup& quark_lookup,
                     const std::vector<RandomIndexCombinationsQ2>& ric_lookup) {

  std::cout << "\tcomputing C40C:";
  clock_t time = clock();

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
    for(auto& corr : correlator)
      corr /= norm/Lt;
    // write data to file
    write_correlators(correlator, c_look);
  }
  time = clock() - time;
  std::cout << "\t\tSUCCESS - " << ((float) time) / CLOCKS_PER_SEC 
            << " seconds" << std::endl;
}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void LapH::Correlators::build_C40B(const Quarklines& quarklines,
                     const std::vector<CorrInfo>& corr_lookup,
                     const QuarklineLookup& quark_lookup,
                     const std::vector<RandomIndexCombinationsQ2>& ric_lookup) {

  std::cout << "\tcomputing C40B:";
  clock_t time = clock();

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
    for(auto& corr : correlator)
      corr /= norm/Lt;
    // write data to file
    write_correlators(correlator, c_look);
  }
  time = clock() - time;
  std::cout << "\t\tSUCCESS - " << ((float) time) / CLOCKS_PER_SEC 
            << " seconds" << std::endl;
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
                     const QuarklineLookup& quark_lookup) {

  // 1. Build all functions which need corrC and free it afterwards.
  build_corrC(perambulators, meson_operator, operator_lookup, 
              corr_lookup.corrC, quark_lookup);
  build_C2c(corr_lookup.C2c);
  build_C4cD(operator_lookup, corr_lookup, quark_lookup);
  build_C4cV(operator_lookup, corr_lookup, quark_lookup);
  // 2. Build all functions which need corr0 and free it afterwards.
  build_corr0(meson_operator, perambulators, corr_lookup.corr0, 
              quark_lookup, operator_lookup);
  // in C3c, also corr0 is build, since this is much faster
  build_C3c(meson_operator, perambulators, operator_lookup, corr_lookup.C3c, 
                                                                 quark_lookup);
  build_C20(corr_lookup.C20);
  build_C40D(operator_lookup, corr_lookup, quark_lookup);
  build_C40V(operator_lookup, corr_lookup, quark_lookup);
  // 3. Build all other correlation functions.
  build_C1(quarklines, corr_lookup.C1, quark_lookup, 
                                                 operator_lookup.ricQ2_lookup);
  build_C4cC(meson_operator, perambulators, operator_lookup, corr_lookup.C4cC, 
                                                                 quark_lookup);
//  build_C4cC(quarklines, meson_operator, operator_lookup, corr_lookup.C4cC, 
//                                                                 quark_lookup);
  build_C4cB(meson_operator, perambulators, operator_lookup, corr_lookup.C4cB, 
                                                                 quark_lookup);
  build_C30(quarklines, corr_lookup.C30, quark_lookup, 
                                                 operator_lookup.ricQ2_lookup);
  build_C40C(quarklines, corr_lookup.C40C, quark_lookup, 
                                                 operator_lookup.ricQ2_lookup);
  build_C40B(quarklines, corr_lookup.C40B, quark_lookup, 
                                                 operator_lookup.ricQ2_lookup);
}






