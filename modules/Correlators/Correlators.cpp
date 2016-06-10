#include "Correlators.h"

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
static void write_correlators(const std::vector<cmplx>& corr, 
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
  // writing the data
  std::ofstream file((corr_info.outpath+corr_info.outfile).c_str(), 
                     std::ios::out | std::ofstream::binary | std::ios::trunc);
  if(file.is_open()){
    file.write(reinterpret_cast<const char*>(&corr[0]), 
               corr.size()*sizeof(cmplx));
    file.close();
  }
  else
    std::cout << "can't open " << (corr_info.outpath+corr_info.outfile).c_str() 
              << std::endl;
}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------





// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
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
  Quarklines quarklines_intern(2, dilT, dilE, nev, quark_lookup, 
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
      corr /= Lt*corr0[c_look.id][0][0].size();
    // write data to file
    write_correlators(correlator, c_look);
  }
}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void LapH::Correlators::build_C40D(const OperatorLookup& operator_lookup, 
                                   const CorrelatorLookup& corr_lookup,
                                   const QuarklineLookup& quark_lookup) {

  for(const auto& c_look : corr_lookup.C40D){
    std::vector<cmplx> correlator(Lt, cmplx(.0,.0));
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

        correlator[t] += 
                   corr0[c_look.lookup[0]][t1][t2].at(&rnd0 - &ric0[0]) *
                   corr0[c_look.lookup[0]][t1][t2].at(&rnd1 - &ric1[0]);
        norm++;
      }
    }}
    // normalisation
    for(auto& corr : correlator)
      corr /= norm/Lt;
    // write data to file
    write_correlators(correlator, c_look);
  }
}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void LapH::Correlators::build_C40V(const OperatorLookup& operator_lookup, 
                                   const CorrelatorLookup& corr_lookup,
                                   const QuarklineLookup& quark_lookup) {

  for(const auto& c_look : corr_lookup.C40V){
    std::vector<cmplx> correlator(Lt, cmplx(.0,.0));
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

        correlator[t] += 
                   corr0[c_look.lookup[0]][t1][t1].at(&rnd0 - &ric0[0]) *
                   corr0[c_look.lookup[0]][t2][t2].at(&rnd1 - &ric1[0]);
        norm++;
      }
    }}
    // normalisation
    for(auto& corr : correlator)
      corr /= norm/Lt;
    // write data to file
    write_correlators(correlator, c_look);
  }
}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------





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
  Quarklines quarklines(2*dilT, dilT, dilE, nev, quark_lookup, 
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
}// omp paralle ends here
  time = clock() - time;
  std::cout << "\t\tSUCCESS - " << ((float) time) / CLOCKS_PER_SEC 
            << " seconds" << std::endl;
}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void LapH::Correlators::build_C2c(const std::vector<CorrInfo>& corr_lookup) {

  for(const auto& c_look : corr_lookup){
    std::vector<cmplx> correlator(Lt, cmplx(.0,.0));
    for(int t1 = 0; t1 < Lt; t1++){
    for(int t2 = 0; t2 < Lt; t2++){
      int t = abs((t2 - t1 - (int)Lt) % (int)Lt);
      for(const auto& corr : corrC[c_look.lookup[0]][t1][t2]){
        correlator[t] += corr;
      }
    }}
    // normalisation
    for(auto& corr : correlator)
      corr /= Lt*corrC[c_look.id][0][0].size();
    // write data to file
    write_correlators(correlator, c_look);
  }
}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void LapH::Correlators::build_C4cD(const OperatorLookup& operator_lookup, 
                                   const CorrelatorLookup& corr_lookup,
                                   const QuarklineLookup& quark_lookup) {

  for(const auto& c_look : corr_lookup.C4cD){
    std::vector<cmplx> correlator(Lt, cmplx(.0,.0));
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

        correlator[t] += 
                   corrC[c_look.lookup[0]][t1][t2].at(&rnd0 - &ric0[0]) *
                   corrC[c_look.lookup[0]][t1][t2].at(&rnd1 - &ric1[0]);
        norm++;
      }
    }}
    // normalisation
    for(auto& corr : correlator)
      corr /= norm/Lt;
    // write data to file
    write_correlators(correlator, c_look);
  }
}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void LapH::Correlators::build_C4cV(const OperatorLookup& operator_lookup, 
                                   const CorrelatorLookup& corr_lookup,
                                   const QuarklineLookup& quark_lookup) {

  for(const auto& c_look : corr_lookup.C4cV){
    std::vector<cmplx> correlator(Lt, cmplx(.0,.0));
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

        correlator[t] += 
                   corrC[c_look.lookup[0]][t1][t1].at(&rnd0 - &ric0[0]) *
                   corrC[c_look.lookup[0]][t2][t2].at(&rnd1 - &ric1[0]);
        norm++;
      }
    }}
    // normalisation
    for(auto& corr : correlator)
      corr /= norm/Lt;
    // write data to file
    write_correlators(correlator, c_look);
  }
}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------






// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void LapH::Correlators::build_C4cC(const Quarklines& quarklines, 
                                   const OperatorsForMesons& meson_operator,
                                   const OperatorLookup& operator_lookup,
                                   const std::vector<CorrInfo>& corr_lookup,
                                   const QuarklineLookup& quark_lookup) {

  std::cout << "\tcomputing C4cC:";
  clock_t time = clock();

  for(const auto& c_look : corr_lookup){
    std::vector<cmplx> correlator(Lt, cmplx(.0,.0));
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
      std::cout << "rnd combinations are not the same in build_corr0" 
                << std::endl;
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
      if(rnd0.first != rnd1.first && rnd0.second == rnd1.second){
        const size_t id0 = &rnd0 - &ric0[0];
        const size_t id1 = &rnd1 - &ric1[0];
        Eigen::MatrixXcd M1 = Eigen::MatrixXcd::Zero(4 * dilE, 4 * dilE);
        for(size_t col = 0; col < 4; col++){

          const cmplx value = quarklines.return_gamma_val(c_look.gamma[0], col);
          const size_t gamma_index = quarklines.return_gamma_row(
                                                          c_look.gamma[0], col);
          M1.block(0, col*dilE, 4*dilE, dilE) = value *
            quarklines.return_Q2V(t1, t2/dilT, c_look.lookup[0], id0).
                               block(0, gamma_index*dilE, 4*dilE, dilE) *
            meson_operator.return_rvdaggervr(c_look.lookup[1], t2, id1).
                                block(gamma_index*dilE, col*dilE, dilE, dilE);
        }
        for(const auto& rnd2 : ric2){
        for(const auto& rnd3 : ric3){
        if(rnd2.first != rnd3.first && rnd2.second == rnd3.second && 
           rnd1.first == rnd2.first && rnd0.first == rnd3.first &&
           rnd0.second != rnd3.second){
          const size_t id2 = &rnd2 - &ric2[0];
          const size_t id3 = &rnd3 - &ric3[0];
          Eigen::MatrixXcd M2 = Eigen::MatrixXcd::Zero(4 * dilE, 4 * dilE);
          for(size_t col = 0; col < 4; col++){
            const cmplx value = quarklines.return_gamma_val(c_look.gamma[1], col);
            const size_t gamma_index = quarklines.return_gamma_row(
                                                            c_look.gamma[1], col);
            M2.block(0, col*dilE, 4*dilE, dilE) = value *
              quarklines.return_Q2V(t1, t2/dilT, c_look.lookup[2], id2).
                                 block(0, gamma_index*dilE, 4*dilE, dilE) *
              meson_operator.return_rvdaggervr(c_look.lookup[3], t2, id3).
                                  block(gamma_index*dilE, col*dilE, dilE, dilE);
          }
          C[t] += (M1*M2).trace();
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
  std::cout << "\t\t\tSUCCESS - " << ((float) time) / CLOCKS_PER_SEC 
            << " seconds" << std::endl;
}
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
  Quarklines quarklines(2*dilT, dilT, dilE, nev, quark_lookup, 
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
         id_Q2L_2 = (dir+1)%dilT*dilT + t2%t2_min;                           
       else                                                                     
         id_Q2L_2 = (dir+1)%dilT*dilT + t2;               
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
  Quarklines quarklines(2*dilT, dilT, dilE, nev, quark_lookup, 
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
         id_Q2L_2 = (bla+1)%dilT*dilT + t2%t2_min;                           
       else                                                                     
         id_Q2L_2 = (bla+1)%dilT*dilT + t2;               
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
      std::cout << "rnd combinations are not the same in build_corr0" 
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
    if(ric0.size() != ric1.size() || ric0.size() != ric2.size() || 
       ric0.size() != ric3.size()){
      std::cout << "rnd combinations are not the same in build_corr0" 
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
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------





// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
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
  build_C20(corr_lookup.C20);
  build_C40D(operator_lookup, corr_lookup, quark_lookup);
  build_C40V(operator_lookup, corr_lookup, quark_lookup);
  // 3. Build all other correlation functions.
  build_C1(quarklines, corr_lookup.C1, quark_lookup, 
                                                 operator_lookup.ricQ2_lookup);
  build_C3c(meson_operator, perambulators, operator_lookup, corr_lookup.C3c, 
                                                                 quark_lookup);
  build_C4cC(quarklines, meson_operator, operator_lookup, corr_lookup.C4cC, 
                                                                 quark_lookup);
  build_C4cB(meson_operator, perambulators, operator_lookup, corr_lookup.C4cB, 
                                                                 quark_lookup);
  build_C30(quarklines, corr_lookup.C30, quark_lookup, 
                                                 operator_lookup.ricQ2_lookup);
  build_C40C(quarklines, corr_lookup.C40C, quark_lookup, 
                                                 operator_lookup.ricQ2_lookup);
  build_C40B(quarklines, corr_lookup.C40B, quark_lookup, 
                                                 operator_lookup.ricQ2_lookup);
}






