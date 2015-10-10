#include "Correlators.h"

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
LapH::Correlators::Correlators (
                        const size_t Lt, const size_t dilT, const size_t dilE, 
                        const size_t nev, const CorrelatorLookup& corr_lookup) :
                                      Lt(Lt), dilT(dilT), dilE(dilE), nev(nev) {      

  C20.resize(boost::extents[corr_lookup.C20.size()][Lt]);
  std::fill(C20.data(), C20.data()+C20.num_elements(), cmplx(.0,.0));
  C2c.resize(boost::extents[corr_lookup.C2c.size()][Lt]);
  std::fill(C2c.data(), C2c.data()+C2c.num_elements(), cmplx(.0,.0));

}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void LapH::Correlators::build_corr0(const Quarklines& quarklines,
                     const std::vector<CorrInfo>& corr_lookup,
                     const QuarklineLookup& quark_lookup,
                     const std::vector<RandomIndexCombinationsQ2>& ric_lookup) {

  std::cout << "\tcomputing corr0:";
  clock_t time = clock();

  corr0.resize(boost::extents[corr_lookup.size()][Lt][Lt]);
  for(const auto& c_look : corr_lookup){
    const auto& ric0 = ric_lookup[quark_lookup.Q1[c_look.lookup[0]].
                                                     id_ric_lookup].rnd_vec_ids;
    const auto& ric1 = ric_lookup[quark_lookup.Q1[c_look.lookup[1]].
                                                     id_ric_lookup].rnd_vec_ids;
    if(ric0.size() != ric1.size()){
      std::cout << "rnd combinations are not the same in build_corr0" 
                << std::endl;
      exit(0);
    }
    for(size_t t1 = 0; t1 < Lt; t1++){
    for(size_t t2 = 0; t2 < Lt; t2++){
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
        //for(size_t block = 0; block < 4; block++)
          corr0[c_look.id][t1][t2][id] += 
                      (quarklines.return_Q1(t1, t2/dilT, c_look.lookup[0], id) *
                       quarklines.return_Q1(t2, t1/dilT, c_look.lookup[1], 
                                                     it1-ric1.begin())).trace();
      }
    }}
  }
  time = clock() - time;
  std::cout << "\t\tSUCCESS - " << ((float) time) / CLOCKS_PER_SEC 
            << " seconds" << std::endl;
}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void LapH::Correlators::build_corrC() {

  std::cout << "\tcomputing corrC:";
  clock_t time = clock();




  time = clock() - time;
  std::cout << "\t\tSUCCESS - " << ((float) time) / CLOCKS_PER_SEC 
            << " seconds" << std::endl;
}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void LapH::Correlators::build_C20() {


  for(int t1 = 0; t1 < Lt; t1++){
  for(int t2 = 0; t2 < Lt; t2++){
    int t = abs((t1 - t2 - (int)Lt) % (int)Lt);
    for(const auto& corr : corr0[0][t1][t2])
      C20[0][t] += corr;
  }}
  for(auto& corr : C20[0]){
    // normalisation
    corr /= Lt*corr0[0][0][0].size();
    std::cout << std::setprecision(5) << &corr - &C20[0][0] << "\t" << corr 
              << std::endl;
  }
}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void LapH::Correlators::build_C2c() {
}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void LapH::Correlators::contract (const Quarklines& quarklines, 
                     const OperatorsForMesons& meson_operator,
                     const CorrelatorLookup& corr_lookup, 
                     const QuarklineLookup& quark_lookup,
                     const std::vector<RandomIndexCombinationsQ2>& ric_lookup) {

  // 1. Build all functions which need corrC and free it afterwards.
  build_corrC();
  build_C2c();
//  corrC.resize(boost::extents[0][0][0][0][0][0]);
  // 2. Build all functions which need corr0 and free it afterwards.
  build_corr0(quarklines, corr_lookup.corr0, quark_lookup, ric_lookup);
  build_C20();
//  corr0.resize(boost::extents[0][0][0][0][0][0]);
  // 3. Build all other correlation functions.
}






