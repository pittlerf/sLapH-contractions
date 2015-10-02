#include "Quarklines.h"

LapH::Quarklines::Quarklines(
                     const size_t Lt, const size_t dilT, const size_t dilE, 
                     const QuarklineLookup& quarkline_lookuptable,
                     const std::vector<RandomIndexCombinationsQ2>& ric_lookup) :
                                    Lt(Lt), dilT(dilT), dilE(dilE) {

  // TODO: just for testing
  size_t nb_rnd = 1;

  Q1.resize(boost::extents[Lt][Lt/dilT][quarkline_lookuptable.Q1.size()]);
  Q2V.resize(boost::extents[Lt][Lt/dilT][quarkline_lookuptable.Q2V.size()]);
  Q2L.resize(boost::extents[Lt][Lt/dilT][quarkline_lookuptable.Q2L.size()]);
  for(size_t t1 = 0; t1 < Lt; t1++){                  
  for(size_t t2 = 0; t2 < Lt/dilT; t2++){
    for(size_t op = 0; op < quarkline_lookuptable.Q1.size(); op++){ 
      Q1[t1][t2][op].resize(nb_rnd);                                          
      for(size_t rnd1 = 0; rnd1 < nb_rnd; rnd1++){                                
        Q1[t1][t2][op][rnd1] = Eigen::MatrixXcd::Zero(4*dilE, 4*dilE); 
      }                                                                           
    }
    for(size_t op = 0; op < quarkline_lookuptable.Q2V.size(); op++){ 
      Q2V[t1][t2][op].resize(nb_rnd);                                          
      for(size_t rnd1 = 0; rnd1 < nb_rnd; rnd1++){                                
        Q2V[t1][t2][op][rnd1] = Eigen::MatrixXcd::Zero(4*dilE, 4*dilE); 
      }                                                                           
    }               
    for(size_t op = 0; op < quarkline_lookuptable.Q2L.size(); op++){ 
      Q2L[t1][t2][op].resize(nb_rnd);                                          
      for(size_t rnd1 = 0; rnd1 < nb_rnd; rnd1++){                                
        Q2L[t1][t2][op][rnd1] = Eigen::MatrixXcd::Zero(4*dilE, 4*dilE); 
      }                                                                           
    }               
  }}               
}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void LapH::Quarklines::build_Q1(const Perambulator& peram,
              const OperatorsForMesons& meson_operator,
              const std::vector<QuarklineQ1Indices>& lookup) {
}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void LapH::Quarklines::build_Q2V(const Perambulator& peram,
               const OperatorsForMesons& meson_operator,
               const std::vector<QuarklineQ2Indices>& lookup){
}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void LapH::Quarklines::build_Q2L(const Perambulator& peram,
               const OperatorsForMesons& meson_operator,
               const std::vector<QuarklineQ2Indices>& lookup){
}
// ------------------------------- INTERFACE -----------------------------------
// -----------------------------------------------------------------------------
void LapH::Quarklines::create_quarklines(const Perambulator& peram, 
                                 const OperatorsForMesons& meson_operator,
                                 const QuarklineLookup& quarkline_lookuptable) {

  build_Q1(peram, meson_operator, quarkline_lookuptable.Q1);
  build_Q2L(peram, meson_operator, quarkline_lookuptable.Q2L);
  build_Q2V(peram, meson_operator, quarkline_lookuptable.Q2V);
}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------






