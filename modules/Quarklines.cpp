#include "Quarklines.h"

namespace { // some internal namespace

static const std::complex<double> I(0.0, 1.0);
// Look-up table for gamma matrices. For every Gamma structure (currently 0-15)
// the four non-zero values are specified.
static void create_gamma (std::vector<LapH::gamma_lookup>& gamma, const int i) {
  try {
    switch(i) {
    case 0: // gamma_0
      gamma[0].row[0] = 2;
      gamma[0].value[0] = 1;
      gamma[0].row[1] = 3;
      gamma[0].value[1] = 1;
      gamma[0].row[2] = 0;
      gamma[0].value[2] = 1;
      gamma[0].row[3] = 1;
      gamma[0].value[3] = 1;
      break;

    case 1: // gamma_1
      gamma[1].row[0] = 3;
      gamma[1].value[0] = -I; 
      gamma[1].row[1] = 2;
      gamma[1].value[1] = -I; 
      gamma[1].row[2] = 1;
      gamma[1].value[2] = I;
      gamma[1].row[3] = 0;
      gamma[1].value[3] = I;
      break;

    case 2: // gamma_2
      gamma[2].row[0] = 3;
      gamma[2].value[0] = -1; 
      gamma[2].row[1] = 2;
      gamma[2].value[1] = 1;
      gamma[2].row[2] = 1;
      gamma[2].value[2] = 1;
      gamma[2].row[3] = 0;
      gamma[2].value[3] = -1; 
      break;

    case 3: // gamma_3
      gamma[3].row[0] = 2;
      gamma[3].value[0] = -I; 
      gamma[3].row[1] = 3;
      gamma[3].value[1] = I;
      gamma[3].row[2] = 0;
      gamma[3].value[2] = I;
      gamma[3].row[3] = 1;
      gamma[3].value[3] = -I; 
      break;

    case 4: // unity
      gamma[4].row[0] = 0;
      gamma[4].value[0] = 1;
      gamma[4].row[1] = 1;
      gamma[4].value[1] = 1;
      gamma[4].row[2] = 2;
      gamma[4].value[2] = 1;
      gamma[4].row[3] = 3;
      gamma[4].value[3] = 1;
      break;

    case 5: // gamma_5
      gamma[5].row[0] = 0;
      gamma[5].value[0] = 1;
      gamma[5].row[1] = 1;
      gamma[5].value[1] = 1;
      gamma[5].row[2] = 2;
      gamma[5].value[2] = -1;
      gamma[5].row[3] = 3;
      gamma[5].value[3] = -1;
      break;

    case 6: // gamma_0 * gamma_5
      gamma[6].row[0] = 2;
      gamma[6].value[0] = -1;
      gamma[6].row[1] = 3;
      gamma[6].value[1] = -1;
      gamma[6].row[2] = 0;
      gamma[6].value[2] = 1;
      gamma[6].row[3] = 1;
      gamma[6].value[3] = 1;
      break;

    case 7: // gamma_1 * gamma_5
      gamma[7].row[0] = 3;
      gamma[7].value[0] = I;
      gamma[7].row[1] = 2;
      gamma[7].value[1] = I;
      gamma[7].row[2] = 1;
      gamma[7].value[2] = I;
      gamma[7].row[3] = 0;
      gamma[7].value[3] = I;
      break;

    case 8: // gamma_2 * gamma_5
      gamma[8].row[0] = 3;
      gamma[8].value[0] = 1;
      gamma[8].row[1] = 2;
      gamma[8].value[1] = -1;
      gamma[8].row[2] = 1;
      gamma[8].value[2] = 1;
      gamma[8].row[3] = 0;
      gamma[8].value[3] = -1;
      break;

    case 9: // gamma_3 * gamma_5
      gamma[9].row[0] = 2;
      gamma[9].value[0] = I;
      gamma[9].row[1] = 3;
      gamma[9].value[1] = -I;
      gamma[9].row[2] = 0;
      gamma[9].value[2] = I;
      gamma[9].row[3] = 1;
      gamma[9].value[3] = -I;
      break;

    case 10: // gamma_0 * gamma_1
      gamma[10].row[0] = 1;
      gamma[10].value[0] = I;
      gamma[10].row[1] = 0;
      gamma[10].value[1] = I;
      gamma[10].row[2] = 3;
      gamma[10].value[2] = -I;
      gamma[10].row[3] = 2;
      gamma[10].value[3] = -I;
      break;

    case 11: // gamma_0 * gamma_2
      gamma[11].row[0] = 1;
      gamma[11].value[0] = 1;
      gamma[11].row[1] = 0;
      gamma[11].value[1] = -1;
      gamma[11].row[2] = 3;
      gamma[11].value[2] = -1;
      gamma[11].row[3] = 2;
      gamma[11].value[3] = 1;
      break;

    case 12: // gamma_0 * gamma_3
      gamma[12].row[0] = 0;
      gamma[12].value[0] = I;
      gamma[12].row[1] = 1;
      gamma[12].value[1] = -I;
      gamma[12].row[2] = 2;
      gamma[12].value[2] = -I;
      gamma[12].row[3] = 3;
      gamma[12].value[3] = I;
      break;

    case 13: // gamma_0 * gamma_1 * gamma_5
      gamma[13].row[0] = 1;
      gamma[13].value[0] = I;
      gamma[13].row[1] = 0;
      gamma[13].value[1] = I;
      gamma[13].row[2] = 3;
      gamma[13].value[2] = I;
      gamma[13].row[3] = 2;
      gamma[13].value[3] = I;
      break;

    case 14: // gamma_0 * gamma_2 * gamma_5
      gamma[14].row[0] = 1;
      gamma[14].value[0] = 1;
      gamma[14].row[1] = 0;
      gamma[14].value[1] = -1;
      gamma[14].row[2] = 3;
      gamma[14].value[2] = 1;
      gamma[14].row[3] = 2;
      gamma[14].value[3] = -1;
      break;

    case 15: // gamma_0 * gamma_3 * gamma_5
      gamma[15].row[0] = 0;
      gamma[15].value[0] = I;
      gamma[15].row[1] = 1;
      gamma[15].value[1] = -I;
      gamma[15].row[2] = 2;
      gamma[15].value[2] = I;
      gamma[15].row[3] = 3;
      gamma[15].value[3] = -I;
      break;

    default:
      printf("Dirac component %d not found in BasicOperator::create_gamma\n",i);
      exit(0);
    }
  return;
  }
  catch(std::exception& e){
    std::cout << e.what() << "in: BasicOperator::create_gamma\n";
    exit(0);
  }
}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
} // internal namespace ends here

LapH::Quarklines::Quarklines(
                     const size_t Lt, const size_t dilT, 
                     const size_t dilE, const size_t nev, 
                     const QuarklineLookup& quarkline_lookuptable,
                     const std::vector<RandomIndexCombinationsQ2>& ric_lookup) :
                                    Lt(Lt), dilT(dilT), dilE(dilE), nev(nev) {

  Q1.resize(boost::extents[Lt][Lt/dilT][quarkline_lookuptable.Q1.size()]);
  Q2V.resize(boost::extents[Lt][Lt/dilT][quarkline_lookuptable.Q2V.size()]);
  Q2L.resize(boost::extents[Lt][Lt/dilT][quarkline_lookuptable.Q2L.size()]);
  for(size_t t1 = 0; t1 < Lt; t1++){                  
  for(size_t t2 = 0; t2 < Lt/dilT; t2++){
    for(size_t op = 0; op < quarkline_lookuptable.Q1.size(); op++){ 
      size_t nb_rnd = ric_lookup[(quarkline_lookuptable.Q1[op]).
                                 id_ric_lookup].rnd_vec_ids.size();
      Q1[t1][t2][op].resize(nb_rnd);                                          
      for(size_t rnd1 = 0; rnd1 < nb_rnd; rnd1++){
        Q1[t1][t2][op][rnd1] = Eigen::MatrixXcd::Zero(4*dilE, 4*dilE); 
      } 
    }
    for(size_t op = 0; op < quarkline_lookuptable.Q2V.size(); op++){ 
      size_t nb_rnd = ric_lookup[(quarkline_lookuptable.Q2V[op]).
                                 id_ric_lookup].rnd_vec_ids.size();
      Q2V[t1][t2][op].resize(nb_rnd);                                          
      for(size_t rnd1 = 0; rnd1 < nb_rnd; rnd1++){
        Q2V[t1][t2][op][rnd1] = Eigen::MatrixXcd::Zero(4*dilE, 4*dilE); 
      } 
    }               
    for(size_t op = 0; op < quarkline_lookuptable.Q2L.size(); op++){ 
      size_t nb_rnd = ric_lookup[(quarkline_lookuptable.Q2L[op]).
                                 id_ric_lookup].rnd_vec_ids.size();
      Q2L[t1][t2][op].resize(nb_rnd);                                          
      for(size_t rnd1 = 0; rnd1 < nb_rnd; rnd1++){
        Q2L[t1][t2][op][rnd1] = Eigen::MatrixXcd::Zero(4*dilE, 4*dilE); 
      } 
    }               
  }}              
  // creating gamma matrices
  gamma.resize(16);
  for(int i = 0; i < 16; ++i)
    create_gamma(gamma, i);
}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void LapH::Quarklines::build_Q1(const Perambulator& peram,
              const OperatorsForMesons& meson_operator,
              const std::vector<QuarklineQ1Indices>& ql_lookup,
              const std::vector<RandomIndexCombinationsQ2>& ric_lookup){

  std::cout << "\tcomputing Q1:";
  clock_t time = clock();

// setting memory to zero
#pragma omp parallel for schedule(dynamic)
  for(size_t t1 = 0; t1 < Lt; t1++){                  
  for(size_t t2 = 0; t2 < Lt/dilT; t2++){
    for(size_t op = 0; op < ql_lookup.size(); op++){ 
      size_t nb_rnd = ric_lookup[(ql_lookup[op]).
                                 id_ric_lookup].rnd_vec_ids.size();
      for(size_t rnd1 = 0; rnd1 < nb_rnd; rnd1++){
        Q1[t1][t2][op][rnd1].setZero(); 
      } 
    }
  }}
#pragma omp parallel for schedule(dynamic)
  for(size_t t1 = 0; t1 < Lt; t1++){                  
  for(size_t t2 = 0; t2 < Lt/dilT; t2++){
    for(const auto& qll : ql_lookup){
      const size_t offset = ric_lookup[qll.id_ric_lookup].offset.first;
      size_t rnd_counter = 0;
      for(const auto& rnd_id : ric_lookup[qll.id_ric_lookup].rnd_vec_ids){
        const size_t rid1 = rnd_id.first - offset; 
        const size_t gamma_id = qll.gamma[0]; // TODO: hard coded! VERY BAD!!!
        for(size_t row = 0; row < 4; row++){
        for(size_t col = 0; col < 4; col++){
          Q1[t1][t2][qll.id][rnd_counter].block(row*dilE, col*dilE, dilE, dilE)=
            gamma[gamma_id].value[row] *  
            meson_operator.return_rvdaggerv(qll.id_rvdaggerv, t1, rid1).
                                                  block(row*dilE, 0, dilE, nev)*
            peram[rnd_id.second].block((t1*4 + gamma[gamma_id].row[row])*nev, 
                                       (t2*4 + col)*dilE, 
                                       nev, dilE);
        }}
        rnd_counter++;

      }
    }
  }}

  time = clock() - time;
  std::cout << "\t\t\tSUCCESS - " << ((float) time) / CLOCKS_PER_SEC 
            << " seconds" << std::endl;
}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void LapH::Quarklines::build_Q2V(const Perambulator& peram,
               const OperatorsForMesons& meson_operator,
               const std::vector<QuarklineQ2Indices>& ql_lookup,
               const std::vector<RandomIndexCombinationsQ2>& ric_lookup){

  std::cout << "\tcomputing Q2V:";
  clock_t time = clock();
#pragma omp parallel 
  {
  Eigen::MatrixXcd M = Eigen::MatrixXcd::Zero(4 * dilE, 4 * nev);

// setting memory to zero
#pragma omp for schedule(dynamic)
  for(size_t t1 = 0; t1 < Lt; t1++){                  
  for(size_t t2 = 0; t2 < Lt/dilT; t2++){
    for(size_t op = 0; op < ql_lookup.size(); op++){ 
      size_t nb_rnd = ric_lookup[(ql_lookup[op]).
                                 id_ric_lookup].rnd_vec_ids.size();
      for(size_t rnd1 = 0; rnd1 < nb_rnd; rnd1++){
        Q2V[t1][t2][op][rnd1].setZero(); 
      } 
    }
  }}               


#pragma omp for schedule(dynamic)
  for(size_t t1 = 0; t1 < Lt; t1++){                  
  for(size_t t2 = 0; t2 < Lt/dilT; t2++){

    for(const auto& qll : ql_lookup){
      size_t rnd_counter = 0;
      int check = -1;
      for(const auto& rnd_id : ric_lookup[qll.id_ric_lookup].rnd_vec_ids){

        if(check != rnd_id.first){ // this avoids recomputation
          for(size_t row = 0; row < 4; row++){
          for(size_t col = 0; col < 4; col++){
            if(!qll.need_vdaggerv_dag)
              M.block(col*dilE, row*nev, dilE, nev) =
                peram[rnd_id.first].block((t1*4 + row)*nev, (t2*4 + col)*dilE, 
                                          nev, dilE).adjoint() *
                meson_operator.return_vdaggerv(qll.id_vdaggerv, t1);
            else
              M.block(col*dilE, row*nev, dilE, nev) =
                peram[rnd_id.first].block((t1*4 + row)*nev, (t2*4 + col)*dilE, 
                                          nev, dilE).adjoint() *
                meson_operator.return_vdaggerv(qll.id_vdaggerv, t1).adjoint();
            // gamma_5 trick
            if( ((row + col) == 3) || (abs(row - col) > 1) )
              M.block(col*dilE, row*nev, dilE, nev) *= -1.;
          }}
        }
        Q2V[t1][t2][qll.id][rnd_counter].Zero(4*dilE, 4*dilE);

        const size_t gamma_id = qll.gamma[0]; 

        for(size_t block_dil = 0; block_dil < 4; block_dil++) {
          const cmplx value = gamma[gamma_id].value[block_dil];
          const size_t gamma_index = gamma[gamma_id].row[block_dil];
          for(size_t row = 0; row < 4; row++){
          for(size_t col = 0; col < 4; col++){

          Q2V[t1][t2][qll.id][rnd_counter].
                        block(row*dilE, col*dilE, dilE, dilE) +=
               value * 
               M.block(row*dilE, block_dil*nev, dilE, nev) *
               peram[rnd_id.second].block(
                          (t1*4 + gamma_index)*nev, 
                          (t2*4 + col)*dilE, nev, dilE);

          }}
        }
        check = rnd_id.first;
        rnd_counter++;
      }
    }
  }}
} // pragma omp ends

  time = clock() - time;
  std::cout << "\t\t\tSUCCESS - " << ((float) time) / CLOCKS_PER_SEC 
            << " seconds" << std::endl;
}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void LapH::Quarklines::build_Q2L(const Perambulator& peram,
               const OperatorsForMesons& meson_operator,
               const std::vector<QuarklineQ2Indices>& ql_lookup,
               const std::vector<RandomIndexCombinationsQ2>& ric_lookup){

  std::cout << "\tcomputing Q2L:";
  clock_t time = clock();
#pragma omp parallel 
  {
  Eigen::MatrixXcd M = Eigen::MatrixXcd::Zero(4 * dilE, 4 * nev);

  for(size_t t1 = 0; t1 < Lt; t1++){                  
  for(size_t t2 = 0; t2 < Lt/dilT; t2++){
    for(size_t op = 0; op < ql_lookup.size(); op++){ 
      size_t nb_rnd = ric_lookup[(ql_lookup[op]).
                                 id_ric_lookup].rnd_vec_ids.size();
      for(size_t rnd1 = 0; rnd1 < nb_rnd; rnd1++){
        Q2L[t1][t2][op][rnd1].setZero(); 
      } 
    }               
  }}              

#pragma omp for schedule(dynamic)
  for(size_t t1 = 0; t1 < Lt; t1++){                  
  for(const auto& qll : ql_lookup){
    size_t rnd_counter = 0;
    int check = -1;
    for(const auto& rnd_id : ric_lookup[qll.id_ric_lookup].rnd_vec_ids){
      if(check != rnd_id.first){ // this avoids recomputation
        for(size_t row = 0; row < 4; row++){
        for(size_t col = 0; col < 4; col++){
          if(!qll.need_vdaggerv_dag)
            M.block(col*dilE, row*nev, dilE, nev) =
              peram[rnd_id.first].block((t1*4 + row)*nev, 
                                        (t1/dilT*4 + col)*dilE, 
                                        nev, dilE).adjoint() *
              meson_operator.return_vdaggerv(qll.id_vdaggerv, t1);
          else
            M.block(col*dilE, row*nev, dilE, nev) =
              peram[rnd_id.first].block((t1*4 + row)*nev, 
                                        (t1/dilT*4 + col)*dilE, 
                                        nev, dilE).adjoint() *
              meson_operator.return_vdaggerv(qll.id_vdaggerv, t1).adjoint();
          // gamma_5 trick
          if( ((row + col) == 3) || (abs(row - col) > 1) )
            M.block(col*dilE, row*nev, dilE, nev) *= -1.;
        }}
      }
      for(size_t t2 = 0; t2 < Lt/dilT; t2++){
        Q2L[t1][t2][qll.id][rnd_counter].Zero(4*dilE, 4*dilE);

        const size_t gamma_id = qll.gamma[0]; 

        for(size_t block_dil = 0; block_dil < 4; block_dil++) {
          const cmplx value = gamma[gamma_id].value[block_dil];
          const size_t gamma_index = gamma[gamma_id].row[block_dil];
          for(size_t row = 0; row < 4; row++){
          for(size_t col = 0; col < 4; col++){

          Q2L[t1][t2][qll.id][rnd_counter].
                        block(row*dilE, col*dilE, dilE, dilE) +=
               value * 
               M.block(row*dilE, block_dil*nev, dilE, nev) *
               peram[rnd_id.second].block(
                          (t1*4 + gamma_index)*nev, 
                          (t2*4 + col)*dilE, nev, dilE);

          }}
        }
      }
      check = rnd_id.first;
      rnd_counter++;
    }
  }}
} // pragma omp ends

  time = clock() - time;
  std::cout << "\t\t\tSUCCESS - " << ((float) time) / CLOCKS_PER_SEC 
            << " seconds" << std::endl;
}
// ------------------------------- INTERFACE -----------------------------------
// -----------------------------------------------------------------------------
void LapH::Quarklines::create_quarklines(const Perambulator& peram, 
                      const OperatorsForMesons& meson_operator,
                      const QuarklineLookup& quarkline_lookuptable,
                      const std::vector<RandomIndexCombinationsQ2>& ric_lookup){

  build_Q1(peram, meson_operator, quarkline_lookuptable.Q1, ric_lookup);
  build_Q2L(peram, meson_operator, quarkline_lookuptable.Q2L, ric_lookup);
  build_Q2V(peram, meson_operator, quarkline_lookuptable.Q2V, ric_lookup);
}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------






