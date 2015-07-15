//============================================================================
// Name        : contract.cpp
// Author      : BK
// Version     :
// Copyright   : Copies are prohibited so far
// Description : Contraction code for Laphs Perambulators
//============================================================================

#include <iostream>

#include "omp.h"

#include "global_data.h"

int main (int ac, char* av[]) {

  // reading in global parameters from input file
  GlobalData* global_data = GlobalData::Instance();
  global_data->read_parameters(ac, av);

  // initialization of OMP paralization
  Eigen::initParallel();
  omp_set_num_threads(global_data->get_nb_omp_threads());
  Eigen::setNbThreads(global_data->get_nb_eigen_threads());

//  // ***************************************************************************
//  // Creating instances of perambulators, random vectors, operators, and 
//  // correlators. The eigenvectors are read from disc in the operator class;
//  peram = perambulator();
//  randomvector = randomvector();
//  operators = operators();  
//  correlators = correlators();
//
//  // ***************************************************************************
//  // Loop over all configurations **********************************************
//  // ***************************************************************************
//  for(size_t config_i  = global_data->get_start_config(); 
//             config_i <= global_data->get_end_config(); 
//             config_i += global_data->get_delta_config()){
//    std::cout << "\nprocessing configuration: " << config_i << "\n\n";
//    // read perambulators
//    peram.read();
//    // read random vectors
//    randomvector.read();
//    // read eigenvectors and build operators
//    operators.build();
//    // doing all the contractions
//    correlators.contract(peram, operators);
//  }
//  // That's all Folks!
//  return 0;
}

