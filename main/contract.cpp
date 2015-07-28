// *****************************************************************************
// Name         : contract.cpp
// Author       : Bastian Knippschild
// Contributors : Christopher Helmes, Christian Jost, Liuming Liu, Markus Werner
// Version      : 0.1
// Copyright    : Copies are prohibited so far
// Description  : Contraction code for Laphs Perambulators
// *****************************************************************************

#include <iostream>

#include "omp.h"

#include "global_data.h"
#include "Perambulator.h"
#include "RandomVector.h"

int main (int ac, char* av[]) {

  // reading global parameters from input file
  GlobalData* global_data = GlobalData::Instance();
  global_data->read_parameters(ac, av);

  // initialization of OMP paralization
  Eigen::initParallel();
  omp_set_num_threads(global_data->get_nb_omp_threads());
  Eigen::setNbThreads(global_data->get_nb_eigen_threads());

  // ***************************************************************************
  // Creating instances of perambulators, random vectors, operators, and 
  // correlators. The eigenvectors are read from disc in the operator class.
  LapH::Perambulator perambulators(
                                 global_data->get_peram_construct().nb_entities,
                                 global_data->get_peram_construct().size_rows,
                                 global_data->get_peram_construct().size_cols);
  LapH::RandomVector randomvectors(
                               global_data->get_rnd_vec_construct().nb_entities,
                               global_data->get_rnd_vec_construct().length);
//  operators = operators();  
//  correlators = correlators();

  // ***************************************************************************
  // Loop over all configurations stated in the infile *************************
  for(size_t config_i  = global_data->get_start_config(); 
             config_i <= global_data->get_end_config(); 
             config_i += global_data->get_delta_config()){
    std::cout << "\nprocessing configuration: " << config_i << "\n\n";
    global_data->build_IO_names(config_i);
    // read perambulators
    perambulators.read_perambulators_from_separate_files(
                              global_data->get_Lt(),
                              global_data->get_number_of_eigen_vec(),
                              global_data->get_quarks(),
                              global_data->get_peram_construct().filename_list);
    // read random vectors
    randomvectors.read_random_vectors_from_separate_files(
                            global_data->get_rnd_vec_construct().filename_list);
//    // read eigenvectors and build operators
//    operators.build();
//    // doing all the contractions
//    correlators.contract(peram, operators);
  }
  // That's all Folks!
  return 0;
}

