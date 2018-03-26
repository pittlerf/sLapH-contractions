/*! @file contract.cpp
 *  Main function of Contraction code for Laphs Perambulators
 *
 *  @author Bastian Knippschild
 *  @author Christopher Helmes
 *  @author Christian Jost
 *  @author Liuming Liu
 *  @author Markus Werner
 *
 *  @version 0.1
 *  @copyright Copies are prohibited so far
 */ 

#include <iostream>

#include "omp.h"

#include "Correlators.h" // contains all other headers
#include "global_data.h"
#include "git.h"


/*! Read parameters from infile and perform the specified contractions
 *
 *  In succession instanciate GlobalData, Perambulator,
 *  RandomVector, OperatorsForMesons and Correlators. 
 *  - Get paths, physical quantum numbers and desired operators from infile 
 *  - Loop over Configuration
 *  - Read perambulators, randomvectors and contract
 *
 *  The flow of this function is depicted in the Flowchart below. The colors
 *  uniquely encode the class a member function belongs to. Beige fields are 
 *  not classmembers
 *  @image html contract.pdf "Flow chart of main()"
 */
int main (int ac, char* av[]) {
  std::cout << "This is sLapH-contractions:\n"
            << "  git branch " << git_refspec << "\n"
            << "  git revision " << git_sha1 << "\n"
            << "  git state " << git_changes << "\n"
            << "  compiled by " << git_user << " on " << git_host << "\n";

  // reading global parameters from input file
  GlobalData* global_data = GlobalData::Instance();
  global_data->read_parameters(ac, av);

  // initialization of OMP paralization
  Eigen::initParallel();
  omp_set_num_threads(global_data->get_nb_omp_threads());
  Eigen::setNbThreads(global_data->get_nb_eigen_threads());

  // ---------------------------------------------------------------------------
  // Creating instances of perambulators, random vectors, operators, and 
  // correlators. The eigenvectors are read from disc in the operator class.
  Perambulator perambulators(global_data->get_peram_construct().nb_entities,
                             global_data->get_peram_construct().size_rows,
                             global_data->get_peram_construct().size_cols);
  RandomVector randomvectors(global_data->get_rnd_vec_construct().nb_entities,
                             global_data->get_rnd_vec_construct().length);

  OperatorFactory meson_operators(
                            global_data->get_Lt(), global_data->get_Lx(),
                            global_data->get_Ly(), global_data->get_Lz(),
                            global_data->get_number_of_eigen_vec(),
                            (global_data->get_quarks())[0].number_of_dilution_E,
                            global_data->get_operator_lookuptable(),
                            global_data->get_handling_vdaggerv(),
                            global_data->get_path_vdaggerv(),
                            global_data->get_path_config());

  // ---------------------------------------------------------------------------
  // Loop over all configurations stated in the infile -------------------------
  for(size_t config_i  = global_data->get_start_config(); 
             config_i <= global_data->get_end_config(); 
             config_i += global_data->get_delta_config()){

    std::cout << "\nprocessing configuration: " << config_i 
              << "\n\n" << std::endl;
    // changes all paths and names which depend on the configuration
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
    // read eigenvectors and build operators
    meson_operators.create_operators(global_data->get_filename_eigenvectors(),
                                                       randomvectors, config_i);

    // doing all the contractions
    contract(global_data->get_Lt(),
             (global_data->get_quarks())[0].number_of_dilution_T,
             (global_data->get_quarks())[0].number_of_dilution_E,
             global_data->get_number_of_eigen_vec(),
             meson_operators,
             randomvectors,
             perambulators,
             global_data->get_operator_lookuptable(),
             global_data->get_correlator_lookuptable(),
             global_data->get_quarkline_lookuptable(),
             global_data->get_output_path(),
             global_data->get_filename_ending_correlators());
  }
}

