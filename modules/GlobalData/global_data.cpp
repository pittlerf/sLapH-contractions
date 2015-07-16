/*
 * GlobalData.cpp
 *
 *  Created on: Mar 28, 2013
 *      Author: knippsch
 */

#include "global_data.h"

namespace po = boost::program_options;

GlobalData* GlobalData::instance_ = 0;

GlobalData* GlobalData::Instance () {

  if(instance_ == 0) instance_ = new GlobalData;

  return instance_;
}
// *****************************************************************************
/// @brief Convenience function for when a 'store_to' value is being provided
///        to typed_value.
///
/// @param store_to The variable that will hold the parsed value upon notify.
///
/// @return Pointer to a type_value.
template<typename T>
boost::program_options::typed_value<T>* make_value (T* store_to) {
  return boost::program_options::value<T>(store_to);
}

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
void GlobalData::read_parameters (int ac, char* av[]) {

  std::string input_file;
  std::string output_file;
  // Variables that will store parsed values for quarks.
  std::vector<std::string> quark_configs;
  // Variables that will store parsed values for operators.
  std::vector<std::string> operator_list_configs;
  // Variables that will store parsed values for correlators.
  std::vector<std::string> correlator_list_configs;

  // Declare a group of options that will be allowed only on command line
  po::options_description generic("Command line options");
  generic.add_options()("help,h", "produce help message")("version,v",
      "print version string")("verbose",
      "does additional tests and prints more details")("input,i",
      po::value<std::string>(&input_file)->default_value("LapHs.in"),
      "name of input file.")("output,o",
      po::value<std::string>(&output_file)->default_value("LapHs.out"),
      "name of output file.");

  // Declare a group of options that will be
  // allowed both on command line and in input file
  po::options_description config("Input file options");
  // parallelisation options
  config.add_options()
    ("nb_omp_threads",  
      po::value<size_t>(&nb_omp_threads)->default_value(1),
      "nb_omp_threads: number of openMP threads")
    ("nb_eigen_threads",
      po::value<size_t>(&nb_eigen_threads)->default_value(1),
      "nb_eigen_threads: number of threads Eigen uses internally");
  // lattice options
  config.add_options()
    ("output_path",
      po::value<std::string>(&path_output)-> 
                                          default_value("../../contractions"),
      "path for output")
    ("config_path",
      po::value<std::string>(&path_config)->default_value("../../configs"),
      "path for configurations")
    ("lattice", 
      po::value<std::string>(&name_lattice)-> default_value("lattice"),
      "Codename of the lattice")
    ("Lt", 
      po::value<int>(&Lt)->default_value(0),
      "Lt: temporal lattice extend")
    ("Lx",
      po::value<int>(&Lx)->default_value(0),
      "Lx: lattice extend in x direction")
    ("Ly",
      po::value<int>(&Ly)->default_value(0),
      "Ly: lattice extend in y direction")
    ("Lz",
      po::value<int>(&Lz)->default_value(0),
      "Lz: lattice extend in z direction");
  // eigenvector options
  config.add_options()
    ("number_of_eigen_vec",
      po::value<int>(&number_of_eigen_vec)->default_value(0),
      "Number of eigen vectors")
    ("path_eigenvectors",
      po::value<std::string>(&path_eigenvectors)->default_value("."),
      "directory of eigenvectors")
    ("name_eigenvectors",
      po::value<std::string>(&name_eigenvectors)->
                                                 default_value("eigenvector"),
      "name of eigenvectors\nThe full name is internally created to:\n"
      "\"name_of_eigenvectors.eigenvector\n. time slice.configuration\"");
  // quark options
  config.add_options()
    ("quarks.quark", make_value(&quark_configs),
     "quark input, must be of type:\n"
     "quark = \n type:number of rnd. vec.:\n"
     " dil type time:number of dil time:\n"
     " dil type ev:number of dil ev:\n"
     " dil type Dirac:number of dil Dirac:\n"
     " path of the perambulators for these quarks");
  // operator list options
  config.add_options()
    ("operator_lists.operator_list", make_value(&operator_list_configs),
     "operator input is rather complicated - see documentation!!");
  // correlator list options
  config.add_options()
    ("correlator_lists.correlator_list", make_value(&correlator_list_configs),
     "correlator input is rather complicated - see documentation!!");
  // configuration options
  config.add_options()
    ("start_config",
      po::value<int>(&start_config)->default_value(-1), 
      "First configuration")
    ("end_config", 
      po::value<int>(&end_config)->default_value(0),
      "Last configuration")
    ("delta_config",
      po::value<int>(&delta_config)->default_value(0),
      "Stepsize between two configurations");

  po::options_description cmdline_options;
  cmdline_options.add(generic).add(config);

  po::options_description input_file_options;
  input_file_options.add(config);

  po::options_description visible("Allowed options");
  visible.add(generic).add(config);
  po::positional_options_description p;
  p.add("input-file", -1);

  po::variables_map vm;
  po::store(po::command_line_parser(ac, av).options(cmdline_options).
                                            positional(p).run(), vm);
  po::notify(vm);

  // *************************************************************************
  // command line options ****************************************************
  if(vm.count("help")){
    std::cout << visible << "\n";
    exit(0);
  }
  if(vm.count("verbose")){
    verbose = 1;
  }
  else verbose = 0;
  if(vm.count("version")){
    std::cout << "Contraction code for LapHs perambulators. Version 0.1. \n";
    exit(0);
  }
  std::ifstream ifs(input_file.c_str());
  if(!ifs){
    std::cout << "CANNOT open input file: " << input_file << "\n";
    exit(0);
  }
  else{
    po::store(parse_config_file(ifs, input_file_options), vm);
    po::notify(vm);
  }
  ifs.close();

  // reading input file options
  input_handling(quark_configs, operator_list_configs, 
                 correlator_list_configs);

  // setting the lookup tables for all needed quantum numbers to calculate
  // the wanted correlators
  init_lookup_tables();

  // TODO: should be put in a separate function
  // setting the sizes and numbers of random vectors and perambulators
  rnd_vec_construct.nb_entities = 0;
  for(const auto& q: quarks)
    rnd_vec_construct.nb_entities += q.number_of_rnd_vec;
  rnd_vec_construct.length = Lt * 4 * number_of_eigen_vec;
  peram_construct.nb_entities = rnd_vec_construct.nb_entities;
  peram_construct.size_rows = rnd_vec_construct.length;
  // TODO: At this point it is not possible to have different dilution 
  //       schemes for different quarks! It is all the same.
  peram_construct.size_cols = (Lt / quarks[0].number_of_dilution_T) * 
                               quarks[0].number_of_dilution_E * 
                               quarks[0].number_of_dilution_D;
}
