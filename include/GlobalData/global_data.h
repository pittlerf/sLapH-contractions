/*
 * GlobalData.h
 *
 *  Created on: Mar 28, 2013
 *      Author: knippsch
 */

#ifndef GLOBALDATA_H_
#define GLOBALDATA_H_

#include <array>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <iterator>
#include <string>
#include <vector>

#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/SparseCore>

#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/program_options.hpp>

#include <sys/stat.h>

#include "global_data_typedefs.h"
#include "typedefs.h"

class GlobalData {

private:
  //! A pointer on the class itself
  static GlobalData* instance_;

  //! globally accessible data - naming should be clear and understandable
  int Lx, Ly, Lz, Lt;
  int dim_row, V_TS, V_for_lime;
  int number_of_eigen_vec;
  int number_of_rnd_vec;
  int number_of_inversions;
  int start_config, end_config, delta_config;
  int verbose;
  size_t nb_omp_threads, nb_eigen_threads;
  std::string path_eigenvectors;
  std::string name_eigenvectors;
  std::string filename_eigenvectors;
  std::string path_perambulators;
  std::string name_perambulators;
  std::string name_lattice;
  std::string path_output;
  std::string overwrite;
  std::string path_config;

  RandomVectorConstruction rnd_vec_construct;
  PerambulatorConstruction peram_construct;
 
  // The function input handling creates quarks, operator_lists and
  // correlator_lists from the strings created from the infile!
  // in: vectors of strings of: quark_configs
  //                            operator_list_configs
  //                            correlator_list_configs
  // out: there is no output but the vectors of: quarks
  //                                             operator_list
  //                                             correlator_list
  //      are filled in this function. These lists are needed to build the
  //      lookup tables in init_lookup_tables.
  std::vector<quark> quarks;
  std::vector<Operator_list> operator_list;
  Correlator_list correlator_list;
  void input_handling(const std::vector<std::string>& quark_configs,
                      const std::vector<std::string>& operator_list_configs,
                      const std::vector<std::string>& correlator_list_configs);

  // The function init_lookup_tables creates the lookup tables from quarks, 
  // operator_list and correlator_list.
  OperatorLookup operator_lookuptable;
  QuarklineLookup quarkline_lookuptable;
  CorrelatorLookup correlator_lookuptable;
  void init_lookup_tables();

public:
  static GlobalData* Instance ();

  // reading the input parameters from the infile in the main routine
  void read_parameters(int ac, char* av[]);

  // This fills the random vector and prambulator structs with the paths and
  // file names to read the data in depending on the config.
  void build_IO_names(const size_t config);

  // return functions for parameters
  inline std::string get_name_lattice() {
    return name_lattice;
  }
  inline std::string get_output_path() {
    return path_output;
  }
  inline std::string get_config_path() {
    return path_config;
  }
  inline size_t get_nb_omp_threads() {
    return nb_omp_threads;
  }
  inline size_t get_nb_eigen_threads() {
    return nb_eigen_threads;
  }
  inline int get_Lx () {
    return Lx;
  }
  inline int get_Ly () {
    return Ly;
  }
  inline int get_Lz () {
    return Lz;
  }
  inline int get_Lt () {
    return Lt;
  }
  inline int get_dim_row () {
    return dim_row;
  }
  inline int get_V_TS () {
    return V_TS;
  }
  inline int get_V_for_lime () {
    return V_for_lime;
  }
  inline int get_number_of_inversions () {
    return number_of_inversions;
  }
  inline int get_number_of_rnd_vec () {
    return number_of_rnd_vec;
  }
  inline int get_start_config () {
    return start_config;
  }
  inline int get_end_config () {
    return end_config;
  }
  inline int get_delta_config () {
    return delta_config;
  }
  inline int get_number_of_eigen_vec() {
    return number_of_eigen_vec;
  }
  inline int get_verbose() {
    return verbose;
  }
  inline RandomVectorConstruction get_rnd_vec_construct(){
    return rnd_vec_construct;
  }
  inline PerambulatorConstruction get_peram_construct() {
    return peram_construct;
  }
  inline std::string get_path_eigenvectors() {
    return path_eigenvectors;
  }
  inline std::string get_name_eigenvectors() {
    return name_eigenvectors;
  }
  inline std::string get_filename_eigenvectors () {
    return filename_eigenvectors;
  }
  inline std::string get_path_perambulators() {
    return path_perambulators;
  }
  inline std::string get_name_perambulators() {
    return name_perambulators;
  }
  inline std::vector<quark> get_quarks() {
    return quarks;
  }
  inline std::vector<Operator_list>& get_operator_list() {
    return operator_list;
  }
  inline Correlator_list& get_correlator_list() {
    return correlator_list;
  }
  inline const OperatorLookup get_operator_lookuptable(){
    return operator_lookuptable;
  }
  inline const QuarklineLookup get_quarkline_lookuptable(){
    return quarkline_lookuptable;
  }
  //! All con/de-structors are protected to assure that only one instance exists
  //! at once. DO NOT CHANGE!!
protected:
  GlobalData () {}
  GlobalData (const GlobalData& other) {}
  virtual ~GlobalData () {}

};

#endif /* GLOBALDATA_H_ */
