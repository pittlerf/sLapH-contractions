/*! @file global_data.h
 */

#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <sstream>
#include <string>
#include <vector>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/SparseCore>

#include <sys/stat.h>

#include "global_data_typedefs.h"
#include "typedefs.h"

/*! Class containing all metadata for contractions and functions to set them
 *  from infile
 *
 *  Metadata roughly characterized by either
 *  - physical parameters
 *  - flags
 *  - paths
 */
class GlobalData {
 public:
  /*! Checks whether GlobalData was initialized before. If yes, returns the
   *   existing object. If not, allocates new object on the stack.
   */
  static GlobalData *Instance();

  /*! Reading the input parameters from the infile in the main routine and
   *  initializing GlobalData
   */
  void read_parameters(int ac, char *av[]);

  /*! Fills the random vector and perambulator structs with the paths and
   *  file names to read the data
   */
  void build_IO_names(const ssize_t config);

  inline std::string get_name_lattice() { return name_lattice; }
  inline std::string get_output_path() { return path_output; }
  inline std::string get_filename_ending_correlators() {
    return filename_ending_correlators;
  }
  inline std::string get_path_config() { return path_config; }
  inline ssize_t get_nb_omp_threads() { return nb_omp_threads; }
  inline ssize_t get_nb_eigen_threads() { return nb_eigen_threads; }
  inline int get_Lx() { return Lx; }
  inline int get_Ly() { return Ly; }
  inline int get_Lz() { return Lz; }
  inline int get_Lt() { return Lt; }
  inline int get_dim_row() { return dim_row; }
  inline int get_V_TS() { return V_TS; }
  inline int get_V_for_lime() { return V_for_lime; }
  inline int get_number_of_inversions() { return number_of_inversions; }
  inline int get_number_of_rnd_vec() { return number_of_rnd_vec; }
  inline int get_start_config() { return start_config; }
  inline int get_end_config() { return end_config; }
  inline int get_delta_config() { return delta_config; }
  inline int get_number_of_eigen_vec() { return number_of_eigen_vec; }
  inline int get_verbose() { return verbose; }
  /*! Returns all information needed to build and read random vectors */
  inline RandomVectorConstruction get_rnd_vec_construct() { return rnd_vec_construct; }
  /*! Returns all information needed to build and read perambulators */
  inline PerambulatorConstruction get_peram_construct() { return peram_construct; }
  inline std::string get_path_eigenvectors() { return path_eigenvectors; }
  inline std::string get_name_eigenvectors() { return name_eigenvectors; }
  inline std::string get_filename_eigenvectors() { return filename_eigenvectors; }
  inline std::string get_path_perambulators() { return path_perambulators; }
  inline std::string get_name_perambulators() { return name_perambulators; }
  /*! Return the mode how to treat VdaggerV: Reading, or constructing anew.
   *
   *  @see LapH::OperatorsForMesons::create_operators()
   */
  inline std::string get_handling_vdaggerv() { return handling_vdaggerv; }
  inline std::string get_path_vdaggerv() { return path_vdaggerv; }

  /*! Return munged list of quarks as specified in the infile
   *
   * @todo make that private and delete get function
   */
  inline std::vector<quark> get_quarks() { return quarks; }

  /*! Return munged list of field operators as read from infile
   *
   * @todo make that private and delete get function
   */
  inline Operator_list &get_operator_list() { return operator_list; }
  /*! Return munged list of correlators as read from infile
   *
   * @todo make that private and delete get function
   */
  inline Correlator_list &get_correlator_list() { return correlator_list; }
  /*! Return vector of unique sLapH operators */
  inline const OperatorLookup get_operator_lookuptable() { return operator_lookuptable; }
  /*! Return vector of unique physical contents needed for quarklines */
  inline const DilutedFactorIndicesCollection get_quarkline_lookuptable() {
    return quarkline_lookuptable;
  }
  /*! Return vector of unique physical contents needed for correlators */
  inline const DiagramIndicesCollection get_correlator_lookuptable() {
    return correlator_lookuptable;
  }

  int momentum_cutoff_0;
    
 private:
  GlobalData() {}
  GlobalData(const GlobalData &other) {}
  //! @cond
  //! globally accessible data - naming should be clear and understandable
  int Lx, Ly, Lz, Lt;
  int dim_row, V_TS, V_for_lime;
  int number_of_eigen_vec;
  int number_of_rnd_vec;
  int number_of_inversions;
  int start_config, end_config, delta_config;
  int verbose;
  ssize_t nb_omp_threads, nb_eigen_threads;
  std::string path_eigenvectors;
  std::string name_eigenvectors;
  std::string filename_eigenvectors;
  std::string path_perambulators;
  std::string name_perambulators;
  std::string name_lattice;
  std::string filename_ending_correlators;
  std::string path_output;
  std::string path_config;
  std::string handling_vdaggerv;
  std::string path_vdaggerv;
  //! @endcond

  RandomVectorConstruction rnd_vec_construct;
  PerambulatorConstruction peram_construct;

  std::vector<quark> quarks;
  Operator_list operator_list;
  Correlator_list correlator_list;

  /*! Check and print lattice, config and path data and fills quarks,
   *  operator_lists and correlator_lists with the strings created from the
   *  infile
   */
  void input_handling(const std::vector<std::string> &quark_configs,
                      const std::vector<std::string> &operator_list_configs,
                      const std::vector<std::string> &correlator_list_configs);

  DilutedFactorIndicesCollection quarkline_lookuptable;
  OperatorLookup operator_lookuptable;
  DiagramIndicesCollection correlator_lookuptable;
  /*! Creates the lookup tables for quarklines, operators and correlators */
  void init_lookup_tables();
};
