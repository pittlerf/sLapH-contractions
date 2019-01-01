#pragma once

#include "global_data_typedefs.hpp"
#include "typedefs.hpp"

#include <boost/serialization/array.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/vector.hpp>

#include <sys/stat.h>
#include <iosfwd>

typedef std::map<std::string, std::vector<DiagramIndex>> DiagramIndicesCollection;
typedef std::map<std::string, std::vector<Indices>> TraceIndicesCollection;

enum class Location { source, sink };

struct TraceRequest {
  std::string tr_name;
  ssize_t tr_id;
  std::vector<Location> locations;

  bool operator==(TraceRequest const &other) const {
    return tr_name == other.tr_name && tr_id == other.tr_id &&
           locations == other.locations;
  }

  template <class Archive>
  void serialize(Archive &ar, const unsigned int version) {
    ar &tr_name;
    ar &tr_id;
    ar &locations;
  }
};

struct CorrelatorRequest {
  std::vector<TraceRequest> trace_requests;
  std::string hdf5_dataset_name;

  bool operator==(CorrelatorRequest const &other) const {
    return trace_requests == other.trace_requests;
  }

  template <class Archive>
  void serialize(Archive &ar, const unsigned int version) {
    ar &trace_requests;
    ar &hdf5_dataset_name;
  }
};

typedef std::map<std::string, std::vector<CorrelatorRequest>> CorrelatorRequestsMap;

/**
 * Class containing all metadata for contractions and functions to set them
 * from infile
 *
 *  Metadata roughly characterized by either
 *
 *  - physical parameters
 *  - flags
 *  - paths
 */
struct GlobalData {
  GlobalData();

  int Lx, Ly, Lz, Lt;
  int dim_row, V_TS, V_for_lime;
  int number_of_eigen_vec;
  int number_of_inversions;
  int start_config, end_config, delta_config;
  int verbose;
  ssize_t nb_eigen_threads;
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

  RandomVectorConstruction rnd_vec_construct;
  PerambulatorConstruction peram_construct;

  std::vector<quark> quarks;
  Operator_list operator_list;
  Correlator_list correlator_list;
  DilutedFactorIndicesCollection quarkline_lookuptable;
  OperatorLookup operator_lookuptable;
  TraceIndicesCollection trace_indices_map;
  CorrelatorRequestsMap correlator_requests_map;

  /**
   * Cutoff for the sum of individual momenta.
   *
   * The index is the total momentum squared, @f$ |P|^2 @f$, the value the sum
   * of the individual momenta squared, @f$ |p_1|^2 + |p_2|^2 @f$. The values
   * are chosen by hand from a feeling about the signal quality. When building
   * momentum combinations to compute, the condition @f$ |p_1|^2 + |p_2|^2 \le
   * c(|P|^2) @f$ will be enforced.
   *
   * Signal quality gets worse with large individual momenta, therefore it does
   * not make sense to include a very large @f$ p_1 @f$ in the @f$ |P|^2 = 0
   * @f$ case (with @f$ p_2 = - p_1 @f$). The current cutoff of 4 means that
   * only individual momenta up to @f$ (0, 0, 2) @f$ are computed.
   */
  std::map<int, int> momentum_cutoff;

  HypPars hyp_parameters;

  template <class Archive>
  void serialize(Archive &ar, const unsigned int version) {
    ar &Lx;
    ar &Ly;
    ar &Lz;
    ar &Lt;
    ar &dim_row;
    ar &V_TS;
    ar &V_for_lime;
    ar &number_of_eigen_vec;
    ar &number_of_inversions;
    ar &start_config;
    ar &end_config;
    ar &delta_config;
    ar &verbose;
    ar &nb_eigen_threads;
    ar &path_eigenvectors;
    ar &name_eigenvectors;
    ar &filename_eigenvectors;
    ar &path_perambulators;
    ar &name_perambulators;
    ar &name_lattice;
    ar &filename_ending_correlators;
    ar &path_output;
    ar &path_config;
    ar &handling_vdaggerv;
    ar &path_vdaggerv;
    ar &rnd_vec_construct;
    ar &peram_construct;
    ar &quarks;
    ar &operator_list;
    ar &correlator_list;
    ar &quarkline_lookuptable;
    ar &operator_lookuptable;
    ar &trace_indices_map;
    ar &correlator_requests_map;
  }
};

/**
 * Reading the input parameters from the infile in the main routine and
 * initializing GlobalData.
 */
void read_parameters(GlobalData &gd, int ac, char *av[]);

std::ostream &operator<<(std::ostream &os, GlobalData const &gd);
