#ifndef _PERAMBULATOR_H_
#define _PERAMBULATOR_H_

#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

#include "Eigen/Dense"

#include "typedefs.h"
#include "global_data_typedefs.h"

namespace LapH {

/*! Memory allocation and reading routines for perambulators
 */
class Perambulator {

private:
  std::vector<Eigen::MatrixXcd> peram;

public:  
  // just a default constructor
  Perambulator () : peram(0, Eigen::MatrixXcd(0,0)) {};
  // constructor which creates the correct sizes of the matrices
  // input: nb_entitites -> how many perambulators are there
  //        size_rows    -> number of rows for each perambulator
  //        size_cols    -> number of columns for each prambulator
  Perambulator(const size_t nb_entities, 
               const std::vector<size_t>& size_rows, 
               const std::vector<size_t>& size_cols):
                                    peram(nb_entities, Eigen::MatrixXcd(0, 0)) {
    // TODO: Think about putting this in initialisation list (via lambda?)
    for(size_t i = 0; i < nb_entities; i++)
      peram[i].resize(size_rows[i], size_cols[i]);

    std::cout << "\tPerambulators initialised" << std::endl;
  }
  // the default deconstructor - std::vector and Eigen should handle everything
  ~Perambulator() {};

  // [] operator to directly access the elements of perambulator
  inline const Eigen::MatrixXcd& operator[](const size_t entity) const {
    return peram.at(entity);
  }

  // rading one perambulators from some file
  // input: entity       -> the entry where this peram will be stores
  //        Lt           -> total number of timeslices - for each peram the same
  //        nb_eigen_vec -> total number of eigen vecs - for each peram the same
  //        q            -> contains information about dilution scheme and size
  //        filename     -> just the file name
  void read_perambulator(const size_t entity, const size_t Lt, 
                         const size_t nb_eigen_vec, const quark& q,
                         const std::string& filename);
  // rading perambulators where each perambulator is stored in a different file
  // input: Lt           -> total number of timeslices - for each peram the same
  //        nb_eigen_vec -> total number of eigen vecs - for each peram the same
  //        quark        -> contains information about dilution scheme and size
  //        filename_list -> vector which contains all file names
  void read_perambulators_from_separate_files(
                                 const size_t Lt, const size_t nb_eigen_vec,
                                 const std::vector<quark>& quark,
                                 const std::vector<std::string>& filename_list);

};
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

} // end of namespace

#endif // _PERAMBULATOR_H__
