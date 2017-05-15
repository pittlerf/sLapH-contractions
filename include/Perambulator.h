/*! @file Perambulator.h
 *  Class decleration of LapH::Perambulator
 *
 *  @author Bastian Knippschild
 *  @author Markus Werner
 */

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

/*! Memory allocation and reading routines for perambulators */
class Perambulator {

private:
  std::vector<Eigen::MatrixXcd> peram;

public:  

  /*! Default constructor */
  Perambulator () : peram(0, Eigen::MatrixXcd(0,0)) {};

  /*! Constructor which creates the correct sizes of the matrices
   *
   *  @param nb_entitites Number of perambulators
   *  @param size_rows    Number of rows for each perambulator
   *  @param size_cols    Number of columns for each perambulator
   */
  Perambulator(const size_t nb_entities, 
               const std::vector<size_t>& size_rows, 
               const std::vector<size_t>& size_cols):
                                    peram(nb_entities, Eigen::MatrixXcd(0, 0)) {
    // TODO: Think about putting this in initialisation list (via lambda?)
    for(size_t i = 0; i < nb_entities; i++)
      peram[i].resize(size_rows[i], size_cols[i]);

    std::cout << "\tPerambulators initialised" << std::endl;
  }

  /*! Default deconstructor - std::vector and Eigen should handle everything */
  ~Perambulator() {};

  /*! Overloading [] operator for Perambulator objects */
  inline const Eigen::MatrixXcd& operator[](const size_t entity) const {
    return peram.at(entity);
  }

  /*! Reading one perambulators from a single file */
  void read_perambulator(const size_t entity, const size_t Lt, 
                         const size_t nb_eigen_vec, const quark& quark,
                         const std::string& filename);

  /*! Reading perambulators where each perambulator is stored in a different 
   *  file 
   */
  void read_perambulators_from_separate_files(
                                 const size_t Lt, const size_t nb_eigen_vec,
                                 const std::vector<quark>& quark,
                                 const std::vector<std::string>& filename_list);

};

} // end of namespace

#endif // _PERAMBULATOR_H__
