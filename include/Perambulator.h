#ifndef _PERAMBULATOR_H_
#define _PERAMBULATOR_H_

#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

#include "Eigen/Dense"

#include "typedefs.h"

namespace LapH {

class Perambulator {

private:
  std::vector<Eigen::MatrixXcd> peram;
  size_t nb_entities;
  size_t size_rows, size_cols;

public:  
  // standard ctor and dtor are enough - all else is handled by std::vector
  Perambulator () : peram(0, Eigen::MatrixXcd(0,0)), nb_entities(0), 
                    size_rows(0), size_cols(0) {};
  // NOTE: IT is assumed that all perambulators have the same size!
  Perambulator(const size_t nb_entities, 
               const size_t size_rows, 
               const size_t size_cols):
                 peram(nb_entities, Eigen::MatrixXcd(size_rows, size_cols)),
                 nb_entities(nb_entities), size_rows(size_rows), 
                 size_cols(size_cols) {
    std::cout << "\tPerambulators initialised" << std::endl;
  }
  ~Perambulator() {};

  // [] operator to directly access the elements of perambulator
  inline const Eigen::MatrixXcd& operator[](const size_t entity) const {
    return peram[entity];
  }

  // rading one perambulators from some file
  // input: entity   -> the entry where this peram will be stores
  //        filename -> just the file name
  void read_perambulator(const size_t entity, const std::string& filename);
  // rading perambulators where each perambulator is stored in a different file
  // input: filename_list -> vector which contains all file names
  void read_perambulators_from_separate_files(
                                 const std::vector<std::string>& filename_list);

};

} // end of namespace

#endif // _PERAMBULATOR_H__
