#ifndef RANDOMVECTOR_H_
#define RANDOMVECTOR_H_

#include <complex>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "ranlxs.h"

namespace LapH {

typedef std::complex<double> cmplx;

/*! Memory allocation and IO routines for perambulators
 */
class RandomVector {

private:
  // the random vector
  std::vector<cmplx> vec;
  const size_t nb_entities;
  const size_t length;

public:
  // standard ctor and dtor are enough - all else is handled by std::vector
  RandomVector() : vec(0, cmplx(0.0,0.0)), nb_entities(0), length(0) {
    std::cout << "The standard ctor of random vectors is useless!" << std::endl;
    exit(0);
  };
  // NOTE: It is assumed that all random vector have the same length!
  RandomVector(const size_t nb_entities, const size_t length) : 
                                     vec(nb_entities*length, cmplx(0.0,0.0)), 
                                     nb_entities(nb_entities), length(length) {
    std::cout << "\tRandom vectors initialised" << std::endl;
  };
  ~RandomVector() {};

  // () operator to directly access the elements of vec
  inline cmplx operator()(const size_t entity, const size_t entry) const {
    return vec.at(entity*length+entry);
  }

  // computes the random vectors for the sources
  // input: entity -> the random vector you want to set
  //        seed   -> seed for the random vector
  //        length -> length of the random vector
  void set(const size_t entity, const int seed);
  // computes the random vectors for the sources and stores them
  // input: entity   -> the random vector you want to set
  //        seed     -> seed for the random vector
  //        length   -> length of the random vector
  //        filename -> the random vector is direlty stored in this file
  void set(const size_t entity, const int seed, const std::string& filename);
  // writing all random vector to some file
  // input: filename -> the filename
  void write_random_vector(const std::string& filename) const;
  // writing one random vector (entity) to some file
  // input: entity   -> the random vector which will be stored on disk
  //        filename -> the filename
  void write_random_vector(const size_t entity, 
                           const std::string& filename) const;
  // reading all random vectors from some file
  // input: filename -> the filename
  void read_random_vector(const std::string& filename);
  // reading one random vector from some file
  // input: entity   -> the place where this random vector will be stored
  //        filename -> the filename
  void read_random_vector(const size_t entity, const std::string& filename);
  // rading random vectors where each vector is stored in a different file
  // input: filename_list -> vector which contains all file names
  void read_random_vectors_from_separate_files(
                                 const std::vector<std::string>& filename_list);

};
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

} // end of namespace

#endif // RANDOMVECTOR_H_ 
