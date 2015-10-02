#ifndef OPERATORSFORMESONS_H_
#define OPERATORSFORMESONS_H_

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

#include "boost/multi_array.hpp"
#include "Eigen/Dense"

#include "EigenVector.h"
#include "RandomVector.h"
#include "typedefs.h"

namespace LapH {

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
class OperatorsForMesons {

private:
  // contains for operators which are accessible from outside
  array_Xcd_d2_eigen vdaggerv;
  Xcd_d3_eigen rvdaggerv;
  Xcd_d3_eigen rvdaggervr;
  // internal indices etc.
  array_cd_d2 momentum;
  const OperatorLookup operator_lookuptable;
  const size_t Lt, Lx, Ly, Lz, nb_ev, dilE;
  bool is_vdaggerv_set = false;

  // internal functions to build individual operators --> The interface to these
  // functions is 'create_Operators'
  void build_vdaggerv(const std::string& filename);
  void build_rvdaggerv(const LapH::RandomVector& rnd_vec);
  void build_rvdaggervr(const LapH::RandomVector& rnd_vec);

public:
  // standard ctor makes no sence
  OperatorsForMesons () : vdaggerv(), momentum(), operator_lookuptable(),
                          Lt(0), Lx(0), Ly(0), Lz(0), nb_ev(0), dilE(0) {
    std::cout << "Standard constructor for Operators makes no sence" 
              << std::endl;
    exit(0);
  };
  // ctor which builds up all the memory
  OperatorsForMesons (const size_t Lt, const size_t Lx, const size_t Ly, 
                      const size_t Lz, const size_t nb_ev, const size_t dilE,
                      const OperatorLookup& operator_lookuptable);
  // standard dtor - eveything should be handled by Eigen, std::vector, and
  //                 Boost::Mulidimensional arrays.
  ~OperatorsForMesons () {};


  // -------------- INTERFACE FOR BUILDING ALL OPERATORS -----------------------
  // ---------------------------------------------------------------------------
  void create_operators(const std::string& filename,
                        const LapH::RandomVector& rnd_vec);
  // memory of vdaggerv can be freed when it's not needed, eg after building Q2
  void free_memory_vdaggerv();
  // memory of rvdaggerv can be freed when it's not needed, eg after building Q1
  void free_memory_rvdaggerv();


  // -------------- RETURN FUNCTIONS -------------------------------------------
  // ---------------------------------------------------------------------------
  inline const Eigen::MatrixXcd& return_vdaggerv(const size_t index,
                                                 const size_t t) const {
    return vdaggerv[index][t];
  }
  // ---------------------------------------------------------------------------
  inline const Eigen::MatrixXcd& return_rvdaggerv(const size_t index, 
                                                  const size_t t, 
                                                  const size_t rnd_id) const {
    return rvdaggerv[index][t][rnd_id];
  }
  // ---------------------------------------------------------------------------
  inline const Eigen::MatrixXcd& return_rvdaggervr(const size_t index, 
                                                   const size_t t, 
                                                   const size_t rnd_id) const {
    return rvdaggervr[index][t][rnd_id];
  }

};
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

} // end of namespace

#endif // OPERATORSFORMESONS_H_ 


