/*! @file OperatorsForMesons.h
 *  Class declaration of LapH::OperatorsForMesons
 *
 *  @author Bastian Knippschild
 *  @author Markus Werner
 */

#ifndef OPERATORSFORMESONS_H_
#define OPERATORSFORMESONS_H_

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

#include "boost/multi_array.hpp"
#include "boost/filesystem.hpp"
#include "Eigen/Dense"

#include "EigenVector.h"
#include "GaugeField.h"
#include "RandomVector.h"
#include "typedefs.h"

namespace LapH {

/*! Calculates operators as they emerge in correlation functions using the 
 *  stochastic estimates from the stochastic Laplacian Heaviside method
 *
 *  Basically this calculates every operator of the form 
 *  - VdaggerV    : @f$ V^\dagger exp(ipx) V @f$
 *  - rVdaggerV   : @f$ (P^(b)\rho V)^\dagger exp(ipx) V @f$
 *  - rVdaggerVr  : @f$ (P^(b)\rho V)^\dagger exp(ipx) (P^(b)\rho V) @f$
 *
 *  It is straightforward to generalize the operators for Displacement, but 
 *  not implemented
 *
 *  @TODO Implement derivate operators
 */
class OperatorsForMesons {

private:

  // Containers for operators which are accessible from outside
  array_Xcd_d2_eigen vdaggerv;
  Xcd_d3_eigen rvdaggerv;
  Xcd_d3_eigen rvdaggervr;
  /*! @cond
   *  internal indices etc.
   */
  array_cd_d2 momentum;
  /*! @endcond */
  
  /****************************************************************************/
  /*! @TODO comment private members */
  const OperatorLookup operator_lookuptable;
  const size_t Lt, Lx, Ly, Lz;
  const size_t nb_ev, dilE;
  bool is_vdaggerv_set = false;
  std::string handling_vdaggerv;
  std::string path_vdaggerv;
  std::string path_gaugefields;

  // Internal functions to build individual operators --> The interface to these
  // functions is 'create_Operators'
  // input -> filename: name and path of eigenvectors
  void build_vdaggerv(const std::string& filename, const int config);
  void read_vdaggerv(const int config);
  void read_vdaggerv_liuming(const int config);
  void build_rvdaggerv(const LapH::RandomVector& rnd_vec);
  void build_rvdaggervr(const LapH::RandomVector& rnd_vec);

public:
  /*! Constructor which allocates memory for all operators */
  OperatorsForMesons(const size_t Lt, const size_t Lx, const size_t Ly, 
                     const size_t Lz, const size_t nb_ev, const size_t dilE,
                     const OperatorLookup& operator_lookuptable,
                     const std::string& handling_vdaggerv,
                     const std::string& path_vdaggerv,
                     const std::string& path_gaugefields);
  /*! Standard Destructor
   *
   *  Everything should be handled by Eigen, std::vector, and boost::multi_array
   */
  ~OperatorsForMesons () {};

  /****************************************************************************/
  /**************** INTERFACE FOR BUILDING ALL OPERATORS **********************/

  /*! Builds or reads @f$ V^\dagger exp(ipx) V @f$ and performs dilution i.e. 
   *  calculates rvdaggerv and rvdaggervr
   */
  void create_operators(const std::string& filename,
                        const LapH::RandomVector& rnd_vec, const int config);
  /*! Free memory of vdaggerv */
  void free_memory_vdaggerv();
  /*! Free memory of rvdaggerv */
  void free_memory_rvdaggerv();

  inline const Eigen::MatrixXcd& return_vdaggerv(const size_t index,
                                                 const size_t t) const {
    return vdaggerv[index][t];
  }

  inline const Eigen::MatrixXcd& return_rvdaggerv(const size_t index, 
                                                  const size_t t, 
                                                  const size_t rnd_id) const {
    return rvdaggerv.at(index).at(t).at(rnd_id);
  }

  inline const Eigen::MatrixXcd& return_rvdaggervr(const size_t index, 
                                                   const size_t t, 
                                                   const size_t rnd_id) const {
    return rvdaggervr.at(index).at(t).at(rnd_id);
  }

};

} // end of namespace

#endif // OPERATORSFORMESONS_H_ 


