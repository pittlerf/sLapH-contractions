/*! @file OperatorsForMesons.h
 *  Class declaration of OperatorsForMesons
 *
 *  @author Bastian Knippschild
 *  @author Markus Werner
 */

#pragma once

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <Eigen/Dense>

#include "boost/filesystem.hpp"
#include "boost/format.hpp"
#include "boost/multi_array.hpp"


#include "EigenVector.h"
#include "global_data_typedefs.h"
#include "GaugeField.h"
#include "RandomVector.h"
#include "typedefs.h"

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
 *  @note (MW 3.12.2017) rVdaggerVr moved to QuarkLineBlock
 *  @note (MW 21.12.2017) rVdaggerV moved to QuarkLineBlock
 */
class OperatorFactory {
 public:
  /*! Constructor which allocates memory for all operators */
  OperatorFactory(const size_t Lt,
                  const size_t Lx,
                  const size_t Ly,
                  const size_t Lz,
                  const size_t nb_ev,
                  const size_t dilE,
                  const OperatorLookup &operator_lookuptable,
                  const std::string &handling_vdaggerv,
                  const std::string &path_vdaggerv,
                  const std::string &path_config);
  /*! Standard Destructor
   *
   *  Everything should be handled by Eigen, std::vector, and boost::multi_array
   */
  ~OperatorFactory(){};

  /****************************************************************************/
  /**************** INTERFACE FOR BUILDING ALL OPERATORS **********************/

  void create_operators(const std::string &filename,
                        const RandomVector &rnd_vec,
                        const int config);
  /*! Free memory of vdaggerv */
  void free_memory_vdaggerv();

  /*! @todo check of vdaggerv is already build */
  inline const Eigen::MatrixXcd &return_vdaggerv(const size_t index,
                                                 const size_t t) const {
    return vdaggerv[index][t];
  }

private:
  // Containers for operators which are accessible from outside
  array_Xcd_d2_eigen vdaggerv;
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
  std::string path_config;

  // Internal functions to build individual operators --> The interface to these
  // functions is 'create_Operators'
  // input -> filename: name and path of eigenvectors
  void build_vdaggerv(const std::string &filename, const int config);
  void read_vdaggerv(const int config);
  void read_vdaggerv_liuming(const int config);


};
