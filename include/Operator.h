/*! @file Operator.h
 *  Declaration of Operator multiplication
 *
 *  @author Martin Ueding
 *  @author Markus Werner
 *
 */

#pragma once

#include <vector>
#include "Eigen/Dense"

#include "OperatorsForMesons.h"
#include "QuarkLineBlock.h"
#include "typedefs.h"

void Q2VxrVdaggerVr(std::vector<Eigen::MatrixXcd> &result, 
                    QuarkLineBlock<Q2V> const &quarklines,
                    OperatorsForMesons const &meson_operator,
                    int const b2,
                    int const t1,
                    int const t2,
                    std::array<size_t, 3> const look,
                    OperatorLookup &operator_lookup,
                    QuarklineLookup &quark_lookup,
                    size_t const dilE,
                    size_t const dilD);
