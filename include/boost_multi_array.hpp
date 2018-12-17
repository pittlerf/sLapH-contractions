#pragma once

#include "ComplexProduct.hpp"

#include <boost/multi_array.hpp>
#include <Eigen/Dense>

/** Data type for momentum */
typedef boost::multi_array<Complex, 2> array_cd_d2;
typedef boost::multi_array<Complex, 3> array_cd_d3;
typedef boost::multi_array<Complex, 4> array_cd_d4;
typedef boost::multi_array<Complex, 5> array_cd_d5;
typedef boost::multi_array<Complex, 6> array_cd_d6;
typedef boost::multi_array<Complex, 7> array_cd_d7;
typedef boost::multi_array<Complex, 8> array_cd_d8;
typedef boost::multi_array<Complex, 9> array_cd_d9;
typedef boost::multi_array<Complex, 10> array_cd_d10;

/** Special type for Quarklines */
typedef boost::multi_array<std::vector<Eigen::MatrixXcd>, 3> array_quarkline;

// Eigen typedefs
typedef std::vector<Eigen::MatrixXcd> vec_Xcd_eigen;
/** Data type for rvdaggerv and rvdaggervr */
typedef std::vector<std::vector<std::vector<Eigen::MatrixXcd>>> Xcd_d3_eigen;
/** Data type for vdaggerv */
typedef boost::multi_array<Eigen::MatrixXcd, 2> array_Xcd_d2_eigen;
typedef boost::multi_array<Eigen::MatrixXcd, 3> array_Xcd_d3_eigen;
typedef boost::multi_array<Eigen::MatrixXcd, 4> array_Xcd_d4_eigen;
typedef boost::multi_array<Eigen::MatrixXcd, 5> array_Xcd_d5_eigen;
typedef boost::multi_array<Eigen::MatrixXcd, 6> array_Xcd_d6_eigen;
typedef boost::multi_array<Eigen::MatrixXcd, 7> array_Xcd_d7_eigen;
typedef boost::multi_array<Eigen::MatrixXcd, 8> array_Xcd_d8_eigen;
typedef boost::multi_array<Eigen::MatrixXcd, 9> array_Xcd_d9_eigen;
typedef boost::multi_array<Eigen::MatrixXcd, 10> array_Xcd_d10_eigen;

/** Special type for Correlators */
typedef boost::multi_array<std::vector<Complex>, 3> array_corr;

