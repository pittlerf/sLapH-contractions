/*
 * typedefs.h
 *
 *  Created on: Aug 28, 2014
 *      Author: C. Jost
 */

#ifndef _TYPEDEFS_H_
#define _TYPEDEFS_H_

#include <algorithm>
#include <complex>
#include <vector>
#include <list>

#include "Eigen/Core"
#include "Eigen/Dense"
#include "boost/multi_array.hpp"

// cpp standard library typedefs
typedef std::complex<double> cmplx;
typedef std::vector<cmplx> vec;
typedef boost::multi_array<cmplx, 2> array_cd_d2;
typedef boost::multi_array<cmplx, 3> array_cd_d3;
typedef boost::multi_array<cmplx, 4> array_cd_d4;
typedef boost::multi_array<cmplx, 5> array_cd_d5;
typedef boost::multi_array<cmplx, 6> array_cd_d6;
typedef boost::multi_array<cmplx, 7> array_cd_d7;
typedef boost::multi_array<cmplx, 8> array_cd_d8;
typedef boost::multi_array<cmplx, 9> array_cd_d9;
typedef boost::multi_array<cmplx, 10> array_cd_d10;

// special type for Corr and Q2
typedef boost::multi_array<std::vector<std::vector<cmplx> >, 4> array_corr;
typedef boost::multi_array<std::vector<std::vector<Eigen::MatrixXcd> >, 4> 
                                                                       array_q2;

// Eigen typedefs
typedef std::vector<Eigen::MatrixXcd> vec_Xcd_eigen;
typedef std::vector<std::vector<std::vector<Eigen::MatrixXcd> > > Xcd_d3_eigen;
typedef boost::multi_array<Eigen::MatrixXcd, 2> array_Xcd_d2_eigen;
typedef boost::multi_array<Eigen::MatrixXcd, 3> array_Xcd_d3_eigen;
typedef boost::multi_array<Eigen::MatrixXcd, 4> array_Xcd_d4_eigen;
typedef boost::multi_array<Eigen::MatrixXcd, 5> array_Xcd_d5_eigen;
typedef boost::multi_array<Eigen::MatrixXcd, 6> array_Xcd_d6_eigen;
typedef boost::multi_array<Eigen::MatrixXcd, 7> array_Xcd_d7_eigen;
typedef boost::multi_array<Eigen::MatrixXcd, 8> array_Xcd_d8_eigen;
typedef boost::multi_array<Eigen::MatrixXcd, 9> array_Xcd_d9_eigen;
typedef boost::multi_array<Eigen::MatrixXcd, 10> array_Xcd_d10_eigen;

// index typedefs
typedef std::list<size_t> indexlist_1;
typedef std::list<std::pair<size_t, size_t> > indexlist_2;
typedef std::list<std::array<size_t, 3> > indexlist_3;
typedef std::list<std::array<size_t, 4> > indexlist_4;


// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
  struct VdaggerVQuantumNumbers{ 
    size_t id;
    std::array<int, 3> momentum;
    std::array<int, 3> displacement;
    VdaggerVQuantumNumbers(const size_t id, const std::array<int, 3>& momentum, 
             const std::array<int, 3>& displacement) :
             id(id), momentum(momentum), displacement(displacement) {};
  };
  // ---------------------------------------------------------------------------
  // ---------------------------------------------------------------------------
  struct RandomIndexCombinationsQ2{
    size_t id;
    size_t id_q1, id_q2;
    std::vector<std::pair<size_t, size_t> > rnd_vec_ids;
    // just a small constructor to ensure easy filling of its vector form
    RandomIndexCombinationsQ2(const size_t id, 
                   const size_t id_q1, const size_t id_q2, 
                   const std::vector<std::pair<size_t, size_t> >& rnd_vec_ids) :
                         id(id), id_q1(id_q1), id_q2(id_q2), 
                         rnd_vec_ids(rnd_vec_ids) {};
  };  
  // ---------------------------------------------------------------------------
  // ---------------------------------------------------------------------------
  struct RandomIndexCombinationsQ1{
    size_t id;
    size_t id_q1;
    std::vector<size_t> rnd_vec_ids;
    // just a small constructor to ensure easy filling of its vector form
    RandomIndexCombinationsQ1(const size_t id, const size_t id_q1,
                          const std::vector<size_t>& rnd_vec_ids) :
                              id(id), id_q1(id_q1), rnd_vec_ids(rnd_vec_ids) {};
  };  
  // ---------------------------------------------------------------------------
  // ---------------------------------------------------------------------------
  struct VdaggerVRandomLookup{
    size_t id;
    size_t id_vdaggerv;
    size_t id_ricQ_lookup;
    bool need_vdaggerv_daggering;
    // just a small constructor to ensure easy filling of its vector form
    VdaggerVRandomLookup(const size_t id, const size_t id_vdaggerv,
              const size_t id_ricQ_lookup, const bool need_vdaggerv_daggering) :
               id(id), id_vdaggerv(id_vdaggerv), id_ricQ_lookup(id_ricQ_lookup),
               need_vdaggerv_daggering(need_vdaggerv_daggering) {};
  };
  // ---------------------------------------------------------------------------
  // ---------------------------------------------------------------------------
  struct OperatorLookup{

    std::vector<VdaggerVQuantumNumbers> vdaggerv_lookup;

    std::vector<RandomIndexCombinationsQ1> ricQ1_lookup;
    std::vector<RandomIndexCombinationsQ2> ricQ2_lookup;

    std::vector<VdaggerVRandomLookup> rvdaggerv_lookuptable;  
    std::vector<VdaggerVRandomLookup> rvdaggervr_lookuptable;  

    size_t index_of_unity;
  };
  // ---------------------------------------------------------------------------
  // ---------------------------------------------------------------------------
  struct QuarklineQ1Indices {
    size_t id, id_rvdaggerv, id_peram;
    std::vector<int> gamma; 
    // just a small constructor to ensure easy filling of its vector form
    QuarklineQ1Indices(const size_t id, const size_t id_rvdaggerv, 
                       const size_t id_peram, const std::vector<int>& gamma) :
                     id(id), id_rvdaggerv(id_rvdaggerv), id_peram(id_peram),
                     gamma(gamma) {};
  };
  // ---------------------------------------------------------------------------
  // ---------------------------------------------------------------------------
  struct QuarklineQ2Indices {
    size_t id, id_vdaggerv, id_peram1, id_peram2;
    bool need_vdaggerv_dag;
    std::vector<int> gamma; 
    // just a small constructor to ensure easy filling of its vector form
    QuarklineQ2Indices(const size_t id, const size_t id_vdaggerv, 
                       const size_t id_peram1, const size_t id_peram2, 
                       const bool need_vdaggerv_dag,
                       const std::vector<int>& gamma) :
                     id(id), id_vdaggerv(id_vdaggerv), id_peram1(id_peram1),
                     id_peram2(id_peram2), need_vdaggerv_dag(need_vdaggerv_dag),
                     gamma(gamma) {};
  };
  // ---------------------------------------------------------------------------
  // ---------------------------------------------------------------------------
  struct QuarklineLookup{
    std::vector<QuarklineQ1Indices> Q1;
    std::vector<QuarklineQ2Indices> Q2V;
    std::vector<QuarklineQ2Indices> Q2L;
  };
  // ---------------------------------------------------------------------------
  // ---------------------------------------------------------------------------
  struct CorrInfo{
    size_t id;
    std::string outfilename;
    std::vector<size_t> lookup;
    // just a small constructor to ensure easy filling of its vector form
    CorrInfo(const size_t id, const std::string& outfilename, 
             const std::vector<size_t>& lookup) :
                     id(id), outfilename(outfilename), lookup(lookup) {};
  };
  // ---------------------------------------------------------------------------
  // ---------------------------------------------------------------------------
  struct CorrelatorLookup{
    std::vector<CorrInfo> C20;
    std::vector<CorrInfo> C2c;
  };
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------



//// struct which contains all id's from pdg for which VdaggerV is not 
//// redundant. Only half the momenta and no gamma structure is contained
//// id    - number to reference VdaggerV
//// index - id of pdg
//  struct pd{
//    size_t id;
//    size_t index;
//  };
//
//// struct which contains all id's from pdg for which rVdaggerVr is not 
//// redundant. No gamma structure is contained. Half the momenta can be 
//// calculated by adjoining their negative counterpart, which will be done if
//// the adjoint-flag is set to 1
//// id    - number to reference VdaggerV
//// index - id of pdg and id for the adjoint pdg
//  struct pd_r{
//    size_t id;
//    size_t id_adjoint;
//    size_t index;
//    bool adjoint;
//  };
//
//  struct index_2pt {
//    size_t id;
//    size_t index_Q2;
//    size_t index_Corr;
//  };
//
//  struct index_4pt {
//    size_t id;
//    size_t index_Q2[2];
//    size_t index_Corr[2];
//  };
//
//  struct index_IO_1 {
//    size_t id;
//    indexlist_1 index_pt;
//  };
//
//  struct index_IO_2 {
//    size_t id;
//    indexlist_2 index_pt;
//  };
//
//typedef std::vector<pdg> vec_pdg_Corr;  
//typedef std::vector<pd> vec_pd_VdaggerV;
//typedef std::vector<pd_r> vec_pd_rVdaggerVr;
//typedef std::vector<index_2pt> vec_index_2pt;
//typedef std::vector<index_4pt> vec_index_4pt;
//typedef std::vector<index_IO_1> vec_index_IO_1;
//typedef std::vector<index_IO_2> vec_index_IO_2;

#endif // _TYPEDEFS_H_
