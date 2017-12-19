/*! @file typedefs.h
 *  Custom file definitions and abbreviations
 *
 *  @author Christian Jost
 *
 *  @date Aug 28, 2014
 *
 *  @TODO {seems like a lot of data types are deprecated. Comment these out?}
 */

#pragma once

#include <algorithm>
#include <array>
#include <complex>
#include <vector>
#include <list>

#include "Eigen/Core"
#include "Eigen/Dense"
#include "boost/multi_array.hpp"
#include "boost/container/static_vector.hpp"

#include "SmallVector.h"

int constexpr max_rnd_ids = 10;

using RndId = int8_t;
//using SmallVectorRndId = SmallVector<RndId, max_rnd_ids>;
using SmallVectorRndId = boost::container::static_vector<RndId, max_rnd_ids>;

/*! @{ Abbreviation for complex data types */
typedef std::complex<double> cmplx;
typedef std::vector<cmplx> vec;
/*! @} */

struct DilutedFactor {
  using Data = Eigen::MatrixXcd;

  Data data;
  std::pair<RndId, RndId> ric;
  SmallVectorRndId used_rnd_ids;
};

struct DilutedTrace {
  using Data = cmplx;

  Data data;
  SmallVectorRndId used_rnd_ids;
};

/*! Data type for momentum */
typedef boost::multi_array<cmplx, 2> array_cd_d2;
typedef boost::multi_array<cmplx, 3> array_cd_d3;
typedef boost::multi_array<cmplx, 4> array_cd_d4;
typedef boost::multi_array<cmplx, 5> array_cd_d5;
typedef boost::multi_array<cmplx, 6> array_cd_d6;
typedef boost::multi_array<cmplx, 7> array_cd_d7;
typedef boost::multi_array<cmplx, 8> array_cd_d8;
typedef boost::multi_array<cmplx, 9> array_cd_d9;
typedef boost::multi_array<cmplx, 10> array_cd_d10;

/*! Special type for Correlators */
typedef boost::multi_array<std::vector<cmplx>, 3> array_corr;
typedef boost::multi_array<std::vector<DilutedTrace>, 3> DilutedTraceCollection;
typedef boost::multi_array<std::vector<DilutedTrace>, 2> DilutedTraceCollection2;
/*! Special type for Quarklines */
typedef boost::multi_array<std::vector<Eigen::MatrixXcd>, 3> array_quarkline;

//typedef boost::multi_array<std::vector<std::vector<cmplx> >, 4> array_corr;
/*! @TODO {Is that deprecated?} */
typedef boost::multi_array<std::vector<cmplx>, 2> array_C1;

// Eigen typedefs
typedef std::vector<Eigen::MatrixXcd> vec_Xcd_eigen;
/*! Data type for rvdaggerv and rvdaggervr */
typedef std::vector<std::vector<std::vector<Eigen::MatrixXcd> > > Xcd_d3_eigen;
/*! Data type for vdaggerv */
typedef boost::multi_array<Eigen::MatrixXcd, 2> array_Xcd_d2_eigen;
typedef boost::multi_array<Eigen::MatrixXcd, 3> array_Xcd_d3_eigen;
typedef boost::multi_array<Eigen::MatrixXcd, 4> array_Xcd_d4_eigen;
typedef boost::multi_array<Eigen::MatrixXcd, 5> array_Xcd_d5_eigen;
typedef boost::multi_array<Eigen::MatrixXcd, 6> array_Xcd_d6_eigen;
typedef boost::multi_array<Eigen::MatrixXcd, 7> array_Xcd_d7_eigen;
typedef boost::multi_array<Eigen::MatrixXcd, 8> array_Xcd_d8_eigen;
typedef boost::multi_array<Eigen::MatrixXcd, 9> array_Xcd_d9_eigen;
typedef boost::multi_array<Eigen::MatrixXcd, 10> array_Xcd_d10_eigen;

/******************************************************************************/ 
/*! This is just a workaround for complex numbers to get it running for hdf5
 *
 *  @todo Change namespace into something more specific
 */
typedef struct { 
  double re; 
  double im; 
} complex_t; 

// This is the datatype to write 4pt functions and that alike directly
struct compcomp_t { 
  double rere;   
  double reim;
  double imre;
  double imim;   
  compcomp_t(const double rere, const double reim, 
             const double imre, const double imim) : 
                              rere(rere), reim(reim), imre(imre), imim(imim) {};
  compcomp_t& operator+=(compcomp_t const &other){
    rere += other.rere;   
    reim += other.reim;
    imre += other.imre;
    imim += other.imim;   

    return *this;
  }
};

/******************************************************************************/
/*! Struct to uniquely identify a sLapH operator
 *  @f$ V^\dagger exp(i(p + d/2) x) V @f$
 *
 *  In contrast to the field operator the Dirac structure is factored out
 *
 *  @todo check the expression for displacement
 *  @deprecated only displacement 0 is implement atm (26.3.17)
 */
struct VdaggerVQuantumNumbers{ 
  size_t id;
  std::array<int, 3> momentum;      /*!< The -momentum as 3-vector */
  std::array<int, 3> displacement;  /*!< The displacement as 3-vector */

  /*! Constructor */
  VdaggerVQuantumNumbers(const size_t id, const std::array<int, 3>& momentum, 
           const std::array<int, 3>& displacement) :
           id(id), momentum(momentum), displacement(displacement) {};
};

/******************************************************************************/
/*! Struct containing random index combinations for quarklines depending on two
 *  indices.
 *
 *  The random indices are uniquely identifying quark and random vector. Thus
 *  There are @f$ \sum_i q_i nb_rnd(q_i) @f$ random indices.
 *
 *  @todo is the offset still needed?
 */
struct RandomIndexCombinationsQ2{
  size_t id;
  /*! @{ 
   *  Numbers identifying quarks in quark_list
   */
  size_t id_q1, id_q2;
  /*! @} */
  /*! The entries of the pair correspond to the first and second random index.
   *  Contains the first random index corresponding to @em id_q1 and @em id_q2
   *  respectively, so that rnd_vec_id - offset is the actual random seed number
   */
  std::pair<size_t, size_t> offset;
  /*! The entries of the pair correspond to the first and second random index.
   *  List of all possible combinations of random vector indices for quarks
   *  specified by @em id_q1 and @em id_q2
   */
  std::vector<std::pair<size_t, size_t> > rnd_vec_ids;

  /*! Just a small constructor to ensure easy filling of its vector form */
  RandomIndexCombinationsQ2(const size_t id, 
                 const size_t id_q1, const size_t id_q2, 
                 const std::pair<size_t, size_t> offset, 
                 const std::vector<std::pair<size_t, size_t> >& rnd_vec_ids) :
                       id(id), id_q1(id_q1), id_q2(id_q2), offset(offset),
                       rnd_vec_ids(rnd_vec_ids) {};
};  

/******************************************************************************/
/*! Struct containing random index combinations for quarklines depending on one
 *  index.
 *
 *  The random indices are uniquely identifying quark and random vector. Thus
 *  There are @f$ \sum_i q_i nb_rnd(q_i) @f$ random indices.
 *
 *  @todo why is here no offset like in RandomIndexCombinationsQ2? 
 */
struct RandomIndexCombinationsQ1{
  size_t id;

  /*!  Identifier of quark in quark_list */
  size_t id_q1;
  /*! List of all possible random vector indices for a quark specified by 
   *  @em id_q1
   */
  std::vector<size_t> rnd_vec_ids;

  /*! just a small constructor to ensure easy filling of its vector form */
  RandomIndexCombinationsQ1(const size_t id, const size_t id_q1,
                        const std::vector<size_t>& rnd_vec_ids) :
                            id(id), id_q1(id_q1), rnd_vec_ids(rnd_vec_ids) {};
};  

/******************************************************************************/
/*! Struct that holds all information on which VdaggerV must be diluted with 
 *  which random vector.
 *
 *  For rVdaggerV and rVdaggerVr the VdaggerV-operators are additionaly 
 *  multiplied with random vectors. For both, VdaggerV and random index
 *  combinations there are lookuptables in OperatorLookup. This struct contains
 *  the id's of vdaggerv and ric which belong together.
 */
struct VdaggerVRandomLookup{
  size_t id;
  /*! id of vdaggerv_lookup in OperatorLookup */
  size_t id_vdaggerv;           
  /*! id of ricQ1_lookup for rvdaggerv or ricQ2_lookup for rvdaggervr in 
   *  OperatorLookup 
   */
  size_t id_ric_lookup;
  /*! Flag that indicates whether VdaggerV must be daggered (prior to 
   *  multiplication with random vectors) to get the correct quantum numbers
   */
  bool need_vdaggerv_daggering;

  /*! Just a small constructor to ensure easy filling of its vector form */
  VdaggerVRandomLookup(const size_t id, const size_t id_vdaggerv,
            const size_t id_ric_lookup, const bool need_vdaggerv_daggering) :
             id(id), id_vdaggerv(id_vdaggerv), id_ric_lookup(id_ric_lookup),
             need_vdaggerv_daggering(need_vdaggerv_daggering) {};
};

/******************************************************************************/
/*! Struct that contains all information for a sLapH operator 
 *
 *  @todo confusing because rvdaggerv_lookuptable and rvdaggervr_lookuptable
 *        contain indices of vdaggerv_lookup and ric_lookup. Multiple 
 *        hierarchy levels in a single struct
 *  @todo confusing because either Q1 and rvdaggerv OR Q2 and rvdaggervr are
 *        used, but in both cases half the members are spurious
 */
struct OperatorLookup{

  /*! Specifies physical content of quark field operator (with Dirac structure
   *  factored out)
   */
  std::vector<VdaggerVQuantumNumbers> vdaggerv_lookup;

  /*! Specifies the random vector index in case of rvdaggerv */
  std::vector<RandomIndexCombinationsQ1> ricQ1_lookup;
  /*! Specifies the random vector indices in case of rvdaggervr */
  std::vector<RandomIndexCombinationsQ2> ricQ2_lookup;

  /*! @{ 
   *  Specifies which entries of @em vdaggerv_lookup and @em ric_lookup shall
   *  be combined
   */
  std::vector<VdaggerVRandomLookup> rvdaggerv_lookuptable;  
  std::vector<VdaggerVRandomLookup> rvdaggervr_lookuptable;  
  /*! @} */

  /*! For @f$ \vec{p} = 0 @f$, @f$ V^\dagger exp(ipx) V = \mathbb{1} @f$. If 
   *  applicable This contains the index of @em vdaggerv_lookup where it can 
   *  be replaced by a unit matrix. If @f$ \vec{p} = 0 @f$ is not needed, 
   *  @em index_of_unity is set to -1
   */
  int index_of_unity;

};

/******************************************************************************/
/*! Indices needed to uniquely identify Q1 objects
 *
 *  Because vdaggerv is diagonal in dirac space, gamma may be factored out and 
 *  it proves useful to calculate and reuse 
 *
 *    Q1 = rvdaggerv * gamma * peram
 */
struct QuarklineQ1Indices {
  size_t id;
  /*! Identifies physical content and random index of rvdaggerv */
  size_t id_rvdaggerv;    
  /*! Identifies physical content and random index of peram */
  size_t id_peram;        
  /*! Pair of random vectors: 
   *  - First used for rvdaggerv 
   *  - Second used for peram
   */
  size_t id_ric_lookup; 
  std::vector<int> gamma; /*!< List of necessarry gamma combinations */
  /*! Just a small constructor to ensure easy filling of its vector form */
  QuarklineQ1Indices(const size_t id,
                     const size_t id_rvdaggerv,
                     const size_t id_peram,
                     const size_t id_ric_lookup,
                     const std::vector<int> &gamma)
      : id(id),
        id_rvdaggerv(id_rvdaggerv),
        id_peram(id_peram),
        id_ric_lookup(id_ric_lookup),
        gamma(gamma){};
};

/******************************************************************************/
/*! Indices needed to uniquely identify Q2 objects
 *
 *  Because vdaggerv is diagonal in dirac space, gamma may be factored out. 
 *
 *  If the @f$ \gamma_5 @f$-trick is used it proves useful to calculate and 
 *  reuse 
 *
 *    Q2 = @f$ \gamma_5 @f$ peram1@f$ ^\dagger \gamma_5 @f$ * vdaggerv * 
 *          gamma * peram2
 */
struct QuarklineQ2Indices {
  size_t id;
  /*! Identifies physical content of vdaggerv */
  size_t id_vdaggerv;     
  /*! Identifies perambulator with @f$\gamma_5@f$-trick 
   *
   *  @todo Is that deprecated?
   */
  size_t id_peram1;       
  /*! Identifies perambulator without @f$\gamma_5@f$-trick 
   *
   *  @todo Is that deprecated?
   */
  size_t id_peram2; 
  /*! id specifying the random vector indicies for peram1 and peram2 */
  size_t id_ric_lookup;
  /*! Flag that indicates whether VdaggerV must be daggered (prior to 
   *  multiplication with random vectors) to get the correct physical content 
   */
  bool need_vdaggerv_dag;
  std::vector<int> gamma; /*!< List of necessarry gamma combinations */

  /*! Just a small constructor to ensure easy filling of its vector form */
  QuarklineQ2Indices(const size_t id, const size_t id_vdaggerv, 
                     const size_t id_peram1, const size_t id_peram2, 
                     const size_t id_ric_lookup, const bool need_vdaggerv_dag,
                     const std::vector<int>& gamma) :
                   id(id), id_vdaggerv(id_vdaggerv), id_peram1(id_peram1),
                   id_peram2(id_peram2), id_ric_lookup(id_ric_lookup), 
                   need_vdaggerv_dag(need_vdaggerv_dag), gamma(gamma) {};
};

/******************************************************************************/
/*! Maps index from CorrelatorLookup to QuarklineQ1Indicies or 
 *  QuarklineQ2Indicies, depending on the quarkline needed in Correlator
 *
 *  @deprecated No longer necessary in memory_optimised branch as the quarklines
 *              are built on the fly
 */
struct QuarklineLookup{
  std::vector<QuarklineQ1Indices> Q1;
  std::vector<QuarklineQ2Indices> Q2V;
  std::vector<QuarklineQ2Indices> Q2L;
};


// Q0 formerly called rVdaggerVr
enum class QuarkLineType {Q0, Q1, Q2L, Q2V};

/******************************************************************************/
/*! All information needed to build and write the correlator given the 
 *  quarklines were calculated beforehand
 * 
 *  - id
 *  - Indices for the quarklines and gamma structure
 *  - Paths and information for IO
 */
struct CorrInfo{
  size_t id;
  std::string hdf5_dataset_name;
  std::vector<size_t> lookup;
  std::vector<int> gamma;
  /*! Just a small constructor to ensure easy filling of its vector form */
  CorrInfo(const size_t id,
           const std::string &hdf5_dataset_name,
           const std::vector<size_t> &lookup,
           const std::vector<int> &gamma)
      : id(id), hdf5_dataset_name(hdf5_dataset_name), lookup(lookup), gamma(gamma){};
};

/******************************************************************************/
/*! Contains information on all correlators
 *  
 *  @todo modular programming looks different
 */
struct CorrelatorLookup{
  std::vector<CorrInfo> C1;
  std::vector<CorrInfo> C1T;

  std::vector<CorrInfo> corr0;
  std::vector<CorrInfo> C20;
  std::vector<CorrInfo> C40D;
  std::vector<CorrInfo> C40V;

  std::vector<CorrInfo> corrC;
  std::vector<CorrInfo> C2c;
  std::vector<CorrInfo> C4cD;
  std::vector<CorrInfo> C4cV;

  std::vector<CorrInfo> C30;
  std::vector<CorrInfo> C3c;

  std::vector<CorrInfo> C40C;
  std::vector<CorrInfo> C4cC;
  std::vector<CorrInfo> C40B;
  std::vector<CorrInfo> C4cB;
};

/*! typetrait class which allows to use QuarklineQ1Indices for Q1 and
 *  QuarklineQ2Indices for Q2L and Q2V
 */
template <QuarkLineType qlt>
struct QuarkLineIndices {};

/*! @todo QuarkLineType is a bad name in this case. That's a proxy for
 *        CorrInfo.lookup
 */
template <>
struct QuarkLineIndices<QuarkLineType::Q0> {
  typedef std::vector<VdaggerVRandomLookup> type;
};

template <>
struct QuarkLineIndices<QuarkLineType::Q1> {
  typedef std::vector<QuarklineQ1Indices> type;
};

template <>
struct QuarkLineIndices<QuarkLineType::Q2L> {
  typedef std::vector<QuarklineQ2Indices> type;
};

template <>
struct QuarkLineIndices<QuarkLineType::Q2V> {
  typedef std::vector<QuarklineQ2Indices> type;
};
