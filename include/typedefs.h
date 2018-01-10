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
#include <list>
#include <vector>

#include "Eigen/Core"
#include "Eigen/Dense"
#include "boost/container/static_vector.hpp"
#include "boost/multi_array.hpp"

int constexpr max_rnd_ids = 10;

using RndId = int8_t;

typedef std::complex<double> cmplx;

template <size_t rvecs>
using SmallVectorRndId = boost::container::static_vector<RndId, rvecs>;

template <size_t rvecs>
struct DilutedFactor {
  using Data = Eigen::MatrixXcd;

  Data data;
  std::pair<RndId, RndId> ric;
  boost::container::static_vector<RndId, rvecs> used_rnd_ids;
};

template <size_t rvecs>
struct DilutedTrace {
  using Data = cmplx;

  Data data;
  boost::container::static_vector<RndId, rvecs> used_rnd_ids;
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

// TODO (Martin Ueding 2017-12-21): Rename this to DilutedTraceCollection3.
template <size_t rvecs>
using DilutedTraceCollection = boost::multi_array<std::vector<DilutedTrace<rvecs>>, 3>;

template <size_t rvecs>
using DilutedTraceCollection2 = boost::multi_array<std::vector<DilutedTrace<rvecs>>, 2>;

/*! Special type for Quarklines */
typedef boost::multi_array<std::vector<Eigen::MatrixXcd>, 3> array_quarkline;

// typedef boost::multi_array<std::vector<std::vector<cmplx> >, 4> array_corr;
/*! @TODO {Is that deprecated?} */
typedef boost::multi_array<std::vector<cmplx>, 2> array_C1;

// Eigen typedefs
typedef std::vector<Eigen::MatrixXcd> vec_Xcd_eigen;
/*! Data type for rvdaggerv and rvdaggervr */
typedef std::vector<std::vector<std::vector<Eigen::MatrixXcd>>> Xcd_d3_eigen;
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

  compcomp_t() : rere(0.0), reim(0.0), imre(0.0), imim(0.0) {}

  compcomp_t(const double rere, const double reim, const double imre, const double imim)
      : rere(rere), reim(reim), imre(imre), imim(imim) {}

  compcomp_t &operator+=(compcomp_t const &other) {
    rere += other.rere;
    reim += other.reim;
    imre += other.imre;
    imim += other.imim;

    return *this;
  }

  compcomp_t &operator/=(double const &right) {
    rere /= right;
    reim /= right;
    imre /= right;
    imim /= right;

    return *this;
  }

  compcomp_t operator/(double const &right) const {
    auto res = *this;
    res /= right;
    return res;
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
struct VdaggerVQuantumNumbers {
  size_t id;
  std::array<int, 3> momentum;     /*!< The -momentum as 3-vector */
  std::array<int, 3> displacement; /*!< The displacement as 3-vector */

  /*! Constructor */
  VdaggerVQuantumNumbers(const size_t id,
                         const std::array<int, 3> &momentum,
                         const std::array<int, 3> &displacement)
      : id(id), momentum(momentum), displacement(displacement){};
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
struct RandomIndexCombinationsQ2 {
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
  std::vector<std::pair<size_t, size_t>> rnd_vec_ids;

  /*! Just a small constructor to ensure easy filling of its vector form */
  RandomIndexCombinationsQ2(const size_t id,
                            const size_t id_q1,
                            const size_t id_q2,
                            const std::pair<size_t, size_t> offset,
                            const std::vector<std::pair<size_t, size_t>> &rnd_vec_ids)
      : id(id), id_q1(id_q1), id_q2(id_q2), offset(offset), rnd_vec_ids(rnd_vec_ids){};
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
struct OperatorLookup {
  /*! Specifies physical content of quark field operator (with Dirac structure
   *  factored out)
   */
  std::vector<VdaggerVQuantumNumbers> vdaggerv_lookup;

  /*! Specifies the random vector indices in case of rvdaggervr */
  std::vector<RandomIndexCombinationsQ2> ricQ2_lookup;

  /*! For @f$ \vec{p} = 0 @f$, @f$ V^\dagger exp(ipx) V = \mathbb{1} @f$. If
   *  applicable This contains the index of @em vdaggerv_lookup where it can
   *  be replaced by a unit matrix. If @f$ \vec{p} = 0 @f$ is not needed,
   *  @em index_of_unity is set to -1
   */
  int index_of_unity;
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
/*! Indices needed to uniquely identify Q1 objects
 *
 *  Because vdaggerv is diagonal in dirac space, gamma may be factored out and
 *  it proves useful to calculate and reuse
 *
 *    Q1 = rvdaggerv * gamma * peram
 */
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
struct QuarklineIndices {
  /*! Identifies physical content of @f$ V^dagger V @f$ */
  size_t id_vdaggerv;
  /*! Flag that indicates whether VdaggerV must be daggered (prior to
   *  multiplication with random vectors) to get the correct quantum numbers
   */
  bool need_vdaggerv_daggering;

  /*! List of necessarry gamma combinations */
  std::vector<int> gamma;

  /*! The entries of the pair correspond to the first and second random index.
   *  List of all possible combinations of random vector indices for quarks
   *  specified by @em id_q1 and @em id_q2
   */
  std::vector<std::pair<size_t, size_t>> rnd_vec_ids;
};

inline bool operator==(QuarklineIndices const &first, QuarklineIndices const second) {
  if ((first.id_vdaggerv == second.id_vdaggerv) &&
      (first.need_vdaggerv_daggering == second.need_vdaggerv_daggering) &&
      (first.gamma == second.gamma) && (first.rnd_vec_ids == second.rnd_vec_ids))
    return true;
  else
    return false;
}
/******************************************************************************/
/*! Maps index from CorrelatorLookup to QuarklineQ1Indicies or
 *  QuarklineQ2Indicies, depending on the quarkline needed in Correlator
 *
 *  @deprecated No longer necessary in memory_optimised branch as the quarklines
 *              are built on the fly
 */
struct QuarklineLookup {
  std::vector<QuarklineIndices> Q0;
  std::vector<QuarklineIndices> Q1;
  std::vector<QuarklineIndices> Q2V;
  std::vector<QuarklineIndices> Q2L;
};

// Q0 formerly called rVdaggerVr
enum class QuarkLineType { Q0, Q1, Q2, Q2L, Q2V };

/******************************************************************************/
/*! All information needed to build and write the correlator given the
 *  quarklines were calculated beforehand
 *
 *  - id
 *  - Indices for the quarklines and gamma structure
 *  - Paths and information for IO
 */
struct CorrInfo {
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
struct CorrelatorLookup {
  std::vector<CorrInfo> C1;
  std::vector<CorrInfo> C1T;
  std::vector<CorrInfo> C20V;

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

template <QuarkLineType qlt>
struct QuarkLineIndices {};

/*! @todo QuarkLineType is a bad name in this case. That's a proxy for
 *        CorrInfo.lookup
 */
template <>
struct QuarkLineIndices<QuarkLineType::Q0> {
  typedef std::vector<QuarklineIndices> type;
  static size_t constexpr num_times = 1;
};

template <>
struct QuarkLineIndices<QuarkLineType::Q1> {
  typedef std::vector<QuarklineIndices> type;
  static size_t constexpr num_times = 2;
};

template <>
struct QuarkLineIndices<QuarkLineType::Q2> {
  typedef std::vector<QuarklineIndices> type;
  static size_t constexpr num_times = 3;
};

template <>
struct QuarkLineIndices<QuarkLineType::Q2L> {
  typedef std::vector<QuarklineIndices> type;
  static size_t constexpr num_times = 3;
};

template <>
struct QuarkLineIndices<QuarkLineType::Q2V> {
  typedef std::vector<QuarklineIndices> type;
  static size_t constexpr num_times = 3;
};

#define MU_DEBUG(x) std::cout << std::setw(20) << #x << ": " << (x) << std::endl;
