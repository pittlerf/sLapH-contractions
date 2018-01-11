#include "h5-wrapper.h"

#include "DilutedFactor.h"

#include "boost/filesystem.hpp"

void create_folder(std::string const &path) {
  if (access(path.c_str(), 0) != 0) {
    std::cout << "\tdirectory " << path << " does not exist and will be created";
    boost::filesystem::path dir(path.c_str());
    if (boost::filesystem::create_directories(dir))
      std::cout << "\tSuccess" << std::endl;
    else
      std::cout << "\tFailure" << std::endl;
  }
}

struct Complex1Times {
  int t;
  double re;
  double im;
};

template <>
H5::CompType make_comp_type<cmplx>() {
  H5::CompType cmplx_w(2 * sizeof(double));
  auto type = H5::PredType::NATIVE_DOUBLE;
  cmplx_w.insertMember("re", HOFFSET(complex_t, re), type);
  cmplx_w.insertMember("im", HOFFSET(complex_t, im), type);

  return cmplx_w;
}

template <>
H5::CompType make_comp_type<compcomp_t>() {
  H5::CompType cmplxcmplx_w(4 * sizeof(double));
  auto type = H5::PredType::NATIVE_DOUBLE;
  cmplxcmplx_w.insertMember("rere", HOFFSET(compcomp_t, rere), type);
  cmplxcmplx_w.insertMember("reim", HOFFSET(compcomp_t, reim), type);
  cmplxcmplx_w.insertMember("imre", HOFFSET(compcomp_t, imre), type);
  cmplxcmplx_w.insertMember("imim", HOFFSET(compcomp_t, imim), type);

  return cmplxcmplx_w;
}

template <>
H5::CompType make_comp_type<Complex1Times>() {
  H5::CompType comp_type(sizeof(Complex1Times));
  auto const type_int = H5::PredType::NATIVE_INT;
  auto const type = H5::PredType::NATIVE_DOUBLE;
  comp_type.insertMember("t", HOFFSET(Complex1Times, t), type_int);
  comp_type.insertMember("re", HOFFSET(Complex1Times, re), type);
  comp_type.insertMember("im", HOFFSET(Complex1Times, im), type);

  return comp_type;
}

template <>
void write_homogenious(H5::DataSet &data_set, std::vector<cmplx> const &payload) {
  data_set.write(payload.data(), make_comp_type<cmplx>());
}

template <>
void write_homogenious(H5::DataSet &data_set, std::vector<compcomp_t> const &payload) {
  data_set.write(payload.data(), make_comp_type<compcomp_t>());
}

/*! Creates compound datatype to write complex numbers from complex_t
 *  vectors to HDF5 file
 *
 *  @Returns cmplx_w   HDF5 compound datatype for complex numbers
 */
template <>
H5::CompType comp_type_factory<cmplx>() {
  H5::CompType cmplx_w(2 * sizeof(double));
  auto type = H5::PredType::NATIVE_DOUBLE;
  cmplx_w.insertMember("re", HOFFSET(complex_t, re), type);
  cmplx_w.insertMember("im", HOFFSET(complex_t, im), type);

  return cmplx_w;
}

/*! Creates compound datatype to write complex numbers from compcomp_t
 *  vectors to HDF5 file
 *
 *  @Returns cmplx_w   HDF5 compound datatype for structs of four doubles
 */
template <>
H5::CompType comp_type_factory<compcomp_t>() {
  H5::CompType cmplxcmplx_w(4 * sizeof(double));
  auto type = H5::PredType::NATIVE_DOUBLE;
  cmplxcmplx_w.insertMember("rere", HOFFSET(compcomp_t, rere), type);
  cmplxcmplx_w.insertMember("reim", HOFFSET(compcomp_t, reim), type);
  cmplxcmplx_w.insertMember("imre", HOFFSET(compcomp_t, imre), type);
  cmplxcmplx_w.insertMember("imim", HOFFSET(compcomp_t, imim), type);

  return cmplxcmplx_w;
}
