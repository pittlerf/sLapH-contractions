#include "h5-wrapper.h"

#include "boost/filesystem.hpp"

#include <iostream>
#include <map>
#include <sstream>
#include <string>

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

#if 0
template <>
void write_heterogenious(H5::Group &group, std::vector<DilutedTrace> const &payload) {
  std::map<std::string, std::vector<cmplx>> data;

  for (auto const &elem : payload) {
    std::ostringstream oss;
    bool first = false;
    for (auto const &id : elem.used_rnd_ids) {
      if (!first) {
        first = true;
        oss << "_";
      }
      oss << "rnd" << id;
    }
    std::string const key = oss.str();

    data[key].push_back(elem.data);
  }

  for (auto const &pair : data) {
    hsize_t dim(pair.second.size());
    H5::DataSpace dspace(1, &dim);
    auto data_set =
        group.createDataSet(pair.first.c_str(), make_comp_type<cmplx>(), dspace);
    data_set.write(pair.second.data(), make_comp_type<cmplx>());
  }
}

template <int rvecs>
void write_heterogenious(
    H5::Group &group,
    boost::detail::multi_array::sub_array<std::vector<DilutedTrace<rvecs>>, 1ul> const
        &payload) {
  std::map<std::string, std::vector<Complex1Times>> data;

  for (int t = 0; t < payload.size(); ++t) {
    for (auto const &elem : payload[t]) {
      std::ostringstream oss;
      bool first = true;
      for (auto const &id : elem.used_rnd_ids) {
        if (!first) {
          first = false;
          oss << "_";
        }
        oss << "rnd";
        // The `static_cast<int>` is needed because otherwise the `uint8_t` is interpreted
        // as `signed char` and then this will define the string representation. We want
        // it interpreted as an integer.
        oss << static_cast<int>(id);
      }
      std::string const key = oss.str();

      data[key].push_back(Complex1Times{t, elem.data.real(), elem.data.imag()});
    }
  }

  auto const &comp_type = make_comp_type<Complex1Times>();

  for (auto const &pair : data) {
    hsize_t dim(pair.second.size());
    H5::DataSpace dspace(1, &dim);
    auto data_set = group.createDataSet(pair.first.c_str(), comp_type, dspace);
    data_set.write(pair.second.data(), comp_type);
  }
}
#endif
