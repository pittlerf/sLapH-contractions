#include "h5-wrapper.hpp"

#include "DilutedFactor.hpp"

#include <boost/filesystem.hpp>

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

template <>
H5::CompType make_comp_type<Complex>() {
  H5::CompType cmplxcmplx_w(2 * sizeof(double));
  auto type = H5::PredType::NATIVE_DOUBLE;
  cmplxcmplx_w.insertMember("re", 0 * sizeof(double), type);
  cmplxcmplx_w.insertMember("im", 1 * sizeof(double), type);

  return cmplxcmplx_w;
}

template <>
void write_homogenious(H5::DataSet &data_set, std::vector<Complex> const &payload) {
  data_set.write(payload.data(), make_comp_type<Complex>());
}
