#include "h5-wrapper.h"

#include "boost/filesystem.hpp"

#include <iostream>

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
void write(HDF5Handle &handle, std::vector<cmplx> const &payload) {
}
