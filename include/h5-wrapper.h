#pragma once

#include "typedefs.h"

#include "H5Cpp.h"


void create_folder(std::string const &path);

class HDF5Handle {
 public:
  HDF5Handle(std::string const &output_path,
             std::string const &diagram,
             std::string const &output_filename,
             H5::CompType const &_comp_type) {

    H5::Exception::dontPrint();

    create_folder(output_path);
    H5std_string const file_name((output_path + "/" + diagram + output_filename).c_str());
    file_ = H5::H5File(file_name, H5F_ACC_TRUNC);
  }

 private:
  H5::H5File file_;
};

template <typename Payload>
void write(HDF5Handle &handle, Payload const &payload);
