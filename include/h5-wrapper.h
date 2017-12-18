#pragma once

#include "typedefs.h"

#include "H5Cpp.h"


void create_folder(std::string const &path);

class HDF5Handle {
 public:
  HDF5Handle(std::string const &output_path,
             std::string const &diagram,
             std::string const &output_filename) {
    H5::Exception::dontPrint();

    create_folder(output_path);
    H5std_string const file_name((output_path + "/" + diagram + output_filename).c_str());
    file_ = H5::H5File(file_name, H5F_ACC_TRUNC);
  }

  H5::DataSet create_data_set(std::string const &dataset_name,
                              H5::CompType const &comp_type,
                              H5::DataSpace const &data_space) {
    return file_.createDataSet(dataset_name.c_str(), comp_type, data_space);
  }

  H5::Group create_group(std::string const &name) {
    return file_.createGroup(name.c_str());
  }

 private:
  H5::H5File file_;
};

template <typename Payload>
void write_homogenious(H5::DataSet &handle, Payload const &payload);

template <typename Payload>
void write_heterogenious(H5::Group &handle, Payload const &payload);

template <typename T>
H5::CompType make_comp_type();
