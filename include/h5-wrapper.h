#pragma once

#include "typedefs.h"

#include "H5Cpp.h"

#include <iostream>
#include <map>
#include <sstream>
#include <string>

#define PRINT(x) std::cout << #x << ": " << (x) << std::endl;

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

template <size_t rvecs>
void write_heterogenious(
    H5::Group &group,
    boost::detail::multi_array::sub_array<std::vector<DilutedTrace<rvecs>>, 1ul> const
        &payload) {
  std::map<std::string, std::vector<cmplx>> data;

  for (int t = 0; t < payload.size(); ++t) {
    PRINT(payload[t].size());
    auto const &subs = sub_accumulate(payload[t]);
    PRINT(subs.size());

    for (auto const &elem : subs) {
      std::ostringstream oss;
      bool first = true;
      for (auto const &id : elem.first) {
        if (first) {
          first = false;
        } else {
          oss << "_";
        }

        oss << "rnds";
        // The `static_cast<int>` is needed because otherwise the `uint8_t` is interpreted
        // as `signed char` and then this will define the string representation. We want
        // it interpreted as an integer.
        oss << static_cast<int>(id);
      }
      std::string const key = oss.str();

      std::cout << key << " " << elem.second << std::endl;

      data[key].push_back(elem.second);
    }
  }

  auto const &comp_type = make_comp_type<cmplx>();

  for (auto const &pair : data) {
    hsize_t dim(pair.second.size());
    H5::DataSpace dspace(1, &dim);
    auto data_set = group.createDataSet(pair.first.c_str(), comp_type, dspace);
    data_set.write(pair.second.data(), comp_type);
  }
}
