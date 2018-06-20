#include <typedefs.h>

#include "boost/multi_array.hpp"

template <size_t rvecs>
struct DilutedTraceCollection3{
  DilutedTraceCollection3(size_t const size0,
                          ssize_t const Lt){

  tr.resize(boost::extents[size0][Lt][Lt]);
  }

//  Value const &operator[](Key const &key) {
//    if (Ql.count(key) == 0) {
//      build(key);
//    }
//
//    return Ql.at(key);
//  }

//  void build(Key const &time_key);

  void clear(){
    return;
  }

  boost::multi_array<std::vector<DilutedTrace<rvecs>>, 3> tr;
};

template <size_t rvecs>
struct DilutedTraceCollection2{
  DilutedTraceCollection2(size_t const size0,
                          ssize_t const Lt){

    tr.resize(boost::extents[size0][Lt]);
  }

  void clear(){
    return;
  }

  boost::multi_array<std::vector<DilutedTrace<rvecs>>, 2> tr;
};


