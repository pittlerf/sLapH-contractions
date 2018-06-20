#include "boost/multi_array.hpp"

#include "DilutedFactorY.h"

template <DilutedFactorType qlt1, DilutedFactorType qlt2, size_t rvecs>
struct DilutedTraceCollection{
  public: 

  /*! num_times is the sum of times of contained factors -1 for each continuity 
   *  condition of the quarkline diagram 
   */
  static constexpr int num_times = DilutedFactorTypeTraits<qlt1>::num_times + 
    DilutedFactorTypeTraits<qlt2>::num_times - 2;

  using Key = std::array<int, num_times>;
  using Value = std::vector<DilutedTrace<rvecs>>;

  DilutedTraceCollection(size_t const size0, ssize_t const Lt){
    tr.resize(boost::extents[Lt][Lt][size0]);
  }

  Value const &at(Key const &key, size_t const size0) {
//    if (Ql.count(key) == 0) {
//      build(key);
//    }

    return tr[key[0]][key[1]][size0];
  }

  void build(DilutedFactorFactory<qlt1> &q,
      DilutedFactorFactory<qlt2> &q_bla,
      DiagramIndex const &c_look,
                 int const t1,
                 int const t2,
                 int const b1,
                 int const b2);

  void clear(){
    return;
  }

  private:
  boost::multi_array<std::vector<DilutedTrace<rvecs>>, num_times + 1> tr;

};

template <DilutedFactorType qlt, size_t rvecs>
struct DilutedTraceCollection2{

  static constexpr int num_times = DilutedFactorTypeTraits<qlt>::num_times - 1;

  using Key = std::array<int, num_times>;
  using Value = DilutedTrace<rvecs>;

  DilutedTraceCollection2(size_t const size0,
                          ssize_t const Lt){

    tr.resize(boost::extents[size0][Lt]);
  }

  void clear(){
    return;
  }

  boost::multi_array<std::vector<DilutedTrace<rvecs>>, 1 + num_times> tr;
};


