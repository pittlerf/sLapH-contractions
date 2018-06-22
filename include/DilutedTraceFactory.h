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
  using Value = DilutedTraces<rvecs>;

  DilutedTraceCollection(DilutedFactorFactory<qlt1> &_df1,
       DilutedFactorFactory<qlt2> &_df2) :
    df1(_df1),
    df2(_df2){
  }

  Value const &at(Key const &key) {
//    if (Ql.count(key) == 0) {
//      build(key);
//    }

    return tr[{key[0],key[1]}];
  }

  void build(DiagramIndex const &c_look,
                 int const t1,
                 int const t2,
                 int const b1,
                 int const b2);

  void clear(){
    return;
  }

  private:
  std::map<Key, Value> tr;
  DilutedFactorFactory<qlt1> df1;
  DilutedFactorFactory<qlt2> df2;
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


