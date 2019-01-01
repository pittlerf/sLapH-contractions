#include "Correlators.hpp"
#include "git.hpp"
#include "global_data.hpp"
#include "global_data_build_IO_names.hpp"

#include <omp.h>
#include <boost/archive/binary_oarchive.hpp>

#include <iostream>

int main(int ac, char *av[]) {
  std::cout << "This is sLapH-contractions:\n"
            << "  git branch " << git_refspec << "\n"
            << "  git revision " << git_sha1 << "\n"
            << "  git state " << git_changes << "\n"
            << "  compiled by " << git_user << " on " << git_host << "\n"
            << "  running with up to " << omp_get_max_threads() << " OpenMP threads\n";

  // reading global parameters from input file
  GlobalData gd;

  read_parameters(gd, ac, av);

  std::ofstream ofs("global_data.bin", std::ios::binary);
  boost::archive::binary_oarchive archive(ofs);
  archive << gd;
}
