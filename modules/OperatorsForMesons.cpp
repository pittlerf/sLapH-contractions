/*! @file OperatorsForMesons.cpp
 *  Class definition of OperatorsForMesons
 *
 *  @author Bastian Knippschild
 *  @author Markus Werner
 */

#include "OperatorsForMesons.h"

#include <boost/format.hpp>

#include <iomanip>

namespace {

/*! Creates a two-dimensional vector containing the momenta for the operators
 *
 *  @param[in] Lx, Ly, Lz      Lattice extent in spatial directions
 *  @param[in] vdaggerv_lookup Contains the momenta
 *  @param[in,out] momentum    Two dimensional array where the momenta are
 *                             stored
 */
void create_momenta(const size_t Lx,
                    const size_t Ly,
                    const size_t Lz,
                    const std::vector<VdaggerVQuantumNumbers> &vdaggerv_lookup,
                    array_cd_d2 &momentum) {
  static const std::complex<double> I(0.0, 1.0);

  /*! To calculate Vdagger exp(i*p*x) V only the momenta corresponding to the
   *  quantum number id in op_VdaggerV will be used. The rest can be obtained
   *  by adjoining
   */
  for (const auto &op : vdaggerv_lookup) {
    // op_VdaggerV contains the index of one (redundancy) op_Corr which
    // allows to deduce the quantum numbers (momentum)
    const double ipx = op.momentum[0] * 2. * M_PI / (double)Lx;
    const double ipy = op.momentum[1] * 2. * M_PI / (double)Ly;
    const double ipz = op.momentum[2] * 2. * M_PI / (double)Lz;
    // calculate \vec{p} \cdot \vec{x} for all \vec{x} on the lattice
    for (int x = 0; x < Lx; ++x) {
      const int xH = x * Ly * Lz;   // helper variable
      const double ipxH = ipx * x;  // helper variable
      for (int y = 0; y < Ly; ++y) {
        const int xHyH = xH + y * Lz;            // helper variable
        const double ipxHipyH = ipxH + ipy * y;  // helper variable
        for (int z = 0; z < Lz; ++z) {
          // multiply \vec{p} \cdot \vec{x} with complex unit and exponentiate
          momentum[op.id][xHyH + z] = exp(-I * (ipxHipyH + ipz * z));
        }
      }
    }  // loops over spatial vectors end here
  }    // loop over redundant quantum numbers ends here
}

/******************************************************************************/
void write_vdaggerv(const std::string &pathname,
                    const std::string &filename,
                    const Eigen::MatrixXcd &Vt) {
  // writing the data
  std::ofstream file((pathname + filename).c_str(),
                     std::ofstream::binary | std::ofstream::trunc);

  if (file.is_open()) {
    std::cout << "\twriting VdaggerV to file:" << pathname + filename << std::endl;
    // buffer for writing
    std::vector<Complex> eigen_vec(Vt.size());
    for (size_t ncol = 0; ncol < Vt.cols(); ncol++) {
      for (size_t nrow = 0; nrow < Vt.rows(); nrow++) {
        eigen_vec.at(ncol * Vt.rows() + nrow) = (Vt)(nrow, ncol);
      }
    }
    file.write(reinterpret_cast<const char *>(&eigen_vec[0]),
               Vt.size() * sizeof(Complex));
    if (!file.good())
      std::cout << "Problems while write to " << (pathname + filename).c_str()
                << std::endl;
    file.close();
  } else
    std::cout << "can't open " << (pathname + filename).c_str() << std::endl;
}

}  // namespace

/******************************************************************************/
/*!
 * @param Lt, Lx, Ly, Lz  Temporal and spatial lattice extent
 * @param nb_ev           Number of eigenvectors
 * @param dilE            Number of diluted blocks in eigenvector space
 * @param operator_lookuptable ?
 * @param handling_vdaggerv
 * @param path_vdaggerv
 *
 * The initialization of the container attributes of OperatorsForMesons
 * is done in the member initializer list of the constructor. The allocation
 * of heap memory is delegated to boost::multi_array::resize
 */
OperatorFactory::OperatorFactory(const size_t Lt,
                                 const size_t Lx,
                                 const size_t Ly,
                                 const size_t Lz,
                                 const size_t nb_ev,
                                 const size_t dilE,
                                 const OperatorLookup &operator_lookuptable,
                                 const std::string &handling_vdaggerv,
                                 const std::string &path_vdaggerv,
                                 const std::string &path_config)
    : vdaggerv(),
      momentum(),
      operator_lookuptable(operator_lookuptable),
      Lt(Lt),
      Lx(Lx),
      Ly(Ly),
      Lz(Lz),
      nb_ev(nb_ev),
      dilE(dilE),
      handling_vdaggerv(handling_vdaggerv),
      path_vdaggerv(path_vdaggerv),
      path_config(path_config) {

  // resizing containers to their correct size
  vdaggerv.resize(boost::extents[operator_lookuptable.vdaggerv_lookup.size()][Lt]);

  // the momenta only need to be calculated for a subset of quantum numbers
  // (see VdaggerV::build_vdaggerv)
  momentum.resize(
      boost::extents[operator_lookuptable.vdaggerv_lookup.size()][Lx * Ly * Lz]);
  create_momenta(Lx, Ly, Lz, operator_lookuptable.vdaggerv_lookup, momentum);

  std::cout << "\tMeson operators initialised" << std::endl;
}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void OperatorFactory::build_vdaggerv(const std::string &filename, const int config) {
  clock_t t2 = clock();
  const size_t dim_row = 3 * Lx * Ly * Lz;
  const int id_unity = operator_lookuptable.index_of_unity;

  // prepare full path for writing
  std::string const full_path =
      (boost::format("/%s/cnfg%04d/") % path_vdaggerv % config).str();

  // check if directory exists
  if (handling_vdaggerv == "write" && access(full_path.c_str(), 0) != 0) {
    std::cout << "\tdirectory " << full_path.c_str()
              << " does not exist and will be created";
    boost::filesystem::path dir(full_path.c_str());
    if (!boost::filesystem::create_directories(dir))
      std::cout << "\tSuccess" << std::endl;
    else
      std::cout << "\tFailure" << std::endl;
  }

  // resizing each matrix in vdaggerv
  // TODO: check if it is better to use for_each and resize instead of std::fill
  std::fill(vdaggerv.origin(),
            vdaggerv.origin() + vdaggerv.num_elements(),
            Eigen::MatrixXcd::Zero(nb_ev, nb_ev));

#pragma omp parallel
  {
    // Read in gauge field from lime configuration, need all directions
    //GaugeField gauge = GaugeField(Lt, Lx, Ly, Lz, PATH_GAUGE_IN, 0, Lt-1, 4);
    //gauge.read_gauge_field(CONFIG,0,Lt-1);
    Eigen::VectorXcd mom = Eigen::VectorXcd::Zero(dim_row);
    Eigen::VectorXcd dis = Eigen::VectorXcd::Zero(dim_row);
    // Check if gauges are needed
    bool need_gauge = false;
    for (auto const& op : operator_lookuptable.vdaggerv_lookup)
        if(!op.displacement.empty()) need_gauge=true;
    GaugeField gauge = GaugeField(Lt, Lx, Ly, Lz, path_config,0,Lt-1,4);
    if(need_gauge) gauge.read_gauge_field(config,0,Lt-1); 
    EigenVector V_t(1, dim_row, nb_ev);  // each thread needs its own copy
#pragma omp for schedule(dynamic)
    for (ssize_t t = 0; t < Lt; ++t) {
      // creating full filename for eigenvectors and reading them in
      if (!(operator_lookuptable.vdaggerv_lookup.size() == 1 &&
            operator_lookuptable.vdaggerv_lookup[0].id == id_unity)) {
        auto const inter_name = (boost::format("%s%03d") % filename % t).str();
        V_t.read_eigen_vector(inter_name.c_str(), 0, 0);  // reading eigenvectors
      }

      // VdaggerV is independent of the gamma structure and momenta connected by
      // sign flip are related by adjoining VdaggerV. Thus the expensive
      // calculation must only be performed for a subset of quantum numbers given
      // in op_VdaggerV.
      for (const auto &op : operator_lookuptable.vdaggerv_lookup) {
        // For zero momentum and displacement VdaggerV is the unit matrix, thus
        // the calculation is not performed
        if (op.id != id_unity) {
          // Forward derivative
          Eigen::MatrixXcd W_t;
          if (!op.displacement.empty()){
            W_t = gauge.displace_eigenvectors(V_t[0],t,op.displacement,1);
            vdaggerv[op.id][t] = V_t[0].adjoint() * W_t;
          }
          else {
            W_t = V_t[0];
          }
          // momentum vector contains exp(-i p x). Divisor 3 for colour index.
          // All three colours on same lattice site get the same momentum.
          for (size_t x = 0; x < dim_row; ++x) {
            mom(x) = momentum[op.id][x / 3];
          }
          vdaggerv[op.id][t] = V_t[0].adjoint() * mom.asDiagonal() * W_t;
          // writing vdaggerv to disk
          if (handling_vdaggerv == "write") {
            auto const dummy2 = (boost::format("operators.%04d.p_") % config).str();
            std::string dummy = dummy2 + std::to_string(op.momentum[0]) +
                                std::to_string(op.momentum[1]) +
                                std::to_string(op.momentum[2]);

            auto const outfile = (boost::format("%s_.t_%03d") % dummy % t).str();
            write_vdaggerv(full_path, outfile, vdaggerv[op.id][t]);
          }
        } else  // zero momentum
          vdaggerv[op.id][t] = Eigen::MatrixXcd::Identity(nb_ev, nb_ev);
      }
    }  // loop over time
  }    // pragma omp parallel ends here

  t2 = clock() - t2;
  std::cout << std::setprecision(1) << "\t\t\tSUCCESS - " << std::fixed
            << ((float)t2) / CLOCKS_PER_SEC << " seconds" << std::endl;
  is_vdaggerv_set = true;
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void OperatorFactory::read_vdaggerv(const int config) {
  clock_t t2 = clock();
  const size_t dim_row = 3 * Lx * Ly * Lz;
  const int id_unity = operator_lookuptable.index_of_unity;

  // prepare full path for reading
  auto const full_path =
      (boost::format("/%s/cnfg%04d/operators.%04d") % path_vdaggerv % config % config)
          .str();

  // resizing each matrix in vdaggerv
  std::fill(vdaggerv.origin(),
            vdaggerv.origin() + vdaggerv.num_elements(),
            Eigen::MatrixXcd::Zero(nb_ev, nb_ev));

#pragma omp parallel
  {
#pragma omp for schedule(dynamic)
    for (ssize_t t = 0; t < Lt; ++t) {
      for (const auto &op : operator_lookuptable.vdaggerv_lookup) {
        // For zero momentum and displacement VdaggerV is the unit matrix, thus
        // the calculation is not performed
        if (op.id != id_unity) {
          // creating full filename for vdaggerv and reading them in
          std::string dummy = full_path + ".p_" + std::to_string(op.momentum[0]) +
                              std::to_string(op.momentum[1]) +
                              std::to_string(op.momentum[2]);

          auto const infile = (boost::format("%s_.t_%03d") % dummy % t).str();

          // writing the data
          std::ifstream file(infile, std::ifstream::binary);

          if (file.is_open()) {
            std::cout << "\treading VdaggerV from file:" << infile << std::endl;

            // buffer for reading
            std::vector<Complex> eigen_vec(vdaggerv[op.id][t].size());
            file.read(reinterpret_cast<char *>(&eigen_vec[0]),
                      vdaggerv[op.id][t].size() * sizeof(Complex));
            for (size_t ncol = 0; ncol < vdaggerv[op.id][t].cols(); ncol++) {
              for (size_t nrow = 0; nrow < vdaggerv[op.id][t].rows(); nrow++) {
                (vdaggerv[op.id][t])(nrow, ncol) =
                    eigen_vec.at(ncol * vdaggerv[op.id][t].rows() + nrow);
              }
            }
            if (!file.good()) {
              std::cout << "Problems while reading from " << infile << std::endl;
              exit(0);
            }
            file.close();
          } else {
            std::cout << "can't open " << infile << std::endl;
            exit(0);
          }
        } else  // zero momentum
          vdaggerv[op.id][t] = Eigen::MatrixXcd::Identity(nb_ev, nb_ev);
      }
    }  // loop over time
  }    // pragma omp parallel ends here

  t2 = clock() - t2;
  std::cout << std::setprecision(1) << "\t\t\tSUCCESS - " << std::fixed
            << ((float)t2) / CLOCKS_PER_SEC << " seconds" << std::endl;
  is_vdaggerv_set = true;
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void OperatorFactory::read_vdaggerv_liuming(const int config) {
  clock_t t2 = clock();
  const size_t dim_row = 3 * Lx * Ly * Lz;
  const int id_unity = operator_lookuptable.index_of_unity;

  // prepare full path for reading
  auto const full_path = (boost::format("/%s/VdaggerV.") % path_vdaggerv).str();

  // resizing each matrix in vdaggerv
  std::fill(vdaggerv.origin(),
            vdaggerv.origin() + vdaggerv.num_elements(),
            Eigen::MatrixXcd::Zero(nb_ev, nb_ev));

#pragma omp parallel
  {
  //  #pragma omp for schedule(dynamic)
  //    for(const auto& op : operator_lookuptable.vdaggerv_lookup){
#pragma omp for schedule(dynamic)
    for (size_t i = 0; i < operator_lookuptable.vdaggerv_lookup.size(); ++i) {
      const auto op = (operator_lookuptable.vdaggerv_lookup[i]);
      // For zero momentum and displacement VdaggerV is the unit matrix, thus
      // the calculation is not performed
      if (op.id != id_unity) {
        // creating full filename for vdaggerv and reading them in
        // both possibilities must be checked
        std::string dummy1 = full_path + "p" + std::to_string(-op.momentum[0]) + "p" +
                             std::to_string(-op.momentum[1]) + "p" +
                             std::to_string(-op.momentum[2]) + ".conf";
        auto const infile1 = (boost::format("%s%04d") % dummy1 % config).str();
        std::ifstream file1(infile1, std::ifstream::binary);

        // second possibility for a name
        std::string dummy2 = full_path + "p" + std::to_string(op.momentum[0]) + "p" +
                             std::to_string(op.momentum[1]) + "p" +
                             std::to_string(op.momentum[2]) + ".conf";
        auto const infile2 = (boost::format("%s%04d") % dummy2 % config).str();
        std::ifstream file2(infile2, std::ifstream::binary);

        if (file1.is_open()) {
          std::cout << "\treading VdaggerV from file:" << infile1 << std::endl;
          for (size_t t = 0; t < Lt; ++t) {
            // buffer for reading
            std::vector<Complex> eigen_vec(vdaggerv[op.id][t].size());
            file1.read(reinterpret_cast<char *>(&eigen_vec[0]),
                       vdaggerv[op.id][t].size() * sizeof(Complex));
            for (size_t ncol = 0; ncol < vdaggerv[op.id][t].cols(); ncol++) {
              for (size_t nrow = 0; nrow < vdaggerv[op.id][t].rows(); nrow++) {
                (vdaggerv[op.id][t])(nrow, ncol) =
                    eigen_vec.at(nrow * vdaggerv[op.id][t].cols() + ncol);
              }
            }
            vdaggerv[op.id][t].adjointInPlace();
            if (!file1.good()) {
              std::cout << "Problems while reading from " << infile1 << std::endl;
              exit(0);
            }
          }  // loop over time
          file1.close();
        } else if (file2.is_open()) {
          std::cout << "\treading VdaggerV from file:" << infile2 << std::endl;
          for (size_t t = 0; t < Lt; ++t) {
            // buffer for reading
            std::vector<Complex> eigen_vec(vdaggerv[op.id][t].size());
            file2.read(reinterpret_cast<char *>(&eigen_vec[0]),
                       vdaggerv[op.id][t].size() * sizeof(Complex));
            for (size_t ncol = 0; ncol < vdaggerv[op.id][t].cols(); ncol++) {
              for (size_t nrow = 0; nrow < vdaggerv[op.id][t].rows(); nrow++) {
                (vdaggerv[op.id][t])(nrow, ncol) =
                    eigen_vec.at(nrow * vdaggerv[op.id][t].cols() + ncol);
              }
            }
//            // @todo Check whether that must be adjoint before file 1 is closed
//            // (master branch)
//            vdaggerv[op.id][t].adjointInPlace();
            if (!file2.good()) {
              std::cout << "Problems while reading from " << infile2 << std::endl;
              exit(0);
            }
          }  // loop over time
          file2.close();
        } else {
          std::cout << "can't open " << infile1 << " NOR " << infile2 << std::endl;
          exit(0);
        }
      } else  // zero momentum
        for (size_t t = 0; t < Lt; ++t)
          vdaggerv[op.id][t] = Eigen::MatrixXcd::Identity(nb_ev, nb_ev);
    }
  }  // pragma omp parallel ends here

  t2 = clock() - t2;
  std::cout << std::setprecision(1) << "\t\t\tSUCCESS - " << std::fixed
            << ((float)t2) / CLOCKS_PER_SEC << " seconds" << std::endl;
  is_vdaggerv_set = true;
}

// ------------------------ INTERFACE ------------------------------------------
// -----------------------------------------------------------------------------

/*!
 *  @param filename The name to write to / read from the V^\dagger V operators
 *  @param rnd_vec  The random vector
 *  @param config   The configuration number to be read. Unused if
 *                  handling_vdaggerv is "build" or "write"
 *
 *  Behavior of this function depends on handling_vdaggerv flag.
 *  - "read" | "liuming" The operators are read in the corresponding format.
 *  - "build"            The operators are constructed from the eigenvectors
 *  - "write"            The operators are constructed and additionaly written
 *                       out.
 */
void OperatorFactory::create_operators(const std::string &filename,
                                       const RandomVector &rnd_vec,
                                       const int config) {
  is_vdaggerv_set = false;
  if (handling_vdaggerv == "write" || handling_vdaggerv == "build")
    build_vdaggerv(filename, config);
  else if (handling_vdaggerv == "read")
    read_vdaggerv(config);
  else if (handling_vdaggerv == "liuming")
    read_vdaggerv_liuming(config);
  else {
    std::cout << "\n\tThe flag handling_vdaggerv in input file is wrong!!\n\n"
              << std::endl;
    exit(0);
  }
}

/******************************************************************************/
/*!
 *  E.g. after building Quarkline Q2, vdaggerv is no longer needed and can be
 *  deleted to free up space
 *
 *  Resizes vdaggerv to 0
 */
void OperatorFactory::free_memory_vdaggerv() {
  std::for_each(vdaggerv.origin(),
                vdaggerv.origin() + vdaggerv.num_elements(),
                [](Eigen::MatrixXcd m) { m.resize(0, 0); });
}
