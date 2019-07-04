#include "OperatorsForMesons.hpp"
#include "StopWatch.hpp"

#include <boost/format.hpp>

#include <omp.h>
#include <iomanip>

static ssize_t map_char_to_dir(const char dir) {
  ssize_t integer_dir;
  if (dir == 'x')
    integer_dir = 0;
  if (dir == 'y')
    integer_dir = 1;
  if (dir == 'z')
    integer_dir = 2;
  return integer_dir;
}

/** Displace Matrix of eigenvectors by looping over displacement
 vectors

 v is the Eigensystem at a specific timeslice
 t is the timeslice to take into account
 disp is a vector of displacement pairs
 verbose is an integer controlling the debug level
 The eigenvectors of timeslice t are displaced using the vector of
 displacement pairs. The First entry of each pair specifies the direction
 of the right acting derivative, ">" meaning forward and "<" meaning backward

 @todo Swap order of loops?
 */
Eigen::MatrixXcd displace_eigenvectors(Eigen::MatrixXcd const &v,
                                       GaugeField const &gauge,
                                       ssize_t const t,
                                       DisplacementDirection const disp,
                                       ssize_t const verbose) {
  auto res = v;
  // iterate over displacement vector
  for (const auto &d : disp) {
    // Loop over all eigenvectors in v
    for (int ev = 0; ev < v.cols(); ++ev) {
      // Multiply eigenvector with according gauge matrix
      for (int spatial_ind = 0; spatial_ind < v.rows() / 3; ++spatial_ind) {
        res.col(ev).segment(3 * spatial_ind, 3) =
            (d.first == '>')
                ? gauge.forward_uv(res.col(ev), t, spatial_ind, map_char_to_dir(d.second))
                : gauge.backward_uv(
                      res.col(ev), t, spatial_ind, map_char_to_dir(d.second));
      }  // End spatial loop

    }  // end eigenvector loop
  }
  return res;
}

/** Creates a two-dimensional vector containing the momenta for the operators
 *
 *  @param[in] Lx, Ly, Lz      Lattice extent in spatial directions
 *  @param[in] vdaggerv_lookup Contains the momenta
 *  @param[in,out] momentum    Two dimensional array where the momenta are
 *                             stored
 */
void create_momenta(ssize_t const Lx,
                    ssize_t const Ly,
                    ssize_t const Lz,
                    std::vector<VdaggerVQuantumNumbers> const &vdaggerv_lookup,
                    array_cd_d2 &momentum) {
  static const std::complex<double> I(0.0, 1.0);

  /** To calculate Vdagger exp(i*p*x) V only the momenta corresponding to the
   *  quantum number id in op_VdaggerV will be used. The rest can be obtained
   *  by adjoining
   */
  for (const auto &op : vdaggerv_lookup) {
    // op_VdaggerV contains the index of one (redundancy) op_Corr which
    // allows to deduce the quantum numbers (momentum)
    const double ipx = op.momentum[0] * 2. * M_PI / static_cast<double>(Lx);
    const double ipy = op.momentum[1] * 2. * M_PI / static_cast<double>(Ly);
    const double ipz = op.momentum[2] * 2. * M_PI / static_cast<double>(Lz);
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

std::vector<Complex> roots_of_unity(int momentum, ssize_t size) {
  Complex const I(0.0, 1.0);

  double const two_pi_p_L = momentum * 2.0 * M_PI / static_cast<double>(size);

  std::vector<Complex> roots(size);
  for (ssize_t i = 0; i < size; ++i) {
    roots[i] = exp(-I * two_pi_p_L * static_cast<double>(i));
  }

  return roots;
}

void create_momenta_new(ssize_t const Lx,
                        ssize_t const Ly,
                        ssize_t const Lz,
                        std::vector<VdaggerVQuantumNumbers> const &vdaggerv_lookup,
                        array_cd_d2 &momentum) {
  /** To calculate Vdagger exp(i*p*x) V only the momenta corresponding to the
   *  quantum number id in op_VdaggerV will be used. The rest can be obtained
   *  by adjoining
   */
  for (auto const &op : vdaggerv_lookup) {
    auto const phases_x = roots_of_unity(op.momentum[0], Lx);
    auto const phases_y = roots_of_unity(op.momentum[1], Ly);
    auto const phases_z = roots_of_unity(op.momentum[2], Lz);

    // Calculate $\vec{p} \cdot \vec{x}$ for all $\vec{x}$ on the lattice.
    for (ssize_t x = 0; x < Lx; ++x) {
      auto const xH = x * Ly * Lz;
      auto const phase_x = phases_x[x];

      for (ssize_t y = 0; y < Ly; ++y) {
        auto const xHyH = xH + y * Lz;
        auto const phase_xy = phase_x * phases_y[y];

        for (ssize_t z = 0; z < Lz; ++z) {
          auto const phase_xyz = phase_xy * phases_z[z];
          momentum[op.id][xHyH + z] = phase_xyz;
        }
      }
    }
  }
}

namespace {

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
    for (ssize_t ncol = 0; ncol < Vt.cols(); ncol++) {
      for (ssize_t nrow = 0; nrow < Vt.rows(); nrow++) {
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

/**
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
OperatorFactory::OperatorFactory(const ssize_t Lt,
                                 const ssize_t Lx,
                                 const ssize_t Ly,
                                 const ssize_t Lz,
                                 const ssize_t nb_ev,
                                 const ssize_t dilE,
                                 const OperatorLookup &operator_lookuptable,
                                 const std::string &handling_vdaggerv,
                                 const std::string &path_vdaggerv,
                                 const std::string &path_config,
                                 const HypPars &hyp_parameters)
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
      path_config(path_config),
      hyp_parameters(hyp_parameters) {
  // resizing containers to their correct size
  vdaggerv.resize(boost::extents[operator_lookuptable.vdaggerv_lookup.size()][Lt]);

  // the momenta only need to be calculated for a subset of quantum numbers
  // (see VdaggerV::build_vdaggerv)
  momentum.resize(
      boost::extents[operator_lookuptable.vdaggerv_lookup.size()][Lx * Ly * Lz]);
  create_momenta(Lx, Ly, Lz, operator_lookuptable.vdaggerv_lookup, momentum);

  std::cout << "\tMeson operators initialised" << std::endl;
}

void OperatorFactory::build_vdaggerv(const std::string &filename, const int config) {
  const ssize_t dim_row = 3 * Lx * Ly * Lz;
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

  StopWatch swatch("Eigenvector and Gauge I/O");
  // resizing each matrix in vdaggerv
  // TODO: check if it is better to use for_each and resize instead of std::fill
  std::fill(vdaggerv.origin(),
            vdaggerv.origin() + vdaggerv.num_elements(),
            Eigen::MatrixXcd::Zero(nb_ev, nb_ev));

  // Read gauge field only if it is needed.
  /** @todo might be useful to parallelize */
  std::unique_ptr<GaugeField> gauge(nullptr);
  if (operator_lookuptable.need_gaugefield) {
    // If parameters for smearing are set, smear operator
    gauge.reset(new GaugeField(Lt, Lx, Ly, Lz, path_config, 0, Lt - 1, 4));
    gauge->read_gauge_field(config, 0, Lt - 1);
    if (hyp_parameters.iterations > 0) {
      const double alpha1 = hyp_parameters.alpha1;
      const double alpha2 = hyp_parameters.alpha2;
      const size_t iter = hyp_parameters.iterations;
      for (ssize_t t = 0; t < Lt; ++t) {
        gauge->smearing_hyp(t, alpha1, alpha2, iter);
      }
    }
  }

  omp_set_nested(1);
  Eigen::setNbThreads(4);

//  int nb_omp=omp_get_num_threads();
  omp_set_num_threads(8);



#pragma omp parallel
  {

    swatch.start();
    Eigen::VectorXcd mom = Eigen::VectorXcd::Zero(dim_row);

    EigenVector V_t(1, dim_row, nb_ev);  // each thread needs its own copy
    Eigen::MatrixXcd W_t;
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
          if (!op.displacement.empty()) {
            W_t = displace_eigenvectors(V_t[0], *gauge, t, op.displacement, 1);
          } else {
            // momentum vector contains exp(-i p x). Divisor 3 for colour index.
            // All three colours on same lattice site get the same momentum.
            for (ssize_t x = 0; x < dim_row; ++x) {
              mom(x) = momentum[op.id][x / 3];
            }
            W_t = mom.asDiagonal() * V_t[0];
          }

          vdaggerv[op.id][t] = V_t[0].adjoint() * W_t;

          // writing vdaggerv to disk
          if (handling_vdaggerv == "write") {
            std::string momentum_string = std::to_string(op.momentum[0]) +
                                          std::to_string(op.momentum[1]) +
                                          std::to_string(op.momentum[2]);
            std::string displacement_string = to_string(op.displacement);
            std::string outfile =
                (boost::format("operators.%04d.p_%s.d_%s.t_%03d") % config %
                 momentum_string % displacement_string % (int)t)
                    .str();
            write_vdaggerv(full_path, std::string(outfile), vdaggerv[op.id][t]);
          }

        } else {
          // zero momentum and no displacement
          vdaggerv[op.id][t] = Eigen::MatrixXcd::Identity(nb_ev, nb_ev);
        }
      }
    }  // loop over time
    swatch.stop();
  }  // pragma omp parallel ends here

  swatch.print();
  is_vdaggerv_set = true;
  exit(1);
}

void OperatorFactory::read_vdaggerv(const int config) {
  const int id_unity = operator_lookuptable.index_of_unity;

  // prepare full path for reading
  auto const full_path =
      (boost::format("/%s/cnfg%04d/operators.%04d") % path_vdaggerv % config % config)
          .str();

  // resizing each matrix in vdaggerv
  std::fill(vdaggerv.origin(),
            vdaggerv.origin() + vdaggerv.num_elements(),
            Eigen::MatrixXcd::Zero(nb_ev, nb_ev));
  StopWatch swatch("VdaggerV I/O");

#pragma omp parallel
  {
    swatch.start();
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
            for (ssize_t ncol = 0; ncol < vdaggerv[op.id][t].cols(); ncol++) {
              for (ssize_t nrow = 0; nrow < vdaggerv[op.id][t].rows(); nrow++) {
                (vdaggerv[op.id][t])(nrow, ncol) =
                    eigen_vec.at(ncol * vdaggerv[op.id][t].rows() + nrow);
              }
            }
            if (!file.good()) {
              std::ostringstream oss;
              oss << "Problems while reading from " << infile;
              std::runtime_error(oss.str());
            }
            file.close();
          } else {
            std::ostringstream oss;
            oss << "Can't open " << infile;
            std::runtime_error(oss.str());
          }
        } else  // zero momentum
          vdaggerv[op.id][t] = Eigen::MatrixXcd::Identity(nb_ev, nb_ev);
      }
    }  // loop over time
    swatch.stop();
  }  // pragma omp parallel ends here

  swatch.print();
  is_vdaggerv_set = true;
}

void OperatorFactory::read_vdaggerv_liuming(const int config) {
  const int id_unity = operator_lookuptable.index_of_unity;

  // prepare full path for reading
  auto const full_path = (boost::format("/%s/VdaggerV.") % path_vdaggerv).str();

  // resizing each matrix in vdaggerv
  std::fill(vdaggerv.origin(),
            vdaggerv.origin() + vdaggerv.num_elements(),
            Eigen::MatrixXcd::Zero(nb_ev, nb_ev));

  StopWatch swatch("Liuming VdaggerV I/O");
#pragma omp parallel
  {
    swatch.start();
    //  #pragma omp for schedule(dynamic)
    //    for(const auto& op : operator_lookuptable.vdaggerv_lookup){
#pragma omp for schedule(dynamic)
    for (ssize_t i = 0; i < ssize(operator_lookuptable.vdaggerv_lookup); ++i) {
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
          for (ssize_t t = 0; t < Lt; ++t) {
            // buffer for reading
            std::vector<Complex> eigen_vec(vdaggerv[op.id][t].size());
            file1.read(reinterpret_cast<char *>(&eigen_vec[0]),
                       vdaggerv[op.id][t].size() * sizeof(Complex));
            for (ssize_t ncol = 0; ncol < vdaggerv[op.id][t].cols(); ncol++) {
              for (ssize_t nrow = 0; nrow < vdaggerv[op.id][t].rows(); nrow++) {
                (vdaggerv[op.id][t])(nrow, ncol) =
                    eigen_vec.at(nrow * vdaggerv[op.id][t].cols() + ncol);
              }
            }
            vdaggerv[op.id][t].adjointInPlace();
            if (!file1.good()) {
              std::ostringstream oss;
              oss << "Problems while reading from " << infile1;
              std::runtime_error(oss.str());
            }
          }  // loop over time
          file1.close();
        } else if (file2.is_open()) {
          std::cout << "\treading VdaggerV from file:" << infile2 << std::endl;
          for (ssize_t t = 0; t < Lt; ++t) {
            // buffer for reading
            std::vector<Complex> eigen_vec(vdaggerv[op.id][t].size());
            file2.read(reinterpret_cast<char *>(&eigen_vec[0]),
                       vdaggerv[op.id][t].size() * sizeof(Complex));
            for (ssize_t ncol = 0; ncol < vdaggerv[op.id][t].cols(); ncol++) {
              for (ssize_t nrow = 0; nrow < vdaggerv[op.id][t].rows(); nrow++) {
                (vdaggerv[op.id][t])(nrow, ncol) =
                    eigen_vec.at(nrow * vdaggerv[op.id][t].cols() + ncol);
              }
            }
            //            // @todo Check whether that must be adjoint before file 1 is
            //            closed
            //            // (master branch)
            //            vdaggerv[op.id][t].adjointInPlace();
            if (!file2.good()) {
              std::ostringstream oss;
              oss << "Problems while reading from " << infile2;
              std::runtime_error(oss.str());
            }
          }  // loop over time
          file2.close();
        } else {
          std::ostringstream oss;
          oss << "can't open " << infile1 << " NOR " << infile2;
          std::runtime_error(oss.str());
        }
      } else  // zero momentum
        for (ssize_t t = 0; t < Lt; ++t)
          vdaggerv[op.id][t] = Eigen::MatrixXcd::Identity(nb_ev, nb_ev);
    }
    swatch.stop();
  }  // pragma omp parallel ends here

  swatch.print();
  is_vdaggerv_set = true;
}

/**
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
    throw std::runtime_error("The flag handling_vdaggerv in input file is wrong!");
  }
}

/**
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
