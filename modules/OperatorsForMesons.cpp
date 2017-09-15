/*! @file OperatorsForMesons.cpp
 *  Class definition of LapH::OperatorsForMesons
 *
 *  @author Bastian Knippschild
 *  @author Markus Werner
 */

#include "OperatorsForMesons.h"

namespace {

/*! Creates a two-dimensional vector containing the momenta for the operators 
 *
 *  @param[in] Lx, Ly, Lz      Lattice extent in spatial directions
 *  @param[in] vdaggerv_lookup Contains the momenta
 *  @param[in,out] momentum    Two dimensional array where the momenta are 
 *                             stored
 */
void create_momenta(const size_t Lx, const size_t Ly, const size_t Lz, 
                    const std::vector<VdaggerVQuantumNumbers>& vdaggerv_lookup, 
                    array_cd_d2& momentum){
  static const std::complex<double> I(0.0, 1.0);

  /*! To calculate Vdagger exp(i*p*x) V only the momenta corresponding to the
   *  quantum number id in op_VdaggerV will be used. The rest can be obtained
   *  by adjoining
   */
  for(const auto& op : vdaggerv_lookup){
    // op_VdaggerV contains the index of one (redundancy) op_Corr which
    // allows to deduce the quantum numbers (momentum)
    const double ipx = op.momentum[0] * 2. * M_PI / (double) Lx; 
    const double ipy = op.momentum[1] * 2. * M_PI / (double) Ly;
    const double ipz = op.momentum[2] * 2. * M_PI / (double) Lz;
    // calculate \vec{p} \cdot \vec{x} for all \vec{x} on the lattice
    for(int x = 0; x < Lx; ++x){
      const int xH = x * Ly * Lz; // helper variable
      const double ipxH = ipx * x; // helper variable
      for(int y = 0; y < Ly; ++y){
        const int xHyH = xH + y * Lz; // helper variable
        const double ipxHipyH = ipxH + ipy * y; // helper variable
        for(int z = 0; z < Lz; ++z){
          // multiply \vec{p} \cdot \vec{x} with complex unit and exponentiate
          momentum[op.id][xHyH + z] = exp(-I * (ipxHipyH + ipz * z));
    }}}//loops over spatial vectors end here
  }//loop over redundant quantum numbers ends here
}

/******************************************************************************/
void write_vdaggerv(const std::string& pathname, const std::string& filename, 
                    const Eigen::MatrixXcd& Vt){

  // writing the data
  std::ofstream file((pathname+filename).c_str(), 
                     std::ofstream::binary | std::ofstream::trunc);

  if(file.is_open()){
    std::cout << "\twriting VdaggerV to file:" << pathname+filename 
              << std::endl;
    // buffer for writing
    vec eigen_vec(Vt.size());
    for (size_t ncol = 0; ncol < Vt.cols(); ncol++) {
      for(size_t nrow = 0; nrow < Vt.rows(); nrow++){
        eigen_vec.at(ncol*Vt.rows() + nrow) = (Vt)(nrow, ncol);
      }
    }
    file.write(reinterpret_cast<const char*>(&eigen_vec[0]), 
               Vt.size()*sizeof(cmplx));
    if(!file.good())
      std::cout << "Problems while write to " << (pathname+filename).c_str() 
                << std::endl;
    file.close();
  }
  else
    std::cout << "can't open " << (pathname+filename).c_str() 
              << std::endl;
}

} // internal namespace ends here

/******************************************************************************/
/*!
 * @param Lt, Lx, Ly, Lz  Temporal and spatial lattice extent
 * @param nb_ev           Number of eigenvectors
 * @param dilE            Number of diluted blocks in eigenvector space
 * @param operator_lookuptable ?
 * @param handling_vdaggerv
 * @param path_vdaggerv
 *
 * The initialization of the container attributes of LapH::OperatorsForMesons
 * is done in the member initializer list of the constructor. The allocation
 * of heap memory is delegated to boost::multi_array::resize
 */
LapH::OperatorsForMesons::OperatorsForMesons
                        (const size_t Lt, const size_t Lx, const size_t Ly, 
                         const size_t Lz, const size_t nb_ev, const size_t dilE,
                         const OperatorLookup& operator_lookuptable,
                         const std::string& handling_vdaggerv,
                         const std::string& path_vdaggerv,
                         const std::string& path_gaugefields) : 
                               vdaggerv(), momentum(), 
                               operator_lookuptable(operator_lookuptable),
                               Lt(Lt), Lx(Lx), Ly(Ly), Lz(Lz), nb_ev(nb_ev), 
                               dilE(dilE), handling_vdaggerv(handling_vdaggerv),
                               path_vdaggerv(path_vdaggerv),
                               path_gaugefields(path_gaugefields){
  // if path_gaugefields is not empty
  // resizing containers to their correct size
  vdaggerv.resize(boost::extents[
                             operator_lookuptable.vdaggerv_lookup.size()][Lt]);

  rvdaggerv.resize(operator_lookuptable.rvdaggerv_lookuptable.size());
  size_t counter = 0;
  for(auto& rvdv_level1 : rvdaggerv){
    rvdv_level1.resize(Lt);
    size_t nb_rnd_combinations = 
        operator_lookuptable.ricQ1_lookup[
                           operator_lookuptable.rvdaggerv_lookuptable[counter].
                    id_ricQ_lookup].rnd_vec_ids.size();
    counter++;
    for(auto& rvdv_level2 : rvdv_level1)
      rvdv_level2.resize(nb_rnd_combinations);
  }

  rvdaggervr.resize(operator_lookuptable.rvdaggervr_lookuptable.size());
  counter = 0;
  for(auto& rvdvr_level1 : rvdaggervr){
    rvdvr_level1.resize(Lt);
    size_t nb_rnd_combinations = 
        operator_lookuptable.ricQ2_lookup[
                           operator_lookuptable.rvdaggervr_lookuptable[counter].
                    id_ricQ_lookup].rnd_vec_ids.size();
    counter++;
    for(auto& rvdvr_level2 : rvdvr_level1)
      rvdvr_level2.resize(nb_rnd_combinations);
  }

  // the momenta only need to be calculated for a subset of quantum numbers
  // (see VdaggerV::build_vdaggerv)
  momentum.resize(boost::extents[
                       operator_lookuptable.vdaggerv_lookup.size()][Lx*Ly*Lz]);
  create_momenta(Lx, Ly, Lz, operator_lookuptable.vdaggerv_lookup, momentum);

  std::cout << "\tMeson operators initialised" << std::endl;
}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void LapH::OperatorsForMesons::build_vdaggerv(const std::string& filename,
                                              const int config) {

  clock_t t2 = clock();
  const size_t dim_row = 3*Lx*Ly*Lz;
  const int id_unity = operator_lookuptable.index_of_unity;
  const bool need_gaugefields = operator_lookuptable.need_gaugefields;

  // prepare full path for writing
  char dummy_path[200];
  sprintf(dummy_path, "/%s/cnfg%04d/", path_vdaggerv.c_str(), config);
  const std::string full_path(dummy_path);
  // check if directory exists
  if(handling_vdaggerv == "write" && access(full_path.c_str(), 0 ) != 0) {
    std::cout << "\tdirectory " << full_path.c_str() 
              << " does not exist and will be created";
    boost::filesystem::path dir(full_path.c_str());
    if(!boost::filesystem::create_directories(dir))
      std::cout << "\tSuccess" << std::endl;
    else
      std::cout << "\tFailure" << std::endl;
  }

  // resizing each matrix in vdaggerv
  // TODO: check if it is better to use for_each and resize instead of std::fill
  std::fill(vdaggerv.origin(), vdaggerv.origin() + vdaggerv.num_elements(), 
            Eigen::MatrixXcd::Zero(nb_ev, nb_ev));

  std::cout << "Input to gauge fields: Lt, Lx, Ly, Lz" << Lt << " " << Lx <<
    " " << Ly << " " << Lz << std::endl;
  GaugeField gauge = GaugeField(Lt, Lx, Ly, Lz, 
             path_gaugefields, size_t(0), 
             size_t(Lt-1), size_t(4));
  // TODO: enable read in of one timeslice
  if(need_gaugefields){
    gauge.read_gauge_field(config,size_t(0),size_t(Lt-1));
  } 

#pragma omp parallel shared(gauge)
{
  Eigen::VectorXcd mom = Eigen::VectorXcd::Zero(dim_row);
  LapH::EigenVector V_t(1, dim_row, nb_ev);// each thread needs its own copy
  // Loop over timeslices
  #pragma omp for schedule(dynamic)
  for(size_t t = 0; t < Lt; ++t){

    // creating full filename for eigenvectors and reading them in
    if(!(operator_lookuptable.vdaggerv_lookup.size() == 1 &&
       operator_lookuptable.vdaggerv_lookup[0].id == id_unity)){
      char inter_name[200];
      sprintf(inter_name, "%s%03d", filename.c_str(), (int) t);
      V_t.read_eigen_vector(inter_name, 0, 0); // reading eigenvectors
    }
    // VdaggerV is independent of the gamma structure and momenta connected by
    // sign flip are related by adjoining VdaggerV. Thus the expensive 
    // calculation must only be performed for a subset of quantum numbers given
    // in op_VdaggerV.
    for(const auto& op : operator_lookuptable.vdaggerv_lookup){
      // For zero momentum and displacement VdaggerV is the unit matrix, thus
      // the calculation is not performed
      std::cout << "In build_vdaggerv: id_unity is:" << id_unity << std::endl;
      std::cout << "In build_vdaggerv: displacement is:";
      for (auto& e : op.displacement) std::cout << e << " ";
      std::cout << std::endl;
      if(op.id != id_unity){
        // check whether displacement is wanted and determine the direction
        // (parallel to gamma)
        size_t dir = 0;
        //TODO: Order of displacements matters
        //TODO: At the moment only support for d > 0!!!!
        Eigen::MatrixXcd W_t = V_t[0];
        //for(auto& d : op_Corr[op.index].dis3){
        //op.displacement is a 3-vector of (x,y,z) displacements
        for(auto& d : op.displacement){ 
          if(d > 0){
            // displace d times in direction dir
            for(size_t nb_derv_one_dir = 0; nb_derv_one_dir < d; nb_derv_one_dir++){ 
              // LapH::EigenVector W_t(1,dim_row, nb_ev);
              //if(nb_derv_one_dir == 0)
              //  W_t = gauge.disp(V_t[0], t, dir, false);
              //else
              W_t = gauge.disp(W_t, t, dir, true);
            }
          }
          dir++;
        }
        vdaggerv[op.id][t] = V_t[0].adjoint() * W_t;
       // Eigen::MatrixXcd Trash = vdaggerv[op.id][t].adjoint();
       // vdaggerv[op.id][t] -= Trash; 

        // momentum vector contains exp(-i p x). Divisor 3 for colour index. 
        // All three colours on same lattice site get the same momentum.
        for(size_t x = 0; x < dim_row; ++x) {
          mom(x) = momentum[op.id][x/3];
        }
        //vdaggerv[op.id][t] = V_t[0].adjoint() * mom.asDiagonal() * W_t;
        // writing vdaggerv to disk
        if(handling_vdaggerv == "write"){
          char dummy2[200];
          sprintf(dummy2, "operators.%04d.p_", config);
          std::string dummy = std::string(dummy2) + 
                              std::to_string(op.momentum[0]) + 
                              std::to_string(op.momentum[1]) + 
                              std::to_string(op.momentum[2]);
          char outfile[200];
          sprintf(outfile, "%s_.t_%03d", dummy.c_str(), (int) t);
          write_vdaggerv(full_path, std::string(outfile), vdaggerv[op.id][t]);
        }
      }
      else // zero momentum
        vdaggerv[op.id][t] = Eigen::MatrixXcd::Identity(nb_ev, nb_ev);
    }
  } // loop over time
}// pragma omp parallel ends here

  t2 = clock() - t2;
  std::cout << std::setprecision(1) << "\t\t\tSUCCESS - " << std::fixed 
    << ((float) t2)/CLOCKS_PER_SEC << " seconds" << std::endl;
  is_vdaggerv_set = true;
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void LapH::OperatorsForMesons::read_vdaggerv(const int config){

  clock_t t2 = clock();
  const size_t dim_row = 3*Lx*Ly*Lz;
  const int id_unity = operator_lookuptable.index_of_unity;

  // prepare full path for reading
  char dummy_path[200];
  sprintf(dummy_path, "/%s/cnfg%04d/operators.%04d", path_vdaggerv.c_str(), 
                                                                config, config);
  std::string full_path(dummy_path);

  // resizing each matrix in vdaggerv
  std::fill(vdaggerv.origin(), vdaggerv.origin() + vdaggerv.num_elements(), 
            Eigen::MatrixXcd::Zero(nb_ev, nb_ev));

#pragma omp parallel
{
  #pragma omp for schedule(dynamic)
  for(size_t t = 0; t < Lt; ++t){
    for(const auto& op : operator_lookuptable.vdaggerv_lookup){
      // For zero momentum and displacement VdaggerV is the unit matrix, thus
      // the calculation is not performed
      if(op.id != id_unity){

        // creating full filename for vdaggerv and reading them in
        std::string dummy = full_path + ".p_" + 
                            std::to_string(op.momentum[0]) + 
                            std::to_string(op.momentum[1]) + 
                            std::to_string(op.momentum[2]);

        char infile[200];
        sprintf(infile, "%s_.t_%03d", dummy.c_str(), (int) t);

        // writing the data
        std::ifstream file(infile, std::ifstream::binary);
      
        if(file.is_open()){
          std::cout << "\treading VdaggerV from file:" << infile << std::endl;

          // buffer for reading
          vec eigen_vec(vdaggerv[op.id][t].size());
          file.read(reinterpret_cast<char*>(&eigen_vec[0]), 
                    vdaggerv[op.id][t].size()*sizeof(cmplx));
          for (size_t ncol = 0; ncol < vdaggerv[op.id][t].cols(); ncol++) {
            for(size_t nrow = 0; nrow < vdaggerv[op.id][t].rows(); nrow++){
               (vdaggerv[op.id][t])(nrow, ncol) = 
                            eigen_vec.at(ncol*vdaggerv[op.id][t].rows() + nrow);
            }
          }
          if(!file.good()){
            std::cout << "Problems while reading from " << infile << std::endl;
            exit(0);
          }
          file.close();
        }
        else{
          std::cout << "can't open " << infile << std::endl;
          exit(0);
        }
      }
      else // zero momentum
        vdaggerv[op.id][t] = Eigen::MatrixXcd::Identity(nb_ev, nb_ev);
    }
  } // loop over time
}// pragma omp parallel ends here

  t2 = clock() - t2;
  std::cout << std::setprecision(1) << "\t\t\tSUCCESS - " << std::fixed 
    << ((float) t2)/CLOCKS_PER_SEC << " seconds" << std::endl;
  is_vdaggerv_set = true;
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void LapH::OperatorsForMesons::read_vdaggerv_liuming(const int config){

  clock_t t2 = clock();
  const size_t dim_row = 3*Lx*Ly*Lz;
  const int id_unity = operator_lookuptable.index_of_unity;

  // prepare full path for reading
  char dummy_path[200];
  sprintf(dummy_path, "/%s/VdaggerV.", path_vdaggerv.c_str());
  std::string full_path(dummy_path);

  // resizing each matrix in vdaggerv
  std::fill(vdaggerv.origin(), vdaggerv.origin() + vdaggerv.num_elements(), 
            Eigen::MatrixXcd::Zero(nb_ev, nb_ev));

#pragma omp parallel
{
//  #pragma omp for schedule(dynamic)
//    for(const auto& op : operator_lookuptable.vdaggerv_lookup){
  #pragma omp for schedule(dynamic)
  for(size_t i = 0; i < operator_lookuptable.vdaggerv_lookup.size(); ++i){
    const auto op = (operator_lookuptable.vdaggerv_lookup[i]);
    // For zero momentum and displacement VdaggerV is the unit matrix, thus
    // the calculation is not performed
    if(op.id != id_unity){

      // creating full filename for vdaggerv and reading them in
      // both possibilities must be checked
      std::string dummy1 = full_path + "p" + 
                          std::to_string(-op.momentum[0]) + "p" +
                          std::to_string(-op.momentum[1]) + "p" +
                          std::to_string(-op.momentum[2]) + ".conf";
      char infile1[200];
      sprintf(infile1, "%s%04d", dummy1.c_str(), config);
      std::ifstream file1(infile1, std::ifstream::binary);
      // second possibility for a name
      std::string dummy2 = full_path + "p" + 
                          std::to_string(op.momentum[0]) + "p" +
                          std::to_string(op.momentum[1]) + "p" +
                          std::to_string(op.momentum[2]) + ".conf";
      char infile2[200];
      sprintf(infile2, "%s%04d", dummy2.c_str(), config);
      std::ifstream file2(infile2, std::ifstream::binary);
    
      if(file1.is_open()){
        std::cout << "\treading VdaggerV from file:" << infile1 << std::endl;
        for(size_t t = 0; t < Lt; ++t){
          // buffer for reading
          vec eigen_vec(vdaggerv[op.id][t].size());
          file1.read(reinterpret_cast<char*>(&eigen_vec[0]), 
                    vdaggerv[op.id][t].size()*sizeof(cmplx));
          for (size_t ncol = 0; ncol < vdaggerv[op.id][t].cols(); ncol++) {
            for(size_t nrow = 0; nrow < vdaggerv[op.id][t].rows(); nrow++){
               (vdaggerv[op.id][t])(nrow, ncol) = 
                            eigen_vec.at(nrow*vdaggerv[op.id][t].cols() + ncol);
            }
          }
          vdaggerv[op.id][t].adjointInPlace();
          if(!file1.good()){
            std::cout << "Problems while reading from " << infile1 << std::endl;
            exit(0);
          }
        } // loop over time
        file1.close();
      }
      else if(file2.is_open()){
        std::cout << "\treading VdaggerV from file:" << infile2 << std::endl;
        for(size_t t = 0; t < Lt; ++t){
          // buffer for reading
          vec eigen_vec(vdaggerv[op.id][t].size());
          file2.read(reinterpret_cast<char*>(&eigen_vec[0]), 
                    vdaggerv[op.id][t].size()*sizeof(cmplx));
          for (size_t ncol = 0; ncol < vdaggerv[op.id][t].cols(); ncol++) {
            for(size_t nrow = 0; nrow < vdaggerv[op.id][t].rows(); nrow++){
               (vdaggerv[op.id][t])(nrow, ncol) = 
                            eigen_vec.at(nrow*vdaggerv[op.id][t].cols() + ncol);
            }
          }
          // @todo Check whether that must be adjoint before file 1 is closed 
          // (master branch)
          vdaggerv[op.id][t].adjointInPlace();
          if(!file2.good()){
            std::cout << "Problems while reading from " << infile2 << std::endl;
            exit(0);
          }
        } // loop over time
        file2.close();
      }
      else{
        std::cout << "can't open " << infile1 << " NOR " << infile2 
                  << std::endl;
        exit(0);
      }
    }
    else // zero momentum
      for(size_t t = 0; t < Lt; ++t)
        vdaggerv[op.id][t] = Eigen::MatrixXcd::Identity(nb_ev, nb_ev);
    }
}// pragma omp parallel ends here


  t2 = clock() - t2;
  std::cout << std::setprecision(1) << "\t\t\tSUCCESS - " << std::fixed 
    << ((float) t2)/CLOCKS_PER_SEC << " seconds" << std::endl;
  is_vdaggerv_set = true;
  
}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void LapH::OperatorsForMesons::build_rvdaggerv(
                                            const LapH::RandomVector& rnd_vec) {

  // check if vdaggerv is already build
  if(not is_vdaggerv_set){
    std::cout << "\n\n\tCaution: vdaggerv is not set and rvdaggerv cannot be" 
              << " computed\n\n" << std::endl;
    exit(0);
  }

  clock_t t2 = clock();
  std::cout << "\tbuild rvdaggerv:";

  for(auto& rvdv_level1 : rvdaggerv)
    for(auto& rvdv_level2 : rvdv_level1)
      for(auto& rvdv_level3 : rvdv_level2)
        rvdv_level3 = Eigen::MatrixXcd::Zero(4*dilE, nb_ev);

#pragma omp parallel for schedule(dynamic)
  for(size_t t = 0; t < Lt; t++){

  // rvdaggerv is calculated by multiplying vdaggerv with the same quantum
  // numbers with random vectors from the left.
  for(const auto& op : operator_lookuptable.rvdaggerv_lookuptable){

    Eigen::MatrixXcd vdv;
    if(op.need_vdaggerv_daggering == false)
      vdv = vdaggerv[op.id_vdaggerv][t];
    else
      vdv = vdaggerv[op.id_vdaggerv][t].adjoint();

    size_t rid = 0;
    for(const auto& rnd_id : 
              operator_lookuptable.ricQ1_lookup[op.id_ricQ_lookup].rnd_vec_ids){

      for(size_t block = 0; block < 4; block++){
      for(size_t vec_i = 0; vec_i < nb_ev; ++vec_i) {
        size_t blk =  block + vec_i * 4 + 4 * nb_ev * t;
        
        rvdaggerv[op.id][t][rid].block(vec_i%dilE + dilE*block, 0, 1, nb_ev) += 
             vdv.row(vec_i) * std::conj(rnd_vec(rnd_id, blk));
      }}
      rid++;
    }
  }}// time and operator loops end here

  t2 = clock() - t2;
  std::cout << std::setprecision(1) << "\t\tSUCCESS - " << std::fixed 
    << ((float) t2)/CLOCKS_PER_SEC << " seconds" << std::endl;
}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void LapH::OperatorsForMesons::build_rvdaggervr(
                                            const LapH::RandomVector& rnd_vec) {

  // check of vdaggerv is already build
  if(not is_vdaggerv_set){
    std::cout << "\n\n\tCaution: vdaggerv is not set and rvdaggervr cannot be" 
              << " computed\n\n" << std::endl;
    exit(0);
  }

  clock_t t2 = clock();
  std::cout << "\tbuild rvdaggervr:";

  for(auto& rvdvr_level1 : rvdaggervr)
    for(auto& rvdvr_level2 : rvdvr_level1)
      for(auto& rvdvr_level3 : rvdvr_level2)
        rvdvr_level3 = Eigen::MatrixXcd::Zero(4*dilE, 4*dilE);

#pragma omp parallel for schedule(dynamic)
  for(size_t t = 0; t < Lt; t++){

  // rvdaggervr is calculated by multiplying vdaggerv with the same quantum
  // numbers with random vectors from right and left.
  for(const auto& op : operator_lookuptable.rvdaggervr_lookuptable){

    Eigen::MatrixXcd vdv;
    if(op.need_vdaggerv_daggering == false)
      vdv = vdaggerv[op.id_vdaggerv][t];
    else
      vdv = vdaggerv[op.id_vdaggerv][t].adjoint();

    size_t rid = 0;
    int check = -1;
    Eigen::MatrixXcd M; // Intermediate memory
    for(const auto& rnd_id : 
              operator_lookuptable.ricQ2_lookup[op.id_ricQ_lookup].rnd_vec_ids){

      if(check != rnd_id.first){ // this avoids recomputation
        M = Eigen::MatrixXcd::Zero(nb_ev, 4*dilE);
        for(size_t block = 0; block < 4; block++){
        for(size_t vec_i = 0; vec_i < nb_ev; vec_i++) {
          size_t blk =  block + (vec_i + nb_ev * t) * 4;
          M.block(0, vec_i%dilE + dilE*block, nb_ev, 1) += 
               vdv.col(vec_i) * rnd_vec(rnd_id.first, blk);
        }}
      }
      for(size_t block_x = 0; block_x < 4; block_x++){
      for(size_t block_y = 0; block_y < 4; block_y++){
      for(size_t vec_y = 0; vec_y < nb_ev; ++vec_y) {
        size_t blk =  block_y + (vec_y + nb_ev * t) * 4;
        rvdaggervr[op.id][t][rid].block(
                            dilE*block_y + vec_y%dilE, dilE*block_x, 1, dilE) +=
                M.block(vec_y, dilE*block_x, 1, dilE) * 
                std::conj(rnd_vec(rnd_id.second, blk));
      }}} 
      check = rnd_id.first;
      rid++;
    }
  }}// time and operator loops end here

  t2 = clock() - t2;
  std::cout << std::setprecision(1) << "\t\tSUCCESS - " << std::fixed 
    << ((float) t2)/CLOCKS_PER_SEC << " seconds" << std::endl;
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
void LapH::OperatorsForMesons::create_operators(const std::string& filename, 
                                            const LapH::RandomVector& rnd_vec,
                                            const int config) {
  is_vdaggerv_set = false;
  if(handling_vdaggerv == "write" || handling_vdaggerv == "build")
    build_vdaggerv(filename, config);
  else if(handling_vdaggerv == "read")
    read_vdaggerv(config);
  else if(handling_vdaggerv == "liuming")
    read_vdaggerv_liuming(config);
  else{
    std::cout << "\n\tThe flag handling_vdaggerv in input file is wrong!!\n\n"
              << std::endl;
    exit(0);
  }
  build_rvdaggerv(rnd_vec);
  build_rvdaggervr(rnd_vec);
}

/******************************************************************************/
/*!
 *  E.g. after building Quarkline Q1, vdaggerv is no longer needed and can be 
 *  deleted to free up space
 *
 *  Resizes rvdaggerv to 0
 */void LapH::OperatorsForMesons::free_memory_rvdaggerv(){
  for(auto& rvdv_level1 : rvdaggerv)
    for(auto& rvdv_level2 : rvdv_level1)
      for(auto& rvdv_level3 : rvdv_level2)
        rvdv_level3.resize(0, 0);
}

/******************************************************************************/
/*!
 *  E.g. after building Quarkline Q2, vdaggerv is no longer needed and can be 
 *  deleted to free up space
 *
 *  Resizes vdaggerv to 0
 */
void LapH::OperatorsForMesons::free_memory_vdaggerv(){
  std::for_each(vdaggerv.origin(), vdaggerv.origin() + vdaggerv.num_elements(), 
                [](Eigen::MatrixXcd m){m.resize(0, 0);});
}














