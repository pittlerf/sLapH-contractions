#include "OperatorsForMesons.h"

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
namespace {
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// creates a two-dimensional vector containing the momenta for the operators
// input: Lx, Ly, Lz      -> lattice extend in x, y, and z direction
//        vdaggerv_lookup -> contains the momenta
// output: momentum -> two dimensional array where the momenta are stored
void create_momenta(const size_t Lx, const size_t Ly, const size_t Lz, 
                    const std::vector<VdaggerVQuantumNumbers>& vdaggerv_lookup, 
                    array_cd_d2& momentum){

  static const std::complex<double> I(0.0, 1.0);

  // To calculate Vdagger exp(i*p*x) V only the momenta corresponding to the
  // quantum number id in op_VdaggerV will be used. The rest can be obtained
  // by adjoining
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
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
} // internal namespace ends here
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
LapH::OperatorsForMesons::OperatorsForMesons
                        (const size_t Lt, const size_t Lx, const size_t Ly, 
                         const size_t Lz, const size_t nb_ev, const size_t dilE,
                         const OperatorLookup& operator_lookuptable) : 
                               vdaggerv(), momentum(), 
                               operator_lookuptable(operator_lookuptable),
                               Lt(Lt), Lx(Lx), Ly(Ly), Lz(Lz), nb_ev(nb_ev), 
                               dilE(dilE) {

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
    for(auto& rvdv_level2 : rvdv_level1){
      rvdv_level2.resize(nb_rnd_combinations);
  }}

  rvdaggervr.resize(operator_lookuptable.rvdaggervr_lookuptable.size());
  counter = 0;
  for(auto& rvdvr_level1 : rvdaggervr){
    rvdvr_level1.resize(Lt);
    size_t nb_rnd_combinations = 
        operator_lookuptable.ricQ2_lookup[
                           operator_lookuptable.rvdaggervr_lookuptable[counter].
                    id_ricQ_lookup].rnd_vec_ids.size();
    counter++;
    for(auto& rvdvr_level2 : rvdvr_level1){
      rvdvr_level2.resize(nb_rnd_combinations);
  }}

  // the momenta only need to be calculated for a subset of quantum numbers
  // (see VdaggerV::build_vdaggerv)
  momentum.resize(boost::extents[
                       operator_lookuptable.vdaggerv_lookup.size()][Lx*Ly*Lz]);
  create_momenta(Lx, Ly, Lz, operator_lookuptable.vdaggerv_lookup, momentum);

}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void LapH::OperatorsForMesons::build_vdaggerv(const std::string& filename) {

  clock_t t2 = clock();
  const size_t dim_row = 3*Lx*Ly*Lz;
  const size_t id_unity = operator_lookuptable.index_of_unity;

  // resizing each matrix in vdaggerv
  // TODO: check if it is better to use for_each and resize instead of std::fill
  std::fill(vdaggerv.origin(), vdaggerv.origin() + vdaggerv.num_elements(), 
            Eigen::MatrixXcd::Zero(nb_ev, nb_ev));

#pragma omp parallel
{
  Eigen::VectorXcd mom = Eigen::VectorXcd::Zero(dim_row);
  LapH::EigenVector V_t(1, dim_row, nb_ev);// each thread needs its own copy
  #pragma omp for schedule(dynamic)
  for(size_t t = 0; t < Lt; ++t){

    // creating full filename for eigenvectors and reading them in
    if(!((operator_lookuptable.vdaggerv_lookup.size() == 1) &&
         (operator_lookuptable.vdaggerv_lookup[0].id == id_unity))){
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
      if(op.id != id_unity){
        // momentum vector contains exp(-i p x). Divisor 3 for colour index. 
        // All three colours on same lattice site get the same momentum.
        for(size_t x = 0; x < dim_row; ++x) {
          mom(x) = momentum[op.id][x/3];
        }
        vdaggerv[op.id][t] = V_t[0].adjoint() * mom.asDiagonal() * V_t[0];
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
void LapH::OperatorsForMesons::build_rvdaggerv(
                                            const LapH::RandomVector& rnd_vec) {

  // check if vdaggerv is already build
  if(not is_vdaggerv_set){
    std::cout << "\n\n\tCaution: vdaggerv is not set and rvdaggervr cannot be" 
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

  // rvdaggervr is calculated by multiplying vdaggerv with the same quantum
  // numbers with random vectors from right and left.
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
// -----------------------------------------------------------------------------
void LapH::OperatorsForMesons::create_operators(const std::string& filename, 
                                            const LapH::RandomVector& rnd_vec) {
  is_vdaggerv_set = false;
  build_vdaggerv(filename);
  build_rvdaggerv(rnd_vec);
  build_rvdaggervr(rnd_vec);
}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void LapH::OperatorsForMesons::free_memory_rvdaggerv(){
  for(auto& rvdv_level1 : rvdaggerv)
    for(auto& rvdv_level2 : rvdv_level1)
      for(auto& rvdv_level3 : rvdv_level2)
        rvdv_level3.resize(0, 0);
}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void LapH::OperatorsForMesons::free_memory_vdaggerv(){
  std::for_each(vdaggerv.origin(), vdaggerv.origin() + vdaggerv.num_elements(), 
                [](Eigen::MatrixXcd m){m.resize(0, 0);});
}














