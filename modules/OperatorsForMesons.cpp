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
                            const size_t Lz, const size_t nb_ev,
                            const OperatorLookup& operator_lookuptable) : 
                                  vdaggerv(), momentum(), 
                                  operator_lookuptable(operator_lookuptable),
                                  Lt(Lt), Lx(Lx), Ly(Ly), Lz(Lz), nb_ev(nb_ev) {

//  const std::vector<quark> quarks = global_data->get_quarks();
//  const size_t nb_rnd = quarks[0].number_of_rnd_vec;
//  const vec_pd_rVdaggerVr op_rVdaggerVr = global_data->get_lookup_rVdaggerVr();
//  const size_t nb_rVdaggerVr = op_rVdaggerVr.size();
//  rvdaggervr.resize(boost::extents[nb_rVdaggerVr][Lt][nb_rnd][nb_rnd]);

  // only half of the array is stored to save memory. But be careful, it 
  // must be mapped correctly from outside by addressing the momentum
  // correctly and daggering
  vdaggerv.resize(boost::extents[
                             operator_lookuptable.vdaggerv_lookup.size()][Lt]);


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
  std::fill(vdaggerv.origin(), vdaggerv.origin() + vdaggerv.num_elements(), 
            Eigen::MatrixXcd::Zero(nb_ev, nb_ev));

#pragma omp parallel
{
  Eigen::VectorXcd mom = Eigen::VectorXcd::Zero(dim_row);
  LapH::EigenVector V_t(1, dim_row, nb_ev);// each thread needs its own copy
  #pragma omp for schedule(dynamic)
  for(size_t t = 0; t < Lt; ++t){

    // creating full filename for eigenvectors and reading them in
    char inter_name[200];
    sprintf(inter_name, "%s%03d", filename.c_str(), (int) t);
    V_t.read_eigen_vector(inter_name, 0, 0); // reading eigenvectors

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

}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
//void LapH::Operators::build_rvdaggervr(const int config_i,
//                     const std::vector<std::vector<LapH::RandomVector> >& rnd_vec) {
//
//  // check of vdaggerv is already build
//  if(not is_vdaggerv_set){
//    std::cout << "\n\n\tCaution: vdaggerv is not set and rvdaggervr cannot be" 
//              << " computed\n\n" << std::endl;
//    exit(0);
//  }
//
//  clock_t t2 = clock();
//  std::cout << "\tbuild rvdaggervr:";
//
//  const size_t Lt = global_data->get_Lt();
//  const size_t nb_ev = global_data->get_number_of_eigen_vec();
//  const std::vector<quark> quarks = global_data->get_quarks();
//  const size_t dilE = quarks[0].number_of_dilution_E;
//  const size_t nb_rnd = quarks[0].number_of_rnd_vec;
//
//  const vec_pd_rVdaggerVr op_rVdaggerVr = global_data->get_lookup_rVdaggerVr();
//  const vec_pdg_Corr op_Corr = global_data->get_lookup_corr();
//
//  std::fill(rvdaggervr.data(), rvdaggervr.data() + rvdaggervr.num_elements(), 
//            Eigen::MatrixXcd::Zero(dilE, 4*dilE));
//
//  // TODO: just a workaround
//  // can be changed to op by running over p = op/nb_dg, but dis currently
//  // not supported.
//
//  #pragma omp parallel for schedule(dynamic)
//  for(size_t t = 0; t < Lt; t++){
//
//  // rvdaggervr is calculated by multiplying the vdaggerv with the same quantum
//  // numbers with random vectors from right and left. rvdaggervr for momenta 
//  // with opposing sign are related by adjoining. Thus it suffices to calculate
//  // it for half of the indices for which the flag op.adjoint < 0.
//  for(const auto& op : op_rVdaggerVr){
//    if(op.adjoint == false){
//
//      size_t id_VdaggerV = op_Corr[op.index].id_vdv;
//
//      for(size_t rnd_i = 0; rnd_i < nb_rnd; ++rnd_i) {
//        Eigen::MatrixXcd M = Eigen::MatrixXcd::Zero(nb_ev, 4*dilE);
//        // dilution from left
//        for(size_t block= 0; block < 4; block++){
//        for(size_t vec_i = 0; vec_i < nb_ev; ++vec_i) {
//          size_t blk_i =  block + vec_i * 4 + 4 * nb_ev * t;
//          
//          M.block(0, vec_i%dilE + dilE*block, nb_ev, 1) += 
//               vdaggerv[id_VdaggerV][t].col(vec_i) * 
//               rnd_vec[1][rnd_i][blk_i];
//        }}// end of dilution
//        for(size_t rnd_j = 0; rnd_j < nb_rnd; ++rnd_j){
//        if(rnd_i != rnd_j){
//          // dilution from right
//          for(size_t block = 0; block < 4; block++){
//          for(size_t vec_j = 0; vec_j < nb_ev; ++vec_j) {
//            size_t blk_j =  block + vec_j * 4 + 4 * nb_ev * t;
//            rvdaggervr[op.id][t][rnd_j][rnd_i]
//                          .block(vec_j%dilE, dilE*block , 1, dilE) +=
//                M.block(vec_j, dilE*block, 1, dilE) * 
//                std::conj(rnd_vec[1][rnd_j][blk_j]);
//          }}// end of dilution
//        }}// rnd_j loop ends here
//      }// rnd_i loop ends here
//    }
//  }
//
//  // rvdaggervr for momenta with opposing sign are related by adjoining. Thus
//  // for half of the indices, the calculation reduces to adjoining the
//  // corresponding rvdaggervr and swapping the random vectors (as the order
//  // of multiplication is reversed). The index of corresponding quantum numbers
//  // is op.adjoint. It serves as flag for adjoining simultaneously, as it is
//  // positive if and only if it shall be adjoined.
//  // Need to loop twice as the corresponding rvdaggervr must all be calculated
//  // already.
//  for(const auto& op : op_rVdaggerVr){
//    if(op.adjoint == true){
//
//      for(size_t rnd_i = 0; rnd_i < nb_rnd; ++rnd_i) {
//      for(size_t rnd_j = 0; rnd_j < nb_rnd; ++rnd_j){
//      if(rnd_i != rnd_j){
//
//        // rvdaggervr is a blockdiagonal 4*dilE x 4*dilE matrix. To save memory,
//        // only the diagonal blocks are saved and it is written as a column 
//        // vector of blocks. To reproduce the correct behavior under adjoining, 
//        // the blocks have to be adjoined seperately.
//        // is .adjoint().transpose() faster?
//        for(size_t block = 0; block < 4; block++){
//          rvdaggervr[op.id][t][rnd_j][rnd_i]
//                              .block(0, block*dilE, dilE, dilE) =
//            (rvdaggervr[op.id_adjoint][t][rnd_i][rnd_j]
//                              .block(0, block*dilE, dilE, dilE)).adjoint();
//        }
//      }}}// loops over rnd vecs
//
//    }
//  }
//
//  }// time, momemtum and displacement loop ends here
//
//  t2 = clock() - t2;
//  std::cout << std::setprecision(1) << "\t\tSUCCESS - " << std::fixed 
//    << ((float) t2)/CLOCKS_PER_SEC << " seconds" << std::endl;
//
//}

// ------------------------ INTERFACE ------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void LapH::OperatorsForMesons::create_operators(const std::string& filename) {

  build_vdaggerv(filename);

}
