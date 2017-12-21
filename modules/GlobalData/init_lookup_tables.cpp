/*! @file
 *
 *  Functions translating lists build from the infile into lookup_tables.
 *
 *  @author Bastian Knippschild
 *  @author Markus Werner
 *  
 *  The lookup tables only contain unique quantum number combinations and 
 *  from the lists and latter are replaced by indexlists referring to those 
 *  lookup tables.
 *
 *  This automatically avoids recalculation of operators or even complete 
 *  correlators 
 *
 *  @todo Why are most functions static and not simply in an unnamed namespace?
 */

#include "global_data.h"
#include "global_data_utils.h"

namespace {

using Vector = QuantumNumbers::VectorData;

/******************************************************************************/

/*! 
 *  Check whether a momentum is within the set of desired momenta specified 
 *  in the correlator_list
 *
 *  @param[in]   p_tot Total 3-momentum at source or sink
 *  @param[out]  P     The specfied momenta of the reference frame
 *
 *  @returns     true if @f$|\mathtt{p_tot}|^2 \in \mathtt{P}@f$, false otherwise     
 */
bool desired_total_momentum(Vector const &p_tot,
                            std::vector<Vector> const &P){

  /*! If no total momentum is specified, there is no selection.
   *  @todo It is better to force the user to specify P. Catch that in 
   *        global_data_input_handling
   */
  if( P.empty() ){
    return true;
  }

  if( std::find(P.begin(), P.end(), p_tot) == P.end() ){
    return false;
  }
  else{
    return true;
  }
}

/*!
 *  For multi-meson operators check whether the sum of momenta at the same 
 *  time slice is above a given cutoff
 *
 *  @param[in]  p1  Momentum of the first meson
 *  @param[in]  p2  Momentum of the second meson
 *
 *  @returns        true if @f$|\mathtt{p1}|^2 + |\mathtt{p2}|^2 < \mathtt{cutoff} @f$, 
 *                  false otherwise
 *
 *  @warning Cutoff is hardcoded in this function!
 *  @warning No cutoff for P > 4 is implemented! 
 */
static bool momenta_below_cutoff(Vector const &p1,
                                 Vector const &p2) {

  std::map<int,int> cutoff;
  // default value is the default initializer of type int: cutoff[default] = 0;
  cutoff[0] = 4;
  cutoff[1] = 5;
  cutoff[2] = 6;
  cutoff[3] = 7;
  cutoff[4] = 4;

  const Vector p_tot = p1+p2;

  if(p_tot.squaredNorm() > 4){
    std::cout << "In momenta_below_cutoff(): WARNING! No cutoff for P > 4"
              << " implemented" << std::endl;
  }

  if( p1.squaredNorm() + p2.squaredNorm() > cutoff[p_tot.squaredNorm()] ) {
    return false;
  }
  else {
    return true;
  }
}

} // end of unnamed namespace

/******************************************************************************/
/*! Build an array with all the quantum numbers needed for a particular 
 *  correlation function respecting physical conservation laws
 *
 *  @param[in]  correlator      A single correlator specified in the infile 
 *                              and processed into the Correlators struct
 *  @param[in]  operator_list   List of all operators specified in the infile
 *                              and processed into Operators struct
 *  @param[out] quantum_numbers A list of all physical quantum numbers as 
 *                              specified in the QuantumNumbers struct that are
 *                              possible for @em correlator
 *
 *  @em correlator contains multiple operator_numbers. From combinatorics a 
 *  large number of combinations arise. In general only a subset of them are 
 *  physically allowed or necessary to calculate. 
 *  In this function momentum conservation is enforced and multiple cutoffs
 *  introduced.
 */
void build_quantum_numbers_from_correlator_list(const Correlators_2& correlator, 
                    const Operator_list& operator_list,
                    std::vector<std::vector<QuantumNumbers> >& quantum_numbers){


  std::vector<Operators> qn_op;
  for(const auto& op_number : correlator.operator_numbers){
    qn_op.emplace_back(operator_list[op_number]);
  }

  /*! Restriction to what shall actually be computed is done in if statements 
   *  for each diagram because it depends on the number of quarks.
   *
   *  @todo Think about a way to avoid these if conditions. 
   */
  if (correlator.type == "C1" || correlator.type == "C1T") {
    for(const auto& op0 : qn_op[0])
      quantum_numbers.emplace_back(std::vector<QuantumNumbers>({op0}));
  }
  else if (correlator.type == "C2+" || correlator.type == "C20" || 
           correlator.type == "C20V" || correlator.type == "Check") {

    // Build all combinations of operators and impose momentum conservation 
    // and cutoffs
    for(const auto& op0 : qn_op[0]){
      Vector p_so = op0.momentum;

      if( desired_total_momentum(p_so, correlator.tot_mom) ){

        for(const auto& op1 : qn_op[1]){ 
          Vector p_si = op1.momentum;
  
          // momentum at source and sink must always be the same for 2pt fcts.
          if(p_so == -p_si){
  
            std::vector<QuantumNumbers> single_vec_qn = {op0, op1};
            quantum_numbers.emplace_back(single_vec_qn);
          }
        } 
      }
    }
  }

  else if (correlator.type == "C3+" || correlator.type == "C30") {
    std::cout << "Constructing momentum combinations for C3"  << std::endl;
    
    std::map<int, int> counter; /*! initialized with zero */

    for(const auto& op0 : qn_op[0]){
    for(const auto& op2 : qn_op[2]){ 
      Vector p_so_1 = op0.momentum;
      Vector p_so_2 = op2.momentum;
      Vector p_so = p_so_1 + p_so_2;

      if( desired_total_momentum(p_so, correlator.tot_mom) &&
          momenta_below_cutoff(p_so_1, p_so_2) ){

        for(const auto& op1 : qn_op[1]){ 
          Vector p_si = op1.momentum;

          if( desired_total_momentum(p_si, correlator.tot_mom) ){

            if(p_so == -p_si){

              const int p_tot = p_si.squaredNorm();
              counter[p_tot]++;
  
              quantum_numbers.emplace_back(
                  std::vector<QuantumNumbers>{op0, op1, op2});
            }
          }
        }
      }
    }}

    int total_number_of_combinations = 0;
    for(const auto c : counter){
      std::cout << "\tCombinations for P = " << c.first << ": " << c.second 
                << std::endl;
      total_number_of_combinations += c.second;
    }
    std::cout << "\tTest finished - Combinations: " 
              << total_number_of_combinations << std::endl;
  }

  else if (correlator.type == "C4+D") {
    std::cout << "Constructing momentum combinations for C4+D"  << std::endl;
    
    std::map<int, int> counter; /*! initialized with zero */

    for(const auto& op0 : qn_op[0]){
    for(const auto& op2 : qn_op[2]){ 
      Vector p_so_1 = op0.momentum;
      Vector p_so_2 = op2.momentum;
      Vector p_so = p_so_1 + p_so_2;

      if( desired_total_momentum(p_so, correlator.tot_mom) &&
          momenta_below_cutoff(p_so_1, p_so_2) ){

        for(const auto& op1 : qn_op[1]){ 
        for(const auto& op3 : qn_op[3]){
          Vector p_si_1 = op1.momentum;
          Vector p_si_2 = op3.momentum;
          Vector p_si = p_si_1 + p_si_2;
    
          if( desired_total_momentum(p_si, correlator.tot_mom) &&
              momenta_below_cutoff(p_si_1, p_si_2) ){
  
            if(p_so == -p_si){

              const int p_tot = p_si.squaredNorm();
              counter[p_tot]++;
  
              // create combinations
              quantum_numbers.emplace_back(
                  std::vector<QuantumNumbers>{op0, op1, op2, op3});
            }
          }
        }}
      }
    }}

    int total_number_of_combinations = 0;
    for(const auto c : counter){
      std::cout << "\tCombinations for P = " << c.first << ": " << c.second 
                << std::endl;
      total_number_of_combinations += c.second;
    }
    std::cout << "\tTest finished - Combinations: " 
              << total_number_of_combinations << std::endl;
  }

  else if (correlator.type == "C4+B") {
    std::cout << "Constructing momentum combinations for C4+B"  << std::endl;
    
    std::map<int, int> counter; /*! initialized with zero */

    for(const auto& op0 : qn_op[0]){
    for(const auto& op3 : qn_op[3]){
      Vector p_so_1 = op0.momentum;
      Vector p_so_2 = op3.momentum;
      Vector p_so = p_so_1 + p_so_2;

      if( desired_total_momentum(p_so, correlator.tot_mom) &&
          momenta_below_cutoff(p_so_1, p_so_2) ){

        for(const auto& op1 : qn_op[1]){ 
        for(const auto& op2 : qn_op[2]){
          Vector p_si_1 = op1.momentum;
          Vector p_si_2 = op2.momentum;
          Vector p_si = p_si_1 + p_si_2;
    
          if( desired_total_momentum(p_si, correlator.tot_mom) &&
              momenta_below_cutoff(p_si_1, p_si_2) ){
  
            if(p_so == -p_si){

              const int p_tot = p_si.squaredNorm();
              counter[p_tot]++;
  
              // create combinations
              quantum_numbers.emplace_back(
                  std::vector<QuantumNumbers>{op0, op1, op2, op3});
            }
          }
        }}
      }
    }}

    int total_number_of_combinations = 0;
    for(const auto c : counter){
      std::cout << "\tCombinations for P = " << c.first << ": " << c.second 
                << std::endl;
      total_number_of_combinations += c.second;
    }
    std::cout << "\tTest finished - Combinations: " 
              << total_number_of_combinations << std::endl;

  }

  /*! @todo Check whether that is identical to C4+D */
  else if (correlator.type == "C4+C") {
    std::cout << "Constructing momentum combinations for C4+C"  << std::endl;
    
    std::map<int, int> counter; /*! initialized with zero */

    for(const auto& op0 : qn_op[0]){
    for(const auto& op2 : qn_op[2]){ 
      Vector p_so_1 = op0.momentum;
      Vector p_so_2 = op2.momentum;
      Vector p_so = p_so_1 + p_so_2;

      if( desired_total_momentum(p_so, correlator.tot_mom) &&
          momenta_below_cutoff(p_so_1, p_so_2) ){

        for(const auto& op1 : qn_op[1]){ 
        for(const auto& op3 : qn_op[3]){
          Vector p_si_1 = op1.momentum;
          Vector p_si_2 = op3.momentum;
          Vector p_si = p_si_1 + p_si_2;
    
          if( desired_total_momentum(p_si, correlator.tot_mom) &&
              momenta_below_cutoff(p_si_1, p_si_2) ){
  
            if(p_so == -p_si){

              const int p_tot = p_si.squaredNorm();
              counter[p_tot]++;
  
              // create combinations
              quantum_numbers.emplace_back(
                  std::vector<QuantumNumbers>{op0, op1, op2, op3});
            }
          }
        }}
      }
    }}

    int total_number_of_combinations = 0;
    for(const auto c : counter){
      std::cout << "\tCombinations for P = " << c.first << ": " << c.second 
                << std::endl;
      total_number_of_combinations += c.second;
    }
    std::cout << "\tTest finished - Combinations: " 
              << total_number_of_combinations << std::endl;
  }

  /*! @todo: For C40D, C40B, C40V, C40C, C4+V, C4+C still all combinations
   *         are built. 
   *         This must be changed later if GEVP should be used!!!!!!!!!!!!!!!
   */
  else if (correlator.type == "C40D" || correlator.type == "C40V" ||
           correlator.type == "C40B" || correlator.type == "C40C" ||
           correlator.type == "C4+V") {
    for(const auto& op0 : qn_op[0]){
    for(const auto& op1 : qn_op[1]){ 
    for(const auto& op2 : qn_op[2]){ 
    for(const auto& op3 : qn_op[3]){ // all combinations of operators
      quantum_numbers.emplace_back(
          std::vector<QuantumNumbers>{op0, op1, op2, op3});
    }}}}
  }
}

/******************************************************************************/
/*! Create the names for output files and hdf5 datasets.
 *
 *  @param[in]  corr_type {C1,C2+,C20,C20V,C3+,C30,C4+D,C4+V,C4+C,C4+B,C40D,
 *                         C40V,C40C,C40B} :
 *  @param[in]  cnfg :            Number of first gauge configuration
 *  @param[in]  outpath           Output path from the infile.
 *  @param[in]  overwrite {yes,no} : deprecated
 *  @param[in]  quark_types       Flavor of the quarks
 *  @param[in]  quantum_numbers   Physical quantum numbers 
 *  @param[out] hdf5_dataset_name Names for the datasets in one-to-one
 *                                correspondence to @em quantum_numbers
 *
 *  The output path is constructed by appending a "/" to @em outpath.
 *  The output filename is built from @em corr_type and @em cnfg.
 *  The dataset name is built from @em corr_type, a letter for each 
 *  @em quark_type, and the quantum numbers.
 *
 *  @todo Why don't we just build the complete path here already?
 */
static void build_correlator_names(const std::string& corr_type, int cnfg,  
               const std::string& outpath, const std::string& overwrite,
               const std::vector<std::string>& quark_types, 
               const std::vector<std::vector<QuantumNumbers> >& quantum_numbers,
               std::vector<std::string>& hdf5_dataset_name){
  
  for(const auto& qn_row : quantum_numbers){
    std::string filename =  corr_type + "_";
    for(const auto& qt : quark_types) // adding quark content
      filename += qt;
    size_t id = 0;
    for(const auto& qn : qn_row){ // adding quantum numbers
      filename += std::string("_p") + to_string(qn.momentum);
      filename += std::string(".d") + to_string(qn.displacement);
      filename += std::string(".g") + to_string(qn.gamma);
    }
    hdf5_dataset_name.emplace_back(filename);
  }
}

/******************************************************************************/
/*! Translate list of QuantumNumbers into lookuptable for VdaggerV 
 *
 *  @param[in]  quantum_numbers List of all quantum numbers operators are needed
 *                              for
 *  @param[out] vdaggerv_lookup Unique list of all VdaggerV operators needed.
 *                              Corresponds to @em quantum_numbers, but in 
 *                              contrast does not contain Dirac structure.
 *                              Part of GlobalData::operator_lookup
 *  @param[out] vdv_indices     Indexlist referring to @em vdaggerv_lookup
 *                              to replace @em quantum_numbers
 *                              The first index is the id of VdaggerV, the 
 *                              second tells us if VdaggerV must be daggered to 
 *                              get the desired quantum numbers.
 */
void build_VdaggerV_lookup(
             const std::vector<std::vector<QuantumNumbers> >& quantum_numbers,
             std::vector<VdaggerVQuantumNumbers>& vdaggerv_lookup,
             std::vector<std::vector<std::pair<size_t, bool> > >& vdv_indices) {

  for(const auto& qn_vec : quantum_numbers){
    std::vector<std::pair<size_t, bool> > vdv_indices_row;
    for(const auto& qn : qn_vec){
      bool dagger;
      // checking if the combination of quantum numbers already exists in 
      // vdaggerv. The position is stored in the iterator it.
      auto it = std::find_if(vdaggerv_lookup.begin(), vdaggerv_lookup.end(),
                           [&qn, &dagger](VdaggerVQuantumNumbers vdv_qn)
                           {
                             auto c1 = (Vector(vdv_qn.displacement.data()) == qn.displacement);
                             auto c2 = (Vector(vdv_qn.momentum.data()) == qn.momentum);
                             // also negative momentum is checked
                             const Vector pm = (-1)*qn.momentum;
                             auto c3 = (Vector(vdv_qn.momentum.data()) == pm);
                             // TODO: Think about the daggering!!
                             const Vector zero(0,0,0);
                             if (c1 and c2){
                               dagger = false;
                               return true;
                             } 
                             else if ((c1 and c3) and 
                                      (qn.displacement == zero)){
                               dagger = true;
                               return true;
                             }
                             else
                               return false;
                           });
      // If the quantum number combination already exists only the id is needed
      // otherwise a new element is created at the end of the lookuptable. 
      if(it != vdaggerv_lookup.end()) {
        vdv_indices_row.emplace_back((*it).id, dagger);
      }
      else {

        vdaggerv_lookup.emplace_back(VdaggerVQuantumNumbers(
            vdaggerv_lookup.size(), 
            {qn.momentum[0], qn.momentum[1], qn.momentum[2]}, 
            {qn.displacement[0], qn.displacement[1], qn.displacement[2]}));
        vdv_indices_row.emplace_back(vdaggerv_lookup.back().id, false);
      }
    }
    vdv_indices.emplace_back(vdv_indices_row);
  }
}

/******************************************************************************/
/*! @brief  Obtain index combinations of random vectors for charged correlator 
 *          i.e. correlator utilizing @f$ \gamma_5 @f$-trick
 *
 *  @param[in]  quarks      Quarks as read from the infile and processed into 
 *                          quark struct
 *  @param[in]  id_q1       Specifies which quark the first random index 
 *                          belongs to
 *  @param[in]  id_q2       Specifies which quark the second random index 
 *                          belongs to
 *  @param[in]  C1          Flag distinguishing whether the indexcombinations
 *                          are for C1 or not. 
 *  @param[out] rnd_vec_ids The random index combinations possible. If 
 *                          operator_list.rnd_vec_ids already contained these
 *                          combinations, just return the corresponding index
 *  @return                 index to operator_list.rnd_vec_ids pointing to the
 *                          random vector combinations built
 *
 *  For every quark propagator a statistical 1 in the form 
 *  @f$ ( P^{(b)} \rho) \cdot (P^{(b)} \rho)^\dagger @f$
 *  is introduced. 
 *  
 *  As explained in GlobalData, when factorizing the correlators this ones 
 *  are always split. To reconstruct the correct random index combinations, 
 *  this function constructs all allowed combinations of random indices for
 *  a quarkline with two random indices. To avoid bias, two different random
 *  vectors must always have different seed and thus different indices.
 *
 *  The random indices are uniquely identifying quark and random vector. Thus
 *  There are @f$ \sum_i q_i N_\text{rnd}(q_i) @f$ random indices.
 *  
 */
static size_t set_rnd_vec_charged(const std::vector<quark>& quarks, 
                          const size_t id_q1, const size_t id_q2, const bool C1,
                          std::vector<RandomIndexCombinationsQ2>& rnd_vec_ids) {

  // set start and end points of rnd numbers
  auto rndq1_start = 0;
  for(auto i = 0; i < id_q1; i++)
    rndq1_start += quarks[i].number_of_rnd_vec;
  auto rndq2_start = 0;
  for(auto i = 0; i < id_q2; i++)
    rndq2_start += quarks[i].number_of_rnd_vec;

  auto rndq1_end = rndq1_start + quarks[id_q1].number_of_rnd_vec;
  auto rndq2_end = rndq2_start + quarks[id_q2].number_of_rnd_vec;

  // check if there are enough random vectors
  if( (quarks[id_q1].number_of_rnd_vec < 1) || 
      (quarks[id_q2].number_of_rnd_vec < 1) ||
      (id_q1 == id_q2 && quarks[id_q1].number_of_rnd_vec < 2) ){
      std::cerr << "There are not enough random vectors for charged correlators"
                << std::endl;
      exit(-1);
  }

  // finally filling the array
  std::pair<size_t, size_t> offset = std::make_pair(rndq1_start, rndq2_start);
  std::vector<std::pair<size_t, size_t> > rnd_vec_comb;
  if(!C1){
    for(size_t i = rndq1_start; i < rndq1_end; ++i)
      for(size_t j = rndq2_start; j < rndq2_end; ++j)
        // To avoid bias, different random vectors must have different indices.
        if(i != j)
          rnd_vec_comb.emplace_back(i, j);
  }
  else {
    for(size_t i = rndq1_start; i < rndq1_end; ++i)
      // if C1 == True there is only one random vector and thus only same index
      // combinations are possible
      rnd_vec_comb.emplace_back(i, i);
  }

  // check if the random vector index combinations already exist
  const auto it = std::find_if(rnd_vec_ids.begin(), rnd_vec_ids.end(), 
                               [&](RandomIndexCombinationsQ2& vec){
                                 return (vec.rnd_vec_ids == rnd_vec_comb);
                               });
  if(it == rnd_vec_ids.end()) {
    rnd_vec_ids.emplace_back(RandomIndexCombinationsQ2(rnd_vec_ids.size(), id_q1, 
                                                  id_q2, offset, rnd_vec_comb));
    return rnd_vec_ids.back().id;
  }
  else
    return (*it).id;
}

/******************************************************************************/
/*! @brief  Obtain index combinations of random vectors for an uncharged 
 *          correlator i.e. correlator utilizing @f$ \gamma_5 @f$-trick
 *
 *  @param[in]  quarks      Quarks as read from the infile and processed into 
 *                          quark struct
 *  @param[in]  id_q1       Specifies which quark the first random index 
 *                          belongs to
 *  @param[out] rnd_vec_ids The random index combinations possible. If 
 *                          operator_list.rnd_vec_ids already contained these
 *                          combinations, just return the corresponding index
 *  @return                 index to operator_list.rnd_vec_ids pointing to the
 *                          random vector combinations built
 *
 *  For every quark propagator a statistical 1 in the form 
 *  @f$ ( P^{(b)} \rho) \cdot (P^{(b)} \rho)^\dagger @f$
 *  is introduced. 
 *  
 *  As explained in GlobalData, when factorizing the correlators this ones 
 *  are always split. To reconstruct the correct random index combinations, 
 *  this function constructs all allowed combinations of random indices for
 *  a quarkline with two random indices
 *
 */
static size_t set_rnd_vec_uncharged(const std::vector<quark>& quarks, 
                          const size_t id_q1, 
                          std::vector<RandomIndexCombinationsQ1>& rnd_vec_ids) {

  // First, check if the random vector indices already exists
  for(const auto& r_id : rnd_vec_ids)
    if(r_id.id_q1 == id_q1) 
      return r_id.id;

  // Set start and end points of rnd numbers
  auto rndq1_start = 0;
  for(auto i = 0; i < id_q1; i++)
    rndq1_start += quarks[i].number_of_rnd_vec;
  auto rndq1_end = rndq1_start + quarks[id_q1].number_of_rnd_vec;

  if(quarks[id_q1].number_of_rnd_vec < 2){
    std::cerr << "There are not enough random vectors for uncharged correlators"
              << std::endl;
    exit(-1);
  }

  // finally filling the array
  std::vector<size_t> rnd_vec_comb;
  for(size_t i = rndq1_start; i < rndq1_end; ++i)
    rnd_vec_comb.emplace_back(i);
  rnd_vec_ids.emplace_back(RandomIndexCombinationsQ1(rnd_vec_ids.size(), id_q1, 
                                                                 rnd_vec_comb));
  return rnd_vec_ids.back().id;
}

/******************************************************************************/
/*! For each correlator create lookuptable where to find the VdaggerV-operators
 *  and random index combinations necessary to build rVdaggerVr
 *  
 *  @param[in]  rnd_vec_id        Indices of random index combinations wanted 
 *                                for rVdaggerVr
 *  @param[in]  vdv_indices       Indices of QuantumNumbers  wanted for 
 *                                rVdaggerVr
 *  @param[out] rvdaggervr_lookup Unique list of combinations of @em rnd_vec_id
 *                                and @em vdv_indices needed
 *  @param[out] rvdvr_indices     List of indices referring to 
 *                                @em rvdaggervr_lookup. Entries correspond to
 *                                entries of @em vdv_indices
 */
static void build_rVdaggerVr_lookup(const std::vector<size_t>& rnd_vec_id, 
         const std::vector<std::vector<std::pair<size_t, bool> > >& vdv_indices,
         std::vector<VdaggerVRandomLookup>& rvdaggervr_lookup,
         std::vector<std::vector<size_t> >& rvdvr_indices) {

  // vdv_row contains one operator for each quark
  /*! Every correlator is treated seperately. */
  for(const auto& vdv_row : vdv_indices){
    std::vector<size_t> rvdvr_indices_row;

    for(size_t vdv_id = 0; vdv_id < vdv_row.size(); vdv_id++){

      /*! In case of multi-meson correlation functions multiple quarklines and 
       *  VdaggerV-operators must be combined into a single correlation 
       *  function. In these cases rnd_vec_id contains multiple entries. The 
       *  first two vdv_id's correspond to the first ricQ2lookup referenced by
       *  rnd_index, the next two to the second.
       */
      size_t rnd_index = 0; // construct to get correct random number indices
      if(vdv_id == 2 || vdv_id == 3)
        rnd_index = 1;

      const auto vdv = vdv_row[vdv_id];
     
      /*! Checks if the vdv and random vector index combination already are in 
       *  rvdaggervr_lookup
       */
      auto it = std::find_if(rvdaggervr_lookup.begin(), rvdaggervr_lookup.end(),
                             [&](VdaggerVRandomLookup vdv_qn)
                             {
                               auto c1 = (vdv_qn.id_ric_lookup == 
                                          rnd_vec_id[rnd_index]);
                               auto c2 = (vdv_qn.id_vdaggerv == vdv.first);
                               auto c3 = (vdv_qn.need_vdaggerv_daggering == 
                                          vdv.second); 
                               return (c1 && c2 && c3);
                             });

      /*! If yes, the corresponding index are saved to rvdvr_indices */
      if(it != rvdaggervr_lookup.end()) {
        rvdvr_indices_row.emplace_back((*it).id);
      }
      /*! If not, a new entry in rvdaggervr_lookup is created and then the 
       *  corresponding index is saved to rvdvr_indices 
       */
      else {
        rvdaggervr_lookup.emplace_back(VdaggerVRandomLookup(
                                       rvdaggervr_lookup.size(), vdv.first, 
                                       rnd_vec_id[rnd_index], vdv.second));
        rvdvr_indices_row.emplace_back(rvdaggervr_lookup.back().id);
      }
    }
    rvdvr_indices.emplace_back(rvdvr_indices_row);
  }
}

/******************************************************************************/
/*! For each correlator create lookuptable where to find the VdaggerV-operators
 *  and random index combinations necessary to build rVdaggerV
 *  
 *  @param[in]  rnd_vec_id        Indices of random index combinations wanted 
 *                                for rVdaggerV
 *  @param[in]  vdv_indices       Indices of QuantumNumbers  wanted for 
 *                                rVdaggerV
 *  @param[out] rvdaggervr_lookup Unique list of combinations of @em rnd_vec_id
 *                                and @em vdv_indices needed
 *  @param[out] rvdv_indices      List of indices referring to 
 *                                @em rvdaggerv_lookup. Entries correspond to
 *                                entries of @em vdv_indices
 *
 *  For detailed summary @see build_rVdaggerVr_lookup(). In contrast 
 *  @em rnd_vec_id refers to RandomIndexCombinationsQ1 and thus there is 
 *  exactly the same number of entries as for @vdv_indices
 *
 *  @todo Think about merging this and build_rVdaggerVr_lookup into one
 *        function. For n-point functions also the other one needs a vector
 *        of rnd_vec_id!
 */
static void build_rVdaggerV_lookup(const std::vector<size_t> rnd_vec_id, 
         const std::vector<std::vector<std::pair<size_t, bool> > >& vdv_indices,
         std::vector<VdaggerVRandomLookup>& rvdaggerv_lookup,
         std::vector<std::vector<size_t> >& rvdv_indices) {

  for(const auto& vdv_row : vdv_indices){

    std::vector<size_t> rvdv_indices_row;
    for(size_t i = 0; i < vdv_row.size(); i++){
      const auto& vdv = vdv_row.at(i);
      const auto& rnd = rnd_vec_id.at(i);

      auto it = std::find_if(rvdaggerv_lookup.begin(), rvdaggerv_lookup.end(),
                           [&](VdaggerVRandomLookup vdv_qn)
                           {
                             auto c1 = (vdv_qn.id_ric_lookup == rnd);
                             auto c2 = (vdv_qn.id_vdaggerv == vdv.first);
                             auto c3 = (vdv_qn.need_vdaggerv_daggering == 
                                        vdv.second); 
                             return (c1 && c2 && c3);
                           });

      if(it != rvdaggerv_lookup.end()) {
        rvdv_indices_row.emplace_back((*it).id);
      }
      else {
        rvdaggerv_lookup.emplace_back(VdaggerVRandomLookup(
                         rvdaggerv_lookup.size(), vdv.first, rnd, vdv.second));
        rvdv_indices_row.emplace_back(rvdaggerv_lookup.back().id);
      }
    }
    rvdv_indices.emplace_back(rvdv_indices_row);
  }
}

/******************************************************************************/
/******************************************************************************/
/*! Create lookuptable where to find the perambulators, VdaggerV-operators
 *  and random index combinations necessary to build Q2 quarkline
 *
 *  @param[in]  id_quark1       Identifies perambulator with gamma_5 trick
 *  @param[in]  id_quark2       Identifies perambulator without gamma_5 trick
 *  @param[in]  operator_id     Identifies needed operator from @em vdv_indices
 *  @param[in]  quantum_numbers Physical field operators for all correlators 
 *                              with Dirac structure factored out
 *  @param[in]  quarks          Quarks as specified in the infile munged into 
 *                              quark struct
 *  @param[in]  vdv_indices     Indices identifying VdaggerV operators
 *  @param[in]  ric_lookup      Indices identifying random vector index 
 *                              combinations
 *  @param[out] Q2V             Lookuptable containing unique combinations of 
 *                              peram-, vdv-, and ric-indices needed to built Q2
 *  @param[out] Q2_indices      List of indices refering to @em Q2V. Entries' 
 *                              outer vector corresponds to entries of 
 *                              @em quantum_numbers, the inner is specified by 
 *                              @em operator_id
 *
 *  The Quarkline with to quarks is given by
 *
 *    Q2 = @f$ \gamma_5 @f$ peram1@f$ ^\dagger \gamma_5 @f$ * vdaggerv * 
 *          gamma * peram2
 *
 *  This function creates a Lookup table of all unique index combination Q2 is 
 *  needed for and a list referring to the lookup table which gives relates
 *  the entries of the lookup table to each correlator specified in the infile
 *
 *  @todo Are id_quark1 and id_quark2 deprecated?
 */
static void build_Q2_lookup(const size_t id_quark1, const size_t id_quark2,
         const size_t operator_id,
         const std::vector<std::vector<QuantumNumbers> >& quantum_numbers, 
         const std::vector<quark>& quarks,
         const std::vector<std::vector<std::pair<size_t, bool> > >& vdv_indices,
         std::vector<RandomIndexCombinationsQ2>& ric_lookup,
         std::vector<QuarklineQ2Indices>& Q2V,
         std::vector<std::vector<size_t> >& Q2_indices){

  for(size_t row = 0; row < quantum_numbers.size(); row++){
    const auto qn = quantum_numbers[row][operator_id];
    const auto vdv = vdv_indices[row][operator_id];

    // Find out which index of ric_lookup contains the desired quantum 
    // numbers and use that for Q2V
    size_t rnd_index = set_rnd_vec_charged(quarks, id_quark1, id_quark2, 
                                           false, ric_lookup);
    // If Q2V already contains the particular row and physical content, just 
    // set the index to the existing QuarklineQ2Indicies, otherwise generate
    // it and set the index to the new one.
    QuarklineQ2Indices const candidate{Q2V.size(), vdv.first, id_quark1, id_quark2, 
          rnd_index, vdv.second, qn.gamma};
    auto it = std::find(Q2V.begin(), Q2V.end(), candidate);

    if(it != Q2V.end()) {
      Q2_indices[row][operator_id] = (*it).id;
    }
    else {
      Q2V.emplace_back(candidate);
      Q2_indices[row][operator_id] = Q2V.back().id;
    }
  }
}

/******************************************************************************/
/*! Create lookuptable where to find the perambulators, VdaggerV-operators
 *  and random index combinations necessary to build Q1 quarkline with field 
 *  operator @em operator_id
 *
 *  @param[in]  id_quark1       Specifies which quark the first random index 
 *                              belongs to    
 *  @param[in]  id_quark_connected Specifies which quark the second random index 
 *                              belongs to
 *  @param[in]  operator_id     Identifies needed operator from @em vdv_indices
 *  @param[in]  const bool C1   Flag distinguishing whether the 
 *                              indexcombinations are for C1 or not. 
 *  @param[in]  quantum_numbers Physical field operators for all correlators 
 *                              with Dirac structure factored out
 *  @param[in]  quarks          Quarks as specified in the infile munged into 
 *                              quark struct
 *  @param[in]  rvdv_indices    Indices identifying rVdaggerV operators
 *  @param[in]  ric_lookup      Indices identifying random vector index 
 *                              combinations
 *  @param[out] Q1              Lookuptable containing unique combinations of 
 *                              peram-, vdv-, and ric-indices needed to built Q1
 *  @param[out] Q1_indices      List of indices refering to @em Q1. Entries' 
 *                              outer vector corresponds to entries of 
 *                              @em quantum_numbers, the inner is specified by 
 *                              @em operator_id
 *
 *  The Quarkline with one quark is given by
 *
 *    Q1 =  rvdaggerv * gamma * peram
 *
 *  This function creates a Lookup table of all unique index combination Q2 is 
 *  needed for and a list referring to the lookup table which gives relates
 *  the entries of the lookup table to each correlator specified in the infile
 *
 *  @todo Is id_quark_used deprecated? Maybe it is only necessary with 
 *        different flavors supported
 *  @bug I think there is a bug in the order of quark ids. The first one 
 *        should be connected as it belongs to rVdV, while the second one 
 *        should belong to the quark as it is part of the perambulator. 
 *        Currently it is the other way round (MW, 27.3.17)
 */
static void build_Q1_lookup(const size_t id_quark_used, 
         const size_t id_quark_connected, const size_t operator_id, 
         const bool C1,
         const std::vector<std::vector<QuantumNumbers> >& quantum_numbers, 
         const std::vector<quark>& quarks,
         const std::vector<std::vector<size_t> >& rvdv_indices,
         std::vector<RandomIndexCombinationsQ2>& ric_lookup,
         std::vector<QuarklineQ1Indices>& Q1,
         std::vector<std::vector<size_t> >& Q1_indices){

  for(size_t row = 0; row < quantum_numbers.size(); row++){
    const auto qn = quantum_numbers[row][operator_id];
    const auto rvdv = rvdv_indices[row][operator_id];
    // Find out which index of ric_lookup contains the desired quantum 
    // numbers and use that for Q2V. Bug with order?
    const size_t rnd_index = set_rnd_vec_charged(quarks, id_quark_used,
        // must be id_quark_connected, id_quark used after bugfix in 
        // init_lookup_tables()
                                           id_quark_connected, C1, ric_lookup);

    // If Q2V already contains the particular row and physical content, just 
    // set the index to the existing QuarklineQ2Indicies, otherwise generate
    // it and set the index to the new one.
    auto it = std::find_if(Q1.begin(), Q1.end(),
                         [&](QuarklineQ1Indices q1)
                         {
                           auto c1 = (q1.id_peram == id_quark_used);
                           auto c2 = (q1.gamma == qn.gamma);
                           auto c3 = (q1.id_rvdaggerv == rvdv);
                           auto c4 = (q1.rnd_vec_ids == ric_lookup[rnd_index].rnd_vec_ids);
                           return c1 && c2 && c3 && c4;
                         });
    if(it != Q1.end()) {
        Q1_indices[row][operator_id] = (*it).id;
    }
    else {
      Q1.emplace_back(QuarklineQ1Indices{Q1.size(), 
                                         rvdv, 
                                         id_quark_used, 
                                         ric_lookup[rnd_index].offset, 
                                         ric_lookup[rnd_index].rnd_vec_ids, 
                                         qn.gamma});
      Q1_indices[row][operator_id] = Q1.back().id;
    }
  }
}

/******************************************************************************/
/******************************************************************************/
/*! Create lookuptable where to find the quarkline to build C1
 *  
 *  @param[in]  quantum_numbers   A list of all physical quantum numbers as 
 *                                specified in the QuantumNumbers struct that 
 *                                are possible for @em correlator
 *  @param[in]  hdf5_dataset_name Names for the datasets in one-to-one
 *                                correspondence to @em quantum_numbers
 *  @param[in]  Q1_indices        List of indices referring to lookup table
 *                                for Q1
 *  @param[out] corr_lookup       Lookup table containing lookup tables for 
 *                                all correlators this code can calculate.
 *                                This function sets corr_lookup.C1
 */

static void build_C1_lookup(
      const std::vector<std::vector<QuantumNumbers> >& quantum_numbers, 
      const std::vector<std::string>& hdf5_dataset_name,
      const std::vector<std::vector<size_t> >& Q1_indices, 
      CorrelatorLookup& corr_lookup){

  size_t row = 0;
  for(const auto& Q1 : Q1_indices){

    /*! C1 can be obtained from Q1 by just closing the trace. The index for C1 
     *  is a vector containing only one element for the Q1 quarkline      
     */
    std::vector<size_t> indices = {Q1[0]};

    auto it_C1 = std::find_if(corr_lookup.C1.begin(), corr_lookup.C1.end(),
                          [&](CorrInfo corr)
                          {
                            return (corr.hdf5_dataset_name == 
                                    hdf5_dataset_name[row]); 
                          });
    if(it_C1 == corr_lookup.C1.end()){
      corr_lookup.C1.emplace_back(CorrInfo(corr_lookup.C1.size(), 
                      hdf5_dataset_name[row], indices, 
                      quantum_numbers[row][0].gamma));
    }
    row++;
  }  
}

/******************************************************************************/
/******************************************************************************/
/*! Create lookuptable where to find the quarklines and rVdaggerVr-operators,
 *  to build C2c. Also sets corrC.
 *  
 *  @param[in]  quantum_numbers   A list of all physical quantum numbers as 
 *                                specified in the QuantumNumbers struct that 
 *                                are possible for @em correlator
 *  @param[in]  hdf5_dataset_name Names for the datasets in one-to-one
 *                                correspondence to @em quantum_numbers
 *  @param[in]  rvdvr_indices     List of indices referring to lookup table
 *                                for rvdaggerv. 
 *  @param[in]  Q2_indices        List of indices referring to lookup table
 *                                for Q2
 *  @param[out] corr_lookup       Lookup table containing lookup tables for 
 *                                all correlators this code can calculate.
 *                                This function sets corr_lookup.C2c and adds
 *                                all used unique combinations to 
 *                                corr_lookup.corrC
 *
 *  C2c is a part of C4cD and C4cV. To reuse C2c, they all contain indices
 *  of corrC which in turn contains the indices for rVdaggerVr and Q2.
 *
 */
static void build_C2c_lookup( 
      const std::vector<std::vector<QuantumNumbers> >& quantum_numbers, 
      const std::vector<std::string>& hdf5_dataset_name,
      const std::vector<std::vector<size_t> >& rvdvr_indices,
      const std::vector<std::vector<size_t> >& Q2_indices, 
      CorrelatorLookup& corr_lookup){

  /*! Loop over all distinct physical quantum numbers desired for C2c */
  for(size_t row = 0; row < hdf5_dataset_name.size(); row++){

    /*! Index for corrC is set to vector containing index for Quarkline first 
     *  and index for rVdaggerVr second
     *  Explicitly builds Op0:Op1
     */
    std::vector<size_t> indices = {Q2_indices[row][0], rvdvr_indices[row][1]};


    /*! If the quantum numbers are not yet specified in corr_lookup.C2c
     *  add them to the list 
     */
    auto it_C2c = std::find_if(corr_lookup.C2c.begin(), corr_lookup.C2c.end(),
                         [&](CorrInfo corr)
                         {
                           return (corr.hdf5_dataset_name == 
                                   hdf5_dataset_name[row]); 
                         });
    if(it_C2c == corr_lookup.C2c.end()){
      /*! If they refer to an existing corrC, just set the index, otherwise
       *  also add them to correalator_list.corrC. Note, that the Dirac 
       *  structure has been factored out from corrC.
       */
      auto it = std::find_if(corr_lookup.corrC.begin(), corr_lookup.corrC.end(),
                             [&](CorrInfo corr)
                             {
                               return (corr.lookup == indices); 
                             });
      if(it != corr_lookup.corrC.end()){
        corr_lookup.C2c.emplace_back(CorrInfo(corr_lookup.C2c.size(), 
                      hdf5_dataset_name[row], std::vector<size_t>({(*it).id}),
                      std::vector<int>({})));
      }
      else {
        corr_lookup.corrC.emplace_back(CorrInfo(corr_lookup.corrC.size(), 
                           "", indices, quantum_numbers[row][1].gamma));
        corr_lookup.C2c.emplace_back(CorrInfo(corr_lookup.C2c.size(), 
                           hdf5_dataset_name[row],
                           std::vector<size_t>({corr_lookup.corrC.back().id}),
                           std::vector<int>({})));
      }
    }
  }  
}

/******************************************************************************/
/*! Create lookuptable where to find the quarklines and rVdaggerVr-operators,
 *  to build C4cD. Also sets corrC.
 *  
 *  @param[in]  quantum_numbers   A list of all physical quantum numbers as 
 *                                specified in the QuantumNumbers struct that 
 *                                are possible for @em correlator
 *  @param[in]  hdf5_dataset_name Names for the datasets in one-to-one
 *                                correspondence to @em quantum_numbers
 *  @param[in]  rvdvr_indices     List of indices referring to lookup table
 *                                for rvdaggerv. 
 *  @param[in]  Q2_indices        List of indices referring to lookup table
 *                                for Q2
 *  @param[out] corr_lookup       Lookup table containing lookup tables for 
 *                                all correlators this code can calculate.
 *                                This function sets corr_lookup.C4cD and adds
 *                                all used unique combinations to 
 *                                corr_lookup.corrC
 *
 *  C4cD like C4cV contains C2c. To reuse C2c, they all contain indices
 *  of corrC which in turn contains the indices for rVdaggerVr and Q2.
 *
 */
static void build_C4cD_lookup( 
      const std::vector<std::vector<QuantumNumbers> >& quantum_numbers, 
      const std::vector<std::string>& hdf5_dataset_name,
      const std::vector<std::vector<size_t> >& rvdvr_indices,
      const std::vector<std::vector<size_t> >& Q2_indices, 
      CorrelatorLookup& corr_lookup){

  /*! Loop over all distinct physical quantum numbers desired for C4cD */
  for(size_t row = 0; row < hdf5_dataset_name.size(); row++){
    /*! Index for corrC is set to vector containing index for Quarkline first 
     *  and index for rVdaggerVr second. 
     *  Explicitly builds Op0:Op1 and Op2:Op3
     */
    std::vector<size_t> indices1 = {Q2_indices[row][0], rvdvr_indices[row][1]};
    std::vector<size_t> indices2 = {Q2_indices[row][2], rvdvr_indices[row][3]};

    auto it_C4cD = std::find_if(corr_lookup.C4cD.begin(), 
                                corr_lookup.C4cD.end(),
                          [&](CorrInfo corr)
                          {
                            return (corr.hdf5_dataset_name == 
                                    hdf5_dataset_name[row]); 
                          });

    if(it_C4cD == corr_lookup.C4cD.end()){
      size_t id1, id2;
      auto it1 = std::find_if(corr_lookup.corrC.begin(), 
                              corr_lookup.corrC.end(),
                              [&](CorrInfo corr)
                              {
                                return (corr.lookup == indices1); 
                              });
      if(it1 == corr_lookup.corrC.end()){
        corr_lookup.corrC.emplace_back(CorrInfo(corr_lookup.corrC.size(), 
                          "", indices1, quantum_numbers[row][1].gamma));
        id1 = corr_lookup.corrC.back().id;
      }
      else 
        id1 = (*it1).id;

      auto it2 = std::find_if(corr_lookup.corrC.begin(), 
                              corr_lookup.corrC.end(),
                              [&](CorrInfo corr)
                              {
                                return (corr.lookup == indices2); 
                              });
      if(it2 == corr_lookup.corrC.end()){
        corr_lookup.corrC.emplace_back(CorrInfo(corr_lookup.corrC.size(), 
                          "", indices2, quantum_numbers[row][3].gamma));
        id2 = corr_lookup.corrC.back().id;
      }
      else 
        id2 = (*it2).id;

      corr_lookup.C4cD.emplace_back(CorrInfo(corr_lookup.C4cD.size(), 
                      hdf5_dataset_name[row], std::vector<size_t>({id1, id2}), 
                      std::vector<int>({})));
    }
  }  
}

/******************************************************************************/
/*! Create lookuptable where to find the quarklines and rVdaggerVr-operators,
 *  to build C4cV Also sets corrC.
 *  
 *  @param[in]  quantum_numbers   A list of all physical quantum numbers as 
 *                                specified in the QuantumNumbers struct that 
 *                                are possible for @em correlator
 *  @param[in]  hdf5_dataset_name Names for the datasets in one-to-one
 *                                correspondence to @em quantum_numbers
 *  @param[in]  rvdvr_indices     List of indices referring to lookup table
 *                                for rvdaggerv. 
 *  @param[in]  Q2_indices        List of indices referring to lookup table
 *                                for Q2
 *  @param[out] corr_lookup       Lookup table containing lookup tables for 
 *                                all correlators this code can calculate.
 *                                This function sets corr_lookup.C4cV and adds
 *                                all used unique combinations to 
 *                                corr_lookup.corrC
 *
 *  C4cV like C4cC contains C2c. To reuse C2c, they all contain indices
 *  of corrC which in turn contains the indices for rVdaggerVr and Q2.
 *
 *  @todo The gamma attribut is set for corrC but not for C2c, C4cD, C4cV.
 *        Vice versa the corrlator names and hdf5_dataset_names are set for
 *        the latter but not corrC. Thats bs. (MW)
 */
static void build_C4cV_lookup( 
      const std::vector<std::vector<QuantumNumbers> >& quantum_numbers, 
      const std::vector<std::string>& hdf5_dataset_name,
      const std::vector<std::vector<size_t> >& rvdvr_indices,
      const std::vector<std::vector<size_t> >& Q2_indices, 
      CorrelatorLookup& corr_lookup){

  for(size_t row = 0; row < hdf5_dataset_name.size(); row++){
    /*! Index for corrC is set to vector containing index for Quarkline first 
     *  and index for rVdaggerVr second. 
     *  Explicitly builds Op0:Op1 and Op2:Op3
     */
    std::vector<size_t> indices1 = {Q2_indices[row][0], rvdvr_indices[row][1]};
    std::vector<size_t> indices2 = {Q2_indices[row][2], rvdvr_indices[row][3]};

    auto it_C4cV = std::find_if(corr_lookup.C4cV.begin(), 
                                corr_lookup.C4cV.end(),
                          [&](CorrInfo corr)
                          {
                            return (corr.hdf5_dataset_name == 
                                    hdf5_dataset_name[row]); 
                          });

    if(it_C4cV == corr_lookup.C4cV.end()){
      size_t id1, id2;
      auto it1 = std::find_if(corr_lookup.corrC.begin(), 
                              corr_lookup.corrC.end(),
                              [&](CorrInfo corr)
                              {
                                return (corr.lookup == indices1); 
                              });
      if(it1 == corr_lookup.corrC.end()){
        corr_lookup.corrC.emplace_back(CorrInfo(corr_lookup.corrC.size(), 
                         "", indices1, quantum_numbers[row][1].gamma));
        id1 = corr_lookup.corrC.back().id;
      }
      else 
        id1 = (*it1).id;

      auto it2 = std::find_if(corr_lookup.corrC.begin(), 
                              corr_lookup.corrC.end(),
                              [&](CorrInfo corr)
                              {
                                return (corr.lookup == indices2); 
                              });
      if(it2 == corr_lookup.corrC.end()){
        corr_lookup.corrC.emplace_back(CorrInfo(corr_lookup.corrC.size(), 
                         "", indices2, quantum_numbers[row][3].gamma));
        id2 = corr_lookup.corrC.back().id;
      }
      else 
        id2 = (*it2).id;

      corr_lookup.C4cV.emplace_back(CorrInfo(corr_lookup.C4cV.size(), 
                      hdf5_dataset_name[row], std::vector<size_t>({id1, id2}), 
                      std::vector<int>({})));
    }
  }  
}

/******************************************************************************/
/*! Create lookuptable where to find the quarklines and rVdaggerVr-operators,
 *  to build C4cC 
 *  
 *  @param[in]  quantum_numbers   A list of all physical quantum numbers as 
 *                                specified in the QuantumNumbers struct that 
 *                                are possible for @em correlator
 *  @param[in]  hdf5_dataset_name Names for the datasets in one-to-one
 *                                correspondence to @em quantum_numbers
 *  @param[in]  rvdvr_indices     List of indices referring to lookup table
 *                                for rvdaggerv. 
 *  @param[in]  Q2_indices        List of indices referring to lookup table
 *                                for Q2
 *  @param[out] corr_lookup       Lookup table containing lookup tables for 
 *                                all correlators this code can calculate.
 *                                This function sets corr_lookup.C4cC
 */
static void build_C4cC_lookup( 
      const std::vector<std::vector<QuantumNumbers> >& quantum_numbers, 
      const std::vector<std::string>& hdf5_dataset_name,
      const std::vector<std::vector<size_t> >& rvdvr_indices,
      const std::vector<std::vector<size_t> >& Q2_indices, 
      CorrelatorLookup& corr_lookup){

  for(size_t row = 0; row < hdf5_dataset_name.size(); row++){
    /*! Index for C4cC is set to vector containing index for Quarkline first
     *  and third and index for rVdaggerVr second and forth.
     *  Explicitly builds Op0:Op1:Op2:Op3
     */
    std::vector<size_t> indices = {Q2_indices[row][0], rvdvr_indices[row][1], 
                                   Q2_indices[row][2], rvdvr_indices[row][3]};
    auto it_C4cC = std::find_if(corr_lookup.C4cC.begin(), 
                                corr_lookup.C4cC.end(),
                          [&](CorrInfo corr)
                          {
                            return (corr.hdf5_dataset_name == 
                                    hdf5_dataset_name[row]); 
                          });

    if(it_C4cC == corr_lookup.C4cC.end()){
      std::vector<int> gammas = {{quantum_numbers[row][1].gamma[0], 
                                  quantum_numbers[row][3].gamma[0]}};
      corr_lookup.C4cC.emplace_back(CorrInfo(corr_lookup.C4cC.size(), 
                      hdf5_dataset_name[row], indices, gammas));
    }
  }  
}

/******************************************************************************/
/*! Create lookuptable where to find the quarklines and rVdaggerVr-operators,
 *  to build C4cB
 *  
 *  @param[in]  quantum_numbers   A list of all physical quantum numbers as 
 *                                specified in the QuantumNumbers struct that 
 *                                are possible for @em correlator
 *  @param[in]  hdf5_dataset_name Names for the datasets in one-to-one
 *                                correspondence to @em quantum_numbers
 *  @param[in]  rvdvr_indices     List of indices referring to lookup table
 *                                for rvdaggerv. 
 *  @param[in]  Q2_indices        List of indices referring to lookup table
 *                                for Q2
 *  @param[out] corr_lookup       Lookup table containing lookup tables for 
 *                                all correlators this code can calculate.
 *                                This function sets corr_lookup.C4cB
 */
static void build_C4cB_lookup( 
      const std::vector<std::vector<QuantumNumbers> >& quantum_numbers, 
      const std::vector<std::string>& hdf5_dataset_name,
      const std::vector<std::vector<size_t> >& rvdvr_indices,
      const std::vector<std::vector<size_t> >& Q2_indices, 
      CorrelatorLookup& corr_lookup){

  for(size_t row = 0; row < hdf5_dataset_name.size(); row++){
    /*! Index for C4cB is set to vector containing index for Quarkline first
     *  and third and index for rVdaggerVr second and forth.
     *  Explicitly builds Op0:Op1:Op2:Op3
     */
    std::vector<size_t> indices = {Q2_indices[row][0], rvdvr_indices[row][1], 
                                   Q2_indices[row][2], rvdvr_indices[row][3]};

    auto it_C4cB = std::find_if(corr_lookup.C4cB.begin(), 
                                corr_lookup.C4cB.end(),
                          [&](CorrInfo corr)
                          {
                            return (corr.hdf5_dataset_name == 
                                    hdf5_dataset_name[row]); 
                          });

    if(it_C4cB == corr_lookup.C4cB.end()){
      std::vector<int> gammas = {{quantum_numbers[row][1].gamma[0], 
                                  quantum_numbers[row][3].gamma[0]}};
      corr_lookup.C4cB.emplace_back(CorrInfo(corr_lookup.C4cB.size(), 
                      hdf5_dataset_name[row], indices, gammas));
    }
  } 
}

/******************************************************************************/
/*! Create lookuptable where to find the quarklines and rVdaggerVr-operators,
 *  to build C3c
 *  
 *  @param[in]  quantum_numbers   A list of all physical quantum numbers as 
 *                                specified in the QuantumNumbers struct that 
 *                                are possible for @em correlator
 *  @param[in]  hdf5_dataset_name Names for the datasets in one-to-one
 *                                correspondence to @em quantum_numbers
 *  @param[in]  rvdvr_indices     List of indices referring to lookup table
 *                                for rvdaggerv. 
 *  @param[in]  Q1_indices        List of indices referring to lookup table
 *                                for Q1
 *  @param[in]  Q2_indices        List of indices referring to lookup table
 *                                for Q2
 *  @param[out] corr_lookup       Lookup table containing lookup tables for 
 *                                all correlators this code can calculate.
 *                                This function sets corr_lookup.C3c
 *
 *  @bug I am fairly certain that the quarks are mixed up. It is 
 *        also wrong in init_lookup_tables()
 */
static void build_C3c_lookup( 
      const std::vector<std::vector<QuantumNumbers> >& quantum_numbers, 
      const std::vector<std::string>& hdf5_dataset_name,
      const std::vector<std::vector<size_t> >& rvdvr_indices,
      const std::vector<std::vector<size_t> >& Q1_indices, 
      const std::vector<std::vector<size_t> >& Q2_indices, 
      CorrelatorLookup& corr_lookup){

  for(size_t row = 0; row < hdf5_dataset_name.size(); row++){
    /*! Index for C4cB is set to vector containing index for Quarkline Q2 first,
     *  the index for Quarkline Q1 second and the index for rVdaggerVr third
     *  This is correct at least.
     *  Explicitly builds Op0:Op1:Op2
     */
    std::vector<size_t> indices = {Q2_indices[row][0], Q1_indices[row][1],  
                                   rvdvr_indices[row][2]};
    auto it_C3c = std::find_if(corr_lookup.C3c.begin(), 
                               corr_lookup.C3c.end(),
                          [&](CorrInfo corr)
                          {
                            return (corr.hdf5_dataset_name == 
                                    hdf5_dataset_name[row]); 
                          });
    if(it_C3c == corr_lookup.C3c.end()){
      std::vector<int> gammas = {{quantum_numbers[row][2].gamma[0]}};
      corr_lookup.C3c.emplace_back(CorrInfo(corr_lookup.C3c.size(), 
                      hdf5_dataset_name[row], indices, gammas));
    }
  } 
}

/******************************************************************************/
/******************************************************************************/
/*! Create lookuptable where to find the quarklines to build C20. Also sets
 *  C1.
 *  
 *  @param[in]  quantum_numbers   A list of all physical quantum numbers as 
 *                                specified in the QuantumNumbers struct that 
 *                                are possible for @em correlator
 *  @param[in]  hdf5_dataset_name Names for the datasets in one-to-one
 *                                correspondence to @em quantum_numbers
 *  @param[in]  Q1_indices        List of indices referring to lookup table
 *                                for Q1
 *  @param[out] corr_lookup       Lookup table containing lookup tables for 
 *                                all correlators this code can calculate.
 *                                This function sets corr_lookup.C20
 *
 *  C20V contains C1. To reuse C1, they all contain indices
 *  of C1 which in turn contains the indices for Q1.

 *  @bug I am fairly certain that the quarks are mixed up. It is 
 *        also wrong in init_lookup_tables() (MW 27.3.17)
 */

static void build_C20V_lookup(
      const std::vector<std::string>& hdf5_dataset_name,
      const std::vector<std::vector<size_t> >& Q1_indices, 
      CorrelatorLookup& corr_lookup){

  size_t row = 0;
  for(const auto& Q1 : Q1_indices){
    std::vector<size_t> indices1 = {Q1[0]};
    std::vector<size_t> indices2 = {Q1[1]};
    auto it_C20V = std::find_if(corr_lookup.C20V.begin(), 
                                corr_lookup.C20V.end(),
                          [&](CorrInfo corr)
                          {
                            return (corr.hdf5_dataset_name == 
                                    hdf5_dataset_name[row]); 
                          });
    if(it_C20V == corr_lookup.C20V.end()){
      size_t id1, id2;
      auto it1 = std::find_if(corr_lookup.C1.begin(), 
                              corr_lookup.C1.end(),
                              [&](CorrInfo corr)
                              {
                                return (corr.lookup == indices1); 
                              });
      if((it1 == corr_lookup.C1.end())){
        corr_lookup.C1.emplace_back(CorrInfo(corr_lookup.C1.size(), 
                                  "", indices1, std::vector<int>({})));
        id1 = corr_lookup.C1.back().id;
      }
      else
        id1 = (*it1).id;
      auto it2 = std::find_if(corr_lookup.C1.begin(), 
                              corr_lookup.C1.end(),
                              [&](CorrInfo corr)
                              {
                                return (corr.lookup == indices2); 
                              });
      if((it2 == corr_lookup.C1.end())){
        corr_lookup.C1.emplace_back(CorrInfo(corr_lookup.C1.size(), 
                                  "", indices2, std::vector<int>({})));
        id2 = corr_lookup.C1.back().id;
      }
      else
        id2 = (*it2).id;

      corr_lookup.C20V.emplace_back(CorrInfo(corr_lookup.C20V.size(), 
                      hdf5_dataset_name[row], std::vector<size_t>({id1, id2}), 
                      std::vector<int>({})));
    }
    row++;
  }  
}

/******************************************************************************/
/******************************************************************************/
/*! Create lookuptable where to find the quarklines to build C20. Also sets
 *  corr0.
 *  
 *  @param[in]  quantum_numbers   A list of all physical quantum numbers as 
 *                                specified in the QuantumNumbers struct that 
 *                                are possible for @em correlator
 *  @param[in]  hdf5_dataset_name Names for the datasets in one-to-one
 *                                correspondence to @em quantum_numbers
 *  @param[in]  Q1_indices        List of indices referring to lookup table
 *                                for Q1
 *  @param[out] corr_lookup       Lookup table containing lookup tables for 
 *                                all correlators this code can calculate.
 *                                This function sets corr_lookup.C20
 *
 *  C20 is a part of C40D and C40V. To reuse C20, they all contain indices
 *  of corr0 which in turn contains the indices for Q1.

 *  @bug I am fairly certain that the quarks are mixed up. It is 
 *        also wrong in init_lookup_tables() (MW 27.3.17)
 */
static void build_C20_lookup(
      const std::vector<std::string>& hdf5_dataset_name,
      const std::vector<std::vector<size_t> >& Q1_indices, 
      CorrelatorLookup& corr_lookup){

  size_t row = 0;
  for(const auto& Q1 : Q1_indices){
    std::vector<size_t> indices = {Q1[0], Q1[1]};
    auto it_C20 = std::find_if(corr_lookup.C20.begin(), corr_lookup.C20.end(),
                          [&](CorrInfo corr)
                          {
                            return (corr.hdf5_dataset_name == 
                                    hdf5_dataset_name[row]); 
                          });
    if(it_C20 == corr_lookup.C20.end()){
      auto it = std::find_if(corr_lookup.corr0.begin(), corr_lookup.corr0.end(),
                             [&](CorrInfo corr)
                             {
                               return (corr.lookup == indices); 
                             });
      if(it != corr_lookup.corr0.end()){
        corr_lookup.C20.emplace_back(CorrInfo(corr_lookup.C20.size(), 
                      hdf5_dataset_name[row], std::vector<size_t>({(*it).id}), 
                      std::vector<int>({})));
      }
      else {
        corr_lookup.corr0.emplace_back(CorrInfo(corr_lookup.corr0.size(), 
                                   "", indices, std::vector<int>({})));
        corr_lookup.C20.emplace_back(CorrInfo(corr_lookup.C20.size(), 
                           hdf5_dataset_name[row],  
                           std::vector<size_t>({corr_lookup.corr0.back().id}),
                           std::vector<int>({})));
      }
    }
    row++;
  }  
}

/******************************************************************************/
/*! Create lookuptable where to find the quarklines to build C30.
 *  
 *  @param[in]  quantum_numbers   A list of all physical quantum numbers as 
 *                                specified in the QuantumNumbers struct that 
 *                                are possible for @em correlator
 *  @param[in]  hdf5_dataset_name Names for the datasets in one-to-one
 *                                correspondence to @em quantum_numbers
 *  @param[in]  Q1_indices        List of indices referring to lookup table
 *                                for Q1
 *  @param[out] corr_lookup       Lookup table containing lookup tables for 
 *                                all correlators this code can calculate.
 *                                This function sets corr_lookup.C30
 *
 *  @bug I am fairly certain that the quarks are mixed up. It is 
 *        also wrong in init_lookup_tables() (MW 27.3.17)
 */
static void build_C30_lookup(
      const std::vector<std::string>& hdf5_dataset_name,
      const std::vector<std::vector<size_t> >& Q1_indices, 
      CorrelatorLookup& corr_lookup){

  size_t row = 0;
  // Indices for correlator_lookup.C30 are one-to-one the indices for Q1_indices
  for(const auto& Q1 : Q1_indices){
    auto it = std::find_if(corr_lookup.C30.begin(), corr_lookup.C30.end(),
                          [&](CorrInfo corr)
                          {
                            return (corr.hdf5_dataset_name == 
                                    hdf5_dataset_name[row]); 
                          });
    if(it == corr_lookup.C30.end()){
      corr_lookup.C30.emplace_back(CorrInfo(corr_lookup.C30.size(), 
                      hdf5_dataset_name[row], Q1, std::vector<int>({})));
    }
    row++;
  }  
}

/******************************************************************************/
/*! Create lookuptable where to find the quarklines to build C40D. Also sets
 *  corr0.
 *  
 *  @param[in]  quantum_numbers   A list of all physical quantum numbers as 
 *                                specified in the QuantumNumbers struct that 
 *                                are possible for @em correlator
 *  @param[in]  hdf5_dataset_name Names for the datasets in one-to-one
 *                                correspondence to @em quantum_numbers
 *  @param[in]  Q1_indices        List of indices referring to lookup table
 *                                for Q1
 *  @param[out] corr_lookup       Lookup table containing lookup tables for 
 *                                all correlators this code can calculate.
 *                                This function sets corr_lookup.C40D
 *
 *  C40D like C40V contains C20. To reuse C20, they all contain indices
 *  of corr0 which in turn contains the indices for Q1.

 *  @bug I am fairly certain that the quarks are mixed up. It is 
 *        also wrong in init_lookup_tables() (MW 27.3.17)
 */
static void build_C40D_lookup(
      const std::vector<std::string>& hdf5_dataset_name,
      const std::vector<std::vector<size_t> >& Q1_indices, 
      CorrelatorLookup& corr_lookup){

  size_t row = 0;
  for(const auto& Q1 : Q1_indices){
    std::vector<size_t> indices1 = {Q1[0], Q1[1]};
    std::vector<size_t> indices2 = {Q1[2], Q1[3]};
    auto it_C40D = std::find_if(corr_lookup.C40D.begin(), 
                                corr_lookup.C40D.end(),
                          [&](CorrInfo corr)
                          {
                            return (corr.hdf5_dataset_name == 
                                    hdf5_dataset_name[row]); 
                          });
    if(it_C40D == corr_lookup.C40D.end()){
      size_t id1, id2;
      auto it1 = std::find_if(corr_lookup.corr0.begin(), 
                              corr_lookup.corr0.end(),
                              [&](CorrInfo corr)
                              {
                                return (corr.lookup == indices1); 
                              });
      if((it1 == corr_lookup.corr0.end())){
        corr_lookup.corr0.emplace_back(CorrInfo(corr_lookup.corr0.size(), 
                                  "", indices1, std::vector<int>({})));
        id1 = corr_lookup.corr0.back().id;
      }
      else
        id1 = (*it1).id;
      auto it2 = std::find_if(corr_lookup.corr0.begin(), 
                              corr_lookup.corr0.end(),
                              [&](CorrInfo corr)
                              {
                                return (corr.lookup == indices2); 
                              });
      if((it2 == corr_lookup.corr0.end())){
        corr_lookup.corr0.emplace_back(CorrInfo(corr_lookup.corr0.size(), 
                                  "", indices2, std::vector<int>({})));
        id2 = corr_lookup.corr0.back().id;
      }
      else
        id2 = (*it2).id;

      corr_lookup.C40D.emplace_back(CorrInfo(corr_lookup.C40D.size(), 
                      hdf5_dataset_name[row], std::vector<size_t>({id1, id2}), 
                      std::vector<int>({})));
    }
    row++;
  }  
}

/******************************************************************************/
/*! Create lookuptable where to find the quarklines to build C40V. Also sets
 *  corr0.
 *  
 *  @param[in]  quantum_numbers   A list of all physical quantum numbers as 
 *                                specified in the QuantumNumbers struct that 
 *                                are possible for @em correlator
 *  @param[in]  hdf5_dataset_name Names for the datasets in one-to-one
 *                                correspondence to @em quantum_numbers
 *  @param[in]  Q1_indices        List of indices referring to lookup table
 *                                for Q1
 *  @param[out] corr_lookup       Lookup table containing lookup tables for 
 *                                all correlators this code can calculate.
 *                                This function sets corr_lookup.C40V
 *
 *  C40V like C40D contains C20. To reuse C20, they all contain indices
 *  of corr0 which in turn contains the indices for Q1.

 *  @bug I am fairly certain that the quarks are mixed up. It is 
 *        also wrong in init_lookup_tables() (MW 27.3.17)
 */
static void build_C40V_lookup(
      const std::vector<std::string>& hdf5_dataset_name,
      const std::vector<std::vector<size_t> >& Q1_indices, 
      CorrelatorLookup& corr_lookup){

  size_t row = 0;
  for(const auto& Q1 : Q1_indices){
    std::vector<size_t> indices1 = {Q1[0], Q1[1]};
    std::vector<size_t> indices2 = {Q1[2], Q1[3]};
    auto it_C40V = std::find_if(corr_lookup.C40V.begin(), 
                                corr_lookup.C40V.end(),
                          [&](CorrInfo corr)
                          {
                            return (corr.hdf5_dataset_name == 
                                    hdf5_dataset_name[row]); 
                          });
    if(it_C40V == corr_lookup.C40V.end()){
      size_t id1, id2;
      auto it1 = std::find_if(corr_lookup.corr0.begin(), 
                              corr_lookup.corr0.end(),
                              [&](CorrInfo corr)
                              {
                                return (corr.lookup == indices1); 
                              });
      if((it1 == corr_lookup.corr0.end())){
        corr_lookup.corr0.emplace_back(CorrInfo(corr_lookup.corr0.size(), 
                                 "", indices1, std::vector<int>({})));
        id1 = corr_lookup.corr0.back().id;
      }
      else
        id1 = (*it1).id;
      auto it2 = std::find_if(corr_lookup.corr0.begin(), 
                              corr_lookup.corr0.end(),
                              [&](CorrInfo corr)
                              {
                                return (corr.lookup == indices2); 
                              });
      if((it2 == corr_lookup.corr0.end())){
        corr_lookup.corr0.emplace_back(CorrInfo(corr_lookup.corr0.size(), 
                                  "", indices2, std::vector<int>({})));
        id2 = corr_lookup.corr0.back().id;
      }
      else
        id2 = (*it2).id;

      corr_lookup.C40V.emplace_back(CorrInfo(corr_lookup.C40V.size(), 
                      hdf5_dataset_name[row], std::vector<size_t>({id1, id2}), 
                      std::vector<int>({})));
    }
    row++;
  }  
}

/******************************************************************************/
/*! Create lookuptable where to find the quarklines to build C40C.
 *  
 *  @param[in]  quantum_numbers   A list of all physical quantum numbers as 
 *                                specified in the QuantumNumbers struct that 
 *                                are possible for @em correlator
 *  @param[in]  hdf5_dataset_name Names for the datasets in one-to-one
 *                                correspondence to @em quantum_numbers
 *  @param[in]  Q1_indices        List of indices referring to lookup table
 *                                for Q1
 *  @param[out] corr_lookup       Lookup table containing lookup tables for 
 *                                all correlators this code can calculate.
 *                                This function sets corr_lookup.C40C
 *
 *  @bug I am fairly certain that the quarks are mixed up. It is 
 *        also wrong in init_lookup_tables() (MW 27.3.17)
 */
static void build_C40C_lookup(
      const std::vector<std::string>& hdf5_dataset_name,
      const std::vector<std::vector<size_t> >& Q1_indices, 
      CorrelatorLookup& corr_lookup){

  size_t row = 0;
  // Indices for correlator_lookup.C40C are one-to-one the indices for Q1_indices
  for(const auto& Q1 : Q1_indices){
    auto it = std::find_if(corr_lookup.C40C.begin(), corr_lookup.C40C.end(),
                          [&](CorrInfo corr)
                          {
                            return (corr.hdf5_dataset_name == 
                                    hdf5_dataset_name[row]); 
                          });
    if(it == corr_lookup.C40C.end()){
      corr_lookup.C40C.emplace_back(CorrInfo(corr_lookup.C40C.size(), 
                      hdf5_dataset_name[row], Q1, std::vector<int>({})));
    }
    row++;
  }  
}

/******************************************************************************/
/*! Create lookuptable where to find the quarklines to build C40B.
 *  
 *  @param[in]  quantum_numbers   A list of all physical quantum numbers as 
 *                                specified in the QuantumNumbers struct that 
 *                                are possible for @em correlator
 *  @param[in]  hdf5_dataset_name Names for the datasets in one-to-one
 *                                correspondence to @em quantum_numbers
 *  @param[in]  Q1_indices        List of indices referring to lookup table
 *                                for Q1
 *  @param[out] corr_lookup       Lookup table containing lookup tables for 
 *                                all correlators this code can calculate.
 *                                This function sets corr_lookup.C40B
 *
 *  @bug I am fairly certain that the quarks are mixed up. It is 
 *        also wrong in init_lookup_tables() (MW 27.3.17)
 */
static void build_C40B_lookup(
      const std::vector<std::string>& hdf5_dataset_name,
      const std::vector<std::vector<size_t> >& Q1_indices, 
      CorrelatorLookup& corr_lookup){

  size_t row = 0;
  // Indices for correlator_lookup.C40B are one-to-one the indices for Q1_indices
  for(const auto& Q1 : Q1_indices){
    auto it = std::find_if(corr_lookup.C40B.begin(), corr_lookup.C40B.end(),
                          [&](CorrInfo corr)
                          {
                            return (corr.hdf5_dataset_name == 
                                    hdf5_dataset_name[row]); 
                          });
    if(it == corr_lookup.C40B.end()){
      corr_lookup.C40B.emplace_back(CorrInfo(corr_lookup.C40B.size(), 
                      hdf5_dataset_name[row], Q1, std::vector<int>({})));
    }
    row++;
  }  
}



/******************************************************************************/
/******************************************************************************/
/*! 
 *  from GlobalData::correlators_list, GlobalData::operator::list and 
 *  GlobalData::quarks
 *
 *  @bug In build_Q1_lookup the order of quarks given is consistently switched.
 */
void GlobalData::init_lookup_tables() {

  for (const auto& correlator : correlator_list){

    /*! 1. Build an array (quantum_numbers) with all the quantum numbers needed 
     *      for this particular correlation function.
     */
    std::vector<std::vector<QuantumNumbers> > quantum_numbers;
    build_quantum_numbers_from_correlator_list(correlator, operator_list, 
                                               quantum_numbers);

    // Build the correlator and dataset names for hdf5 output files
    std::vector<std::string> quark_types; 
    for(const auto& id : correlator.quark_numbers)
      quark_types.emplace_back(quarks[id].type);
    std::vector<std::string> hdf5_dataset_name;
    build_correlator_names(correlator.type, start_config, path_output, 
                     overwrite, quark_types, quantum_numbers,
                     hdf5_dataset_name);

    /*! 2. Build the lookuptable for VdaggerV and return an array of indices
     *      corresponding to @em quantum_numbers computed in step 1. In 
     *      @em vdv_indices the first entry is the id of vdv, the second tells 
     *      us if vdv must be daggered to get the correct quantum numbers.
     */
    std::vector<std::vector<std::pair<size_t, bool> > > vdv_indices;
    build_VdaggerV_lookup(quantum_numbers, operator_lookuptable.vdaggerv_lookup,
                                                                   vdv_indices);
    if (correlator.type == "C1") {
      std::vector<size_t> rnd_vec_id;
      rnd_vec_id.emplace_back( set_rnd_vec_uncharged(quarks, 
                                 correlator.quark_numbers[0], 
                                 operator_lookuptable.ricQ1_lookup) );
      std::vector<std::vector<size_t> > rvdv_indices;
      build_rVdaggerV_lookup(rnd_vec_id, vdv_indices,
                             operator_lookuptable.rvdaggerv_lookuptable,
                             rvdv_indices);
      std::vector<std::vector<size_t> > Q1_indices(rvdv_indices.size(),
                            std::vector<size_t>(rvdv_indices[0].size()));
      build_Q1_lookup(correlator.quark_numbers[0], correlator.quark_numbers[0],
                      0, true, quantum_numbers, quarks, rvdv_indices, 
                      operator_lookuptable.ricQ2_lookup,
                      quarkline_lookuptable.Q1, Q1_indices);
      build_C1_lookup(quantum_numbers, hdf5_dataset_name,
                      Q1_indices, correlator_lookuptable);
    }
    else if (correlator.type == "C2+" || correlator.type == "Check") {
      /*! 3. Build the lookuptable for rVdaggerVr and return an array of indices
       *      corresponding to the 'quantum_numbers' computed in step 1.
       *
       *  @todo There should be a warning if more than 2 entries are to be 
       *        written into rnd_vec_id, as this will break 
       *        build_rVdaggerVr_lookup
       */
      std::vector<size_t> rnd_vec_id;
      rnd_vec_id.emplace_back(set_rnd_vec_charged(quarks, 
                                            correlator.quark_numbers[1], 
                                            correlator.quark_numbers[0], false,
                                            operator_lookuptable.ricQ2_lookup));
      std::vector<std::vector<size_t> > rvdvr_indices;
      build_rVdaggerVr_lookup(rnd_vec_id, vdv_indices,
                              operator_lookuptable.rvdaggervr_lookuptable,
                              rvdvr_indices);
      /*! 4. Build the lookuptable for Q2 and return an array of indices
       *      corresponding to the 'quantum_numbers' computed in step 1.
       */
      std::vector<std::vector<size_t> > Q2_indices(rvdvr_indices.size(),
                                  std::vector<size_t>(rvdvr_indices[0].size()));
      build_Q2_lookup(correlator.quark_numbers[0], correlator.quark_numbers[1],
                      0, quantum_numbers, quarks, vdv_indices, 
                      operator_lookuptable.ricQ2_lookup,
                      quarkline_lookuptable.Q2V, Q2_indices);
      /*! 5. Build the lookuptable for the correlation functions */
      build_C2c_lookup(quantum_numbers, hdf5_dataset_name, 
                       rvdvr_indices, Q2_indices, correlator_lookuptable);
    }
    /*! 6. Repeat steps 1.-5. for all correlators in correlator_list. Where 
     *      applicable rVdaggerVr must be replaced by rVdaggerV in step 3. and 
     *      Q2 by Q1 in step 4.
     */
    else if (correlator.type == "C3+") {
      std::vector<size_t> rnd_vec_id_uncharged;
      rnd_vec_id_uncharged.emplace_back( set_rnd_vec_uncharged(quarks, 
                                 correlator.quark_numbers[1], 
                                 operator_lookuptable.ricQ1_lookup) );
      // these are just dummies
      rnd_vec_id_uncharged.emplace_back( rnd_vec_id_uncharged[0] );
      rnd_vec_id_uncharged.emplace_back( rnd_vec_id_uncharged[0] );
      
      std::vector<std::vector<size_t> > rvdv_indices;
      build_rVdaggerV_lookup(rnd_vec_id_uncharged, vdv_indices,
                             operator_lookuptable.rvdaggerv_lookuptable,
                             rvdv_indices);

      std::vector<size_t> rnd_vec_id_charged;
      for( auto const &op : operator_lookuptable.ricQ2_lookup)
      rnd_vec_id_charged.emplace_back(set_rnd_vec_charged(quarks, 
                                           correlator.quark_numbers[2], 
                                           correlator.quark_numbers[0], false,
                                           operator_lookuptable.ricQ2_lookup));
      // these are just dummies
      rnd_vec_id_charged.emplace_back( rnd_vec_id_uncharged[0] );
      rnd_vec_id_charged.emplace_back( rnd_vec_id_uncharged[0] );

      std::vector<std::vector<size_t> > rvdvr_indices;
      build_rVdaggerVr_lookup(rnd_vec_id_charged, vdv_indices,
                              operator_lookuptable.rvdaggervr_lookuptable,
                              rvdvr_indices);

      std::vector<std::vector<size_t> > Q1_indices(rvdv_indices.size(),
                            std::vector<size_t>(rvdv_indices[0].size()));
      build_Q1_lookup(correlator.quark_numbers[1], correlator.quark_numbers[2],
                      1, false, quantum_numbers, quarks, rvdv_indices, 
                      operator_lookuptable.ricQ2_lookup,
                      quarkline_lookuptable.Q1, Q1_indices);


      std::vector<std::vector<size_t> > Q2_indices(rvdvr_indices.size(),
                                 std::vector<size_t>(rvdvr_indices[0].size()));
      build_Q2_lookup(correlator.quark_numbers[2], correlator.quark_numbers[0],
                      0, quantum_numbers, quarks, vdv_indices, 
                      operator_lookuptable.ricQ2_lookup,
                      quarkline_lookuptable.Q2L, Q2_indices);

      build_C3c_lookup(quantum_numbers, hdf5_dataset_name, 
                       rvdvr_indices, Q1_indices, Q2_indices, 
                       correlator_lookuptable);
    }
    else if (correlator.type == "C4+D") {
      std::vector<size_t> rnd_vec_id;
      rnd_vec_id.emplace_back(set_rnd_vec_charged(quarks, 
                                            correlator.quark_numbers[1], 
                                            correlator.quark_numbers[0], false,
                                            operator_lookuptable.ricQ2_lookup));
      rnd_vec_id.emplace_back(set_rnd_vec_charged(quarks, 
                                            correlator.quark_numbers[3], 
                                            correlator.quark_numbers[2], false,
                                            operator_lookuptable.ricQ2_lookup));
      std::vector<std::vector<size_t> > rvdvr_indices;
      build_rVdaggerVr_lookup(rnd_vec_id, vdv_indices,
                              operator_lookuptable.rvdaggervr_lookuptable,
                              rvdvr_indices);
      std::vector<std::vector<size_t> > Q2_indices(rvdvr_indices.size(),
                                  std::vector<size_t>(rvdvr_indices[0].size()));
      build_Q2_lookup(correlator.quark_numbers[0], correlator.quark_numbers[1],
                      0, quantum_numbers, quarks, vdv_indices, 
                      operator_lookuptable.ricQ2_lookup,
                      quarkline_lookuptable.Q2V, Q2_indices);
      build_Q2_lookup(correlator.quark_numbers[2], correlator.quark_numbers[3],
                      2, quantum_numbers, quarks, vdv_indices, 
                      operator_lookuptable.ricQ2_lookup,
                      quarkline_lookuptable.Q2V, Q2_indices);
      build_C4cD_lookup(quantum_numbers, hdf5_dataset_name, 
                        rvdvr_indices, Q2_indices, correlator_lookuptable);
    }
    else if (correlator.type == "C4+V") {
      std::vector<size_t> rnd_vec_id;
      rnd_vec_id.emplace_back(set_rnd_vec_charged(quarks, 
                                            correlator.quark_numbers[1], 
                                            correlator.quark_numbers[0], false,
                                            operator_lookuptable.ricQ2_lookup));
      rnd_vec_id.emplace_back(set_rnd_vec_charged(quarks, 
                                            correlator.quark_numbers[3], 
                                            correlator.quark_numbers[2], false,
                                            operator_lookuptable.ricQ2_lookup));
      std::vector<std::vector<size_t> > rvdvr_indices;
      build_rVdaggerVr_lookup(rnd_vec_id, vdv_indices,
                              operator_lookuptable.rvdaggervr_lookuptable,
                              rvdvr_indices);
      std::vector<std::vector<size_t> > Q2_indices(rvdvr_indices.size(),
                                  std::vector<size_t>(rvdvr_indices[0].size()));
      build_Q2_lookup(correlator.quark_numbers[0], correlator.quark_numbers[1],
                      0, quantum_numbers, quarks, vdv_indices, 
                      operator_lookuptable.ricQ2_lookup,
                      quarkline_lookuptable.Q2V, Q2_indices);
      build_Q2_lookup(correlator.quark_numbers[2], correlator.quark_numbers[3],
                      2, quantum_numbers, quarks, vdv_indices, 
                      operator_lookuptable.ricQ2_lookup,
                      quarkline_lookuptable.Q2V, Q2_indices);
      build_C4cV_lookup(quantum_numbers, hdf5_dataset_name, 
                        rvdvr_indices, Q2_indices, correlator_lookuptable);
    }
    else if (correlator.type == "C4+C") {
      std::vector<size_t> rnd_vec_id;
      rnd_vec_id.emplace_back(set_rnd_vec_charged(quarks, 
                                            correlator.quark_numbers[0], 
                                            correlator.quark_numbers[1], false,
                                            operator_lookuptable.ricQ2_lookup));
      rnd_vec_id.emplace_back(set_rnd_vec_charged(quarks, 
                                            correlator.quark_numbers[2], 
                                            correlator.quark_numbers[3], false,
                                            operator_lookuptable.ricQ2_lookup));
      std::vector<std::vector<size_t> > rvdvr_indices;
      build_rVdaggerVr_lookup(rnd_vec_id, vdv_indices,
                              operator_lookuptable.rvdaggervr_lookuptable,
                              rvdvr_indices);
      std::vector<std::vector<size_t> > Q2_indices(rvdvr_indices.size(),
                                  std::vector<size_t>(rvdvr_indices[0].size()));
      build_Q2_lookup(correlator.quark_numbers[0], correlator.quark_numbers[1],
                      0, quantum_numbers, quarks, vdv_indices, 
                      operator_lookuptable.ricQ2_lookup,
                      quarkline_lookuptable.Q2V, Q2_indices);
      build_Q2_lookup(correlator.quark_numbers[2], correlator.quark_numbers[3],
                      2, quantum_numbers, quarks, vdv_indices, 
                      operator_lookuptable.ricQ2_lookup,
                      quarkline_lookuptable.Q2V, Q2_indices);
      build_C4cC_lookup(quantum_numbers, hdf5_dataset_name, 
                        rvdvr_indices, Q2_indices, correlator_lookuptable);
    }
    else if (correlator.type == "C4+B") {
      std::vector<size_t> rnd_vec_id;
      rnd_vec_id.emplace_back(set_rnd_vec_charged(quarks, 
                                            correlator.quark_numbers[0], 
                                            correlator.quark_numbers[1], false,
                                            operator_lookuptable.ricQ2_lookup));
      rnd_vec_id.emplace_back(set_rnd_vec_charged(quarks, 
                                            correlator.quark_numbers[2], 
                                            correlator.quark_numbers[3], false,
                                            operator_lookuptable.ricQ2_lookup));
      std::vector<std::vector<size_t> > rvdvr_indices;
      build_rVdaggerVr_lookup(rnd_vec_id, vdv_indices,
                              operator_lookuptable.rvdaggervr_lookuptable,
                              rvdvr_indices);
      std::vector<std::vector<size_t> > Q2_indices(rvdvr_indices.size(),
                                  std::vector<size_t>(rvdvr_indices[0].size()));
      build_Q2_lookup(correlator.quark_numbers[0], correlator.quark_numbers[1],
                      0, quantum_numbers, quarks, vdv_indices, 
                      operator_lookuptable.ricQ2_lookup,
                      quarkline_lookuptable.Q2L, Q2_indices);
      build_Q2_lookup(correlator.quark_numbers[2], correlator.quark_numbers[3],
                      2, quantum_numbers, quarks, vdv_indices, 
                      operator_lookuptable.ricQ2_lookup,
                      quarkline_lookuptable.Q2L, Q2_indices);
      build_C4cB_lookup(quantum_numbers, hdf5_dataset_name, 
                        rvdvr_indices, Q2_indices, correlator_lookuptable);
    }
    else if (correlator.type == "C20V") {
      // 3. build the lookuptable for rVdaggerV and return an array of indices
      //    corresponding to the 'quantum_numbers' computed in step 1.
      std::vector<size_t> rnd_vec_id;
      rnd_vec_id.emplace_back( set_rnd_vec_uncharged(quarks, 
                                 correlator.quark_numbers[0], 
                                 operator_lookuptable.ricQ1_lookup) );
      rnd_vec_id.emplace_back( set_rnd_vec_uncharged(quarks, 
                                 correlator.quark_numbers[1], 
                                 operator_lookuptable.ricQ1_lookup) );

      std::vector<std::vector<size_t> > rvdv_indices;
      build_rVdaggerV_lookup(rnd_vec_id, vdv_indices,
                             operator_lookuptable.rvdaggerv_lookuptable,
                             rvdv_indices);
      // 4. build the lookuptable for Q1 and return an array of indices
      //    corresponding to the 'quantum_numbers' computed in step 1.
      // The size of this lookuptable needs to be known beforehand, because it
      // is build recursevely!
      std::vector<std::vector<size_t> > Q1_indices(rvdv_indices.size(),
                            std::vector<size_t>(rvdv_indices[0].size()));
      build_Q1_lookup(correlator.quark_numbers[0], correlator.quark_numbers[1],
                      0, true, quantum_numbers, quarks, rvdv_indices, 
                      operator_lookuptable.ricQ2_lookup,
                      quarkline_lookuptable.Q1, Q1_indices);
      build_Q1_lookup(correlator.quark_numbers[1], correlator.quark_numbers[0],
                      1, true, quantum_numbers, quarks, rvdv_indices, 
                      operator_lookuptable.ricQ2_lookup,
                      quarkline_lookuptable.Q1, Q1_indices);
      // 5. build the lookuptable for the correlation functions
      build_C20V_lookup(hdf5_dataset_name, Q1_indices, 
                       correlator_lookuptable);
    }
    else if (correlator.type == "C20") {
      // 3. build the lookuptable for rVdaggerV and return an array of indices
      //    corresponding to the 'quantum_numbers' computed in step 1.
      std::vector<size_t> rnd_vec_id;
      rnd_vec_id.emplace_back( set_rnd_vec_uncharged(quarks, 
                                 correlator.quark_numbers[0], 
                                 operator_lookuptable.ricQ1_lookup) );
      rnd_vec_id.emplace_back( set_rnd_vec_uncharged(quarks, 
                                 correlator.quark_numbers[1], 
                                 operator_lookuptable.ricQ1_lookup) );

      std::vector<std::vector<size_t> > rvdv_indices;
      build_rVdaggerV_lookup(rnd_vec_id, vdv_indices,
                             operator_lookuptable.rvdaggerv_lookuptable,
                             rvdv_indices);
      // 4. build the lookuptable for Q1 and return an array of indices
      //    corresponding to the 'quantum_numbers' computed in step 1.
      // The size of this lookuptable needs to be known beforehand, because it
      // is build recursevely!
      std::vector<std::vector<size_t> > Q1_indices(rvdv_indices.size(),
                            std::vector<size_t>(rvdv_indices[0].size()));
      build_Q1_lookup(correlator.quark_numbers[0], correlator.quark_numbers[1],
                      0, false, quantum_numbers, quarks, rvdv_indices, 
                      operator_lookuptable.ricQ2_lookup,
                      quarkline_lookuptable.Q1, Q1_indices);
      build_Q1_lookup(correlator.quark_numbers[1], correlator.quark_numbers[0],
                      1, false, quantum_numbers, quarks, rvdv_indices, 
                      operator_lookuptable.ricQ2_lookup,
                      quarkline_lookuptable.Q1, Q1_indices);
      // 5. build the lookuptable for the correlation functions
      build_C20_lookup(hdf5_dataset_name, Q1_indices, 
                       correlator_lookuptable);
    }
    else if (correlator.type == "C30") {
      std::vector<size_t> rnd_vec_id;
      rnd_vec_id.emplace_back( set_rnd_vec_uncharged(quarks, 
                                 correlator.quark_numbers[0], 
                                 operator_lookuptable.ricQ1_lookup) );
      rnd_vec_id.emplace_back( set_rnd_vec_uncharged(quarks, 
                                 correlator.quark_numbers[1], 
                                 operator_lookuptable.ricQ1_lookup) );
      rnd_vec_id.emplace_back( set_rnd_vec_uncharged(quarks, 
                                 correlator.quark_numbers[2], 
                                 operator_lookuptable.ricQ1_lookup) );
      std::vector<std::vector<size_t> > rvdv_indices;
      build_rVdaggerV_lookup(rnd_vec_id, vdv_indices,
                             operator_lookuptable.rvdaggerv_lookuptable,
                             rvdv_indices);
      // The size of this lookuptable needs to be known beforehand, because it
      // is build recursevely!
      std::vector<std::vector<size_t> > Q1_indices(rvdv_indices.size(),
                            std::vector<size_t>(rvdv_indices[0].size()));
      build_Q1_lookup(correlator.quark_numbers[0], correlator.quark_numbers[1],
      // should be correlator.quark_numbers[1], correlator.quark_numbers[0],
                      0, false, quantum_numbers, quarks, rvdv_indices, 
                      operator_lookuptable.ricQ2_lookup,
                      quarkline_lookuptable.Q1, Q1_indices);
      build_Q1_lookup(correlator.quark_numbers[1], correlator.quark_numbers[2],
      // should be correlator.quark_numbers[2], correlator.quark_numbers[1],
                      1, false, quantum_numbers, quarks, rvdv_indices, 
                      operator_lookuptable.ricQ2_lookup,
                      quarkline_lookuptable.Q1, Q1_indices);
      build_Q1_lookup(correlator.quark_numbers[2], correlator.quark_numbers[0],
      // should be correlator.quark_numbers[0], correlator.quark_numbers[2],
                      2, false, quantum_numbers, quarks, rvdv_indices, 
                      operator_lookuptable.ricQ2_lookup,
                      quarkline_lookuptable.Q1, Q1_indices);
      build_C30_lookup(hdf5_dataset_name, Q1_indices, 
                       correlator_lookuptable);
    }
    else if (correlator.type == "C40D" || correlator.type == "C40V") {

      std::vector<size_t> rnd_vec_id;
      rnd_vec_id.emplace_back( set_rnd_vec_uncharged(quarks, 
                                 correlator.quark_numbers[0], 
                                 operator_lookuptable.ricQ1_lookup) );
      rnd_vec_id.emplace_back( set_rnd_vec_uncharged(quarks, 
                                 correlator.quark_numbers[1], 
                                 operator_lookuptable.ricQ1_lookup) );
      rnd_vec_id.emplace_back( set_rnd_vec_uncharged(quarks, 
                                 correlator.quark_numbers[2], 
                                 operator_lookuptable.ricQ1_lookup) );
      rnd_vec_id.emplace_back( set_rnd_vec_uncharged(quarks, 
                                 correlator.quark_numbers[3], 
                                 operator_lookuptable.ricQ1_lookup) );
      std::vector<std::vector<size_t> > rvdv_indices;
      build_rVdaggerV_lookup(rnd_vec_id, vdv_indices,
                             operator_lookuptable.rvdaggerv_lookuptable,
                             rvdv_indices);
      // The size of this lookuptable needs to be known beforehand, because it
      // is build recursevely!
      std::vector<std::vector<size_t> > Q1_indices(rvdv_indices.size(),
                            std::vector<size_t>(rvdv_indices[0].size()));
      build_Q1_lookup(correlator.quark_numbers[0], correlator.quark_numbers[1],
                      0, false, quantum_numbers, quarks, rvdv_indices, 
                      operator_lookuptable.ricQ2_lookup,
                      quarkline_lookuptable.Q1, Q1_indices);
      build_Q1_lookup(correlator.quark_numbers[1], correlator.quark_numbers[0],
                      1, false, quantum_numbers, quarks, rvdv_indices, 
                      operator_lookuptable.ricQ2_lookup,
                      quarkline_lookuptable.Q1, Q1_indices);
      build_Q1_lookup(correlator.quark_numbers[2], correlator.quark_numbers[3],
                      2, false, quantum_numbers, quarks, rvdv_indices, 
                      operator_lookuptable.ricQ2_lookup,
                      quarkline_lookuptable.Q1, Q1_indices);
      build_Q1_lookup(correlator.quark_numbers[3], correlator.quark_numbers[2],
                      3, false, quantum_numbers, quarks, rvdv_indices, 
                      operator_lookuptable.ricQ2_lookup,
                      quarkline_lookuptable.Q1, Q1_indices);
      if (correlator.type == "C40D")
        build_C40D_lookup(hdf5_dataset_name, Q1_indices, 
                          correlator_lookuptable);
      else
        build_C40V_lookup(hdf5_dataset_name, Q1_indices, 
                          correlator_lookuptable);
    }
    else if (correlator.type == "C40C" || correlator.type == "C40B") {
      std::vector<size_t> rnd_vec_id;
      rnd_vec_id.emplace_back( set_rnd_vec_uncharged(quarks, 
                                 correlator.quark_numbers[0], 
                                 operator_lookuptable.ricQ1_lookup) );
      rnd_vec_id.emplace_back( set_rnd_vec_uncharged(quarks, 
                                 correlator.quark_numbers[1], 
                                 operator_lookuptable.ricQ1_lookup) );
      rnd_vec_id.emplace_back( set_rnd_vec_uncharged(quarks, 
                                 correlator.quark_numbers[2], 
                                 operator_lookuptable.ricQ1_lookup) );
      rnd_vec_id.emplace_back( set_rnd_vec_uncharged(quarks, 
                                 correlator.quark_numbers[3], 
                                 operator_lookuptable.ricQ1_lookup) );
      std::vector<std::vector<size_t> > rvdv_indices;
      build_rVdaggerV_lookup(rnd_vec_id, vdv_indices,
                             operator_lookuptable.rvdaggerv_lookuptable,
                             rvdv_indices);
      // The size of this lookuptable needs to be known beforehand, because it
      // is build recursevely!
      std::vector<std::vector<size_t> > Q1_indices(rvdv_indices.size(),
                            std::vector<size_t>(rvdv_indices[0].size()));
      build_Q1_lookup(correlator.quark_numbers[0], correlator.quark_numbers[1],
                      0, false, quantum_numbers, quarks, rvdv_indices, 
                      operator_lookuptable.ricQ2_lookup,
                      quarkline_lookuptable.Q1, Q1_indices);
      build_Q1_lookup(correlator.quark_numbers[1], correlator.quark_numbers[2],
                      1, false, quantum_numbers, quarks, rvdv_indices, 
                      operator_lookuptable.ricQ2_lookup,
                      quarkline_lookuptable.Q1, Q1_indices);
      build_Q1_lookup(correlator.quark_numbers[2], correlator.quark_numbers[3],
                      2, false, quantum_numbers, quarks, rvdv_indices, 
                      operator_lookuptable.ricQ2_lookup,
                      quarkline_lookuptable.Q1, Q1_indices);
      build_Q1_lookup(correlator.quark_numbers[3], correlator.quark_numbers[0],
                      3, false, quantum_numbers, quarks, rvdv_indices, 
                      operator_lookuptable.ricQ2_lookup,
                      quarkline_lookuptable.Q1, Q1_indices);
      if (correlator.type == "C40C")
        build_C40C_lookup(hdf5_dataset_name, Q1_indices, 
                          correlator_lookuptable);
      else
        build_C40B_lookup(hdf5_dataset_name, Q1_indices, 
                          correlator_lookuptable);
    }
    else {
      std::cout << "Correlator type not known!" << std::endl;
      exit(0);
    }
  }

  /*! Sets index_of_unity to the index of operator_lookuptable.vdaggerv_lookup
   *  where momentum and displacement are both zero, or to -1 if no such entry
   *  is found.
   */
  std::array<int, 3> const zero{0,0,0};
  bool found = false;
  for(const auto& op_vdv : operator_lookuptable.vdaggerv_lookup)
    if( (op_vdv.momentum == zero) && (op_vdv.displacement == zero) ){
      operator_lookuptable.index_of_unity = op_vdv.id;
      found = true;
    }
  if(!found)
    operator_lookuptable.index_of_unity = -1;

}


