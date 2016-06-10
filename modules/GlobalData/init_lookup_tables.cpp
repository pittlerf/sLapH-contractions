#include "global_data.h"
#include "global_data_utils.h"

struct QuantumNumbers{
  std::array<int, 3> momentum;
  std::array<int, 3> displacement;
  std::vector<int> gamma;

  inline void write() const{
    std::cout << "\tmomentum: " << momentum[0] << momentum[1] << momentum[2];
    std::cout << "\n\tdisplacement: " << displacement[0] << displacement[1] 
              << displacement[2] << "\n\tgamma struct: ";
    for(const auto& g : gamma)
      std::cout << g;
    std::cout << "\n" << std::endl;
  }
};
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
static std::array<int, 3> change_sign_array(const std::array<int, 3>& in){
  return {{-in[0], -in[1], -in[2]}};
}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
static int compute_norm_squ(const std::array<int, 3>& in){
  return in[0]*in[0] + in[1]*in[1] + in[2]*in[2];
}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
static int add_momenta_squared(const std::array<int, 3>& in1, 
                               const std::array<int, 3>& in2){
  return (in1[0]+in2[0]) * (in1[0]+in2[0]) + 
         (in1[1]+in2[1]) * (in1[1]+in2[1]) + 
         (in1[2]+in2[2]) * (in1[2]+in2[2]);
}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
static std::array<int, 3> add_momenta(const std::array<int, 3>& in1, 
                                      const std::array<int, 3>& in2){
  return {{in1[0]+in2[0], in1[1]+in2[1], in1[2]+in2[2]}};
}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void build_quantum_numbers_from_correlator_list(const Correlators& correlator, 
                    const std::vector<Operator_list>& operator_list,
                    std::vector<std::vector<QuantumNumbers> >& quantum_numbers){

  // Extracting all possible quantum number combinations for this correlation
  // function. The separation of what is actually computed is done in the
  // if statements below because it depends on the number of quarks.
  std::vector<std::vector<QuantumNumbers> > qn_op;
  QuantumNumbers write;
  for(const auto& op_number : correlator.operator_numbers){
    std::vector<QuantumNumbers> single_vec_qn;
    for(const auto& op : operator_list[op_number]){ 
      write.gamma = op.gammas; // Dirac Matrices
      write.displacement = op.dil_vec; // displacement Vector
      for(const auto& mom_vec : op.mom_vec){ // momenta
        for(auto mom : mom_vec){
          write.momentum = mom;
          single_vec_qn.emplace_back(write);
      }}
    }
    qn_op.emplace_back(single_vec_qn);
  }

  // TODO: think about a way to avoid these if conditions 
  if (correlator.type == "C1" || correlator.type == "C1T") {
    for(const auto& op0 : qn_op[0])
      quantum_numbers.emplace_back(std::vector<QuantumNumbers>({op0}));
  }
  else if (correlator.type == "C2+" || correlator.type == "C20") {
    for(const auto& op0 : qn_op[0]){
      for(const auto& op1 : qn_op[1]){ // all combinations of operators
        std::vector<QuantumNumbers> single_vec_qn;
        // momentum at source and sink must always be the same for 2pt fcts.
        if(op0.momentum == change_sign_array(op1.momentum)){ 
          single_vec_qn.emplace_back(op0); // TODO: might be possible to write
          single_vec_qn.emplace_back(op1); //       more elegantly
          quantum_numbers.emplace_back(single_vec_qn);
        }
    }}
  }
  else if (correlator.type == "C3+" || correlator.type == "C30") {
    size_t counter_test = 0;
    size_t counter_mom0 = 0;
    size_t counter_mom1 = 0;
    size_t counter_mom2 = 0;
    size_t counter_mom3 = 0;
    size_t counter_mom4 = 0;
    for(const auto& op0 : qn_op[0]){
    for(const auto& op2 : qn_op[2]){ 
      const int mom0 = compute_norm_squ(op0.momentum);
      const int mom2 = compute_norm_squ(op2.momentum);
      const int tot_mom_l = add_momenta_squared(op0.momentum, op2.momentum);
      std::array<int, 3> tot_mom_v_l = add_momenta(op0.momentum, op2.momentum);
      
    for(const auto& op1 : qn_op[1]){ // all combinations of operators

      // momenta at source and sink must be equal - sign comes from daggering
      if((tot_mom_v_l[0] != -op1.momentum[0]) ||
         (tot_mom_v_l[1] != -op1.momentum[1]) ||
         (tot_mom_v_l[2] != -op1.momentum[2]))
        continue;

      if(tot_mom_l == 0){
        if(mom0 > 3 || mom0 == 0)
          continue;
        counter_mom0++;
      }
      else if(tot_mom_l == 1){
        if((mom0 + mom2) > 5)
          continue;
        counter_mom1++;
      }
      else if(tot_mom_l == 2){
        if((mom0 + mom2) > 6)
          continue;
        counter_mom2++;
      }
      else if(tot_mom_l == 3){
        if((mom0 + mom2) > 7)
          continue;
        counter_mom3++;
      }
      else if(tot_mom_l == 4){
        if((mom0 + mom2) > 4)
          continue;
        counter_mom4++;
      }
      else
        continue; // maximum momentum is 4

      counter_test++;
      std::vector<QuantumNumbers> single_vec_qn = {op0, op1, op2};
      quantum_numbers.emplace_back(single_vec_qn);
    }}}
    std::cout << "test finished - combinations: " << counter_test << std::endl;
    std::cout << "combination mom0: " << counter_mom0 << std::endl;
    std::cout << "combination mom1: " << counter_mom1 << std::endl;
    std::cout << "combination mom2: " << counter_mom2 << std::endl;
    std::cout << "combination mom3: " << counter_mom3 << std::endl;
    std::cout << "combination mom4: " << counter_mom4 << std::endl;
  }
  else if (correlator.type == "C4+D") {
    // momentum combinations on source side ------------------------------------
    size_t counter_test = 0;
    size_t counter_mom0 = 0;
    size_t counter_mom1 = 0;
    size_t counter_mom2 = 0;
    size_t counter_mom3 = 0;
    size_t counter_mom4 = 0;
    for(const auto& op0 : qn_op[0]){
    for(const auto& op2 : qn_op[2]){
      const int mom0 = compute_norm_squ(op0.momentum);
      const int mom2 = compute_norm_squ(op2.momentum);
      const int tot_mom_l = add_momenta_squared(op0.momentum, op2.momentum);
      std::array<int, 3> tot_mom_v_l = add_momenta(op0.momentum, op2.momentum);
      
      if(tot_mom_l == 0){
        if(mom0 > 3)
          continue;
        counter_mom0++;
      }
      else if(tot_mom_l == 1){
        if((mom0 + mom2) > 5)
          continue;
        counter_mom1++;
      }
      else if(tot_mom_l == 2){
        if((mom0 + mom2) > 6)
          continue;
        counter_mom2++;
      }
      else if(tot_mom_l == 3){
        if((mom0 + mom2) > 7)
          continue;
        counter_mom3++;
      }
      else if(tot_mom_l == 4){
        if((mom0 + mom2) > 4)
          continue;
        counter_mom4++;
      }
      else
        continue; // maximum momentum is 4

    // momentum combinations on sink side --------------------------------------
    for(const auto& op1 : qn_op[1]){ 
    for(const auto& op3 : qn_op[3]){ // all combinations of operators
      const int mom1 = compute_norm_squ(op1.momentum);
      const int mom3 = compute_norm_squ(op3.momentum);
      const int tot_mom_r = add_momenta_squared(op1.momentum, op3.momentum);
      std::array<int, 3> tot_mom_v_r = add_momenta(op1.momentum, op3.momentum);
      if((tot_mom_v_r[0] != -tot_mom_v_l[0]) ||
         (tot_mom_v_r[1] != -tot_mom_v_l[1]) ||
         (tot_mom_v_r[2] != -tot_mom_v_l[2]))
        continue; // both total momenta must be equal

      if(tot_mom_r == 0){
        if(mom1 > 3)
          continue;
      }
      else if(tot_mom_r == 1){
        if((mom1 + mom3) > 5)
          continue;
      }
      else if(tot_mom_r == 2){
        if((mom1 + mom3) > 6)
          continue;
      }
      else if(tot_mom_r == 3){
        if((mom1 + mom3) > 7)
          continue;
      }
      else if(tot_mom_r == 4){
        if((mom1 + mom3) > 4)
          continue;
      }
      else
        continue; // maximum momentum is 4

      // create combinations ---------------------------------------------------
      std::vector<QuantumNumbers> single_vec_qn = {op0, op1, op2, op3};
      quantum_numbers.emplace_back(single_vec_qn);
      counter_test++;
    }}}}
    std::cout << "test finished - combinations: " << counter_test << std::endl;
    std::cout << "combination mom0: " << counter_mom0 << std::endl;
    std::cout << "combination mom1: " << counter_mom1 << std::endl;
    std::cout << "combination mom2: " << counter_mom2 << std::endl;
    std::cout << "combination mom3: " << counter_mom3 << std::endl;
    std::cout << "combination mom4: " << counter_mom4 << std::endl;
  }
  else if (correlator.type == "C4+B") {
    // momentum combinations on source side ------------------------------------
    size_t counter_test = 0;
    size_t counter_mom0 = 0;
    size_t counter_mom1 = 0;
    size_t counter_mom2 = 0;
    size_t counter_mom3 = 0;
    size_t counter_mom4 = 0;
    for(const auto& op0 : qn_op[0]){
    for(const auto& op3 : qn_op[3]){
      const int mom0 = compute_norm_squ(op0.momentum);
      const int mom3 = compute_norm_squ(op3.momentum);
      const int tot_mom_l = add_momenta_squared(op0.momentum, op3.momentum);
      std::array<int, 3> tot_mom_v_l = add_momenta(op0.momentum, op3.momentum);
      
      if(tot_mom_l == 0){
        if(mom0 > 3 || mom0 == 0)
          continue;
        counter_mom0++;
      }
      else if(tot_mom_v_l == std::array<int,3>({{0,0,1}})){
        if((mom0 + mom3) > 5)
          continue;
        counter_mom1++;
      }
      else if(tot_mom_v_l == std::array<int,3>({{0,1,1}})){
        if((mom0 + mom3) > 6)
          continue;
        counter_mom2++;
      }
      else if(tot_mom_v_l == std::array<int,3>({{1,1,1}})){
        if((mom0 + mom3) > 7)
          continue;
        counter_mom3++;
      }
      else if(tot_mom_v_l == std::array<int,3>({{0,0,2}})){
        if((mom0 + mom3) > 4)
          continue;
        counter_mom4++;
      }
      else
        continue; // maximum momentum is 4

    // momentum combinations on sink side --------------------------------------
    for(const auto& op1 : qn_op[1]){ 
    for(const auto& op2 : qn_op[2]){ // all combinations of operators
      const int mom1 = compute_norm_squ(op1.momentum);
      const int mom2 = compute_norm_squ(op2.momentum);
      const int tot_mom_r = add_momenta_squared(op1.momentum, op2.momentum);
      std::array<int, 3> tot_mom_v_r = add_momenta(op1.momentum, op2.momentum);
      if((tot_mom_v_r[0] != -tot_mom_v_l[0]) ||
         (tot_mom_v_r[1] != -tot_mom_v_l[1]) ||
         (tot_mom_v_r[2] != -tot_mom_v_l[2]))
        continue; // both total momenta must be equal

      if(tot_mom_r == 0){
        if(mom1 > 3 || mom1 == 0)
          continue;
      }
      else if(tot_mom_v_r == std::array<int,3>({{0,0,-1}})){
        if((mom1 + mom2) > 5)
          continue;
      }
      else if(tot_mom_v_r == std::array<int,3>({{0,-1,-1}})){
        if((mom1 + mom2) > 6)
          continue;
      }
      else if(tot_mom_v_r == std::array<int,3>({{-1,-1,-1}})){
        if((mom1 + mom2) > 7)
          continue;
      }
      else if(tot_mom_v_r == std::array<int,3>({{0,0,-2}})){
        if((mom1 + mom2) > 4)
          continue;
      }
      else
        continue; // maximum momentum is 4

      // create combinations ---------------------------------------------------
      std::vector<QuantumNumbers> single_vec_qn = {op0, op1, op2, op3};
      quantum_numbers.emplace_back(single_vec_qn);
      counter_test++;
    }}}}
    std::cout << "test finished - combinations: " << counter_test << std::endl;
    std::cout << "combination mom0: " << counter_mom0 << std::endl;
    std::cout << "combination mom1: " << counter_mom1 << std::endl;
    std::cout << "combination mom2: " << counter_mom2 << std::endl;
    std::cout << "combination mom3: " << counter_mom3 << std::endl;
    std::cout << "combination mom4: " << counter_mom4 << std::endl;
  }
  else if (correlator.type == "C40D" || correlator.type == "C40V" ||
           correlator.type == "C40B" || correlator.type == "C40C" ||
           correlator.type == "C4+V" ||
           correlator.type == "C4+C") {
    for(const auto& op0 : qn_op[0]){
    for(const auto& op1 : qn_op[1]){ 
    for(const auto& op2 : qn_op[2]){ 
    for(const auto& op3 : qn_op[3]){ // all combinations of operators
      // TODO: This must be changed later if GEVP should be used!!!!!!!!!!!!!!!
      std::vector<QuantumNumbers> single_vec_qn = {op0, op1, op2, op3};
      quantum_numbers.emplace_back(single_vec_qn);
    }}}}
  }
}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
static void build_correlator_names(const std::string& corr_type, int cnfg,  
               const std::string& outpath, const std::string& overwrite,
               const std::vector<std::string>& quark_types, 
               const std::vector<std::vector<QuantumNumbers> >& quantum_numbers,
               std::vector<std::pair<std::string, std::string> >& corr_names){
  
  for(const auto& qn_row : quantum_numbers){
    std::string pathname = outpath + "/cnfg" + std::to_string(cnfg)  + "/" 
                               + corr_type + "/";
    std::string filename =  corr_type + "_";
    for(const auto& qt : quark_types) // adding quark content
      filename += qt;
    size_t id = 0;
    for(const auto& qn : qn_row){ // adding quantum numbers
      std::stringstream result;
      std::copy(qn.momentum.begin(), qn.momentum.end(), 
                std::ostream_iterator<int>(result, ""));
      if(id == 0)
        pathname += ("first_p_" + result.str() + "/");
      id++;
      filename += ("_p" + result.str());
      result.str("");
      std::copy(qn.displacement.begin(), qn.displacement.end(), 
                std::ostream_iterator<int>(result, ""));
      filename += (".d" + result.str());
      result.str("");
      std::copy(qn.gamma.begin(), qn.gamma.end(), 
                std::ostream_iterator<int>(result, ""));
      filename += (".g" + result.str());
    }
    filename += ".dat";
    // check if the file already exists and terminate program if it should not
    // be overwriten
    if(overwrite == "no"){
      struct stat buffer;
      if((stat ((pathname+filename).c_str(), &buffer) == 0)){
        std::cout << "Program terminated because outfile already exists!" 
                  << std::endl;  
        exit(0);
      }
    }
    corr_names.emplace_back(std::make_pair(pathname, filename));
  }
}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
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
                             auto c1 = (vdv_qn.displacement == qn.displacement);
                             auto c2 = (vdv_qn.momentum == qn.momentum);
                             // also negative momentum is checked
                             const std::array<int, 3> pm = {-qn.momentum[0],
                                                            -qn.momentum[1],
                                                            -qn.momentum[2]};
                             auto c3 = (vdv_qn.momentum == pm);
                             // TODO: Think about the daggering!!
                             const std::array<int, 3> zero = {0,0,0};
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
                vdaggerv_lookup.size(), qn.momentum, qn.displacement));
        vdv_indices_row.emplace_back(vdaggerv_lookup.back().id, false);
      }
    }
    vdv_indices.emplace_back(vdv_indices_row);
  }
}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// function to obtain index combinations of random vectors for charged corr
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
  // TODO: if statements can be joined to make it shorter
  if(id_q1 == id_q2){
    if(quarks[id_q1].number_of_rnd_vec < 2){
      std::cerr << "There are not enough random vectors for charged correlators"
                << std::endl;
      exit(-1);
    }
  }
  else{
    if((quarks[id_q1].number_of_rnd_vec < 1 && 
        quarks[id_q2].number_of_rnd_vec < 1)){
      std::cerr << "There are not enough random vectors for charged correlators"
                << std::endl;
      exit(-1);
    }
  }

  // finally filling the array
  std::pair<size_t, size_t> offset = std::make_pair(rndq1_start, rndq2_start);
  std::vector<std::pair<size_t, size_t> > rnd_vec_comb;
  if(!C1){
    for(size_t i = rndq1_start; i < rndq1_end; ++i)
      for(size_t j = rndq2_start; j < rndq2_end; ++j)
        if(i != j)
          rnd_vec_comb.emplace_back(i, j);
  }
  else {
    for(size_t i = rndq1_start; i < rndq1_end; ++i)
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
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// function to obtain index combinations of random vectors for uncharged corrs
static size_t set_rnd_vec_uncharged(const std::vector<quark>& quarks, 
                          const size_t id_q1, 
                          std::vector<RandomIndexCombinationsQ1>& rnd_vec_ids) {

  // First, check if the random vector indices already exists
  for(const auto& r_id : rnd_vec_ids)
    if(r_id.id_q1 == id_q1) 
      return r_id.id;

  // set start and end points of rnd numbers
  auto rndq1_start = 0;
  for(auto i = 0; i < id_q1; i++)
    rndq1_start =+ quarks[i].number_of_rnd_vec;
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
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
static void build_rVdaggerVr_lookup(const std::vector<size_t>& rnd_vec_id, 
         const std::vector<std::vector<std::pair<size_t, bool> > >& vdv_indices,
         std::vector<VdaggerVRandomLookup>& rvdaggervr_lookup,
         std::vector<std::vector<size_t> >& rvdvr_indices) {

  for(const auto& vdv_row : vdv_indices){
    std::vector<size_t> rvdvr_indices_row;

    for(size_t vdv_id = 0; vdv_id < vdv_row.size(); vdv_id++){

      size_t rnd_index = 0; // construct to get correct random number indices
      if(vdv_id == 2 || vdv_id == 3)
        rnd_index = 1;

      const auto vdv = vdv_row[vdv_id];
     
      auto it = std::find_if(rvdaggervr_lookup.begin(), rvdaggervr_lookup.end(),
                             [&](VdaggerVRandomLookup vdv_qn)
                             {
                               auto c1 = (vdv_qn.id_ricQ_lookup == 
                                          rnd_vec_id[rnd_index]);
                               auto c2 = (vdv_qn.id_vdaggerv == vdv.first);
                               auto c3 = (vdv_qn.need_vdaggerv_daggering == 
                                          vdv.second); 
                               return (c1 && c2 && c3);
                             });

      if(it != rvdaggervr_lookup.end()) {
        rvdvr_indices_row.emplace_back((*it).id);
      }
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
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
static void build_rVdaggerV_lookup(const std::vector<size_t> rnd_vec_id, 
         const std::vector<std::vector<std::pair<size_t, bool> > >& vdv_indices,
         std::vector<VdaggerVRandomLookup>& rvdaggerv_lookup,
         std::vector<std::vector<size_t> >& rvdv_indices) {

  for(const auto& vdv_row : vdv_indices){
    std::vector<size_t> rvdv_indices_row;
    // TODO: Think about merging this and build_rVdaggerVr_lookup into one
    //       function. For n-point functions also the other one needs a vector
    //       of rnd_vec_id!
    for(size_t i = 0; i < vdv_row.size(); i++){
      const auto& vdv = vdv_row.at(i);
      const auto& rnd = rnd_vec_id.at(i);

      auto it = std::find_if(rvdaggerv_lookup.begin(), rvdaggerv_lookup.end(),
                           [&](VdaggerVRandomLookup vdv_qn)
                           {
                             auto c1 = (vdv_qn.id_ricQ_lookup == rnd);
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
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
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

    auto it = std::find_if(Q2V.begin(), Q2V.end(),
                         [&](QuarklineQ2Indices q2)
                         {
                           auto c1 = (q2.id_peram1 == id_quark1);
                           auto c2 = (q2.id_peram2 == id_quark2);
                           auto c3 = (q2.gamma == qn.gamma);
                           auto c4 = (q2.need_vdaggerv_dag == vdv.second);
                           auto c5 = (q2.id_vdaggerv == vdv.first);
                           return c1 && c2 && c3 && c4 && c5;
                         });
    if(it != Q2V.end()) {
      Q2_indices[row][operator_id] = (*it).id;
    }
    else {
      size_t rnd_index = set_rnd_vec_charged(quarks, id_quark1, id_quark2, 
                                             false, ric_lookup);
      Q2V.emplace_back(QuarklineQ2Indices(Q2V.size(), vdv.first, id_quark1, 
                                   id_quark2, rnd_index, vdv.second, qn.gamma));
      Q2_indices[row][operator_id] = Q2V.back().id;
    }
  }
}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
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
    const size_t rnd_index = set_rnd_vec_charged(quarks, id_quark_used,
                                           id_quark_connected, C1, ric_lookup);
    auto it = std::find_if(Q1.begin(), Q1.end(),
                         [&](QuarklineQ1Indices q1)
                         {
                           auto c1 = (q1.id_peram == id_quark_used);
                           auto c2 = (q1.gamma == qn.gamma);
                           auto c3 = (q1.id_rvdaggerv == rvdv);
                           auto c4 = (q1.id_ric_lookup == rnd_index);
                           return c1 && c2 && c3 && c4;
                         });
    if(it != Q1.end()) {
        Q1_indices[row][operator_id] = (*it).id;
    }
    else {
      Q1.emplace_back(QuarklineQ1Indices(Q1.size(), rvdv, id_quark_used, 
                                                         rnd_index, qn.gamma));
      Q1_indices[row][operator_id] = Q1.back().id;
    }
  }
}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
static void build_C2c_lookup( 
      const std::vector<std::vector<QuantumNumbers> >& quantum_numbers, 
      const std::vector<std::pair<std::string, std::string> >& correlator_names,
      const std::vector<std::vector<size_t> >& rvdvr_indices,
      const std::vector<std::vector<size_t> >& Q2_indices, 
      CorrelatorLookup& corr_lookup){

  for(size_t row = 0; row < correlator_names.size(); row++){
    std::vector<size_t> indices = {Q2_indices[row][0], rvdvr_indices[row][1]};
    auto it_C2c = std::find_if(corr_lookup.C2c.begin(), corr_lookup.C2c.end(),
                         [&](CorrInfo corr)
                         {
                           return (corr.outfile==correlator_names[row].second); 
                         });
    if(it_C2c == corr_lookup.C2c.end()){
      auto it = std::find_if(corr_lookup.corrC.begin(), corr_lookup.corrC.end(),
                             [&](CorrInfo corr)
                             {
                               return (corr.lookup == indices); 
                             });
      if(it != corr_lookup.corrC.end()){
        corr_lookup.C2c.emplace_back(CorrInfo(corr_lookup.C2c.size(), 
                      correlator_names[row].first, correlator_names[row].second,
                      std::vector<size_t>({(*it).id}),
                      std::vector<int>({})));
      }
      else {
        corr_lookup.corrC.emplace_back(CorrInfo(corr_lookup.corrC.size(), 
                            "", "", indices, quantum_numbers[row][1].gamma));
        corr_lookup.C2c.emplace_back(CorrInfo(corr_lookup.C2c.size(), 
                           correlator_names[row].first,
                           correlator_names[row].second, 
                           std::vector<size_t>({corr_lookup.corrC.back().id}),
                           std::vector<int>({})));
      }
    }
  }  
}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
static void build_C4cD_lookup( 
      const std::vector<std::vector<QuantumNumbers> >& quantum_numbers, 
      const std::vector<std::pair<std::string, std::string> >& correlator_names,
      const std::vector<std::vector<size_t> >& rvdvr_indices,
      const std::vector<std::vector<size_t> >& Q2_indices, 
      CorrelatorLookup& corr_lookup){

  for(size_t row = 0; row < correlator_names.size(); row++){
    std::vector<size_t> indices1 = {Q2_indices[row][0], rvdvr_indices[row][1]};
    std::vector<size_t> indices2 = {Q2_indices[row][2], rvdvr_indices[row][3]};
    auto it_C4cD = std::find_if(corr_lookup.C4cD.begin(), 
                                corr_lookup.C4cD.end(),
                          [&](CorrInfo corr)
                          {
                            return (corr.outfile==correlator_names[row].second); 
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
                           "", "", indices1, quantum_numbers[row][1].gamma));
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
                           "", "", indices2, quantum_numbers[row][3].gamma));
        id2 = corr_lookup.corrC.back().id;
      }
      else 
        id2 = (*it2).id;

      corr_lookup.C4cD.emplace_back(CorrInfo(corr_lookup.C4cD.size(), 
                      correlator_names[row].first, correlator_names[row].second, 
                      std::vector<size_t>({id1, id2}), std::vector<int>({})));
    }
  }  
}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
static void build_C4cV_lookup( 
      const std::vector<std::vector<QuantumNumbers> >& quantum_numbers, 
      const std::vector<std::pair<std::string, std::string> >& correlator_names,
      const std::vector<std::vector<size_t> >& rvdvr_indices,
      const std::vector<std::vector<size_t> >& Q2_indices, 
      CorrelatorLookup& corr_lookup){

  for(size_t row = 0; row < correlator_names.size(); row++){
    std::vector<size_t> indices1 = {Q2_indices[row][0], rvdvr_indices[row][1]};
    std::vector<size_t> indices2 = {Q2_indices[row][2], rvdvr_indices[row][3]};
    auto it_C4cV = std::find_if(corr_lookup.C4cV.begin(), 
                                corr_lookup.C4cV.end(),
                          [&](CorrInfo corr)
                          {
                            return (corr.outfile==correlator_names[row].second); 
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
                           "", "", indices1, quantum_numbers[row][1].gamma));
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
                           "", "", indices2, quantum_numbers[row][3].gamma));
        id2 = corr_lookup.corrC.back().id;
      }
      else 
        id2 = (*it2).id;

      corr_lookup.C4cV.emplace_back(CorrInfo(corr_lookup.C4cV.size(), 
                      correlator_names[row].first, correlator_names[row].second, 
                      std::vector<size_t>({id1, id2}), std::vector<int>({})));
    }
  }  
}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
static void build_C4cC_lookup( 
      const std::vector<std::vector<QuantumNumbers> >& quantum_numbers, 
      const std::vector<std::pair<std::string, std::string> >& correlator_names,
      const std::vector<std::vector<size_t> >& rvdvr_indices,
      const std::vector<std::vector<size_t> >& Q2_indices, 
      CorrelatorLookup& corr_lookup){

  for(size_t row = 0; row < correlator_names.size(); row++){
    std::vector<size_t> indices = {Q2_indices[row][0], rvdvr_indices[row][1], 
                                   Q2_indices[row][2], rvdvr_indices[row][3]};
    auto it_C4cC = std::find_if(corr_lookup.C4cC.begin(), 
                                corr_lookup.C4cC.end(),
                          [&](CorrInfo corr)
                          {
                            return (corr.outfile==correlator_names[row].second); 
                          });

    if(it_C4cC == corr_lookup.C4cC.end()){
      std::vector<int> gammas = {{quantum_numbers[row][1].gamma[0], 
                                  quantum_numbers[row][3].gamma[0]}};
      corr_lookup.C4cC.emplace_back(CorrInfo(corr_lookup.C4cC.size(), 
                      correlator_names[row].first, correlator_names[row].second, 
                      indices, gammas));
    }
  }  
}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
static void build_C4cB_lookup( 
      const std::vector<std::vector<QuantumNumbers> >& quantum_numbers, 
      const std::vector<std::pair<std::string, std::string> >& correlator_names,
      const std::vector<std::vector<size_t> >& rvdvr_indices,
      const std::vector<std::vector<size_t> >& Q2_indices, 
      CorrelatorLookup& corr_lookup){

  for(size_t row = 0; row < correlator_names.size(); row++){
    std::vector<size_t> indices = {Q2_indices[row][0], rvdvr_indices[row][1], 
                                   Q2_indices[row][2], rvdvr_indices[row][3]};
    auto it_C4cB = std::find_if(corr_lookup.C4cB.begin(), 
                                corr_lookup.C4cB.end(),
                          [&](CorrInfo corr)
                          {
                            return (corr.outfile==correlator_names[row].second); 
                          });

    if(it_C4cB == corr_lookup.C4cB.end()){
      std::vector<int> gammas = {{quantum_numbers[row][1].gamma[0], 
                                  quantum_numbers[row][3].gamma[0]}};
      corr_lookup.C4cB.emplace_back(CorrInfo(corr_lookup.C4cB.size(), 
                      correlator_names[row].first, correlator_names[row].second, 
                      indices, gammas));
    }
  } 
}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
static void build_C3c_lookup( 
      const std::vector<std::vector<QuantumNumbers> >& quantum_numbers, 
      const std::vector<std::pair<std::string, std::string> >& correlator_names,
      const std::vector<std::vector<size_t> >& rvdvr_indices,
      const std::vector<std::vector<size_t> >& Q1_indices, 
      const std::vector<std::vector<size_t> >& Q2_indices, 
      CorrelatorLookup& corr_lookup){

  for(size_t row = 0; row < correlator_names.size(); row++){
    std::vector<size_t> indices = {Q2_indices[row][0], Q1_indices[row][1],  
                                   rvdvr_indices[row][2]};
    auto it_C3c = std::find_if(corr_lookup.C3c.begin(), 
                               corr_lookup.C3c.end(),
                          [&](CorrInfo corr)
                          {
                            return (corr.outfile==correlator_names[row].second); 
                          });
    if(it_C3c == corr_lookup.C3c.end()){
      std::vector<int> gammas = {{quantum_numbers[row][2].gamma[0]}};
      corr_lookup.C3c.emplace_back(CorrInfo(corr_lookup.C3c.size(), 
                      correlator_names[row].first, correlator_names[row].second, 
                      indices, gammas));
    }
  } 
}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------






// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
static void build_C20_lookup(
      const std::vector<std::pair<std::string, std::string> >& correlator_names,
      const std::vector<std::vector<size_t> >& Q1_indices, 
      CorrelatorLookup& corr_lookup){

  size_t row = 0;
  for(const auto& Q1 : Q1_indices){
    std::vector<size_t> indices = {Q1[0], Q1[1]};
    auto it_C20 = std::find_if(corr_lookup.C20.begin(), corr_lookup.C20.end(),
                          [&](CorrInfo corr)
                          {
                            return (corr.outfile==correlator_names[row].second);
                          });
    if(it_C20 == corr_lookup.C20.end()){
      auto it = std::find_if(corr_lookup.corr0.begin(), corr_lookup.corr0.end(),
                             [&](CorrInfo corr)
                             {
                               return (corr.lookup == indices); 
                             });
      if(it != corr_lookup.corr0.end()){
        corr_lookup.C20.emplace_back(CorrInfo(corr_lookup.C20.size(), 
                      correlator_names[row].first, correlator_names[row].second,
                      std::vector<size_t>({(*it).id}), 
                      std::vector<int>({})));
      }
      else {
        corr_lookup.corr0.emplace_back(CorrInfo(corr_lookup.corr0.size(), 
                                     "", "", indices, std::vector<int>({})));
        corr_lookup.C20.emplace_back(CorrInfo(corr_lookup.C20.size(), 
                           correlator_names[row].first,
                           correlator_names[row].second, 
                           std::vector<size_t>({corr_lookup.corr0.back().id}),
                           std::vector<int>({})));
      }
    }
    row++;
  }  
}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
static void build_C30_lookup(
      const std::vector<std::pair<std::string, std::string> >& correlator_names,
      const std::vector<std::vector<size_t> >& Q1_indices, 
      CorrelatorLookup& corr_lookup){

  size_t row = 0;
  for(const auto& Q1 : Q1_indices){
    auto it = std::find_if(corr_lookup.C30.begin(), corr_lookup.C30.end(),
                          [&](CorrInfo corr)
                          {
                            return (corr.outfile==correlator_names[row].second);
                          });
    if(it == corr_lookup.C30.end()){
      corr_lookup.C30.emplace_back(CorrInfo(corr_lookup.C30.size(), 
                      correlator_names[row].first, correlator_names[row].second, 
                      Q1, std::vector<int>({})));
    }
    row++;
  }  
}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
static void build_C40D_lookup(
      const std::vector<std::pair<std::string, std::string> >& correlator_names,
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
                            return (corr.outfile==correlator_names[row].second);
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
                                     "", "", indices1, std::vector<int>({})));
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
                                     "", "", indices2, std::vector<int>({})));
        id2 = corr_lookup.corr0.back().id;
      }
      else
        id2 = (*it2).id;

      corr_lookup.C40D.emplace_back(CorrInfo(corr_lookup.C40D.size(), 
                      correlator_names[row].first, correlator_names[row].second, 
                      std::vector<size_t>({id1, id2}), std::vector<int>({})));
    }
    row++;
  }  
}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
static void build_C40V_lookup(
      const std::vector<std::pair<std::string, std::string> >& correlator_names,
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
                            return (corr.outfile==correlator_names[row].second);
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
                                     "", "", indices1, std::vector<int>({})));
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
                                     "", "", indices2, std::vector<int>({})));
        id2 = corr_lookup.corr0.back().id;
      }
      else
        id2 = (*it2).id;

      corr_lookup.C40V.emplace_back(CorrInfo(corr_lookup.C40V.size(), 
                      correlator_names[row].first, correlator_names[row].second, 
                      std::vector<size_t>({id1, id2}), std::vector<int>({})));
    }
    row++;
  }  
}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
static void build_C40C_lookup(
      const std::vector<std::pair<std::string, std::string> >& correlator_names,
      const std::vector<std::vector<size_t> >& Q1_indices, 
      CorrelatorLookup& corr_lookup){

  size_t row = 0;
  for(const auto& Q1 : Q1_indices){
    auto it = std::find_if(corr_lookup.C40C.begin(), corr_lookup.C40C.end(),
                          [&](CorrInfo corr)
                          {
                            return (corr.outfile==correlator_names[row].second);
                          });
    if(it == corr_lookup.C40C.end()){
      corr_lookup.C40C.emplace_back(CorrInfo(corr_lookup.C40C.size(), 
                      correlator_names[row].first, correlator_names[row].second, 
                      Q1, std::vector<int>({})));
    }
    row++;
  }  
}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
static void build_C40B_lookup(
      const std::vector<std::pair<std::string, std::string> >& correlator_names,
      const std::vector<std::vector<size_t> >& Q1_indices, 
      CorrelatorLookup& corr_lookup){

  size_t row = 0;
  for(const auto& Q1 : Q1_indices){
    auto it = std::find_if(corr_lookup.C40B.begin(), corr_lookup.C40B.end(),
                          [&](CorrInfo corr)
                          {
                            return (corr.outfile==correlator_names[row].second);
                          });
    if(it == corr_lookup.C40B.end()){
      corr_lookup.C40B.emplace_back(CorrInfo(corr_lookup.C40B.size(), 
                      correlator_names[row].first, correlator_names[row].second, 
                      Q1, std::vector<int>({})));
    }
    row++;
  }  
}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
static void build_C1_lookup(
      const std::vector<std::vector<QuantumNumbers> >& quantum_numbers, 
      const std::vector<std::pair<std::string, std::string> >& correlator_names,
      const std::vector<std::vector<size_t> >& Q1_indices, 
      CorrelatorLookup& corr_lookup){

  size_t row = 0;
  for(const auto& Q1 : Q1_indices){
    std::vector<size_t> indices = {Q1[0]};
    auto it_C1 = std::find_if(corr_lookup.C1.begin(), corr_lookup.C1.end(),
                          [&](CorrInfo corr)
                          {
                            return (corr.outfile==correlator_names[row].second);
                          });
    if(it_C1 == corr_lookup.C1.end()){
      corr_lookup.C1.emplace_back(CorrInfo(corr_lookup.C1.size(), 
                      correlator_names[row].first, correlator_names[row].second,
                      indices, quantum_numbers[row][0].gamma));
    }
    row++;
  }  
}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------






// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void GlobalData::init_lookup_tables() {

  for (const auto& correlator : correlator_list){

    // 1. build an array (quantum_numbers) with all the quantum numbers needed 
    //    for this particular correlation function.
    std::vector<std::vector<QuantumNumbers> > quantum_numbers;
    build_quantum_numbers_from_correlator_list(correlator, operator_list, 
                                               quantum_numbers);

    // Build the correlator names
    std::vector<std::string> quark_types; 
    for(const auto& id : correlator.quark_numbers)
      quark_types.emplace_back(quarks[id].type);
    std::vector<std::pair<std::string, std::string> > correlator_names;
    build_correlator_names(correlator.type, start_config, path_output, 
                     overwrite, quark_types, quantum_numbers, correlator_names);

    // 2. build the lookuptable for VdaggerV and return an array of indices
    //    corresponding to the 'quantum_numbers' computed in step 1. In 
    //    'vdv_indices' the first entry is the id of vdv, the second tells us
    //    if vdv must be daggered to get the correct quantum numbers.
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
      build_C1_lookup(quantum_numbers, correlator_names, Q1_indices, 
                                                        correlator_lookuptable);
    }
    else if (correlator.type == "C2+") {
      // 3. build the lookuptable for rVdaggerVr and return an array of indices
      //    corresponding to the 'quantum_numbers' computed in step 1.
      std::vector<size_t> rnd_vec_id;
      rnd_vec_id.emplace_back(set_rnd_vec_charged(quarks, 
                                            correlator.quark_numbers[0], 
                                            correlator.quark_numbers[1], false,
                                            operator_lookuptable.ricQ2_lookup));
      std::vector<std::vector<size_t> > rvdvr_indices;
      build_rVdaggerVr_lookup(rnd_vec_id, vdv_indices,
                              operator_lookuptable.rvdaggervr_lookuptable,
                              rvdvr_indices);
      // 4. build the lookuptable for Q2 and return an array of indices
      //    corresponding to the 'quantum_numbers' computed in step 1.
      std::vector<std::vector<size_t> > Q2_indices(rvdvr_indices.size(),
                                  std::vector<size_t>(rvdvr_indices[0].size()));
      build_Q2_lookup(correlator.quark_numbers[0], correlator.quark_numbers[1],
                      0, quantum_numbers, quarks, vdv_indices, 
                      operator_lookuptable.ricQ2_lookup,
                      quarkline_lookuptable.Q2V, Q2_indices);
      // 5. build the lookuptable for the correlation functions
      build_C2c_lookup(quantum_numbers, correlator_names, rvdvr_indices, 
                       Q2_indices, correlator_lookuptable);
    }
    else if (correlator.type == "C3+") {
      std::vector<size_t> rnd_vec_id;
      rnd_vec_id.emplace_back( set_rnd_vec_uncharged(quarks, 
                                 correlator.quark_numbers[1], 
                                 operator_lookuptable.ricQ1_lookup) );
      rnd_vec_id.emplace_back(set_rnd_vec_charged(quarks, 
                                           correlator.quark_numbers[2], 
                                           correlator.quark_numbers[0], false,
                                           operator_lookuptable.ricQ2_lookup));
      // this is just a dummy
      rnd_vec_id.emplace_back( set_rnd_vec_uncharged(quarks, 
                                 correlator.quark_numbers[0], 
                                 operator_lookuptable.ricQ1_lookup) );

      std::vector<std::vector<size_t> > rvdv_indices;
      build_rVdaggerV_lookup(rnd_vec_id, vdv_indices,
                             operator_lookuptable.rvdaggerv_lookuptable,
                             rvdv_indices);
      std::vector<std::vector<size_t> > rvdvr_indices;
      build_rVdaggerVr_lookup(rnd_vec_id, vdv_indices,
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

      build_C3c_lookup(quantum_numbers, correlator_names, rvdvr_indices, 
                       Q1_indices, Q2_indices, correlator_lookuptable);
    }
    else if (correlator.type == "C4+D") {
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
      build_C4cD_lookup(quantum_numbers, correlator_names, rvdvr_indices, 
                        Q2_indices, correlator_lookuptable);
    }
    else if (correlator.type == "C4+V") {
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
      build_C4cV_lookup(quantum_numbers, correlator_names, rvdvr_indices, 
                        Q2_indices, correlator_lookuptable);
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
      build_C4cC_lookup(quantum_numbers, correlator_names, rvdvr_indices, 
                        Q2_indices, correlator_lookuptable);
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
      build_C4cB_lookup(quantum_numbers, correlator_names, rvdvr_indices, 
                        Q2_indices, correlator_lookuptable);
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
      build_C20_lookup(correlator_names, Q1_indices, correlator_lookuptable);
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
                      0, false, quantum_numbers, quarks, rvdv_indices, 
                      operator_lookuptable.ricQ2_lookup,
                      quarkline_lookuptable.Q1, Q1_indices);
      build_Q1_lookup(correlator.quark_numbers[1], correlator.quark_numbers[2],
                      1, false, quantum_numbers, quarks, rvdv_indices, 
                      operator_lookuptable.ricQ2_lookup,
                      quarkline_lookuptable.Q1, Q1_indices);
      build_Q1_lookup(correlator.quark_numbers[2], correlator.quark_numbers[0],
                      2, false, quantum_numbers, quarks, rvdv_indices, 
                      operator_lookuptable.ricQ2_lookup,
                      quarkline_lookuptable.Q1, Q1_indices);
      build_C30_lookup(correlator_names, Q1_indices, correlator_lookuptable);
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
        build_C40D_lookup(correlator_names, Q1_indices, correlator_lookuptable);
      else
        build_C40V_lookup(correlator_names, Q1_indices, correlator_lookuptable);
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
        build_C40C_lookup(correlator_names, Q1_indices, correlator_lookuptable);
      else
        build_C40B_lookup(correlator_names, Q1_indices, correlator_lookuptable);
    }
    else {
      std::cout << "Correlator type not known!" << std::endl;
      exit(0);
    }
  }
  // finding the index where we have no momentum and no displacement
  const std::array<int, 3> zero = {0,0,0};
  bool found = false;
  for(const auto& op_vdv : operator_lookuptable.vdaggerv_lookup)
    if( (op_vdv.momentum == zero) && (op_vdv.displacement == zero) ){
      operator_lookuptable.index_of_unity = op_vdv.id;
      found = true;
    }
  if(!found)
    operator_lookuptable.index_of_unity = -1;

}

























