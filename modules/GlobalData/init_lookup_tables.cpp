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
  
  }
  else if (correlator.type == "C2+" || correlator.type == "C20") {
    for(const auto& op0 : qn_op[0]){
      for(const auto& op1 : qn_op[1]){ // all combinations of operators
        std::vector<QuantumNumbers> single_vec_qn;
        // momentum at source and sink must always be the same for 2pt fcts.
        if(op0.momentum == op1.momentum){ 
          single_vec_qn.emplace_back(op0); // TODO: might be possible to write
          single_vec_qn.emplace_back(op1); //       more elegantly
          quantum_numbers.emplace_back(single_vec_qn);
        }
    }}
  }
  else if (correlator.type == "C3+" || correlator.type == "C30") {
  
  }
  else if (correlator.type == "C4D" || correlator.type == "C4C") {
  
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
    for(const auto& qn : qn_row){ // adding quantum numbers
      std::stringstream result;
      std::copy(qn.momentum.begin(), qn.momentum.end(), 
                std::ostream_iterator<int>(result, ""));
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
                          const size_t id_q1, const size_t id_q2, 
                          std::vector<RandomIndexCombinationsQ2>& rnd_vec_ids) {

  // First, check if the random vector index combination already exists
  for(const auto& r_id : rnd_vec_ids)
    if(((r_id.id_q1 == id_q1) && (r_id.id_q2 == id_q2)))// ||
//       ((r_id.id_q1 == id_q2) && (r_id.id_q2 == id_q1)) )
      return r_id.id;

  // set start and end points of rnd numbers
  auto rndq1_start = 0;
  for(auto i = 0; i < id_q1; i++)
    rndq1_start =+ quarks[i].number_of_rnd_vec;
  auto rndq2_start = 0;
  for(auto i = 0; i < id_q2; i++)
    rndq2_start =+ quarks[i].number_of_rnd_vec;
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
  for(size_t i = rndq1_start; i < rndq1_end; ++i)
    for(size_t j = rndq2_start; j < rndq2_end; ++j)
      if(i != j) 
        rnd_vec_comb.emplace_back(i, j);
  rnd_vec_ids.emplace_back(RandomIndexCombinationsQ2(rnd_vec_ids.size(), id_q1, 
                                                  id_q2, offset, rnd_vec_comb));
  return rnd_vec_ids.back().id;
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
static void build_rVdaggerVr_lookup(const size_t rnd_vec_id, 
         const std::vector<std::vector<std::pair<size_t, bool> > >& vdv_indices,
         std::vector<VdaggerVRandomLookup>& rvdaggervr_lookup,
         std::vector<std::vector<size_t> >& rvdvr_indices) {

  for(const auto& vdv_row : vdv_indices){
    std::vector<size_t> rvdvr_indices_row;
    for(const auto& vdv : vdv_row){
       
      auto it = std::find_if(rvdaggervr_lookup.begin(), rvdaggervr_lookup.end(),
                           [&vdv, &rnd_vec_id](VdaggerVRandomLookup vdv_qn)
                           {
                             auto c1 = (vdv_qn.id_ricQ_lookup == rnd_vec_id);
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
             rvdaggervr_lookup.size(), vdv.first, rnd_vec_id, vdv.second));
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
    // unelegant, but what else can I co? This is necessary because every
    // operator has its own rnd vector.
    // TODO: Think about merging this and build_rVdaggerVr_lookup into one
    //       function. For n-point functions also the other one needs a vector
    //       of rnd_vec_id!
    for(size_t i = 0; i < vdv_row.size(); i++){
      const auto& vdv = vdv_row[i];
      const auto& rnd = rnd_vec_id[i];
 
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
    std::vector<size_t> Q2_indices_row;
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
      Q2_indices_row.emplace_back((*it).id);
    }
    else {
      size_t rnd_index = set_rnd_vec_charged(quarks, id_quark1, id_quark2, 
                                             ric_lookup);
      Q2V.emplace_back(QuarklineQ2Indices(Q2V.size(), vdv.first, id_quark1, 
                                   id_quark2, rnd_index, vdv.second, qn.gamma));
      Q2_indices_row.emplace_back(Q2V.back().id);
    }
    Q2_indices.emplace_back(Q2_indices_row);
  }
}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
static void build_Q1_lookup(const size_t id_quark_used, 
         const size_t id_quark_connected, const size_t operator_id,
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
                                                id_quark_connected, ric_lookup);
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
                         const std::vector<std::pair<
                                  std::string, std::string> >& correlator_names,
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
                      std::vector<size_t>((*it).id)));
      }
      else {
        corr_lookup.corrC.emplace_back(CorrInfo(corr_lookup.corrC.size(), 
                                       "", "", indices));
        corr_lookup.C2c.emplace_back(CorrInfo(corr_lookup.C2c.size(), 
                                     correlator_names[row].first,
                                     correlator_names[row].second, 
                           std::vector<size_t>(corr_lookup.corrC.back().id)));
      }
    }
  }  
}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
static void build_C20_lookup(
                         const std::vector<std::pair<
                                  std::string, std::string> >& correlator_names,
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
                      std::vector<size_t>((*it).id)));
      }
      else {
        corr_lookup.corr0.emplace_back(CorrInfo(corr_lookup.corr0.size(), 
                                       "", "", indices));
        corr_lookup.C20.emplace_back(CorrInfo(corr_lookup.C20.size(), 
                                     correlator_names[row].first,
                                     correlator_names[row].second, 
                           std::vector<size_t>(corr_lookup.corr0.back().id)));
      }
    }
    row++;
  }  
}
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

    if (correlator.type == "C2+") {
      // 3. build the lookuptable for rVdaggerVr and return an array of indices
      //    corresponding to the 'quantum_numbers' computed in step 1.
      size_t rnd_vec_id = set_rnd_vec_charged(quarks, 
                                  correlator.quark_numbers[0], 
                                  correlator.quark_numbers[1], 
                                  operator_lookuptable.ricQ2_lookup);
      std::vector<std::vector<size_t> > rvdvr_indices;
      build_rVdaggerVr_lookup(rnd_vec_id, vdv_indices,
                              operator_lookuptable.rvdaggervr_lookuptable,
                              rvdvr_indices);
      // 4. build the lookuptable for Q2 and return an array of indices
      //    corresponding to the 'quantum_numbers' computed in step 1.
      std::vector<std::vector<size_t> > Q2_indices;
      build_Q2_lookup(correlator.quark_numbers[0], correlator.quark_numbers[1],
                      0, quantum_numbers, quarks, vdv_indices, 
                      operator_lookuptable.ricQ2_lookup,
                      quarkline_lookuptable.Q2V, Q2_indices);
      // 5. build the lookuptable for the correlation functions
      build_C2c_lookup(correlator_names, rvdvr_indices, Q2_indices, 
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
                      0, quantum_numbers, quarks, rvdv_indices, 
                      operator_lookuptable.ricQ2_lookup,
                      quarkline_lookuptable.Q1, Q1_indices);
      build_Q1_lookup(correlator.quark_numbers[1], correlator.quark_numbers[0],
                      1, quantum_numbers, quarks, rvdv_indices, 
                      operator_lookuptable.ricQ2_lookup,
                      quarkline_lookuptable.Q1, Q1_indices);
      // 5. build the lookuptable for the correlation functions
      build_C20_lookup(correlator_names, Q1_indices, correlator_lookuptable);
    }
  }
  // finding the index where we have no momentum and no displacement
  const std::array<int, 3> zero = {0,0,0};
  for(const auto& op_vdv : operator_lookuptable.vdaggerv_lookup)
    if( (op_vdv.momentum == zero) && (op_vdv.displacement == zero) )
      operator_lookuptable.index_of_unity = op_vdv.id;
}
// TODO: Build a function to change outpahts in the correlator lookup tables for
//       new configuration numbers

























