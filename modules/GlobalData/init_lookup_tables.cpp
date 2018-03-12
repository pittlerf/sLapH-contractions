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
//TODO: Need this header? no function from it called here!
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
 *  @returns     true if @f$|\mathtt{p_tot}|^2 \in \mathtt{P}@f$, false
 *otherwise
 */
bool desired_total_momentum(Vector const &p_tot, std::vector<Vector> const &P) {
  /*! If no total momentum is specified, there is no selection.
   *  @todo It is better to force the user to specify P. Catch that in
   *        global_data_input_handling
   */
  if (P.empty()) {
    return true;
  }

  if (std::find(P.begin(), P.end(), p_tot) == P.end()) {
    return false;
  } else {
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
 *  @returns        true if @f$|\mathtt{p1}|^2 + |\mathtt{p2}|^2 <
 *\mathtt{cutoff} @f$,
 *                  false otherwise
 *
 *  @warning Cutoff is hardcoded in this function!
 *  @warning No cutoff for P > 4 is implemented!
 */
static bool momenta_below_cutoff(Vector const &p1, Vector const &p2) {
  std::map<int, int> cutoff;
  // default value is the default initializer of type int: cutoff[default] = 0;
  cutoff[0] = 4;
  cutoff[1] = 5;
  cutoff[2] = 6;
  cutoff[3] = 7;
  cutoff[4] = 4;

  const Vector p_tot = p1 + p2;

  if (p_tot.squaredNorm() > 4) {
    std::cout << "In momenta_below_cutoff(): WARNING! No cutoff for P > 4"
              << " implemented" << std::endl;
  }

  if (p1.squaredNorm() + p2.squaredNorm() > cutoff[p_tot.squaredNorm()]) {
    return false;
  } else {
    return true;
  }
}

}  // end of unnamed namespace

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
void build_quantum_numbers_from_correlator_list(
    const Correlators_2 &correlator,
    const Operator_list &operator_list,
    std::vector<std::vector<QuantumNumbers>> &quantum_numbers) {
  std::vector<Operators> qn_op;
  for (const auto &op_number : correlator.operator_numbers) {
    qn_op.emplace_back(operator_list[op_number]);
  }

  /*! Restriction to what shall actually be computed is done in if statements
   *  for each diagram because it depends on the number of quarks.
   *
   *  @todo Think about a way to avoid these if conditions.
   */
  if (correlator.type == "C1") {
    for (const auto &op0 : qn_op[0])
      quantum_numbers.emplace_back(std::vector<QuantumNumbers>({op0}));
  } else if (correlator.type == "C2+" || correlator.type == "C20" ||
             correlator.type == "C20V" || correlator.type == "Check") {
    // Build all combinations of operators and impose momentum conservation
    // and cutoffs
    for (const auto &op0 : qn_op[0]) {
      Vector p_so = op0.momentum;

      if (desired_total_momentum(p_so, correlator.tot_mom)) {
        for (const auto &op1 : qn_op[1]) {
          Vector p_si = op1.momentum;

          // momentum at source and sink must always be the same for 2pt fcts.
          if (p_so == -p_si) {
            std::vector<QuantumNumbers> single_vec_qn = {op0, op1};
            quantum_numbers.emplace_back(single_vec_qn);
          }
        }
      }
    }
  }

  else if (correlator.type == "C3+" || correlator.type == "C30") {
    std::cout << "Constructing momentum combinations for C3" << std::endl;

    std::map<int, int> counter; /*! initialized with zero */

    for (const auto &op0 : qn_op[0]) {
      for (const auto &op2 : qn_op[2]) {
        Vector p_so_1 = op0.momentum;
        Vector p_so_2 = op2.momentum;
        Vector p_so = p_so_1 + p_so_2;

        if (desired_total_momentum(p_so, correlator.tot_mom) &&
            momenta_below_cutoff(p_so_1, p_so_2)) {
          for (const auto &op1 : qn_op[1]) {
            Vector p_si = op1.momentum;

            if (desired_total_momentum(p_si, correlator.tot_mom)) {
              if (p_so == -p_si) {
                const int p_tot = p_si.squaredNorm();
                counter[p_tot]++;

                quantum_numbers.emplace_back(std::vector<QuantumNumbers>{op0, op1, op2});
              }
            }
          }
        }
      }
    }

    int total_number_of_combinations = 0;
    for (const auto c : counter) {
      std::cout << "\tCombinations for P = " << c.first << ": " << c.second << std::endl;
      total_number_of_combinations += c.second;
    }
    std::cout << "\tTest finished - Combinations: " << total_number_of_combinations
              << std::endl;
  } else if (correlator.type == "C30V") {
    std::cout << "Constructing momentum combinations for C3" << std::endl;

    std::map<int, int> counter; /*! initialized with zero */

    for (const auto &op0 : qn_op[0]) {
      for (const auto &op1 : qn_op[1]) {
        Vector p_so_1 = op0.momentum;
        Vector p_so_2 = op1.momentum;
        Vector p_so = p_so_1 + p_so_2;

        if (desired_total_momentum(p_so, correlator.tot_mom) &&
            momenta_below_cutoff(p_so_1, p_so_2)) {
          for (const auto &op2 : qn_op[2]) {
            Vector p_si = op2.momentum;

            if (desired_total_momentum(p_si, correlator.tot_mom)) {
              if (p_so == -p_si) {
                const int p_tot = p_si.squaredNorm();
                counter[p_tot]++;

                quantum_numbers.emplace_back(std::vector<QuantumNumbers>{op0, op1, op2});
              }
            }
          }
        }
      }
    }

    int total_number_of_combinations = 0;
    for (const auto c : counter) {
      std::cout << "\tCombinations for P = " << c.first << ": " << c.second << std::endl;
      total_number_of_combinations += c.second;
    }
    std::cout << "\tTest finished - Combinations: " << total_number_of_combinations
              << std::endl;
  }

  else if (correlator.type == "C4+D") {
    std::cout << "Constructing momentum combinations for C4+D" << std::endl;

    std::map<int, int> counter; /*! initialized with zero */

    for (const auto &op0 : qn_op[0]) {
      for (const auto &op2 : qn_op[2]) {
        Vector p_so_1 = op0.momentum;
        Vector p_so_2 = op2.momentum;
        Vector p_so = p_so_1 + p_so_2;

        if (desired_total_momentum(p_so, correlator.tot_mom) &&
            momenta_below_cutoff(p_so_1, p_so_2)) {
          for (const auto &op1 : qn_op[1]) {
            for (const auto &op3 : qn_op[3]) {
              Vector p_si_1 = op1.momentum;
              Vector p_si_2 = op3.momentum;
              Vector p_si = p_si_1 + p_si_2;

              if (desired_total_momentum(p_si, correlator.tot_mom) &&
                  momenta_below_cutoff(p_si_1, p_si_2)) {
                if (p_so == -p_si) {
                  const int p_tot = p_si.squaredNorm();
                  counter[p_tot]++;

                  // create combinations
                  quantum_numbers.emplace_back(
                      std::vector<QuantumNumbers>{op0, op1, op2, op3});
                }
              }
            }
          }
        }
      }
    }

    int total_number_of_combinations = 0;
    for (const auto c : counter) {
      std::cout << "\tCombinations for P = " << c.first << ": " << c.second << std::endl;
      total_number_of_combinations += c.second;
    }
    std::cout << "\tTest finished - Combinations: " << total_number_of_combinations
              << std::endl;
  }

  else if (correlator.type == "C4+B") {
    std::cout << "Constructing momentum combinations for C4+B" << std::endl;

    std::map<int, int> counter; /*! initialized with zero */

    for (const auto &op0 : qn_op[0]) {
      for (const auto &op3 : qn_op[3]) {
        Vector p_so_1 = op0.momentum;
        Vector p_so_2 = op3.momentum;
        Vector p_so = p_so_1 + p_so_2;

        if (desired_total_momentum(p_so, correlator.tot_mom) &&
            momenta_below_cutoff(p_so_1, p_so_2)) {
          for (const auto &op1 : qn_op[1]) {
            for (const auto &op2 : qn_op[2]) {
              Vector p_si_1 = op1.momentum;
              Vector p_si_2 = op2.momentum;
              Vector p_si = p_si_1 + p_si_2;

              if (desired_total_momentum(p_si, correlator.tot_mom) &&
                  momenta_below_cutoff(p_si_1, p_si_2)) {
                if (p_so == -p_si) {
                  const int p_tot = p_si.squaredNorm();
                  counter[p_tot]++;

                  // create combinations
                  quantum_numbers.emplace_back(
                      std::vector<QuantumNumbers>{op0, op1, op2, op3});
                }
              }
            }
          }
        }
      }
    }

    int total_number_of_combinations = 0;
    for (const auto c : counter) {
      std::cout << "\tCombinations for P = " << c.first << ": " << c.second << std::endl;
      total_number_of_combinations += c.second;
    }
    std::cout << "\tTest finished - Combinations: " << total_number_of_combinations
              << std::endl;

  }

  /*! @todo Check whether that is identical to C4+D */
  else if (correlator.type == "C4+C") {
    std::cout << "Constructing momentum combinations for C4+C" << std::endl;

    std::map<int, int> counter; /*! initialized with zero */

    for (const auto &op0 : qn_op[0]) {
      for (const auto &op2 : qn_op[2]) {
        Vector p_so_1 = op0.momentum;
        Vector p_so_2 = op2.momentum;
        Vector p_so = p_so_1 + p_so_2;

        if (desired_total_momentum(p_so, correlator.tot_mom) &&
            momenta_below_cutoff(p_so_1, p_so_2)) {
          for (const auto &op1 : qn_op[1]) {
            for (const auto &op3 : qn_op[3]) {
              Vector p_si_1 = op1.momentum;
              Vector p_si_2 = op3.momentum;
              Vector p_si = p_si_1 + p_si_2;

              if (desired_total_momentum(p_si, correlator.tot_mom) &&
                  momenta_below_cutoff(p_si_1, p_si_2)) {
                if (p_so == -p_si) {
                  const int p_tot = p_si.squaredNorm();
                  counter[p_tot]++;

                  // create combinations
                  quantum_numbers.emplace_back(
                      std::vector<QuantumNumbers>{op0, op1, op2, op3});
                }
              }
            }
          }
        }
      }
    }

    int total_number_of_combinations = 0;
    for (const auto c : counter) {
      std::cout << "\tCombinations for P = " << c.first << ": " << c.second << std::endl;
      total_number_of_combinations += c.second;
    }
    std::cout << "\tTest finished - Combinations: " << total_number_of_combinations
              << std::endl;
  }

  /*! @todo: For C40D, C40B, C40V, C40C, C4+V, C4+C still all combinations
   *         are built.
   *         This must be changed later if GEVP should be used!!!!!!!!!!!!!!!
   */
  else if (correlator.type == "C40D" || correlator.type == "C40V" ||
           correlator.type == "C40B" || correlator.type == "C40C" ||
           correlator.type == "C4+V") {
    for (const auto &op0 : qn_op[0]) {
      for (const auto &op1 : qn_op[1]) {
        for (const auto &op2 : qn_op[2]) {
          for (const auto &op3 : qn_op[3]) {  // all combinations of operators
            quantum_numbers.emplace_back(std::vector<QuantumNumbers>{op0, op1, op2, op3});
          }
        }
      }
    }
  }
}

//TODO: Not the right place, not sure where to place it otherwise
/*! Makes a string object of a displacement vector */
static std::string vector_to_string(const std::vector< std::pair<char,char> > &in){
  std::string out;
  if (in.empty()) out = "000";
  for (auto const& dis : in){ 
    out.push_back(dis.first);
    out.push_back(dis.second);
  }
  return out;
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
static void build_correlator_names(
    const std::string &corr_type,
    int cnfg,
    const std::string &outpath,
    const std::string &overwrite,
    const std::vector<std::string> &quark_types,
    const std::vector<std::vector<QuantumNumbers>> &quantum_numbers,
    std::vector<std::string> &hdf5_dataset_name) {
  for (const auto &qn_row : quantum_numbers) {
    std::string filename = corr_type + "_";
    for (const auto &qt : quark_types)  // adding quark content
      filename += qt;
    size_t id = 0;
    for (const auto &qn : qn_row) {  // adding quantum numbers
      filename += std::string("_p") + to_string(qn.momentum);
      filename += std::string(".d") + vector_to_string(qn.displacement);
      filename += std::string(".g") + to_string(qn.gamma);
    }
    hdf5_dataset_name.emplace_back(filename);
  }
}

static std::string const build_hdf5_dataset_name(
    const std::string &corr_type,
    int cnfg,
    const std::string &outpath,
    const std::string &overwrite,
    const std::vector<std::string> &quark_types,
    const std::vector<QuantumNumbers> &qn) {
  std::string filename = corr_type + "_";
  for (const auto &qt : quark_types)  // adding quark content
    filename += qt;
  size_t id = 0;
  for (const auto &op : qn) {  // adding quantum numbers
    filename += std::string("_p") + to_string(op.momentum);
    filename += std::string(".d") + vector_to_string(op.displacement);
    filename += std::string(".g") + to_string(op.gamma);
  }
  return filename;
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
    const std::vector<std::vector<QuantumNumbers>> &quantum_numbers,
    std::vector<VdaggerVQuantumNumbers> &vdaggerv_lookup,
    std::vector<std::vector<std::pair<size_t, bool>>> &vdv_indices) {
  for (const auto &qn_vec : quantum_numbers) {
    std::vector<std::pair<size_t, bool>> vdv_indices_row;
    for (const auto &qn : qn_vec) {
      bool dagger;
      // checking if the combination of quantum numbers already exists in
      // vdaggerv. The position is stored in the iterator it.
      auto it = std::find_if(
          vdaggerv_lookup.begin(),
          vdaggerv_lookup.end(),
          [&qn, &dagger](VdaggerVQuantumNumbers vdv_qn) {
            auto c1 = (vdv_qn.displacement == qn.displacement);
            auto c2 = (Vector(vdv_qn.momentum.data()) == qn.momentum);
            // also negative momentum is checked
            auto c3 = (Vector(vdv_qn.momentum.data()) == (-1) * qn.momentum);
            /*! @TODO: Think about the daggering!! */
            const Vector zero(0, 0, 0);
            if (c1 and c2) {
              dagger = false;
              return true;
            } else if ((c1 and c3) and (qn.displacement.empty() )) {
              dagger = true;
              return true;
            } else
              return false;
          });
      // If the quantum number combination already exists only the id is needed
      // otherwise a new element is created at the end of the lookuptable.
      if (it != vdaggerv_lookup.end()) {
        vdv_indices_row.emplace_back((*it).id, dagger);
      } else {
        vdaggerv_lookup.emplace_back(VdaggerVQuantumNumbers(
            vdaggerv_lookup.size(),
            {qn.momentum[0], qn.momentum[1], qn.momentum[2]},
            qn.displacement));
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
 *  @return                 pair of unique ids specifying quark flavor and
 *  			    number of random seed
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
static std::vector<std::pair<size_t, size_t>> create_rnd_vec_id(
    const std::vector<quark> &quarks,
    const size_t id_q1,
    const size_t id_q2,
    const bool C1) {
  // set start and end points of rnd numbers
  auto rndq1_start = 0;
  for (auto i = 0; i < id_q1; i++)
    rndq1_start += quarks[i].number_of_rnd_vec;
  auto rndq2_start = 0;
  for (auto i = 0; i < id_q2; i++)
    rndq2_start += quarks[i].number_of_rnd_vec;

  auto rndq1_end = rndq1_start + quarks[id_q1].number_of_rnd_vec;
  auto rndq2_end = rndq2_start + quarks[id_q2].number_of_rnd_vec;

  // check if there are enough random vectors
  if ((quarks[id_q1].number_of_rnd_vec < 1) || (quarks[id_q2].number_of_rnd_vec < 1) ||
      (id_q1 == id_q2 && quarks[id_q1].number_of_rnd_vec < 2)) {
    std::cerr << "There are not enough random vectors for charged correlators"
              << std::endl;
    exit(-1);
  }

  // finally filling the array
  std::pair<size_t, size_t> offset = std::make_pair(rndq1_start, rndq2_start);
  std::vector<std::pair<size_t, size_t>> rnd_vec_comb;
  if (!C1) {
    for (size_t i = rndq1_start; i < rndq1_end; ++i)
      for (size_t j = rndq2_start; j < rndq2_end; ++j)
        // To avoid bias, different random vectors must have different indices.
        if (i != j) {
          rnd_vec_comb.emplace_back(i, j);
        }
  } else {
    for (size_t i = rndq1_start; i < rndq1_end; ++i)
      // if C1 == True there is only one random vector and thus only same index
      // combinations are possible
      rnd_vec_comb.emplace_back(i, i);
  }

  return rnd_vec_comb;
}

/******************************************************************************/
/*! Create lookuptable where to find the perambulators, randomvectors and
 *  VdaggerV-operators necessary to build a Quarkline
 *
 *  @param[in]  operator_id     Identifies physical quantum numbers entering
 *  				Quarkline
 *  @param[in]  quantum_numbers Physical quantum field operators for all
 *  				correlators with Dirac structure factored out
 *  @param[in]  vdv_indices     Indices identifying VdaggerV operators
 *  @param[out] Ql_lookup       Lookuptable containing unique combinations of
 *                              peram-, vdv-, and ric-indices needed to built Q1
 *  @param[out] Q1_lookup_ids	List of indices refering to @em Q1. Entries'
 *                              outer vector corresponds to entries of
 *                              @em quantum_numbers, the inner is specified by
 *                              @em operator_id
 *
 *  The Quarkline with one quark is given by
 *
 *    Q1 =  rvdaggerv * gamma * peram
 *
 *  The Quarkline with two quarks is given by
 *
 *    Q2 = @f$ \gamma_5 @f$ peram1@f$ ^\dagger \gamma_5 @f$ * vdaggerv *
 *          gamma * peram2
 *
 *  This function creates a Lookup table of all unique index combination the
 *  Quarkline is needed for and a list referring to the lookup table which gives
 *  relates the entries of the lookup table to each correlator specified in the
 *  infile
 */
static void build_Quarkline_lookup(
    const size_t operator_id,
    const std::vector<std::vector<QuantumNumbers>> &quantum_numbers,
    std::vector<std::vector<std::pair<size_t, bool>>> const &vdv_indices,
    std::vector<std::pair<size_t, size_t>> const &rnd_vec_ids,
    std::vector<DilutedFactorIndex> &Ql_lookup,
    std::vector<std::vector<size_t>> &Ql_lookup_ids) {
  for (size_t row = 0; row < quantum_numbers.size(); row++) {
    const auto qn = quantum_numbers[row][operator_id];

    const auto id_vdaggerv = vdv_indices[row][operator_id].first;
    const auto need_vdaggerv_daggering = vdv_indices[row][operator_id].second;

    // If Ql_lookup already contains the particular row and physical content,
    // just set the index to the existing QuarklineIndices, otherwise generate
    // it and set the index to the new one.
    DilutedFactorIndex const candidate{
        id_vdaggerv, need_vdaggerv_daggering, qn.gamma, rnd_vec_ids};
    auto it = std::find(Ql_lookup.begin(), Ql_lookup.end(), candidate);

    if (it != Ql_lookup.end()) {
      Ql_lookup_ids[row][operator_id] = it - Ql_lookup.begin();
    } else {
      Ql_lookup_ids[row][operator_id] = Ql_lookup.size();
      Ql_lookup.emplace_back(candidate);
    }
  }
}

static void build_Quarkline_lookup_one_qn(
    size_t const operator_id,
    std::vector<QuantumNumbers> const &quantum_numbers,
    std::vector<std::pair<size_t, bool>> const &vdv_indices,
    std::vector<std::pair<size_t, size_t>> const &rnd_vec_ids,
    std::vector<DilutedFactorIndex> &Ql_lookup,
    std::vector<size_t> &Ql_lookup_ids) {
  const auto qn = quantum_numbers[operator_id];

  const auto id_vdaggerv = vdv_indices[operator_id].first;
  const auto need_vdaggerv_daggering = vdv_indices[operator_id].second;

  // If Ql_lookup already contains the particular row and physical content,
  // just set the index to the existing QuarklineIndices, otherwise generate
  // it and set the index to the new one.
  DilutedFactorIndex const candidate{
      id_vdaggerv, need_vdaggerv_daggering, qn.gamma, rnd_vec_ids};
  auto it = std::find(Ql_lookup.begin(), Ql_lookup.end(), candidate);

  if (it != Ql_lookup.end()) {
    Ql_lookup_ids[operator_id] = it - Ql_lookup.begin();
  } else {
    Ql_lookup_ids[operator_id] = Ql_lookup.size();
    Ql_lookup.emplace_back(candidate);
  }
}

static size_t build_trQ1_lookup(std::vector<size_t> const ql_ids,
                                std::vector<DiagramIndex> &trQ1_lookup) {
  DiagramIndex candidate(trQ1_lookup.size(), "", ql_ids);
  auto it = std::find(trQ1_lookup.begin(), trQ1_lookup.end(), candidate);
  if (it == trQ1_lookup.end()) {
    trQ1_lookup.push_back(candidate);
    return trQ1_lookup.back().id;
  } else
    return (it - trQ1_lookup.begin());
}

/*! @BUG If push_back moves the vector somewhere else, it-begin() might not
 *       give the correct id.
 */
static size_t build_corr0_lookup(std::vector<size_t> const ql_ids,
                                 std::vector<DiagramIndex> &trQ1Q1_lookup) {
  DiagramIndex candidate(trQ1Q1_lookup.size(), "", ql_ids);
  auto it = std::find(trQ1Q1_lookup.begin(), trQ1Q1_lookup.end(), candidate);
  if (it == trQ1Q1_lookup.end()) {
    trQ1Q1_lookup.push_back(candidate);
    return trQ1Q1_lookup.back().id;
  } else
    return (it - trQ1Q1_lookup.begin());
}

static size_t build_corrC_lookup(std::vector<size_t> const ql_ids,
                                 std::vector<DiagramIndex> &trQ0Q2_lookup) {
  DiagramIndex candidate(trQ0Q2_lookup.size(), "", ql_ids);
  auto it = std::find(trQ0Q2_lookup.begin(), trQ0Q2_lookup.end(), candidate);
  if (it == trQ0Q2_lookup.end()) {
    trQ0Q2_lookup.push_back(candidate);
    return trQ0Q2_lookup.back().id;
  } else
    return (it - trQ0Q2_lookup.begin());
}
/******************************************************************************/

/******************************************************************************/
/*! Create lookuptable where to find the quarklines to build C20.
 *
 *  @param[in]  quarks            Quarks as read from the infile and processed
 *                                into quark struct
 *  @param[in]  quark_numbers     List which quarks are specified in the infile
 *  @param[in]  start_config      Number of first gauge configuration
 *  @param[in]  path_output       Output path from the infile.
 *  @param[in]  overwrite {yes,no} : deprecated
 *  @param[in]  quantum_numbers   A list of all physical quantum numbers
 *                                quantum field operators for all correlators
 *                                with Dirac structure factored out that are
 *                                possible for @em correlator
 *  @param[in]  vdv_indices       Indices identifying VdaggerV operators
 *  @param[out] Q1_lookup         Lookuptable containing unique combinations of
 *                                peram-, vdv-, and ric-indices needed to built
 *                                Q1
 *  @param[out] trQ1_lookup       Lookuptable containign unique combinations of
 *                                parts tr(Q1Q1).
 *                                Also known as trQ1Q1
 *  @param[out] c_look            Lookup table for C20
 *
 */
static void build_C1_lookup(
    std::vector<quark> const &quarks,
    std::vector<int> const &quark_numbers,
    int const start_config,
    const std::string &path_output,
    const std::string &overwrite,
    std::vector<std::vector<QuantumNumbers>> const &quantum_numbers,
    std::vector<std::vector<std::pair<size_t, bool>>> const &vdv_indices,
    std::vector<DilutedFactorIndex> &Q1_lookup,
    std::vector<DiagramIndex> &trQ1_lookup,
    std::vector<DiagramIndex> &c_look) {
  std::vector<size_t> Q1_indices(std::vector<size_t>(quantum_numbers[0].size()));
  std::vector<std::pair<size_t, size_t>> ric_ids;

  // Build the correlator and dataset names for hdf5 output files
  std::vector<std::string> quark_types;
  for (const auto &id : quark_numbers)
    quark_types.emplace_back(quarks[id].type);

  for (size_t d = 0; d < quantum_numbers.size(); ++d) {
    /*! Here the outgoing and incoming quarkline are identical. Thus ric_ids
     *  must get true
     */
    ric_ids = create_rnd_vec_id(quarks, quark_numbers[0], quark_numbers[0], true);
    build_Quarkline_lookup_one_qn(
        0, quantum_numbers[d], vdv_indices[d], ric_ids, Q1_lookup, Q1_indices);

    /*! @todo create hdf5_dataset name for trQ1Q1. Must restrict quantum
     *  numbers to 0
     */
    auto id1 = build_trQ1_lookup({Q1_indices[0]}, trQ1_lookup);

    std::string hdf5_dataset_name = build_hdf5_dataset_name(
        "C1", start_config, path_output, overwrite, quark_types, quantum_numbers[d]);

    DiagramIndex candidate{c_look.size(), hdf5_dataset_name, {id1}};

    /*! XXX Better with std::set */
    auto it = std::find(c_look.begin(), c_look.end(), candidate);

    if (it == c_look.end()) {
      c_look.push_back(candidate);
    }
  }
}

/******************************************************************************/
/*! Create lookuptable where to find the quarklines to build C20.
 *
 *  @param[in]  quarks            Quarks as read from the infile and processed
 *                                into quark struct
 *  @param[in]  quark_numbers     List which quarks are specified in the infile
 *  @param[in]  start_config      Number of first gauge configuration
 *  @param[in]  path_output       Output path from the infile.
 *  @param[in]  overwrite {yes,no} : deprecated
 *  @param[in]  quantum_numbers   A list of all physical quantum numbers
 *                                quantum field operators for all correlators
 *                                with Dirac structure factored out that are
 *                                possible for @em correlator
 *  @param[in]  vdv_indices       Indices identifying VdaggerV operators
 *  @param[out] Q1_lookup         Lookuptable containing unique combinations of
 *                                peram-, vdv-, and ric-indices needed to built
 *                                Q1
 *  @param[out] trQ1_lookup       Lookuptable containign unique combinations of
 *                                parts tr(Q1Q1).
 *                                Also known as trQ1Q1
 *  @param[out] c_look            Lookup table for C20
 *
 *  @bug I am fairly certain that the quarks are mixed up. It is
 *        also wrong in init_lookup_tables() (MW 27.3.17)
 */
static void build_C20V_lookup(
    std::vector<quark> const &quarks,
    std::vector<int> const &quark_numbers,
    int const start_config,
    const std::string &path_output,
    const std::string &overwrite,
    std::vector<std::vector<QuantumNumbers>> const &quantum_numbers,
    std::vector<std::vector<std::pair<size_t, bool>>> const &vdv_indices,
    std::vector<DilutedFactorIndex> &Q1_lookup,
    std::vector<DiagramIndex> &trQ1_lookup,
    std::vector<DiagramIndex> &c_look) {
  std::vector<size_t> ql_ids(2);
  std::vector<std::pair<size_t, size_t>> ric_ids;

  // Build the correlator and dataset names for hdf5 output files
  std::vector<std::string> quark_types;
  for (const auto &id : quark_numbers)
    quark_types.emplace_back(quarks[id].type);

  for (size_t d = 0; d < quantum_numbers.size(); ++d) {
    ric_ids = create_rnd_vec_id(quarks, quark_numbers[0], quark_numbers[0], true);
    build_Quarkline_lookup_one_qn(
        0, quantum_numbers[d], vdv_indices[d], ric_ids, Q1_lookup, ql_ids);
    ric_ids = create_rnd_vec_id(quarks, quark_numbers[1], quark_numbers[1], true);
    build_Quarkline_lookup_one_qn(
        1, quantum_numbers[d], vdv_indices[d], ric_ids, Q1_lookup, ql_ids);

    /*! @todo create hdf5_dataset name for trQ1Q1. Must restrict quantum
     *  numbers to 0,1
     */
    auto id1 = build_trQ1_lookup({ql_ids[0]}, trQ1_lookup);
    auto id2 = build_trQ1_lookup({ql_ids[1]}, trQ1_lookup);

    std::string hdf5_dataset_name = build_hdf5_dataset_name(
        "C20V", start_config, path_output, overwrite, quark_types, quantum_numbers[d]);

    DiagramIndex candidate{c_look.size(), hdf5_dataset_name, {id1, id2}};

    /*! XXX Better with std::set */
    auto it = std::find(c_look.begin(), c_look.end(), candidate);

    if (it == c_look.end()) {
      c_look.push_back(candidate);
    }
  }
}

/******************************************************************************/
/*! Create lookuptable where to find the quarklines to build C2c.
 *
 *  @param[in]  quarks            Quarks as read from the infile and processed
 *                                into quark struct
 *  @param[in]  quark_numbers     List which quarks are specified in the infile
 *  @param[in]  start_config      Number of first gauge configuration
 *  @param[in]  path_output       Output path from the infile.
 *  @param[in]  overwrite {yes,no} : deprecated
 *  @param[in]  quantum_numbers   A list of all physical quantum numbers
 *                                quantum field operators for all correlators
 *                                with Dirac structure factored out that are
 *                                possible for @em correlator
 *  @param[in]  vdv_indices       Indices identifying VdaggerV operators
 *  @param[out] Q0_lookup         Lookuptable containing unique combinations of
 *                                peram-, vdv-, and ric-indices needed to built
 *                                Q0
 *  @param[out] Q2V_lookup        Lookuptable containing unique combinations of
 *                                peram-, vdv-, and ric-indices needed to built
 *                                Q2V
 *  @param[out] trQ0Q2_lookup     Lookuptable containign unique combinations of
 *                                parts tr(Q0Q2).
 *                                Also known as trQ0Q2
 *  @param[out] c_look            Lookup table for C2c
 *
 *  @bug I am fairly certain that the quarks are mixed up. It is
 *        also wrong in init_lookup_tables() (MW 27.3.17)
 */
static void build_C2c_lookup(
    std::vector<quark> const &quarks,
    std::vector<int> const &quark_numbers,
    int const start_config,
    const std::string &path_output,
    const std::string &overwrite,
    std::vector<std::vector<QuantumNumbers>> const &quantum_numbers,
    std::vector<std::vector<std::pair<size_t, bool>>> const &vdv_indices,
    std::vector<DilutedFactorIndex> &Q0_lookup,
    std::vector<DilutedFactorIndex> &Q2V_lookup,
    std::vector<DiagramIndex> &trQ0Q2_lookup,
    std::vector<DiagramIndex> &c_look) {
  std::vector<size_t> ql_ids(2);
  std::vector<std::pair<size_t, size_t>> ric_ids;

  // Build the correlator and dataset names for hdf5 output files
  std::vector<std::string> quark_types;
  for (const auto &id : quark_numbers)
    quark_types.emplace_back(quarks[id].type);

  for (size_t d = 0; d < quantum_numbers.size(); ++d) {
    ric_ids = create_rnd_vec_id(quarks, quark_numbers[1], quark_numbers[0], false);
    build_Quarkline_lookup_one_qn(
        0, quantum_numbers[d], vdv_indices[d], ric_ids, Q2V_lookup, ql_ids);
    ric_ids = create_rnd_vec_id(quarks, quark_numbers[0], quark_numbers[1], false);
    build_Quarkline_lookup_one_qn(
        1, quantum_numbers[d], vdv_indices[d], ric_ids, Q0_lookup, ql_ids);

    /*! @todo create hdf5_dataset name for trQ1Q1. Must restrict quantum
     *  numbers to 0,1
     */
    auto id1 = build_corrC_lookup({ql_ids[0], ql_ids[1]}, trQ0Q2_lookup);

    std::string hdf5_dataset_name = build_hdf5_dataset_name(
        "C2+", start_config, path_output, overwrite, quark_types, quantum_numbers[d]);

    DiagramIndex candidate{c_look.size(), hdf5_dataset_name, {id1}};

    /*! XXX Better with std::set */
    auto it = std::find(c_look.begin(), c_look.end(), candidate);

    if (it == c_look.end()) {
      c_look.push_back(candidate);
    }
  }
}
/******************************************************************************/
/*! Create lookuptable where to find the quarklines to build C20.
 *
 *  @param[in]  quarks            Quarks as read from the infile and processed
 *                                into quark struct
 *  @param[in]  quark_numbers     List which quarks are specified in the infile
 *  @param[in]  start_config      Number of first gauge configuration
 *  @param[in]  path_output       Output path from the infile.
 *  @param[in]  overwrite {yes,no} : deprecated
 *  @param[in]  quantum_numbers   A list of all physical quantum numbers
 *                                quantum field operators for all correlators
 *                                with Dirac structure factored out that are
 *                                possible for @em correlator
 *  @param[in]  vdv_indices       Indices identifying VdaggerV operators
 *  @param[out] Q1_lookup         Lookuptable containing unique combinations of
 *                                peram-, vdv-, and ric-indices needed to built
 *                                Q1
 *  @param[out] trQ1Q1_lookup     Lookuptable containign unique combinations of
 *                                parts tr(Q1Q1).
 *                                Also known as trQ1Q1
 *  @param[out] c_look            Lookup table for C20
 *
 *  @bug I am fairly certain that the quarks are mixed up. It is
 *        also wrong in init_lookup_tables() (MW 27.3.17)
 */
static void build_C20_lookup(
    std::vector<quark> const &quarks,
    std::vector<int> const &quark_numbers,
    int const start_config,
    const std::string &path_output,
    const std::string &overwrite,
    std::vector<std::vector<QuantumNumbers>> const &quantum_numbers,
    std::vector<std::vector<std::pair<size_t, bool>>> const &vdv_indices,
    std::vector<DilutedFactorIndex> &Q1_lookup,
    std::vector<DiagramIndex> &trQ1Q1_lookup,
    std::vector<DiagramIndex> &c_look) {
  std::vector<size_t> ql_ids(2);
  std::vector<std::pair<size_t, size_t>> ric_ids;

  // Build the correlator and dataset names for hdf5 output files
  std::vector<std::string> quark_types;
  for (const auto &id : quark_numbers)
    quark_types.emplace_back(quarks[id].type);

  for (size_t d = 0; d < quantum_numbers.size(); ++d) {
    ric_ids = create_rnd_vec_id(quarks, quark_numbers[1], quark_numbers[0], false);
    build_Quarkline_lookup_one_qn(
        0, quantum_numbers[d], vdv_indices[d], ric_ids, Q1_lookup, ql_ids);
    ric_ids = create_rnd_vec_id(quarks, quark_numbers[0], quark_numbers[1], false);
    build_Quarkline_lookup_one_qn(
        1, quantum_numbers[d], vdv_indices[d], ric_ids, Q1_lookup, ql_ids);

    /*! @todo create hdf5_dataset name for trQ1Q1. Must restrict quantum
     *  numbers to 0,1
     */
    auto id1 = build_corr0_lookup({ql_ids[0], ql_ids[1]}, trQ1Q1_lookup);

    std::string hdf5_dataset_name = build_hdf5_dataset_name(
        "C20", start_config, path_output, overwrite, quark_types, quantum_numbers[d]);

    DiagramIndex candidate{c_look.size(), hdf5_dataset_name, {id1}};

    /*! XXX Better with std::set */
    auto it = std::find(c_look.begin(), c_look.end(), candidate);

    if (it == c_look.end()) {
      c_look.push_back(candidate);
    }
  }
}

/******************************************************************************/
/*! Create lookuptable where to find the quarklines to build C30V.
 *
 *  @param[in]  quarks            Quarks as read from the infile and processed
 *                                into quark struct
 *  @param[in]  quark_numbers     List which quarks are specified in the infile
 *  @param[in]  start_config      Number of first gauge configuration
 *  @param[in]  path_output       Output path from the infile.
 *  @param[in]  overwrite {yes,no} : deprecated
 *  @param[in]  quantum_numbers   A list of all physical quantum numbers
 *                                quantum field operators for all correlators
 *                                with Dirac structure factored out that are
 *                                possible for @em correlator
 *  @param[in]  vdv_indices       Indices identifying VdaggerV operators
 *  @param[out] Q1_lookup         Lookuptable containing unique combinations of
 *                                peram-, vdv-, and ric-indices needed to built
 *                                Q1
 *  @param[out] trQ1_lookup       Lookuptable containign unique combinations of
 *                                parts tr(Q1).
 *  @param[out] trQ1Q1_lookup     Lookuptable containign unique combinations of
 *                                parts tr(Q1Q1).
 *                                Also known as trQ1Q1
 *  @param[out] c_look            Lookup table for C30V
 *
 *  @bug I am fairly certain that the quarks are mixed up. It is
 *        also wrong in init_lookup_tables() (MW 27.3.17)
 */
static void build_C30V_lookup(
    std::vector<quark> const &quarks,
    std::vector<int> const &quark_numbers,
    int start_config,
    const std::string &path_output,
    const std::string &overwrite,
    std::vector<std::vector<QuantumNumbers>> const &quantum_numbers,
    std::vector<std::vector<std::pair<size_t, bool>>> const &vdv_indices,
    std::vector<DilutedFactorIndex> &Q1_lookup,
    std::vector<DiagramIndex> &trQ1_lookup,
    std::vector<DiagramIndex> &trQ1Q1_lookup,
    std::vector<DiagramIndex> &c_look) {
  std::vector<size_t> ql_ids(3);
  std::vector<std::pair<size_t, size_t>> ric_ids;

  // Build the correlator and dataset names for hdf5 output files
  std::vector<std::string> quark_types;
  for (const auto &id : quark_numbers) {
    quark_types.emplace_back(quarks[id].type);
  }

  for (size_t d = 0; d < quantum_numbers.size(); ++d) {
    ric_ids = create_rnd_vec_id(quarks, quark_numbers[1], quark_numbers[0], false);
    build_Quarkline_lookup_one_qn(
        0, quantum_numbers[d], vdv_indices[d], ric_ids, Q1_lookup, ql_ids);

    ric_ids = create_rnd_vec_id(quarks, quark_numbers[0], quark_numbers[1], false);
    build_Quarkline_lookup_one_qn(
        1, quantum_numbers[d], vdv_indices[d], ric_ids, Q1_lookup, ql_ids);

    ric_ids = create_rnd_vec_id(quarks, quark_numbers[2], quark_numbers[2], true);
    build_Quarkline_lookup_one_qn(
        2, quantum_numbers[d], vdv_indices[d], ric_ids, Q1_lookup, ql_ids);

    auto const id0 = build_corr0_lookup({ql_ids[0], ql_ids[1]}, trQ1Q1_lookup);
    auto const id1 = build_trQ1_lookup({ql_ids[2]}, trQ1_lookup);

    std::string const hdf5_dataset_name = build_hdf5_dataset_name(
        "C30V", start_config, path_output, overwrite, quark_types, quantum_numbers[d]);

    DiagramIndex const candidate{
        c_look.size(), hdf5_dataset_name, {id0, id1}, std::vector<int>{}};

    /*! XXX Better with std::set */
    auto const it = std::find(c_look.begin(), c_look.end(), candidate);

    if (it == c_look.end()) {
      c_look.push_back(candidate);
    }
  }
}

/******************************************************************************/
/*! Create lookuptable where to find the quarklines to build C3c.
 *
 *  @param[in]  quarks            Quarks as read from the infile and processed
 *                                into quark struct
 *  @param[in]  quark_numbers     List which quarks are specified in the infile
 *  @param[in]  start_config      Number of first gauge configuration
 *  @param[in]  path_output       Output path from the infile.
 *  @param[in]  overwrite {yes,no} : deprecated
 *  @param[in]  quantum_numbers   A list of all physical quantum numbers
 *                                quantum field operators for all correlators
 *                                with Dirac structure factored out that are
 *                                possible for @em correlator
 *  @param[in]  vdv_indices       Indices identifying VdaggerV operators
 *  @param[out] Q0_lookup         Lookuptable containing unique combinations of
 *                                peram-, vdv-, and ric-indices needed to built
 *                                Q0
 *  @param[out] Q1_lookup         Lookuptable containing unique combinations of
 *                                peram-, vdv-, and ric-indices needed to built
 *                                Q1
 *  @param[out] Q2L_lookup        Lookuptable containing unique combinations of
 *                                peram-, vdv-, and ric-indices needed to built
 *                                Q2L
 *  @param[out] c_look            Lookup table for C3c
 */
static void build_C3c_lookup(
    std::vector<quark> const &quarks,
    std::vector<int> const &quark_numbers,
    int start_config,
    const std::string &path_output,
    const std::string &overwrite,
    std::vector<std::vector<QuantumNumbers>> const &quantum_numbers,
    std::vector<std::vector<std::pair<size_t, bool>>> const &vdv_indices,
    std::vector<DilutedFactorIndex> &Q0_lookup,
    std::vector<DilutedFactorIndex> &Q1_lookup,
    std::vector<DilutedFactorIndex> &Q2L_lookup,
    std::vector<DiagramIndex> &c_look) {
  std::vector<size_t> ql_ids(3);
  std::vector<std::pair<size_t, size_t>> ric_ids;

  // Build the correlator and dataset names for hdf5 output files
  std::vector<std::string> quark_types;
  for (const auto &id : quark_numbers)
    quark_types.emplace_back(quarks[id].type);

  for (size_t d = 0; d < quantum_numbers.size(); ++d) {
    ric_ids = create_rnd_vec_id(quarks, quark_numbers[2], quark_numbers[0], false);
    build_Quarkline_lookup_one_qn(
        0, quantum_numbers[d], vdv_indices[d], ric_ids, Q2L_lookup, ql_ids);
    ric_ids = create_rnd_vec_id(quarks, quark_numbers[0], quark_numbers[1], false);
    build_Quarkline_lookup_one_qn(
        1, quantum_numbers[d], vdv_indices[d], ric_ids, Q1_lookup, ql_ids);
    ric_ids = create_rnd_vec_id(quarks, quark_numbers[1], quark_numbers[2], false);
    build_Quarkline_lookup_one_qn(
        2, quantum_numbers[d], vdv_indices[d], ric_ids, Q0_lookup, ql_ids);

    std::string hdf5_dataset_name = build_hdf5_dataset_name(
        "C3+", start_config, path_output, overwrite, quark_types, quantum_numbers[d]);

    DiagramIndex candidate{
        c_look.size(), hdf5_dataset_name, ql_ids, std::vector<int>({})};

    /*! XXX Better with std::set */
    auto it = std::find(c_look.begin(), c_look.end(), candidate);

    if (it == c_look.end()) {
      c_look.push_back(candidate);
    }
  }
}

/******************************************************************************/
/*! Create lookuptable where to find the quarklines to build C30.
 *
 *  @param[in]  quarks            Quarks as read from the infile and processed
 *                                into quark struct
 *  @param[in]  quark_numbers     List which quarks are specified in the infile
 *  @param[in]  start_config      Number of first gauge configuration
 *  @param[in]  path_output       Output path from the infile.
 *  @param[in]  overwrite {yes,no} : deprecated
 *  @param[in]  quantum_numbers   A list of all physical quantum numbers
 *                                quantum field operators for all correlators
 *                                with Dirac structure factored out that are
 *                                possible for @em correlator
 *  @param[in]  vdv_indices       Indices identifying VdaggerV operators
 *  @param[out] Q1_lookup         Lookuptable containing unique combinations of
 *                                peram-, vdv-, and ric-indices needed to built
 *                                Q1
 *  @param[out] c_look            Lookup table for C30
 *
 *  @bug I am fairly certain that the quarks are mixed up. It is
 *        also wrong in init_lookup_tables() (MW 27.3.17)
 */
static void build_C30_lookup(
    std::vector<quark> const &quarks,
    std::vector<int> const &quark_numbers,
    int start_config,
    const std::string &path_output,
    const std::string &overwrite,
    std::vector<std::vector<QuantumNumbers>> const &quantum_numbers,
    std::vector<std::vector<std::pair<size_t, bool>>> const &vdv_indices,
    std::vector<DilutedFactorIndex> &Q1_lookup,
    std::vector<DiagramIndex> &c_look) {
  std::vector<size_t> ql_ids(3);
  std::vector<std::pair<size_t, size_t>> ric_ids;

  // Build the correlator and dataset names for hdf5 output files
  std::vector<std::string> quark_types;
  for (const auto &id : quark_numbers)
    quark_types.emplace_back(quarks[id].type);

  for (size_t d = 0; d < quantum_numbers.size(); ++d) {
    ric_ids = create_rnd_vec_id(quarks, quark_numbers[2], quark_numbers[0], false);
    build_Quarkline_lookup_one_qn(
        0, quantum_numbers[d], vdv_indices[d], ric_ids, Q1_lookup, ql_ids);
    ric_ids = create_rnd_vec_id(quarks, quark_numbers[0], quark_numbers[1], false);
    build_Quarkline_lookup_one_qn(
        1, quantum_numbers[d], vdv_indices[d], ric_ids, Q1_lookup, ql_ids);
    ric_ids = create_rnd_vec_id(quarks, quark_numbers[1], quark_numbers[2], false);
    build_Quarkline_lookup_one_qn(
        2, quantum_numbers[d], vdv_indices[d], ric_ids, Q1_lookup, ql_ids);

    std::string hdf5_dataset_name = build_hdf5_dataset_name(
        "C30", start_config, path_output, overwrite, quark_types, quantum_numbers[d]);

    DiagramIndex candidate{
        c_look.size(), hdf5_dataset_name, ql_ids, std::vector<int>({})};

    /*! XXX Better with std::set */
    auto it = std::find(c_look.begin(), c_look.end(), candidate);

    if (it == c_look.end()) {
      c_look.push_back(candidate);
    }
  }
}

/******************************************************************************/
/*! Create lookuptable where to find the quarklines to build C4cD. Also sets
 *  trQ0Q2.
 *
 *  @param[in]  quarks            Quarks as read from the infile and processed
 *                                into quark struct
 *  @param[in]  quark_numbers     List which quarks are specified in the infile
 *  @param[in]  start_config      Number of first gauge configuration
 *  @param[in]  path_output       Output path from the infile.
 *  @param[in]  overwrite {yes,no} : deprecated
 *  @param[in]  quantum_numbers   A list of all physical quantum numbers
 *                                quantum field operators for all correlators
 *                                with Dirac structure factored out that are
 *                                possible for @em correlator
 *  @param[in]  vdv_indices       Indices identifying VdaggerV operators
 *  @param[out] Q0_lookup         Lookuptable containing unique combinations of
 *                                peram-, vdv-, and ric-indices needed to built
 *                                Q0
 *  @param[out] Q2V_lookup        Lookuptable containing unique combinations of
 *                                peram-, vdv-, and ric-indices needed to built
 *                                Q2V
 *  @param[out] trQ0Q2_lookup     Lookuptable containign unique combinations of
 *                                parts tr(Q0Q2).
 *                                Also known as trQ0Q2
 *  @param[out] c_look            Lookup table for C4cD
 *
 *  C4cD like C4cC contains C2c. To reuse C2c, they all contain indices
 *  of trQ0Q2 which in turn contains the indices for rVdaggerVr and Q2.
 *
 *  @bug I am fairly certain that the quarks are mixed up. It is
 *        also wrong in init_lookup_tables() (MW 27.3.17)
 */
static void build_C4cD_lookup(
    std::vector<quark> const &quarks,
    std::vector<int> const &quark_numbers,
    int start_config,
    const std::string &path_output,
    const std::string &overwrite,
    std::vector<std::vector<QuantumNumbers>> const &quantum_numbers,
    std::vector<std::vector<std::pair<size_t, bool>>> const &vdv_indices,
    std::vector<DilutedFactorIndex> &Q0_lookup,
    std::vector<DilutedFactorIndex> &Q2V_lookup,
    std::vector<DiagramIndex> &trQ0Q2_lookup,
    std::vector<DiagramIndex> &c_look) {
  std::vector<size_t> ql_ids(4);
  std::vector<std::pair<size_t, size_t>> ric_ids;

  // Build the correlator and dataset names for hdf5 output files
  std::vector<std::string> quark_types;
  for (const auto &id : quark_numbers)
    quark_types.emplace_back(quarks[id].type);

  for (size_t d = 0; d < quantum_numbers.size(); ++d) {
    ric_ids = create_rnd_vec_id(quarks, quark_numbers[1], quark_numbers[0], false);
    build_Quarkline_lookup_one_qn(
        0, quantum_numbers[d], vdv_indices[d], ric_ids, Q2V_lookup, ql_ids);
    ric_ids = create_rnd_vec_id(quarks, quark_numbers[0], quark_numbers[1], false);
    build_Quarkline_lookup_one_qn(
        1, quantum_numbers[d], vdv_indices[d], ric_ids, Q0_lookup, ql_ids);
    ric_ids = create_rnd_vec_id(quarks, quark_numbers[3], quark_numbers[2], false);
    build_Quarkline_lookup_one_qn(
        2, quantum_numbers[d], vdv_indices[d], ric_ids, Q2V_lookup, ql_ids);
    ric_ids = create_rnd_vec_id(quarks, quark_numbers[2], quark_numbers[3], false);
    build_Quarkline_lookup_one_qn(
        3, quantum_numbers[d], vdv_indices[d], ric_ids, Q0_lookup, ql_ids);

    /*! @todo create hdf5_dataset name for trQ1Q1. Must restrict quantum
     *  numbers to 0,1 / 2,3
     */
    auto id1 = build_corrC_lookup({ql_ids[0], ql_ids[1]}, trQ0Q2_lookup);
    auto id2 = build_corrC_lookup({ql_ids[2], ql_ids[3]}, trQ0Q2_lookup);

    std::string hdf5_dataset_name = build_hdf5_dataset_name(
        "C4+D", start_config, path_output, overwrite, quark_types, quantum_numbers[d]);

    DiagramIndex candidate{c_look.size(),
                           hdf5_dataset_name,
                           std::vector<size_t>({id1, id2}),
                           std::vector<int>({})};

    /*! XXX Better with std::set */
    auto it = std::find(c_look.begin(), c_look.end(), candidate);

    if (it == c_look.end()) {
      c_look.push_back(candidate);
    }
  }
}

/******************************************************************************/
/*! Create lookuptable where to find the quarklines to build C40D.
 *
 *  @param[in]  quarks            Quarks as read from the infile and processed
 *                                into quark struct
 *  @param[in]  quark_numbers     List which quarks are specified in the infile
 *  @param[in]  start_config      Number of first gauge configuration
 *  @param[in]  path_output       Output path from the infile.
 *  @param[in]  overwrite {yes,no} : deprecated
 *  @param[in]  quantum_numbers   A list of all physical quantum numbers
 *                                quantum field operators for all correlators
 *                                with Dirac structure factored out that are
 *                                possible for @em correlator
 *  @param[in]  vdv_indices       Indices identifying VdaggerV operators
 *  @param[out] Q1_lookup         Lookuptable containing unique combinations of
 *                                peram-, vdv-, and ric-indices needed to built
 *                                Q1
 *  @param[out] trQ1Q1_lookup     Lookuptable containign unique combinations of
 *                                parts tr(Q1Q1).
 *                                Also known as trQ1Q1
 *  @param[out] c_look            Lookup table for C40D
 *
 *  @bug I am fairly certain that the quarks are mixed up. It is
 *        also wrong in init_lookup_tables() (MW 27.3.17)
 */
static void build_C40D_lookup(
    std::vector<quark> const &quarks,
    std::vector<int> const &quark_numbers,
    int const start_config,
    const std::string &path_output,
    const std::string &overwrite,
    std::vector<std::vector<QuantumNumbers>> const &quantum_numbers,
    std::vector<std::vector<std::pair<size_t, bool>>> const &vdv_indices,
    std::vector<DilutedFactorIndex> &Q1_lookup,
    std::vector<DiagramIndex> &trQ1Q1_lookup,
    std::vector<DiagramIndex> &c_look) {
  std::vector<size_t> ql_ids(4);
  std::vector<std::pair<size_t, size_t>> ric_ids;

  // Build the correlator and dataset names for hdf5 output files
  std::vector<std::string> quark_types;
  for (const auto &id : quark_numbers)
    quark_types.emplace_back(quarks[id].type);

  for (size_t d = 0; d < quantum_numbers.size(); ++d) {
    ric_ids = create_rnd_vec_id(quarks, quark_numbers[1], quark_numbers[0], false);
    build_Quarkline_lookup_one_qn(
        0, quantum_numbers[d], vdv_indices[d], ric_ids, Q1_lookup, ql_ids);
    ric_ids = create_rnd_vec_id(quarks, quark_numbers[0], quark_numbers[1], false);
    build_Quarkline_lookup_one_qn(
        1, quantum_numbers[d], vdv_indices[d], ric_ids, Q1_lookup, ql_ids);
    ric_ids = create_rnd_vec_id(quarks, quark_numbers[3], quark_numbers[2], false);
    build_Quarkline_lookup_one_qn(
        2, quantum_numbers[d], vdv_indices[d], ric_ids, Q1_lookup, ql_ids);
    ric_ids = create_rnd_vec_id(quarks, quark_numbers[2], quark_numbers[3], false);
    build_Quarkline_lookup_one_qn(
        3, quantum_numbers[d], vdv_indices[d], ric_ids, Q1_lookup, ql_ids);

    /*! @todo create hdf5_dataset name for trQ1Q1. Must restrict quantum
     *  numbers to 0,1 / 2,3
     */
    auto id1 = build_corr0_lookup({ql_ids[0], ql_ids[1]}, trQ1Q1_lookup);
    auto id2 = build_corr0_lookup({ql_ids[2], ql_ids[3]}, trQ1Q1_lookup);

    std::string hdf5_dataset_name = build_hdf5_dataset_name(
        "C40D", start_config, path_output, overwrite, quark_types, quantum_numbers[d]);

    DiagramIndex candidate{c_look.size(),
                           hdf5_dataset_name,
                           std::vector<size_t>({id1, id2}),
                           std::vector<int>({})};

    /*! XXX Better with std::set */
    auto it = std::find(c_look.begin(), c_look.end(), candidate);

    if (it == c_look.end()) {
      c_look.push_back(candidate);
    }
  }
}

/******************************************************************************/
/*! Create lookuptable where to find the quarklines to build C4cV. Also sets
 *  trQ0Q2.
 *
 *  @param[in]  quarks            Quarks as read from the infile and processed
 *                                into quark struct
 *  @param[in]  quark_numbers     List which quarks are specified in the infile
 *  @param[in]  start_config      Number of first gauge configuration
 *  @param[in]  path_output       Output path from the infile.
 *  @param[in]  overwrite {yes,no} : deprecated
 *  @param[in]  quantum_numbers   A list of all physical quantum numbers
 *                                quantum field operators for all correlators
 *                                with Dirac structure factored out that are
 *                                possible for @em correlator
 *  @param[in]  vdv_indices       Indices identifying VdaggerV operators
 *  @param[out] Q0_lookup         Lookuptable containing unique combinations of
 *                                peram-, vdv-, and ric-indices needed to built
 *                                Q0
 *  @param[out] Q2V_lookup        Lookuptable containing unique combinations of
 *                                peram-, vdv-, and ric-indices needed to built
 *                                Q2V
 *  @param[out] trQ0Q2_lookup     Lookuptable containign unique combinations of
 *                                parts tr(Q0Q2).
 *                                Also known as trQ0Q2
 *  @param[out] c_look            Lookup table for C4cV
 *
 *  C4cV like C4cC contains C2c. To reuse C2c, they all contain indices
 *  of trQ0Q2 which in turn contains the indices for rVdaggerVr and Q2.
 *
 *  @bug I am fairly certain that the quarks are mixed up. It is
 *        also wrong in init_lookup_tables() (MW 27.3.17)
 */
static void build_C4cV_lookup(
    std::vector<quark> const &quarks,
    std::vector<int> const &quark_numbers,
    int start_config,
    const std::string &path_output,
    const std::string &overwrite,
    std::vector<std::vector<QuantumNumbers>> const &quantum_numbers,
    std::vector<std::vector<std::pair<size_t, bool>>> const &vdv_indices,
    std::vector<DilutedFactorIndex> &Q0_lookup,
    std::vector<DilutedFactorIndex> &Q2V_lookup,
    std::vector<DiagramIndex> &trQ0Q2_lookup,
    std::vector<DiagramIndex> &c_look) {
  std::vector<size_t> ql_ids(4);
  std::vector<std::pair<size_t, size_t>> ric_ids;

  // Build the correlator and dataset names for hdf5 output files
  std::vector<std::string> quark_types;
  for (const auto &id : quark_numbers)
    quark_types.emplace_back(quarks[id].type);

  for (size_t d = 0; d < quantum_numbers.size(); ++d) {
    ric_ids = create_rnd_vec_id(quarks, quark_numbers[1], quark_numbers[0], false);
    build_Quarkline_lookup_one_qn(
        0, quantum_numbers[d], vdv_indices[d], ric_ids, Q2V_lookup, ql_ids);
    ric_ids = create_rnd_vec_id(quarks, quark_numbers[0], quark_numbers[1], false);
    build_Quarkline_lookup_one_qn(
        1, quantum_numbers[d], vdv_indices[d], ric_ids, Q0_lookup, ql_ids);
    ric_ids = create_rnd_vec_id(quarks, quark_numbers[3], quark_numbers[2], false);
    build_Quarkline_lookup_one_qn(
        2, quantum_numbers[d], vdv_indices[d], ric_ids, Q2V_lookup, ql_ids);
    ric_ids = create_rnd_vec_id(quarks, quark_numbers[2], quark_numbers[3], false);
    build_Quarkline_lookup_one_qn(
        3, quantum_numbers[d], vdv_indices[d], ric_ids, Q0_lookup, ql_ids);

    /*! @todo create hdf5_dataset name for trQ1Q1. Must restrict quantum
     *  numbers to 0,1 / 2,3
     */
    auto id1 = build_corrC_lookup({ql_ids[0], ql_ids[1]}, trQ0Q2_lookup);
    auto id2 = build_corrC_lookup({ql_ids[2], ql_ids[3]}, trQ0Q2_lookup);

    std::string hdf5_dataset_name = build_hdf5_dataset_name(
        "C4+V", start_config, path_output, overwrite, quark_types, quantum_numbers[d]);

    DiagramIndex candidate{c_look.size(),
                           hdf5_dataset_name,
                           std::vector<size_t>({id1, id2}),
                           std::vector<int>({})};

    /*! XXX Better with std::set */
    auto it = std::find(c_look.begin(), c_look.end(), candidate);

    if (it == c_look.end()) {
      c_look.push_back(candidate);
    }
  }
}

/******************************************************************************/
/*! Create lookuptable where to find the quarklines to build C40V.
 *
 *  @param[in]  quarks            Quarks as read from the infile and processed
 *                                into quark struct
 *  @param[in]  quark_numbers     List which quarks are specified in the infile
 *  @param[in]  start_config      Number of first gauge configuration
 *  @param[in]  path_output       Output path from the infile.
 *  @param[in]  overwrite {yes,no} : deprecated
 *  @param[in]  quantum_numbers   A list of all physical quantum numbers
 *                                quantum field operators for all correlators
 *                                with Dirac structure factored out that are
 *                                possible for @em correlator
 *  @param[in]  vdv_indices       Indices identifying VdaggerV operators
 *  @param[out] Q1_lookup         Lookuptable containing unique combinations of
 *                                peram-, vdv-, and ric-indices needed to built
 *                                Q1
 *  @param[out] trQ1Q1_lookup     Lookuptable containign unique combinations of
 *                                parts tr(Q1Q1).
 *                                Also known as trQ1Q1
 *  @param[out] c_look            Lookup table for C40V
 *
 *  @bug I am fairly certain that the quarks are mixed up. It is
 *        also wrong in init_lookup_tables() (MW 27.3.17)
 */
static void build_C40V_lookup(
    std::vector<quark> const &quarks,
    std::vector<int> const &quark_numbers,
    int start_config,
    const std::string &path_output,
    const std::string &overwrite,
    std::vector<std::vector<QuantumNumbers>> const &quantum_numbers,
    std::vector<std::vector<std::pair<size_t, bool>>> const &vdv_indices,
    std::vector<DilutedFactorIndex> &Q1_lookup,
    std::vector<DiagramIndex> &trQ1Q1_lookup,
    std::vector<DiagramIndex> &c_look) {
  std::vector<size_t> ql_ids(4);
  std::vector<std::pair<size_t, size_t>> ric_ids;

  // Build the correlator and dataset names for hdf5 output files
  std::vector<std::string> quark_types;
  for (const auto &id : quark_numbers)
    quark_types.emplace_back(quarks[id].type);

  for (size_t d = 0; d < quantum_numbers.size(); ++d) {
    ric_ids = create_rnd_vec_id(quarks, quark_numbers[1], quark_numbers[0], false);
    build_Quarkline_lookup_one_qn(
        0, quantum_numbers[d], vdv_indices[d], ric_ids, Q1_lookup, ql_ids);
    ric_ids = create_rnd_vec_id(quarks, quark_numbers[0], quark_numbers[1], false);
    build_Quarkline_lookup_one_qn(
        1, quantum_numbers[d], vdv_indices[d], ric_ids, Q1_lookup, ql_ids);
    ric_ids = create_rnd_vec_id(quarks, quark_numbers[3], quark_numbers[2], false);
    build_Quarkline_lookup_one_qn(
        2, quantum_numbers[d], vdv_indices[d], ric_ids, Q1_lookup, ql_ids);
    ric_ids = create_rnd_vec_id(quarks, quark_numbers[2], quark_numbers[3], false);
    build_Quarkline_lookup_one_qn(
        3, quantum_numbers[d], vdv_indices[d], ric_ids, Q1_lookup, ql_ids);

    /*! @todo create hdf5_dataset name for trQ1Q1. Must restrict quantum
     *  numbers to 0,1 / 2,3
     */
    auto id1 = build_corr0_lookup({ql_ids[0], ql_ids[1]}, trQ1Q1_lookup);
    auto id2 = build_corr0_lookup({ql_ids[2], ql_ids[3]}, trQ1Q1_lookup);

    std::string hdf5_dataset_name = build_hdf5_dataset_name(
        "C40V", start_config, path_output, overwrite, quark_types, quantum_numbers[d]);

    DiagramIndex candidate{c_look.size(),
                           hdf5_dataset_name,
                           std::vector<size_t>({id1, id2}),
                           std::vector<int>({})};

    /*! XXX Better with std::set */
    auto it = std::find(c_look.begin(), c_look.end(), candidate);

    if (it == c_look.end()) {
      c_look.push_back(candidate);
    }
  }
}

/******************************************************************************/
/*! Create lookuptable where to find the quarklines to build C4cC.
 *
 *  @param[in]  quarks            Quarks as read from the infile and processed
 *                                into quark struct
 *  @param[in]  quark_numbers     List which quarks are specified in the infile
 *  @param[in]  start_config      Number of first gauge configuration
 *  @param[in]  path_output       Output path from the infile.
 *  @param[in]  overwrite {yes,no} : deprecated
 *  @param[in]  quantum_numbers   A list of all physical quantum numbers
 *                                quantum field operators for all correlators
 *                                with Dirac structure factored out that are
 *                                possible for @em correlator
 *  @param[in]  vdv_indices       Indices identifying VdaggerV operators
 *  @param[out] Q0_lookup         Lookuptable containing unique combinations of
 *                                peram-, vdv-, and ric-indices needed to built
 *                                Q0
 *  @param[out] Q2L_lookup        Lookuptable containing unique combinations of
 *                                peram-, vdv-, and ric-indices needed to built
 *                                Q2L
 *  @param[out] c_look            Lookup table for C4cC
 */
static void build_C4cC_lookup(
    std::vector<quark> const &quarks,
    std::vector<int> const &quark_numbers,
    int start_config,
    const std::string &path_output,
    const std::string &overwrite,
    std::vector<std::vector<QuantumNumbers>> const &quantum_numbers,
    std::vector<std::vector<std::pair<size_t, bool>>> const &vdv_indices,
    std::vector<DilutedFactorIndex> &Q0_lookup,
    std::vector<DilutedFactorIndex> &Q2V_lookup,
    std::vector<DiagramIndex> &c_look) {
  std::vector<size_t> ql_ids(4);
  std::vector<std::pair<size_t, size_t>> ric_ids;

  // Build the correlator and dataset names for hdf5 output files
  std::vector<std::string> quark_types;
  for (const auto &id : quark_numbers)
    quark_types.emplace_back(quarks[id].type);

  for (size_t d = 0; d < quantum_numbers.size(); ++d) {
    ric_ids = create_rnd_vec_id(quarks, quark_numbers[3], quark_numbers[0], false);
    build_Quarkline_lookup_one_qn(
        0, quantum_numbers[d], vdv_indices[d], ric_ids, Q2V_lookup, ql_ids);
    ric_ids = create_rnd_vec_id(quarks, quark_numbers[0], quark_numbers[1], false);
    build_Quarkline_lookup_one_qn(
        1, quantum_numbers[d], vdv_indices[d], ric_ids, Q0_lookup, ql_ids);
    ric_ids = create_rnd_vec_id(quarks, quark_numbers[1], quark_numbers[2], false);
    build_Quarkline_lookup_one_qn(
        2, quantum_numbers[d], vdv_indices[d], ric_ids, Q2V_lookup, ql_ids);
    ric_ids = create_rnd_vec_id(quarks, quark_numbers[2], quark_numbers[3], false);
    build_Quarkline_lookup_one_qn(
        3, quantum_numbers[d], vdv_indices[d], ric_ids, Q0_lookup, ql_ids);

    std::string hdf5_dataset_name = build_hdf5_dataset_name(
        "C4+C", start_config, path_output, overwrite, quark_types, quantum_numbers[d]);

    DiagramIndex candidate(c_look.size(), hdf5_dataset_name, ql_ids);

    /*! XXX Better with std::set */
    auto it = std::find(c_look.begin(), c_look.end(), candidate);

    if (it == c_look.end()) {
      c_look.push_back(candidate);
    }
  }
}

/******************************************************************************/
/*! Create lookuptable where to find the quarklines to build C40C.
 *
 *  @param[in]  quarks            Quarks as read from the infile and processed
 *                                into quark struct
 *  @param[in]  quark_numbers     List which quarks are specified in the infile
 *  @param[in]  start_config      Number of first gauge configuration
 *  @param[in]  path_output       Output path from the infile.
 *  @param[in]  overwrite {yes,no} : deprecated
 *  @param[in]  quantum_numbers   A list of all physical quantum numbers
 *                                quantum field operators for all correlators
 *                                with Dirac structure factored out that are
 *                                possible for @em correlator
 *  @param[in]  vdv_indices       Indices identifying VdaggerV operators
 *  @param[out] Q1_lookup         Lookuptable containing unique combinations of
 *                                peram-, vdv-, and ric-indices needed to built
 *                                Q1
 *  @param[out] c_look            Lookup table for C40C
 *
 *  @bug I am fairly certain that the quarks are mixed up. It is
 *        also wrong in init_lookup_tables() (MW 27.3.17)
 */
static void build_C40C_lookup(
    std::vector<quark> const &quarks,
    std::vector<int> const &quark_numbers,
    int start_config,
    const std::string &path_output,
    const std::string &overwrite,
    std::vector<std::vector<QuantumNumbers>> const &quantum_numbers,
    std::vector<std::vector<std::pair<size_t, bool>>> const &vdv_indices,
    std::vector<DilutedFactorIndex> &Q1_lookup,
    std::vector<DiagramIndex> &c_look) {
  std::vector<size_t> ql_ids(4);
  std::vector<std::pair<size_t, size_t>> ric_ids;

  // Build the correlator and dataset names for hdf5 output files
  std::vector<std::string> quark_types;
  for (const auto &id : quark_numbers)
    quark_types.emplace_back(quarks[id].type);

  for (size_t d = 0; d < quantum_numbers.size(); ++d) {
    ric_ids = create_rnd_vec_id(quarks, quark_numbers[3], quark_numbers[0], false);
    build_Quarkline_lookup_one_qn(
        0, quantum_numbers[d], vdv_indices[d], ric_ids, Q1_lookup, ql_ids);
    ric_ids = create_rnd_vec_id(quarks, quark_numbers[0], quark_numbers[1], false);
    build_Quarkline_lookup_one_qn(
        1, quantum_numbers[d], vdv_indices[d], ric_ids, Q1_lookup, ql_ids);
    ric_ids = create_rnd_vec_id(quarks, quark_numbers[1], quark_numbers[2], false);
    build_Quarkline_lookup_one_qn(
        2, quantum_numbers[d], vdv_indices[d], ric_ids, Q1_lookup, ql_ids);
    ric_ids = create_rnd_vec_id(quarks, quark_numbers[2], quark_numbers[3], false);
    build_Quarkline_lookup_one_qn(
        3, quantum_numbers[d], vdv_indices[d], ric_ids, Q1_lookup, ql_ids);

    std::string hdf5_dataset_name = build_hdf5_dataset_name(
        "C40C", start_config, path_output, overwrite, quark_types, quantum_numbers[d]);

    DiagramIndex candidate{
        c_look.size(), hdf5_dataset_name, ql_ids, std::vector<int>({})};

    /*! XXX Better with std::set */
    auto it = std::find(c_look.begin(), c_look.end(), candidate);

    if (it == c_look.end()) {
      c_look.push_back(candidate);
    }
  }
}

/******************************************************************************/
/*! Create lookuptable where to find the quarklines to build C4cB.
 *
 *  @param[in]  quarks            Quarks as read from the infile and processed
 *                                into quark struct
 *  @param[in]  quark_numbers     List which quarks are specified in the infile
 *  @param[in]  start_config      Number of first gauge configuration
 *  @param[in]  path_output       Output path from the infile.
 *  @param[in]  overwrite {yes,no} : deprecated
 *  @param[in]  quantum_numbers   A list of all physical quantum numbers
 *                                quantum field operators for all correlators
 *                                with Dirac structure factored out that are
 *                                possible for @em correlator
 *  @param[in]  vdv_indices       Indices identifying VdaggerV operators
 *  @param[out] Q0_lookup         Lookuptable containing unique combinations of
 *                                peram-, vdv-, and ric-indices needed to built
 *                                Q0
 *  @param[out] Q2L_lookup        Lookuptable containing unique combinations of
 *                                peram-, vdv-, and ric-indices needed to built
 *                                Q2L
 *  @param[out] c_look            Lookup table for C4cB
 */
static void build_C4cB_lookup(
    std::vector<quark> const &quarks,
    std::vector<int> const &quark_numbers,
    int start_config,
    const std::string &path_output,
    const std::string &overwrite,
    std::vector<std::vector<QuantumNumbers>> const &quantum_numbers,
    std::vector<std::vector<std::pair<size_t, bool>>> const &vdv_indices,
    std::vector<DilutedFactorIndex> &Q0_lookup,
    std::vector<DilutedFactorIndex> &Q2L_lookup,
    std::vector<DiagramIndex> &c_look) {
  std::vector<size_t> ql_ids(4);
  std::vector<std::pair<size_t, size_t>> ric_ids;

  // Build the correlator and dataset names for hdf5 output files
  std::vector<std::string> quark_types;
  for (const auto &id : quark_numbers)
    quark_types.emplace_back(quarks[id].type);

  for (size_t d = 0; d < quantum_numbers.size(); ++d) {
    ric_ids = create_rnd_vec_id(quarks, quark_numbers[3], quark_numbers[0], false);
    build_Quarkline_lookup_one_qn(
        0, quantum_numbers[d], vdv_indices[d], ric_ids, Q2L_lookup, ql_ids);
    ric_ids = create_rnd_vec_id(quarks, quark_numbers[0], quark_numbers[1], false);
    build_Quarkline_lookup_one_qn(
        1, quantum_numbers[d], vdv_indices[d], ric_ids, Q0_lookup, ql_ids);
    ric_ids = create_rnd_vec_id(quarks, quark_numbers[1], quark_numbers[2], false);
    build_Quarkline_lookup_one_qn(
        2, quantum_numbers[d], vdv_indices[d], ric_ids, Q2L_lookup, ql_ids);
    ric_ids = create_rnd_vec_id(quarks, quark_numbers[2], quark_numbers[3], false);
    build_Quarkline_lookup_one_qn(
        3, quantum_numbers[d], vdv_indices[d], ric_ids, Q0_lookup, ql_ids);

    std::string hdf5_dataset_name = build_hdf5_dataset_name(
        "C4+B", start_config, path_output, overwrite, quark_types, quantum_numbers[d]);

    DiagramIndex candidate(c_look.size(), hdf5_dataset_name, ql_ids);

    /*! XXX Better with std::set */
    auto it = std::find(c_look.begin(), c_look.end(), candidate);

    if (it == c_look.end()) {
      c_look.push_back(candidate);
    }
  }
}

/******************************************************************************/
/*! Create lookuptable where to find the quarklines to build C40B.
 *
 *  @param[in]  quarks            Quarks as read from the infile and processed
 *                                into quark struct
 *  @param[in]  quark_numbers     List which quarks are specified in the infile
 *  @param[in]  start_config      Number of first gauge configuration
 *  @param[in]  path_output       Output path from the infile.
 *  @param[in]  overwrite {yes,no} : deprecated
 *  @param[in]  quantum_numbers   A list of all physical quantum numbers
 *                                quantum field operators for all correlators
 *                                with Dirac structure factored out that are
 *                                possible for @em correlator
 *  @param[in]  vdv_indices       Indices identifying VdaggerV operators
 *  @param[out] Q1_lookup         Lookuptable containing unique combinations of
 *                                peram-, vdv-, and ric-indices needed to built
 *                                Q1
 *  @param[out] c_look            Lookup table for C40B
 *
 *  @bug I am fairly certain that the quarks are mixed up. It is
 *        also wrong in init_lookup_tables() (MW 27.3.17)
 */
static void build_C40B_lookup(
    std::vector<quark> const &quarks,
    std::vector<int> const &quark_numbers,
    int start_config,
    const std::string &path_output,
    const std::string &overwrite,
    std::vector<std::vector<QuantumNumbers>> const &quantum_numbers,
    std::vector<std::vector<std::pair<size_t, bool>>> const &vdv_indices,
    std::vector<DilutedFactorIndex> &Q1_lookup,
    std::vector<DiagramIndex> &c_look) {
  std::vector<size_t> ql_ids(4);
  std::vector<std::pair<size_t, size_t>> ric_ids;

  // Build the correlator and dataset names for hdf5 output files
  std::vector<std::string> quark_types;
  for (const auto &id : quark_numbers)
    quark_types.emplace_back(quarks[id].type);

  for (size_t d = 0; d < quantum_numbers.size(); ++d) {
    ric_ids = create_rnd_vec_id(quarks, quark_numbers[3], quark_numbers[0], false);
    build_Quarkline_lookup_one_qn(
        0, quantum_numbers[d], vdv_indices[d], ric_ids, Q1_lookup, ql_ids);
    ric_ids = create_rnd_vec_id(quarks, quark_numbers[0], quark_numbers[1], false);
    build_Quarkline_lookup_one_qn(
        1, quantum_numbers[d], vdv_indices[d], ric_ids, Q1_lookup, ql_ids);
    ric_ids = create_rnd_vec_id(quarks, quark_numbers[1], quark_numbers[2], false);
    build_Quarkline_lookup_one_qn(
        2, quantum_numbers[d], vdv_indices[d], ric_ids, Q1_lookup, ql_ids);
    ric_ids = create_rnd_vec_id(quarks, quark_numbers[2], quark_numbers[3], false);
    build_Quarkline_lookup_one_qn(
        3, quantum_numbers[d], vdv_indices[d], ric_ids, Q1_lookup, ql_ids);

    std::string hdf5_dataset_name = build_hdf5_dataset_name(
        "C40B", start_config, path_output, overwrite, quark_types, quantum_numbers[d]);

    DiagramIndex candidate{
        c_look.size(), hdf5_dataset_name, ql_ids, std::vector<int>({})};

    /*! XXX Better with std::set */
    auto it = std::find(c_look.begin(), c_look.end(), candidate);

    if (it == c_look.end()) {
      c_look.push_back(candidate);
    }
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
  for (const auto &correlator : correlator_list) {
    /*! 1. Build an array (quantum_numbers) with all the quantum numbers needed
     *      for this particular correlation function.
     */
    std::vector<std::vector<QuantumNumbers>> quantum_numbers;
    build_quantum_numbers_from_correlator_list(
        correlator, operator_list, quantum_numbers);

    // Build the correlator and dataset names for hdf5 output files
    std::vector<std::string> quark_types;
    for (const auto &id : correlator.quark_numbers)
      quark_types.emplace_back(quarks[id].type);
    std::vector<std::string> hdf5_dataset_name;
    build_correlator_names(correlator.type,
                           start_config,
                           path_output,
                           overwrite,
                           quark_types,
                           quantum_numbers,
                           hdf5_dataset_name);

    /*! 2. Build the lookuptable for VdaggerV and return an array of indices
     *      corresponding to @em quantum_numbers computed in step 1. In
     *      @em vdv_indices the first entry is the id of vdv, the second tells
     *      us if vdv must be daggered to get the correct quantum numbers.
     */
    std::vector<std::vector<std::pair<size_t, bool>>> vdv_indices;
    build_VdaggerV_lookup(
        quantum_numbers, operator_lookuptable.vdaggerv_lookup, vdv_indices);
    std::vector<std::pair<size_t, size_t>> rnd_index;
    if (correlator.type == "C1") {
      build_C1_lookup(quarks,
                      correlator.quark_numbers,
                      start_config,
                      path_output,
                      overwrite,
                      quantum_numbers,
                      vdv_indices,
                      quarkline_lookuptable.Q1,
                      correlator_lookuptable.trQ1,
                      correlator_lookuptable.C1);
    } else if (correlator.type == "C2+" || correlator.type == "Check") {
      /*! 3. Build the lookuptable for rVdaggerVr and return an array of indices
       *      corresponding to the 'quantum_numbers' computed in step 1.
       *  4. Build the lookuptable for Q2 and return an array of indices
       *      corresponding to the 'quantum_numbers' computed in step 1.
       *  5. Build the lookuptable for the correlation functions
       */
      build_C2c_lookup(quarks,
                       correlator.quark_numbers,
                       start_config,
                       path_output,
                       overwrite,
                       quantum_numbers,
                       vdv_indices,
                       quarkline_lookuptable.Q0,
                       quarkline_lookuptable.Q2V,
                       correlator_lookuptable.trQ0Q2,
                       correlator_lookuptable.C2c);
    }
    /*! 6. Repeat steps 1.-5. for all correlators in correlator_list. Where
     *      applicable rVdaggerVr must be replaced by rVdaggerV in step 3. and
     *      Q2 by Q1 in step 4.
     */
    else if (correlator.type == "C3+") {
      build_C3c_lookup(quarks,
                       correlator.quark_numbers,
                       start_config,
                       path_output,
                       overwrite,
                       quantum_numbers,
                       vdv_indices,
                       quarkline_lookuptable.Q0,
                       quarkline_lookuptable.Q1,
                       quarkline_lookuptable.Q2L,
                       correlator_lookuptable.C3c);
    } else if (correlator.type == "C4+D") {
      build_C4cD_lookup(quarks,
                        correlator.quark_numbers,
                        start_config,
                        path_output,
                        overwrite,
                        quantum_numbers,
                        vdv_indices,
                        quarkline_lookuptable.Q0,
                        quarkline_lookuptable.Q2V,
                        correlator_lookuptable.trQ0Q2,
                        correlator_lookuptable.C4cD);
    } else if (correlator.type == "C4+V") {
      build_C4cV_lookup(quarks,
                        correlator.quark_numbers,
                        start_config,
                        path_output,
                        overwrite,
                        quantum_numbers,
                        vdv_indices,
                        quarkline_lookuptable.Q0,
                        quarkline_lookuptable.Q2V,
                        correlator_lookuptable.trQ0Q2,
                        correlator_lookuptable.C4cV);
    } else if (correlator.type == "C4+C") {
      build_C4cC_lookup(quarks,
                        correlator.quark_numbers,
                        start_config,
                        path_output,
                        overwrite,
                        quantum_numbers,
                        vdv_indices,
                        quarkline_lookuptable.Q0,
                        quarkline_lookuptable.Q2V,
                        correlator_lookuptable.C4cC);
    } else if (correlator.type == "C4+B") {
      build_C4cB_lookup(quarks,
                        correlator.quark_numbers,
                        start_config,
                        path_output,
                        overwrite,
                        quantum_numbers,
                        vdv_indices,
                        quarkline_lookuptable.Q0,
                        quarkline_lookuptable.Q2L,
                        correlator_lookuptable.C4cB);
    } else if (correlator.type == "C20V") {
      build_C20V_lookup(quarks,
                        correlator.quark_numbers,
                        start_config,
                        path_output,
                        overwrite,
                        quantum_numbers,
                        vdv_indices,
                        quarkline_lookuptable.Q1,
                        correlator_lookuptable.trQ1,
                        correlator_lookuptable.C20V);
    } else if (correlator.type == "C20") {
      build_C20_lookup(quarks,
                       correlator.quark_numbers,
                       start_config,
                       path_output,
                       overwrite,
                       quantum_numbers,
                       vdv_indices,
                       quarkline_lookuptable.Q1,
                       correlator_lookuptable.trQ1Q1,
                       correlator_lookuptable.C20);
    } else if (correlator.type == "C30V") {
      build_C30V_lookup(quarks,
                        correlator.quark_numbers,
                        start_config,
                        path_output,
                        overwrite,
                        quantum_numbers,
                        vdv_indices,
                        quarkline_lookuptable.Q1,
                        correlator_lookuptable.trQ1,
                        correlator_lookuptable.trQ1Q1,
                        correlator_lookuptable.C30V);
    } else if (correlator.type == "C30") {
      build_C30_lookup(quarks,
                       correlator.quark_numbers,
                       start_config,
                       path_output,
                       overwrite,
                       quantum_numbers,
                       vdv_indices,
                       quarkline_lookuptable.Q1,
                       correlator_lookuptable.C30);
    } else if (correlator.type == "C40D") {
      build_C40D_lookup(quarks,
                        correlator.quark_numbers,
                        start_config,
                        path_output,
                        overwrite,
                        quantum_numbers,
                        vdv_indices,
                        quarkline_lookuptable.Q1,
                        correlator_lookuptable.trQ1Q1,
                        correlator_lookuptable.C40D);
    } else if (correlator.type == "C40V") {
      build_C40V_lookup(quarks,
                        correlator.quark_numbers,
                        start_config,
                        path_output,
                        overwrite,
                        quantum_numbers,
                        vdv_indices,
                        quarkline_lookuptable.Q1,
                        correlator_lookuptable.trQ1Q1,
                        correlator_lookuptable.C40V);
    } else if (correlator.type == "C40C") {
      build_C40C_lookup(quarks,
                        correlator.quark_numbers,
                        start_config,
                        path_output,
                        overwrite,
                        quantum_numbers,
                        vdv_indices,
                        quarkline_lookuptable.Q1,
                        correlator_lookuptable.C40C);
    } else if (correlator.type == "C40B") {
      build_C40B_lookup(quarks,
                        correlator.quark_numbers,
                        start_config,
                        path_output,
                        overwrite,
                        quantum_numbers,
                        vdv_indices,
                        quarkline_lookuptable.Q1,
                        correlator_lookuptable.C40B);
    } else {
      std::cout << "Correlator type not known!" << std::endl;
      exit(0);
    }
  }

  /*! Sets index_of_unity to the index of operator_lookuptable.vdaggerv_lookup
   *  where momentum and displacement are both zero, or to -1 if no such entry
   *  is found.
   */
  std::array<int, 3> const zero{0, 0, 0};
  bool found = false;
  for (const auto &op_vdv : operator_lookuptable.vdaggerv_lookup)
    if ((op_vdv.momentum == zero) && (op_vdv.displacement.empty() )) {
      operator_lookuptable.index_of_unity = op_vdv.id;
      found = true;
    }
  if (!found)
    operator_lookuptable.index_of_unity = -1;
}
