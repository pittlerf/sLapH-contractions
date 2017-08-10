/*! @file global_data_typedefs.h
 *  Miscallaneous structs needed in GlobalData
 *
 *  @author Bastian Knippschild
 *  @author Markus Werner
 */

#ifndef OPERATORS_H_
#define OPERATORS_H_

#include <typedefs.h>


/*! Small struct which contains all information to build and read random 
 *  vectors.
 *
 *  How many are needed, their length, and what the input paths are.
 */
struct RandomVectorConstruction {

  size_t nb_entities, length;
  std::vector<std::string> filename_list;

};

/*! Small struct which contains all information to build and read 
 *  Perambulators.
 *
 *  How many are needed, their size, and what the input paths are.
 */
struct PerambulatorConstruction {

  size_t nb_entities;
  std::vector<size_t> size_rows;
  std::vector<size_t> size_cols;
  std::vector<std::string> filename_list;

};

/*! Quark type that contains all quark propagator informations: 
 *
 *  Flavor, number of random vectors, dilution scheme, path and a unique id. 
 */
struct quark {

  std::string type; /*!< Flavor */
  int number_of_rnd_vec;

  std::string dilution_T;
  int number_of_dilution_T;

  std::string dilution_E;
  int number_of_dilution_E;

  std::string dilution_D;
  int number_of_dilution_D;

  size_t id;
  std::string path;

  /*! Constructor */
  quark (std::string type, int number_of_rnd_vec, std::string dilution_T,
         int number_of_dilution_T, std::string dilution_E,
         int number_of_dilution_E, std::string dilution_D,
                           int number_of_dilution_D, int id, std::string path) :
      type(type), number_of_rnd_vec(number_of_rnd_vec), 
      dilution_T(dilution_T), number_of_dilution_T(number_of_dilution_T), 
      dilution_E(dilution_E), number_of_dilution_E(number_of_dilution_E), 
      dilution_D(dilution_D), number_of_dilution_D(number_of_dilution_D),
      id(id), path(path) {}
};

/*! Struct that contains all information operator specifying the field operator 
 *  in the infile:
 *
 *  Dirac structure, displacement and momentum vector as 3-vectors
 */
struct Operators {

public: // TODO: should be changed to private at a later point
  std::vector<int> gammas;
  std::array<int, 3> dil_vec;
  std::vector< std::vector<std::array<int, 3> > > mom_vec;

public:
  /*! Constructor */
  Operators (std::vector<int> gammas, std::array<int, 3> dil_vec, 
                       std::vector<std::vector<std::array<int, 3> > > mom_vec) :
      gammas(gammas), dil_vec(dil_vec), mom_vec(mom_vec) {}

};

/*! Struct that contains all information specifying the correlator in the 
 *  infile:
 *
 *  Internal name of the diagram (@ref LapH::Correlators), and which quark from 
 *  quarks and operator from oeprator_list are to be used
 *
 *  @todo Are GEVP and tot_mom used at the moment? They should support 
 *        building off-diagonal elements as well or resticting to a certain 
 *        CMS momentum 
 */
struct Correlators {

/*! @{ */
/*! @todo Change to private at a later point
 */
public: 
  std::string type;
  std::vector<int> quark_numbers;
  std::vector<int> operator_numbers;
  std::string GEVP;
  std::vector<std::array<int, 3> > tot_mom;
/*! @} */

public:
/*! Constructor */
  Correlators(std::string type, std::vector<int> quark_numbers, 
              std::vector<int> operator_numbers, std::string GEVP, 
              std::vector< std::array<int, 3> > tot_mom) :
      type(type), quark_numbers(quark_numbers), 
      operator_numbers(operator_numbers), GEVP(GEVP), tot_mom(tot_mom) {}

};

typedef std::vector<Operators> Operator_list;
typedef std::vector<Correlators> Correlator_list;

#endif /* OPERATORS_H_ */
