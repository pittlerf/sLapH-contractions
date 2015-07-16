#include "Perambulator.h"

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void LapH::Perambulator::read_perambulator(const size_t entity, 
                                           const std::string& filename) {

  clock_t t = clock();
  FILE *fp = NULL;

//  const int Lt = global_data->get_Lt();
//  const int Lx = global_data->get_Lx();
//  const int Ly = global_data->get_Ly();
//  const int Lz = global_data->get_Lz();
//  const int Vs = Lx * Ly * Lz;
//
//
//  const int number_of_eigen_vec = global_data->get_number_of_eigen_vec();
//  const int number_of_rnd_vec = q.number_of_rnd_vec;
//  const int number_of_inversions = (Lt / q.number_of_dilution_T)
//      * q.number_of_dilution_E * q.number_of_dilution_D;
//  const int size_perambulator_entry = number_of_inversions * Lt * 4
//      * number_of_eigen_vec;


  std::cout << "\tReading perambulator from file:\n\t\t" << filename;

  // reading the data into temporary array
  vec perambulator_read(size_rows*size_cols);
  if((fp = fopen(filename.c_str(), "rb")) == NULL){
    std::cout << "failed to open file to read perambulaots: " 
              << filename << "\n" << std::endl;
    exit(0);
  }
  int blabla = fread(&(perambulator_read[0]), sizeof(std::complex<double>),
                     size_rows*size_cols, fp);
  fclose(fp);

//  // re-sorting and copy into matrix structure 
//  size_t col_i, row_i;
//  for(size_t t1 = 0; t1 < Lt; ++t1)
//    for(size_t ev1 = 0; ev1 < number_of_eigen_vec; ++ev1)
//      for(size_t dirac1 = 0; dirac1 < 4; ++dirac1)
//        for(size_t t2 = 0; t2 < (Lt / q.number_of_dilution_T); ++t2)
//          for(size_t ev2 = 0; ev2 < q.number_of_dilution_E; ++ev2)
//            for(size_t dirac2 = 0; dirac2 < q.number_of_dilution_D; ++dirac2){
//              row_i = 4 * number_of_eigen_vec * t1 + 4 * ev1 + dirac1;
//              col_i = q.number_of_dilution_D * 
//                  q.number_of_dilution_E * t2 + 
//                  q.number_of_dilution_D * ev2 + dirac2;
//              perambulator[entity](4 * number_of_eigen_vec * t1 +
//                  number_of_eigen_vec * dirac1 + ev1, 
//                  q.number_of_dilution_E * q.number_of_dilution_D * t2 + 
//                  q.number_of_dilution_E * dirac2 + ev2) = 
//              perambulator_read[row_i * number_of_inversions + col_i];
//            }

  // writing out how long it took to read the file
  t = clock() - t;
  std::cout << "\n\t\tin: " << std::fixed << std::setprecision(1)
            << ((float) t)/CLOCKS_PER_SEC << " seconds" << std::endl;
}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void LapH::Perambulator::read_perambulators_from_separate_files(
                                 const std::vector<std::string>&filename_list) {

  if(filename_list.size() != nb_entities)
    std::cout << "Problem when reading perambulators: The number of "
              << "perambulators read is not the same as the expected one!" 
              << std::endl;
  for(size_t i = 0; i < filename_list.size(); i++)
    read_perambulator(i, filename_list[i]);

}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
