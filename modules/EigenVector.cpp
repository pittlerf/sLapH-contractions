#include "EigenVector.h"

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void LapH::EigenVector::read_eigen_vector(const std::string& filename, 
                                          const size_t t, const size_t verbose){

  //buffer for read in
  vec eigen_vec(V[t].rows());
  std::cout << "\tReading eigenvectors from files:" << filename << std::endl;

  //setting up file
  std::ifstream infile(filename, std::ifstream::binary); 
  if (infile) {
    for (size_t ncol = 0; ncol < V[t].cols(); ++ncol) {
      infile.read( (char*) &(eigen_vec[0]), 2*V[t].rows()*sizeof(double));
      for(size_t nrow = 0; nrow < V[t].rows(); ++nrow){
        (V[t])(nrow, ncol) = eigen_vec[nrow];
      }
    }
  }
  else {
    std::cout << "eigenvector file does not exist!!!\n" << std::endl;
    exit(0);
  }
  infile.close();

  // small test of trace and sum over the eigen vector matrix!
  if(verbose){
    std::cout << "trace of V^d*V" << ":\t"
        << (V[t].adjoint() * V[t]).trace() << std::endl;
    std::cout << "sum over all entries of V^d*V" << ":\t"
        << (V[t].adjoint() * V[t]).sum() << std::endl;
  }   
}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void LapH::EigenVector::read_eigen_vector(const std::string& filename,
                                          const size_t verbose){

  for(int t = 0; t < V.size(); t++){
    char buff[10];
    sprintf(buff, "%03d", t);
    read_eigen_vector(filename + buff, t, verbose);
  }
}


