#include <gtest/gtest.h>
#include "Gamma.h"
namespace{
  const std::complex<double> I {0.0, 1.0};
}
Eigen::Matrix4cd sparse_to_eigen_dense(const gamma_lookup sparse){
  Eigen::Matrix4cd dense;
  for (size_t col = 0; col < 4; ++col){
    dense(sparse.row[col],col) = sparse.value[col];
  }
  return dense;
}

TEST(gamma,checkBasis){
  Eigen::Matrix4cd tmlqcd_gamma1;
  tmlqcd_gamma1 << 0, 0, 1, 0,
                    0, 0, 0, 1,
                    1, 0, 0, 0,
                    0, 1, 0, 0;
  EXPECT_EQ(sparse_to_eigen_dense(gamma_vec[1]),tmlqcd_gamma1);
}
