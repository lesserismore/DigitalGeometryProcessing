#include "closest_rotation.h"
#include <Eigen/Dense>

using namespace Eigen;

void closest_rotation(
  const Eigen::Matrix3d & M,
  Eigen::Matrix3d & R)
{
  //compute Singular Value Decomposition of M
  JacobiSVD<Matrix3d> svd(M,ComputeFullU | ComputeFullV);
  Matrix3d U = svd.matrixU();
  Matrix3d V = svd.matrixV();
  
  //Set Omega, following in-class derivation
  Matrix3d Omega = Matrix3d::Identity();
  Matrix3d V_T = V.transpose();
  Omega(2,2) = (U*V_T).determinant();

  //Calculate R
  R = U*Omega*V_T;
}
