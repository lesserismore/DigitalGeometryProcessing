#include "hausdorff_lower_bound.h"
#include "random_points_on_mesh.h"
#include "point_mesh_distance.h"

using namespace Eigen;

double hausdorff_lower_bound(
  const Eigen::MatrixXd & VX,
  const Eigen::MatrixXi & FX,
  const Eigen::MatrixXd & VY,
  const Eigen::MatrixXi & FY,
  const int n)
{
  double d;
  MatrixXd X;
  X.resize(n,3);

  //cover (VX,FX) with random sample points
  random_points_on_mesh(n,VX,FX,X);

  VectorXd D;
  MatrixXd P;
  MatrixXd N;
  D.resize(n);
  P.resizeLike(X);
  N.resizeLike(X);

  //calculate closest points on (VY<FY) from points X
  point_mesh_distance(X,VY,FY,D,P,N);
  //find the distance corresponding to furthest of closest points
  d = D.maxCoeff();
  return d;
}
