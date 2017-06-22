#include "icp_single_iteration.h"
#include "random_points_on_mesh.h"
#include "point_mesh_distance.h"
#include "point_to_point_rigid_matching.h"
#include "point_to_plane_rigid_matching.h"

using namespace Eigen;

void icp_single_iteration(
  const Eigen::MatrixXd & VX,
  const Eigen::MatrixXi & FX,
  const Eigen::MatrixXd & VY,
  const Eigen::MatrixXi & FY,
  const int num_samples,
  const ICPMethod method,
  Eigen::Matrix3d & R,
  Eigen::RowVector3d & t)
{
  R = Matrix3d::Identity();
  t = RowVector3d::Zero();

  MatrixXd X;
  X.resize(num_samples,3);

  //cover (VX,FX) with random points, and store them in X
  random_points_on_mesh(num_samples,VX,FX,X);

  MatrixXd P;
  MatrixXd N;
  VectorXd D;
  P.resizeLike(X);
  N.resizeLike(X);
  D.resize(num_samples);
  
  //Project random points onto (VY,FY), store distances D, projections P, normals N
  point_mesh_distance(X,VY,FY,D,P,N);

  if(method == ICP_METHOD_POINT_TO_POINT) {
    point_to_point_rigid_matching(X,P,R,t);
  }
  if(method == ICP_METHOD_POINT_TO_PLANE) {
    point_to_plane_rigid_matching(X,P,N,R,t);
  }
}
