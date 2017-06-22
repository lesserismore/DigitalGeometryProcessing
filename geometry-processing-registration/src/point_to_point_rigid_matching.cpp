#include "point_to_point_rigid_matching.h"
#include "closest_rotation.h"

using namespace Eigen;

void point_to_point_rigid_matching(
  const Eigen::MatrixXd & X,
  const Eigen::MatrixXd & P,
  Eigen::Matrix3d & R,
  Eigen::RowVector3d & t)
{
  R = Eigen::Matrix3d::Identity();
  t = Eigen::RowVector3d::Zero();

  //Find centroids of query points and of projection points
  RowVector3d x_c = (X.colwise().sum())*(1.0/(double(X.rows())));
  RowVector3d p_c = (P.colwise().sum())*(1.0/(double(P.rows())));

  //translate points by centroids
  MatrixXd X_c = X.rowwise() - x_c;
  MatrixXd P_c = P.rowwise() - p_c;

  //Calculate closed form of M
  Matrix3d M = (X_c.transpose())*P_c;
 
  //find energy-minimizing rotation and translation 
  closest_rotation(M,R);
  t = p_c - (R*(x_c.transpose())).transpose();
}

