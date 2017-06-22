#include "point_to_plane_rigid_matching.h"
#include "closest_rotation.h"
#include <iostream>
#include <Eigen/Dense>

using namespace Eigen;

void point_to_plane_rigid_matching(
  const Eigen::MatrixXd & X,
  const Eigen::MatrixXd & P,
  const Eigen::MatrixXd & N,
  Eigen::Matrix3d & R,
  Eigen::RowVector3d & t)
{
  R = Eigen::Matrix3d::Identity();
  t = Eigen::RowVector3d::Zero();

  //get displacement vectors
  VectorXd X0 = P.col(0)-X.col(0);
  VectorXd X1 = P.col(1)-X.col(1);
  VectorXd X2 = P.col(2)-X.col(2);
   
  int n = X.rows();

  //stack displacement vectors vertically
  VectorXd Xs(3*n);
  Xs << X0,
        X1,
        X2;

  VectorXd ZeroV = VectorXd::Zero(n);
  VectorXd OneV = VectorXd::Ones(n);
  MatrixXd A(3*n,6);

  //fill A such that A*u gives transformed X
  A.col(0) << ZeroV,1.0*X.col(2),-1.0*X.col(1);
  A.col(1) << -1.0*X.col(2),ZeroV,1.0*X.col(0);
  A.col(2) << 1.0*X.col(1),-1.0*X.col(0),ZeroV;
  A.col(3) << OneV, ZeroV, ZeroV;
  A.col(4) << ZeroV, OneV, ZeroV;
  A.col(5) << ZeroV, ZeroV, OneV;

  //align normal component data into diagonal matrices
  MatrixXd N0 = N.col(0).asDiagonal();
  MatrixXd N1 = N.col(1).asDiagonal();
  MatrixXd N2 = N.col(2).asDiagonal();

  //stack normal data matrices horizontally
  MatrixXd Ns(n,3*n);
  Ns << N0, N1, N2;

  //calculate tranposes
  MatrixXd Ns_T = Ns.transpose();
  MatrixXd A_T = A.transpose();

  MatrixXd temp = A_T*Ns_T*Ns;
  MatrixXd LHS = temp*A;
  VectorXd RHS = temp*Xs;

  //solve for u in 6x6 linear system A_T*Ns_T*Ns*A*u=A_T*Ns_T*Ns*Xs
  VectorXd u = LHS.colPivHouseholderQr().solve(RHS);

  //Get linearized M from u values, add to 3x3 identity
  Matrix3d M;
  M << 1.0, -1.0*u[2], u[1],
       u[2], 1.0, -1.0*u[0],
       -1.0*u[1], u[0], 1.0;
  //find closest rotation
  closest_rotation(M,R);

  //Get t from u values
  t << u[3], u[4], u[5];
}
