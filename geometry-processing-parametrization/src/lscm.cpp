#include "lscm.h"
#include "vector_area_matrix.h"
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include <igl/repdiag.h>
#include <igl/eigs.h>

using namespace Eigen;
using namespace igl;

void lscm(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::MatrixXd & U)
{
  //get vector area matrix
  SparseMatrix<double> A;
  vector_area_matrix(F,A);

  //get cotangent Laplacian
  SparseMatrix<double> L;
  cotmatrix(V,F,L);

  //get mass matrix
  SparseMatrix<double> M;
  massmatrix(V,F,MASSMATRIX_TYPE_BARYCENTRIC,M);

  //adjust L and M to accomodate vectorized UV coordinates
  SparseMatrix<double> LL;
  repdiag(L,2,LL);
  SparseMatrix<double> B;
  repdiag(M,2,B);

  SparseMatrix<double> Q = 0.5*(LL - A);
  
  //minimize 0.5*(U^TQU) subject to U^TBU=1
  int v = V.rows();
  U.resize(v,2);
  Vector3d sS; //eigenvalues
  MatrixXd sU; //eigenvectors
  sU.resize(Q.rows(),3);
  eigs(Q,B,3,EIGS_TYPE_SM,sU,sS);
  U.col(0) = sU.block(0,2,v,1);
  U.col(1) = sU.block(v,2,v,1);

  //find and apply canonical rotation to UV coordinates
  JacobiSVD<MatrixXd> svd(U,ComputeFullV);
  Matrix2d T = svd.matrixV();
  U = U*T;
}
