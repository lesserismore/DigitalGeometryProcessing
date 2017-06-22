#include "fd_grad.h"
#include "fd_partial_derivative.h"
#include <igl/cat.h>

void fd_grad(
  const int nx,
  const int ny,
  const int nz,
  const double h,
  Eigen::SparseMatrix<double> & G)
{
  ////////////////////////////////////////////////////////////////////////////
  // Add your code here
  ////////////////////////////////////////////////////////////////////////////
  
  //#columns of G
  int C = nx*(ny*nz);

  //Rx+Ry+Rz = #rows of Dx + #rows of Dy + #rows of Dz = #rows of G
  int Rx = (nx-1)*(ny*nz);
  int Ry = nx*((ny-1)*nz);
  int Rz = nx*(ny*(nz-1));

  //intialize Dx, Dy, Dz
  Eigen::SparseMatrix<double> Dx = Eigen::SparseMatrix<double>(Rx,C);
  Eigen::SparseMatrix<double> Dy = Eigen::SparseMatrix<double>(Ry,C);
  Eigen::SparseMatrix<double> Dz = Eigen::SparseMatrix<double>(Rz,C);

  //fill Dx, Dy, Dz
  fd_partial_derivative(nx,ny,nz,h,0,Dx);
  fd_partial_derivative(nx,ny,nz,h,1,Dy);
  fd_partial_derivative(nx,ny,nz,h,2,Dz);

  //concatenate Dx, Dy, Dz to fill G
  Eigen::SparseMatrix<double> temp;
  igl::cat(1,Dx,Dy,temp);
  igl::cat(1,temp,Dz,G);

  return;
}
