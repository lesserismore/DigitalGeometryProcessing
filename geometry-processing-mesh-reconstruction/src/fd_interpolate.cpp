#include "fd_interpolate.h"
#include <iostream>
#include <vector>
#include <cmath>

void fd_interpolate(
  const int nx,
  const int ny,
  const int nz,
  const double h,
  const Eigen::RowVector3d & corner,
  const Eigen::MatrixXd & P,
  Eigen::SparseMatrix<double> & W)
{
  ////////////////////////////////////////////////////////////////////////////
  // Add your code here
  ////////////////////////////////////////////////////////////////////////////

  // r = rows of W, c = columns of W
  int R = P.rows();
  int C = (nx*ny)*nz;
  W = Eigen::SparseMatrix<double>(R,C);

  std::vector<Eigen::Triplet<double>> nodeWeights;
  Eigen::RowVector3d p;
  int n_1,n_2,n_3;
  double alpha,beta,gamma;
  Eigen::Triplet<double> t;

  //voxel volume to scale weights
  double h3 = h*h*h;

  for(int r=0; r<R; r++) {
    //get point p
    p = P.row(r);

    //find x,y,z distances from corner to p
    p = p - corner;
    
    //find x,y,z indices of bottom-left-front node of cell containing p
    n_1 = floor(p[0]/h);
    n_2 = floor(p[1]/h);
    n_3 = floor(p[2]/h);

    //calculate parameters of trilinear interpolation
    alpha = p[0]-n_1*h;
    beta = p[1]-n_2*h;
    gamma = p[2]-n_3*h;

    //construct and insert 8 column index and value pairs for each row
    t = Eigen::Triplet<double>(r,n_1+n_2*nx+n_3*(nx*ny),((h-gamma)*((h-beta)*(h-alpha)))/h3);
    nodeWeights.push_back(t);

    t = Eigen::Triplet<double>(r,(n_1+1)+n_2*nx+n_3*(nx*ny),((h-gamma)*((h-beta)*alpha))/h3);
    nodeWeights.push_back(t);

    t = Eigen::Triplet<double>(r,n_1+(n_2+1)*nx+n_3*(nx*ny),((h-gamma)*(beta*(h-alpha)))/h3);
    nodeWeights.push_back(t);

    t = Eigen::Triplet<double>(r,(n_1+1)+(n_2+1)*nx+n_3*(nx*ny),((h-gamma)*(beta*alpha))/h3);
    nodeWeights.push_back(t);

    t = Eigen::Triplet<double>(r,n_1+n_2*nx+(n_3+1)*(nx*ny),(gamma*((h-beta)*(h-alpha))/h3));
    nodeWeights.push_back(t);

    t = Eigen::Triplet<double>(r,(n_1+1)+n_2*nx+(n_3+1)*(nx*ny),(gamma*((h-beta)*alpha))/h3);
    nodeWeights.push_back(t);

    t = Eigen::Triplet<double>(r,n_1+(n_2+1)*nx+(n_3+1)*(nx*ny),((gamma*beta)*(h-alpha))/h3);
    nodeWeights.push_back(t);

    t = Eigen::Triplet<double>(r,(n_1+1)+(n_2+1)*nx+(n_3+1)*(nx*ny),(gamma*(beta*alpha))/h3);
    nodeWeights.push_back(t);
  } 

  W.setFromTriplets(nodeWeights.begin(), nodeWeights.end());
  return;
}
