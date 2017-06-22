#include "poisson_surface_reconstruction.h"
#include "fd_grad.h"
#include "fd_interpolate.h"
#include <igl/copyleft/marching_cubes.h>
#include <algorithm>
#include <iostream>
#include <fstream>

void poisson_surface_reconstruction(
    const Eigen::MatrixXd & P,
    const Eigen::MatrixXd & N,
    Eigen::MatrixXd & V,
    Eigen::MatrixXi & F)
{
  ////////////////////////////////////////////////////////////////////////////
  // Construct FD grid, CONGRATULATIONS! You get this for free!
  ////////////////////////////////////////////////////////////////////////////
  // number of input points
  const int n = P.rows();
  // Grid dimensions
  int nx, ny, nz;
  // Maximum extent (side length of bounding box) of points
  double max_extent =
    (P.colwise().maxCoeff()-P.colwise().minCoeff()).maxCoeff();
  // padding: number of cells beyond bounding box of input points
  const double pad = 8;
  // choose grid spacing (h) so that shortest side gets 30+2*pad samples
  double h  = max_extent/double(30+2*pad);
  // Place bottom-left-front corner of grid at minimum of points minus padding
  Eigen::RowVector3d corner = P.colwise().minCoeff().array()-pad*h;
  // Grid dimensions should be at least 3 
  nx = std::max((P.col(0).maxCoeff()-P.col(0).minCoeff()+(2.*pad)*h)/h,3.);
  ny = std::max((P.col(1).maxCoeff()-P.col(1).minCoeff()+(2.*pad)*h)/h,3.);
  nz = std::max((P.col(2).maxCoeff()-P.col(2).minCoeff()+(2.*pad)*h)/h,3.);
  // Compute positions of grid nodes
  Eigen::MatrixXd x(nx*ny*nz, 3);
  for(int i = 0; i < nx; i++) 
  {
    for(int j = 0; j < ny; j++)
    {
      for(int k = 0; k < nz; k++)
      {
         // Convert subscript to index
         const auto ind = i + nx*(j + k * ny);
         x.row(ind) = corner + h*Eigen::RowVector3d(i,j,k);
      }
    }
  }
  Eigen::VectorXd g = Eigen::VectorXd::Zero(nx*ny*nz);

  //get gradient matrix G
  Eigen::SparseMatrix<double> G = 
    Eigen::SparseMatrix<double>((nx-1)*(ny*nz)+nx*((ny-1)*nz)+nx*(ny*(nz-1)), nx*ny*nz);
  fd_grad(nx,ny,nz,h,G);

  //Initialize interpolated weight matrices for each staggered grid
  Eigen::SparseMatrix<double> Wx = Eigen::SparseMatrix<double>(n,(nx-1)*(ny*nz));
  Eigen::SparseMatrix<double> Wy = Eigen::SparseMatrix<double>(n,nx*((ny-1)*nz));
  Eigen::SparseMatrix<double> Wz = Eigen::SparseMatrix<double>(n,nx*(ny*(nz-1)));

  //Adjust corners of staggered grids
  double stagger = h/2.;
  Eigen::RowVector3d xcorner;
  Eigen::RowVector3d ycorner;
  Eigen::RowVector3d zcorner;

  xcorner[0] = corner[0]+stagger;
  xcorner[1] = corner[1];
  xcorner[2] = corner[2];

  ycorner[0] = corner[0];
  ycorner[1] = corner[1]+stagger;
  ycorner[2] = corner[2];

  zcorner[0] = corner[0];
  zcorner[1] = corner[1];
  zcorner[2] = corner[2]+stagger;

  //interpolate onto staggered grids
  fd_interpolate(nx-1,ny,nz,h,xcorner,P,Wx);
  fd_interpolate(nx,ny-1,nz,h,ycorner,P,Wy);
  fd_interpolate(nx,ny,nz-1,h,zcorner,P,Wz);
  
  //Intialize vectors to contain directional normal data on staggered grids
  Eigen::VectorXd vx((nx-1)*(ny*nz));
  Eigen::VectorXd vy(nx*((ny-1)*nz));
  Eigen::VectorXd vz(nx*(ny*(nz-1)));

  //distribute normal data
  vx = Wx.transpose()*N.col(0);
  vy = Wy.transpose()*N.col(1);
  vz = Wz.transpose()*N.col(2);

  //concatenate normal data vectors together
  Eigen::VectorXd v((nx-1)*(ny*nz)+nx*((ny-1)*nz)+nx*(ny*(nz-1)));
  v << vx,
       vy,
       vz;

  //scale v by h, since gradient G contains only -1,0, or 1 coeffs
  v *= h;

  //solve for g_off in G.transpose()*G*g_off=G.transpose()*v 
  Eigen::VectorXd g_off = Eigen::VectorXd::Zero(nx*ny*nz);
  Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> solver;

  solver.compute(G.transpose()*G);
  if(solver.info()!=Eigen::Success) {
    // decomposition failed
    return;
  }
  g_off = solver.solve(G.transpose()*v);
  if(solver.info()!=Eigen::Success) {
    // solving failed
    return;
  }
  
  //interpolate points P
  Eigen::SparseMatrix<double> W = Eigen::SparseMatrix<double>(n,nx*(ny*nz));
  fd_interpolate(nx,ny,nz,h,corner,P,W);

  //Calculate iso-value sigma
  double d = (W*g_off).sum();
  float sigma = d/n;

  //pre-shift g-values
  g = g_off - sigma*Eigen::VectorXd::Ones(nx*(ny*nz));

  ////////////////////////////////////////////////////////////////////////////
  // Run black box algorithm to compute mesh from implicit function: this
  // function always extracts g=0, so "pre-shift" your g values by -sigma
  ////////////////////////////////////////////////////////////////////////////
  igl::copyleft::marching_cubes(g, x, nx, ny, nz, V, F);
}
