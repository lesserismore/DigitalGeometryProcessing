#include "fd_partial_derivative.h"
#include <vector>
#include <cmath>

void fd_partial_derivative(
  const int nx,
  const int ny,
  const int nz,
  const double h,
  const int dir,
  Eigen::SparseMatrix<double> & D)
{
  ////////////////////////////////////////////////////////////////////////////
  // Add your code here
  ////////////////////////////////////////////////////////////////////////////
 
  std::vector<Eigen::Triplet<double>> partialTriplets;
  Eigen::Triplet<double> t1,t2;

  //#rows of D
  int R;
  int r; //row index

  //#columns of D
  int C = nx*(ny*nz);
  int c1,c2; //column indices

  //values to insert into D
  double v1 = -1;
  double v2 = 1;

  //grid node indices
  int i,j,k;

  //x-direction
  if(dir == 0){
    int dnx = nx-1;
    R = dnx*(ny*nz);
    D = Eigen::SparseMatrix<double>(R,C);

    //from row index, calculate grid indices
    for(r=0; r<R; r++){
      k = floor(r/(dnx*ny));
      j = floor((r-(k*(dnx*ny)))/dnx);
      i = r-((j*dnx)+(k*(dnx*ny)));

      //from grid indices, calculate column indices
      c1 = i+j*nx+k*(nx*ny);
      c2 = (i+1)+j*nx+k*(nx*ny);

      //form triplets and insert into triplet vector
      t1 = Eigen::Triplet<double>(r,c1,v1);
      t2 = Eigen::Triplet<double>(r,c2,v2);
      partialTriplets.push_back(t1);
      partialTriplets.push_back(t2);
    }

    //build partial derivative matrix from triplets
    D.setFromTriplets(partialTriplets.begin(), partialTriplets.end());
  }

  //similar for y-direction
  if(dir == 1){
    int dny = ny-1;
    R = nx*(dny*nz);
    D = Eigen::SparseMatrix<double>(R,C);
    for(r=0; r<R; r++){
      k = floor(r/(nx*dny));
      j = floor((r-(k*(nx*dny)))/nx);
      i = r-((j*nx)+(k*(nx*dny)));

      c1 = i+j*nx+k*(nx*ny);
      c2 = i+(j+1)*nx+k*(nx*ny);

      t1 = Eigen::Triplet<double>(r,c1,v1);
      t2 = Eigen::Triplet<double>(r,c2,v2);
      partialTriplets.push_back(t1);
      partialTriplets.push_back(t2);
    }
    D.setFromTriplets(partialTriplets.begin(), partialTriplets.end());
  }

  //similar for z-direction
  if(dir == 2){
    int dnz = nz-1;
    R = nx*(ny*dnz);
    D = Eigen::SparseMatrix<double>(R,C);
    for(r=0; r<R; r++){
      k = floor(r/(nx*ny));
      j = floor((r-(k*(nx*ny)))/nx);
      i = r-((j*nx)+(k*(nx*ny)));

      c1 = i+j*nx+k*(nx*ny);
      c2 = i+j*nx+(k+1)*(nx*ny);

      t1 = Eigen::Triplet<double>(r,c1,v1);
      t2 = Eigen::Triplet<double>(r,c2,v2);
      partialTriplets.push_back(t1);
      partialTriplets.push_back(t2);
    } 
    D.setFromTriplets(partialTriplets.begin(), partialTriplets.end());
  }
  return;
}
