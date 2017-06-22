#include "random_points_on_mesh.h"
#include <random>
#include <cmath>
#include <igl/cumsum.h>
#include <igl/doublearea.h>

using namespace Eigen;
using namespace igl;

void random_points_on_mesh(
  const int n,
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::MatrixXd & X)
{
  X.resize(n,3);

  int f = F.rows();

  //get areas of every face
  VectorXd areas = VectorXd(f);
  doublearea(V,F,areas);

  //sum up areas
  VectorXd areaCumSum = VectorXd(f);
  cumsum(areas,1,areaCumSum);

  //Compute each cumulative area as fraction of total area
  double surfaceArea = areaCumSum(f-1);
  areaCumSum *= (1/surfaceArea);

  //choose random number between 0 and 1
  std::random_device rd;
  std::mt19937 mt(rd());
  std::uniform_real_distribution<double> distribution(0.0,1.0);

  //use cumsum to select face
  for(int i = 0; i<n; i++) {
    double urn = distribution(mt);
    int L = 0;
    int R = f-1;
    int idx = 0;
    //binary search
    while(L <= R) {
      idx = std::floor((L+R)/2);
      double face_value = areaCumSum[idx];
      if(face_value < urn) {
        L = idx+1;
      } else if(face_value > urn) {
        R = idx-1;
      } else { //face_value is exactly urn, so idx is correct
        L = R+1;
      }
    }
    //get vertex indices
    int v0 = F(idx,0);
    int v1 = F(idx,1);
    int v2 = F(idx,2);

    //compute 2 more random numbers between 0 and 1
    double urn1 = distribution(mt);
    double urn2 = distribution(mt);

    //select point on chosen face
    if (urn1+urn2 > 1){
      //ensure random point will lie on face
      urn1 = 1-urn1;
      urn2 = 1-urn2;
    }
    //for each vertex in chosen face
    for (int j = 0; j<3; j++) {
      double x0 = V(v0,j);
      double x1 = V(v1,j);
      double x2 = V(v2,j);
      //insert into X
      X(i,j) = x0+urn1*(x1-x0)+urn2*(x2-x1); 
    }
  }
}

