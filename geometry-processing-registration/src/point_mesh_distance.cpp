#include "point_mesh_distance.h"
#include "point_triangle_distance.h"
#include <igl/per_face_normals.h>
#include <cmath>

using namespace Eigen;
using namespace igl;

void point_mesh_distance(
  const Eigen::MatrixXd & X,
  const Eigen::MatrixXd & VY,
  const Eigen::MatrixXi & FY,
  Eigen::VectorXd & D,
  Eigen::MatrixXd & P,
  Eigen::MatrixXd & N)
{
  //Brute Force implementation
  P.resizeLike(X);
  D.resize(X.rows());
  N.resizeLike(X);
 
  //Calculate all normals 
  MatrixXd A;
  per_face_normals(VY,FY,A);
  
  //for each query point x
  for(int i = 0; i < X.rows(); i++) {
    RowVector3d x = X.row(i); //get query point
    RowVector3d p0 = RowVector3d::Zero(); //current closest point
    double d0 = std::numeric_limits<double>::infinity(); //current distance to closest point
    int f_idx = 0;
   
    //for each face
    for(int j = 0; j < FY.rows(); j++) {
      RowVector3i face = FY.row(j); //get face
      RowVector3d a = VY.row(face[0]);
      RowVector3d b = VY.row(face[1]);
      RowVector3d c = VY.row(face[2]);
      double d = 0;
      RowVector3d p = RowVector3d::Zero();
      point_triangle_distance(x,a,b,c,d,p); //find closest point to x on face
      //if p is closer to x than p0
      if (d<d0) {
        d0 = d;
        p0 = p;
        f_idx = j;
      }
    }
    //store the closest point's information
    D[i] = d0;
    P.row(i) = p0;
    N.row(i) = A.row(f_idx);
  }
}
