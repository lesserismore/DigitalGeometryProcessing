#include "point_triangle_distance.h"
#include <Eigen/Dense>

using namespace Eigen;

void point_triangle_distance(
  const Eigen::RowVector3d & x,
  const Eigen::RowVector3d & a,
  const Eigen::RowVector3d & b,
  const Eigen::RowVector3d & c,
  double & d,
  Eigen::RowVector3d & p)
{
  //edges
  RowVector3d e1 = b-a;
  RowVector3d e2 = c-a;
  RowVector3d e3 = c-b;
  
  RowVector3d e1_n = e1.normalized();
  RowVector3d e2_n = e2.normalized();
  RowVector3d e3_n = e3.normalized();
  
  //Projection
  RowVector3d n = e1.cross(e2);
  RowVector3d n_n = n.normalized(); 
  RowVector3d dist = x-a;
  double t = dist.dot(n_n);
  RowVector3d proj = x - t*n_n;

  //Barycentric coordinates of projection relative to face
  double d0 = e1.dot(e1);
  double d1 = e1.dot(e2);
  double d2 = e2.dot(e2);
  double d3 = (proj-a).dot(e1);
  double d4 = (proj-a).dot(e2);
  double denom = d0*d2-d1*d1;
  double u = (d2*d3-d1*d4)/denom;
  double v = (d0*d4-d1*d3)/denom;
  double w = 1.0-u-v;

  //get signs of barycentric coordinates
  int u_sgn = (u > 0);
  int v_sgn = (v > 0);
  int w_sgn = (w > 0);
  int sgn_sum = u_sgn+v_sgn+w_sgn;
  
  if(sgn_sum == 3) {
    //projection lies on face
    p = proj;
  } else if(sgn_sum == 2) {
    //closest point to projection lies on an edge or a vertex
    double e_proj = 0;
    if(u < 0) { 
      RowVector3d y = proj-b;
      e_proj = y.dot(e3_n);
      if(e_proj < 0) {
        p = b;
      } else if (e_proj > e3.norm()) {
        p = c;
      } else {
        //closest point lies on edge e3
        p = e_proj*e3_n+b;
      }
    } else if(v < 0) {
      RowVector3d y = proj-a;
      e_proj = y.dot(e2_n);
      if(e_proj < 0) {
        p = a;
      } else if (e_proj > e2.norm()) {
        p = c;
      } else {
        //closest point lies on edge e2
        p = e_proj*e2_n+a;
      }
    } else {//w < 0
      RowVector3d y = proj-a;
      e_proj = y.dot(e3_n);
      if(e_proj < 0) {
        p = a;
      } else if (e_proj > e1.norm()) {
        p = b;
      } else {
        //closest point lies on edge e1
        p = e_proj*e1_n+a;
      }
    }
  } else { //sgn_sum == 1
    //closest point lies on a vertex
    if(u > 0) {
      p = a;
    } else if(v > 0) {
      p = b;
    } else { //w > 0
      p = c;
    }
  }
  d = (x-p).norm();
}
