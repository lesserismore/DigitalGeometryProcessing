#include "euler_characteristic.h"
#include "edges.h"

int euler_characteristic(const Eigen::MatrixXi &F)
{
  int Chi;
  
  int f = F.rows(); //# of faces
  int e; //# of edges
  int v; //# of vertices

  Eigen::MatrixXi E = edges(F); 
  e = E.rows();

  Eigen::VectorXi V;
  //v = F.colwise().maxCoeff(); //Vector containing maximum index in each column of F
  //v = V.maxCoeff()+1; //(maximum index in F)+1 = #rows in the vertex matrix
  v = F.maxCoeff()+1;
  //famous equation for Euler characteristic, for well-defined polygons 
  Chi = f - e + v;
  return Chi;
}
