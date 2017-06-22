#include "massmatrix.h"
#include <vector>
#include <cmath>

using namespace Eigen;

void massmatrix(
  const Eigen::MatrixXd & l,
  const Eigen::MatrixXi & F,
  Eigen::DiagonalMatrix<double,Eigen::Dynamic> & M)
{
  const int f = F.rows();//# of faces
  int V = F.maxCoeff()+1;//# of vertices

  //Vector to store diagonal entries of M
  VectorXd mass_diagonal = VectorXd::Zero(V);

  //For each face... 
  for(int i = 0; i < f; i++){
    //get the edge lengths
	RowVector3d edge_lengths = l.row(i);
    double l0 = edge_lengths[0];
    double l1 = edge_lengths[1];
    double l2 = edge_lengths[2];

    double s = (0.5)*(l0+l1+l2); //semi-perimeter of face
    double area = sqrt(s*(s-l0)*(s-l1)*(s-l2)); //Heron's area formula

    //get vector of vertex indices
    RowVector3i face = F.row(i);
    //For each vertex in the face...
	for(int v = 0; v < 3; v++){
      int vertex = face[v];
      //add the area of the face to the corresponding entry in mass_diagonal
      mass_diagonal(vertex) += area;
    }
  }

  //scale mass_diagonal by 1/3
  float one_third = 1.0/3.0;
  mass_diagonal *= one_third;

  //Construct dense DiagonalMatrix M with diagonal entries from mass_daigonal
  M = DiagonalMatrix<double,Dynamic>(mass_diagonal);
}
