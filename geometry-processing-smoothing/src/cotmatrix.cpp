#include "cotmatrix.h"
#include <vector>
#include <igl/doublearea.h>
#include <cmath>

using namespace Eigen;
using namespace igl;

void cotmatrix(
  const Eigen::MatrixXd & l,
  const Eigen::MatrixXi & F,
  Eigen::SparseMatrix<double> & L)
{
  int V = F.maxCoeff()+1;//# of vertices
  int f = F.rows();//# of faces
  L.resize(V,V);

  //vector to store triplet representation of L
  std::vector<Triplet<double>> LaplaceTriplets;

  //vector to store diagonal entries of L
  VectorXd LaplaceDiag = VectorXd::Zero(V);

  //For each face...
  for(int i = 0; i < f; i++){
    //get the vertex indices and edge lengths
    RowVector3i face = F.row(i);
	RowVector3d edge_lengths = l.row(i);
    int v0 = face[0];
    int v1 = face[1];
    int v2 = face[2];
    double l0 = edge_lengths[0]; //length of edge [v1,v2] opposite v0
    double l1 = edge_lengths[1]; //"            " [v2,v0] "      " v1
    double l2 = edge_lengths[2]; //"            " [v0,v1] "      " v2

    //Calculate the area
    double s = (0.5)*(l0+l1+l2); //semi-perimeter
    double area = sqrt(s*(s-l0)*(s-l1)*(s-l2));//Heron's area formula
	
    //Calculate half the cotangent of the angle at vertex v0,v1,v2, respectively
    double half_cot0 = (l1*l1+l2*l2-l0*l0)/(6.0*area);
    double half_cot1 = (l0*l0+l2*l2-l1*l1)/(6.0*area);
    double half_cot2 = (l0*l0+l1*l1-l2*l2)/(6.0*area);

    //Triplet t[i] has data half_cot[i] and indices corresponding to edge with length l[i]
    Triplet<double> t0 = Triplet<double>(v1,v2,half_cot0);
    Triplet<double> t1 = Triplet<double>(v2,v0,half_cot1);
    Triplet<double> t2 = Triplet<double>(v0,v1,half_cot2);

    //insert triplets
    LaplaceTriplets.push_back(t0);
    LaplaceTriplets.push_back(t1);
    LaplaceTriplets.push_back(t2);

	//transpose versions of above triplets to ensure symmetry of L
    Triplet<double> t0r = Triplet<double>(v2,v1,half_cot0);
    Triplet<double> t1r = Triplet<double>(v0,v2,half_cot1);
    Triplet<double> t2r = Triplet<double>(v1,v0,half_cot2);

    //insert triplets
    LaplaceTriplets.push_back(t0r);
    LaplaceTriplets.push_back(t1r);
    LaplaceTriplets.push_back(t2r);

    //Balance off-diagonal entries in L with diagonal entries, so sum of entries in each row is 0
    LaplaceDiag(v0) += (-1.0)*(half_cot1+half_cot2);
    LaplaceDiag(v1) += (-1.0)*(half_cot0+half_cot2);
    LaplaceDiag(v2) += (-1.0)*(half_cot0+half_cot1);
  }

  //construct and insert triplet representation of diagonal entries
  for(int i = 0; i < V; i++){
    Triplet<double> t = Triplet<double>(i,i,LaplaceDiag(i));
    LaplaceTriplets.push_back(t);
  }

  //construct sparse L from triplet vector
  L.setFromTriplets(LaplaceTriplets.begin(),LaplaceTriplets.end()); 
}

