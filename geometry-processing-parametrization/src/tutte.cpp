#include "tutte.h"
#include <igl/boundary_loop.h>
#include <igl/map_vertices_to_circle.h>
#include <igl/min_quad_with_fixed.h>
#include <vector>

using namespace igl;
using namespace Eigen;

void tutte(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::MatrixXd & U)
{
  //get boundary vertex indices
  VectorXi bnd;
  boundary_loop(F,bnd);

  //map boundary vertices onto unit disc
  MatrixXd bndUV;
  map_vertices_to_circle(V,bnd,bndUV);

  //Calculate discrete Laplacian L
  SparseMatrix<double> L;
  int v = F.maxCoeff()+1;//# of vertices
  int f = F.rows();//# of faces
  L.resize(v,v);

  //vector to store triplet representation of L
  std::vector<Triplet<double>> Triplets;

  //vector to store diagonal entries of L
  VectorXd Diag = VectorXd::Zero(v);

  //For each face...
  for(int i = 0; i < f; i++){
    //get the vertex indices
    RowVector3i face = F.row(i);
    int v0 = face[0];
    int v1 = face[1];
    int v2 = face[2];

    //Triplet t[i] has indices corresponding to edge opposite vertex v[i]
    Triplet<double> t0 = Triplet<double>(v1,v2,1);
    Triplet<double> t1 = Triplet<double>(v2,v0,1);
    Triplet<double> t2 = Triplet<double>(v0,v1,1);

    //insert triplets
    Triplets.push_back(t0);
    Triplets.push_back(t1);
    Triplets.push_back(t2);

	//transpose versions of above triplets to ensure symmetry of L
    Triplet<double> t0r = Triplet<double>(v2,v1,1);
    Triplet<double> t1r = Triplet<double>(v0,v2,1);
    Triplet<double> t2r = Triplet<double>(v1,v0,1);

    //insert transpose triplets
    Triplets.push_back(t0r);
    Triplets.push_back(t1r);
    Triplets.push_back(t2r);

    //Balance off-diagonal entries in L with diagonal entries, so sum of entries in each row is 0
    Diag(v0) += -2;
    Diag(v1) += -2;
    Diag(v2) += -2;
  }

  //construct and insert triplet representation of diagonal entries
  for(int i = 0; i < v; i++){
    Triplet<double> t = Triplet<double>(i,i,Diag(i));
    Triplets.push_back(t);
  }
  L.setFromTriplets(Triplets.begin(),Triplets.end());
    
  //minimize 0.5*trace((U^T)LU) subject to U(bnd,:) = bndUV
  U.resize(v,2);
  MatrixXd B; //dummy linear term
  B.resize(v,1); //match size of U column
  B.setZero();
  VectorXd Beq = VectorXd::Zero(0); //dummy linear condition term
  SparseMatrix<double> Aeq(0,0);//dummy linear condition term
  min_quad_with_fixed_data<double> data;
  min_quad_with_fixed_precompute(L,bnd,Aeq,false,data);
  MatrixXd U0,U1;
  min_quad_with_fixed_solve(data,B,bndUV.col(0),Beq,U0);
  min_quad_with_fixed_solve(data,B,bndUV.col(1),Beq,U1);
  U << U0,U1;
}
