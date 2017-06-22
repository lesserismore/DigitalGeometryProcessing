#include "vector_area_matrix.h"
#include <igl/boundary_facets.h>
#include <vector>

using namespace Eigen;
using namespace igl;

void vector_area_matrix(
  const Eigen::MatrixXi & F,
  Eigen::SparseMatrix<double>& A)
{
  
  int v = F.maxCoeff()+1; //# of vertices
  A.resize(2*v,2*v); //A sized for vectorized UV coordinates
  std::vector<Triplet<double>> areaTriplets;
  
  //get all boundary edges, including boundaries of interior holes
  MatrixXi bE;
  boundary_facets(F,bE);
  
  //For each boundary edge...
  for(int i = 0; i < bE.rows(); i++) {
    Vector2i edge = bE.row(i);
    
    //construct triplets such that A = B+B^T where U^TBU = 0.5*sum(det(edge[0] edge[1]))
    Triplet<double> t0 = Triplet<double>(edge[0],v+edge[1],0.5);
    Triplet<double> t1 = Triplet<double>(v+edge[1],edge[0],0.5);
    Triplet<double> t2 = Triplet<double>(edge[1],v+edge[0],-0.5);
    Triplet<double> t3 = Triplet<double>(v+edge[0],edge[1],-0.5);
    areaTriplets.push_back(t0);
    areaTriplets.push_back(t1);
    areaTriplets.push_back(t2);
    areaTriplets.push_back(t3);
  }

  A.setFromTriplets(areaTriplets.begin(),areaTriplets.end());
}

