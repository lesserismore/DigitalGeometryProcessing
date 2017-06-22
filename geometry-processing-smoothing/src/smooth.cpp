#include "smooth.h"
#include "cotmatrix.h"
#include "massmatrix.h"
#include <Eigen/Cholesky>
#include <igl/edge_lengths.h>

using namespace Eigen;

void smooth(
    const Eigen::MatrixXd & V,
    const Eigen::MatrixXi & F,
    const Eigen::MatrixXd & G,
    double lambda,
    Eigen::MatrixXd & U)
{
  //get edge lengths of mesh (V,F)
  MatrixXd l;
  igl::edge_lengths(V,F,l);

  //Get the discrete Laplacian L
  SparseMatrix<double> L;
  cotmatrix(l,F,L);

  //Get the mass matrix M
  DiagonalMatrix<double,Dynamic> M;
  massmatrix(l,F,M);

  //Convert M to a sparse matrix sparseM
  VectorXd M_diag = M.diagonal();
  std::vector<Triplet<double>> sparseM_triplets;
  int v = M_diag.rows(); //# of vertices
  for(int i = 0; i < v; i++){
    Triplet<double> mass_i = Triplet<double>(i,i,M_diag(i));
    sparseM_triplets.push_back(mass_i);
  }
  SparseMatrix<double> sparseM;
  sparseM.resize(v,v);
  sparseM.setFromTriplets(sparseM_triplets.begin(),sparseM_triplets.end());

  //solve the linear system (sparseM-lambda*L)*U = sparseM*G 
  SimplicialLDLT<SparseMatrix<double>> solver;
  solver.compute(sparseM-(L*lambda));
  if(solver.info()!=Eigen::Success) {
    // decomposition failed
    return;
  }
  U = solver.solve(M*G); //solve for updated (G.cols())-dim data
  if(solver.info()!=Eigen::Success) {
    // solving failed
    return;
  }
}
