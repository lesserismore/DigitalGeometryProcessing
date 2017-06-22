#include "edges.h"
#include <set> //std::set
#include <utility> //std::pair and make_pair

Eigen::MatrixXi edges(const Eigen::MatrixXi &F)
{
  Eigen::MatrixXi E;
  
  E = Eigen::MatrixXi::Zero(1,2); //intialize edge matrix with 2 columns for vertices
  int row_count = 0; //counts populated rows of E
  std::set<std::pair<int,int>> s; //set for storing edges as pairs of ints
  int r = F.rows(); //# of faces

  //for each face
  for (int i=0; i<r; i++ ) {
    //get the indices of the vertices of the face
    int v0 = F(i,0);
    int v1 = F(i,1);
    int v2 = F(i,2);
   
    //intialize the int pairs that will hold the edge data
    std::pair<int,int> e0;
    std::pair<int,int> e1;
    std::pair<int,int> e2;
 
    //assign adjacent vertices to pairs so that first<second
    //(ensures undirected edges can be reliably compared)
    if (v0<v1) {
      e0 = std::make_pair(v0,v1);
    }
    else {
      e0 = std::make_pair(v1,v0);
    }
    if (v1<v2) {
      e1 = std::make_pair(v1,v2);
    }
    else {
      e1 = std::make_pair(v2,v1);
    }
    if (v2<v0) {
      e2 = std::make_pair(v2,v0);
    }
    else {
      e2 = std::make_pair(v0,v2);
    }

    //check if each new edge is already in the set of edges,
    //and if it's not, add it to the set and insert the data in
    //new row of E
    auto edge_check_0 = s.find(e0);
    if (edge_check_0 == s.end()){
      E.conservativeResize(row_count+1,Eigen::NoChange); //make space for new edge data
      s.insert(e0);
      E(row_count,0) = e0.first;
      E(row_count,1)=e0.second;
      row_count++; //keep count of populated rows
    } 
    auto edge_check_1 = s.find(e1);
    if (edge_check_1 == s.end()){
      E.conservativeResize(row_count+1,Eigen::NoChange);
      s.insert(e1);
      E(row_count,0) = e1.first;
      E(row_count,1)=e1.second;
      row_count++; 
    } 
    auto edge_check_2 = s.find(e2);
    if (edge_check_2 == s.end()){
      E.conservativeResize(row_count+1,Eigen::NoChange);
      s.insert(e2);
      E(row_count,0) = e2.first;
      E(row_count,1)=e2.second;
      row_count++;
    } 
  }

  return E;
}
