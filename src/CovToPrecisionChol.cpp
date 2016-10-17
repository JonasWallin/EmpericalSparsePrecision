#include <Rcpp.h>
#include <Eigen/SparseCore>
#include <Eigen/SparseCholesky>
#include "SuiteSparse_config.h"
#include "cholmod.h"
#include "cholmod_cholesky.h"
#include "cholmod_core.h"
#include <Eigen/Cholesky>
#include "EigenToCholmod.h"
#include <vector>
using namespace Rcpp;


typedef Eigen::Triplet<double> T;
//'
//' @title CovToPrecisionCholcpp
//' @keywords internal
//' @export
// [[Rcpp::export]]
Eigen::SparseMatrix<double,0,int> CovToPrecisionChol_cpp(Eigen::SparseMatrix<double,0,int> & L,
                                                         Eigen::MatrixXd Sigma) {

  int d = Sigma.rows();
  Eigen::VectorXd D(d);
  std::vector<T> tripletList;
  for(int i =0; i < d; i++)
    tripletList.push_back(T(i, i, 1.));

  for (int k=0; k<L.outerSize(); ++k){
    int A_size = L.innerVector(k).nonZeros() - 1;
    Eigen::MatrixXd A(A_size, A_size);
    Eigen::VectorXd b(A_size);
    std::vector<int> index(A_size);
    int count = 0;
    for (Eigen::SparseMatrix<double,0,int>::InnerIterator it(L, k); it; ++it)
    {
      if(it.row() != k)
        index[count++] = it.row();
    }
    for(int i =0; i < A_size; i++)
    {
      A(i, i) = Sigma(index[i], index[i]);
      b(i)    = Sigma(index[i], k);
      for(int ii = 0; ii < i ; ii++)
      {
        A(ii, i) = Sigma(index[i],index[ii]) ;
        A(i, ii) = Sigma(index[i],index[ii]) ;
      }

    }
    Eigen::VectorXd Ainvb =  - A.llt().solve(b);
    D(k) = Sigma(k, k ) + b.dot( Ainvb);
    for(int i =0; i < A_size; i++)
      tripletList.push_back(T(index[i], k, Ainvb[i]));

  }
  D = D.cwiseInverse();
  D = D.array().sqrt();
  Eigen::SparseMatrix<double,0,int> T_mat;
  T_mat.resize(d, d);
  T_mat.setFromTriplets(tripletList.begin(), tripletList.end());

  return T_mat * D.asDiagonal();
}

