#include <Rcpp.h>
#include <vector>
using namespace Rcpp;
#include "EigenToCholmod.h"

//' @title reordercpp
//' @keywords internal
//' @export
// [[Rcpp::export]]
Rcpp::List reorder_cpp(Eigen::SparseMatrix<double,0,int> & K) {
  cholmod_common c ;
  cholmod_start (&c) ; /* start CHOLMOD */
  cholmod_sparse K_cholmod = viewAsCholmod(K);
  K_cholmod.stype = -1;
  cholmod_factor *L = cholmod_analyze (&K_cholmod, &c) ;
  std::vector<int> Perm((int)L->n), IPerm((int)L->n);
  for(int i = 0; i < L->n; i++){
    Perm[i]                     = ((int* )L->Perm)[i] + 1;
    IPerm[((int* )L->Perm)[i]]  = i + 1;
  }

  Rcpp::List z;
  z["Perm"]   = Perm;
  z["IPerm"]  = IPerm;
  return z;
}


