#include <Rcpp.h>
#include <Eigen/SparseCore>
#include <Eigen/SparseCholesky>
#include "SuiteSparse_config.h"
#include "cholmod.h"
#include "cholmod_cholesky.h"
#include "cholmod_core.h"
#include "EigenToCholmod.h"

//' @title symbolicCholeskyFactorcpp
//' @keywords internal
//' @export
// [[Rcpp::export]]
Eigen::SparseMatrix<double,0,int> symbolicCholeskyFactor_cpp(Eigen::SparseMatrix<double,0,int> & K) {


  cholmod_common c ;
  cholmod_start (&c) ; /* start CHOLMOD */
  cholmod_sparse K_cholmod = viewAsCholmod(K);

  K_cholmod.stype = -1;
  int n = K_cholmod.nrow ;
  int   *Parent,
         *Post,
         *ColCount,
         *First,
         *Level;

  Parent   = (int*) cholmod_malloc(n, sizeof(int), &c );
  Post     = (int*) cholmod_malloc(n, sizeof (int), &c) ;
  ColCount = (int*) cholmod_malloc(n, sizeof (int), &c) ;
  First    = (int*) cholmod_malloc(n, sizeof (int), &c) ;
  Level    = (int*) cholmod_malloc(n, sizeof (int), &c) ;

  cholmod_sparse *F = cholmod_transpose (&K_cholmod, 0, &c) ;
  cholmod_etree (F, Parent, &c) ;
  cholmod_postorder (Parent, n, NULL, Post, &c);
  cholmod_rowcolcounts (&K_cholmod, NULL, 0, Parent, Post, NULL, ColCount,
                          First, Level, &c) ;
  /* count the total number of entries in L */
  double lnz = 0 ;
  int j;
  for ( j = 0 ; j < n ; j++)
    lnz += ColCount [j] ;

  cholmod_sparse* L = cholmod_allocate_sparse (n, n, lnz, TRUE, TRUE, 0,
                                 CHOLMOD_PATTERN, &c) ;
  int *Lp = (int*) L->p ;
  int *Li = (int*) L->i ;
  /* initialize column pointers */
  lnz = 0 ;

  for (j = 0 ; j < n ; j++)
  {
    Lp [j] = lnz ;
    lnz += ColCount [j] ;
  }
  Lp [j] = lnz ;
  /* create a copy of the column pointers */
  int *W = First ;
  for (j = 0 ; j < n ; j++)
  {
    W [j] = Lp [j] ;
  }

  /* get workspace for computing one row of L */
  cholmod_sparse *R = cholmod_allocate_sparse (n, 1, n, FALSE, TRUE, 0, CHOLMOD_PATTERN,
                                 &c) ;
  int *Rp = (int *) R->p ;
  int *Ri = (int *) R->i ;
  /* compute L one row at a time */
  for (int k = 0 ; k < n ; k++)
  {
    /* get the kth row of L and store in the columns of L */
    cholmod_row_subtree (F, NULL, k, Parent, R, &c) ;
    for (int p = 0 ; p < Rp [1] ; p++)
      Li [W [Ri [p]]++] = k ;

    /* add the diagonal entry */
    Li [W [k]++] = k ;
  }
  cholmod_free_sparse (&R, &c);
  L->x = cholmod_malloc (lnz, sizeof (double), &c) ;
  double* Lx = (double*) L->x ;
  for (int p = 0 ; p < lnz ; p++)
    Lx [p] = 1 ;

  L->xtype = CHOLMOD_REAL ;

  Eigen::SparseMatrix<double,0,int> K2 = viewAsEigen<double,0,int>(*L);
  cholmod_free (n, sizeof (int), Parent, &c) ;
  cholmod_free (n, sizeof (int), Post, &c) ;
  cholmod_free (n, sizeof (int), ColCount, &c) ;
  cholmod_free (n, sizeof (int), First, &c) ;
  cholmod_free (n, sizeof (int), Level, &c) ;
  cholmod_free_sparse (&F, &c);
  cholmod_finish(&c);
  return K2;
}
