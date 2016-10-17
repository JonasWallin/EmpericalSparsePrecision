
#ifndef ETCHOLMOD
#define ETCHOLMOD

#include <Rcpp.h>
#include <Eigen/SparseCore>
#include <Eigen/SparseCholesky>
#include "SuiteSparse_config.h"
#include "cholmod.h"
#include "cholmod_cholesky.h"
#include "cholmod_core.h"

 namespace Eigen {

 namespace internal {

 template<typename Scalar, typename CholmodType>
 void cholmod_configure_matrix(CholmodType& mat)
 {
   if (internal::is_same<Scalar,float>::value)
   {
     mat.xtype = CHOLMOD_REAL;
     mat.dtype = CHOLMOD_SINGLE;
   }
   else if (internal::is_same<Scalar,double>::value)
   {
     mat.xtype = CHOLMOD_REAL;
     mat.dtype = CHOLMOD_DOUBLE;
   }
   else if (internal::is_same<Scalar,std::complex<float> >::value)
   {
     mat.xtype = CHOLMOD_COMPLEX;
     mat.dtype = CHOLMOD_SINGLE;
   }
   else if (internal::is_same<Scalar,std::complex<double> >::value)
   {
     mat.xtype = CHOLMOD_COMPLEX;
     mat.dtype = CHOLMOD_DOUBLE;
   }
   else
   {
     eigen_assert(false && "Scalar type not supported by CHOLMOD");
   }
 }

#endif
 } // namespace internal

}

 /** Wraps the Eigen sparse matrix \a mat into a Cholmod sparse matrix object.
  * Note that the data are shared.
  */
 template<typename _Scalar, int _Options, typename _Index>
 cholmod_sparse viewAsCholmod(Eigen::SparseMatrix<_Scalar,_Options,_Index>& mat)
 {
   cholmod_sparse res;
   res.nzmax   = mat.nonZeros();
   res.nrow    = mat.rows();;
   res.ncol    = mat.cols();
   res.p       = mat.outerIndexPtr();
   res.i       = mat.innerIndexPtr();
   res.x       = mat.valuePtr();
   res.z       = 0;
   res.sorted  = 1;
   if(mat.isCompressed())
   {
     res.packed  = 1;
     res.nz = 0;
   }
   else
   {
     res.packed  = 0;
     res.nz = mat.innerNonZeroPtr();
   }

   res.dtype   = 0;
   res.stype   = -1;

   if (Eigen::internal::is_same<_Index,int>::value)
   {
     res.itype = CHOLMOD_INT;
   }
   else if (Eigen::internal::is_same<_Index,UF_long>::value)
   {
     res.itype = CHOLMOD_LONG;
   }
   else
   {
     eigen_assert(false && "Index type not supported yet");
   }

   // setup res.xtype
   Eigen::internal::cholmod_configure_matrix<_Scalar>(res);

   res.stype = 0;

   return res;
 }

 template<typename _Scalar, int _Options, typename _Index>
 const cholmod_sparse viewAsCholmod(const Eigen::SparseMatrix<_Scalar,_Options,_Index>& mat)
 {
   cholmod_sparse res = viewAsCholmod(mat.const_cast_derived());
   return res;
 }

 /** Returns a view of the Eigen sparse matrix \a mat as Cholmod sparse matrix.
 * The data are not copied but shared. */
 template<typename _Scalar, int _Options, typename _Index, unsigned int UpLo>
 cholmod_sparse viewAsCholmod(const Eigen::SparseSelfAdjointView<Eigen::SparseMatrix<_Scalar,_Options,_Index>, UpLo>& mat)
 {
   cholmod_sparse res = viewAsCholmod(mat.matrix().const_cast_derived());

   if(UpLo==Eigen::Upper) res.stype =  1;
   if(UpLo==Eigen::Lower) res.stype = -1;

   return res;
 }


 /** Returns a view of the Cholmod sparse matrix \a cm as an Eigen sparse matrix.
  * The data are not copied but shared. */
template<typename Scalar, int Flags, typename Index>
Eigen::MappedSparseMatrix<Scalar,Flags,Index> viewAsEigen(cholmod_sparse& cm)
{
   return Eigen::MappedSparseMatrix<Scalar,Flags,Index>
   (cm.nrow, cm.ncol, static_cast<Index*>(cm.p)[cm.ncol],
    static_cast<Index*>(cm.p), static_cast<Index*>(cm.i),static_cast<Scalar*>(cm.x) );
}
