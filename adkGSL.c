/*
/
/ adkGSL.c
/
/ homebrewed addons to gsl */

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_cblas.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
//#include "../extern/include/f2c.h"
//#include "../extern/include/blaswrap.h"
//#include <clapack.h>
//#include <vecLib/cblas.h>
//#include <accelerate/accelerate.h>
#include "adkGSL.h"

#define SUM_LOG_THRESHOLD  -10


//gsl_matrix_covariance-- computes sample covariance matrix and stores it in matrix cov
void gsl_matrix_covariance(gsl_matrix *data, gsl_matrix *cov){
	gsl_vector_view a, b;
	size_t i, j;
	double v;
	for (i = 0; i < data->size2; i++) {
 		for (j = 0; j < data->size2; j++) {
 			a = gsl_matrix_column (data, i);
 			b = gsl_matrix_column (data, j);
 			v = gsl_stats_covariance (a.vector.data, a.vector.stride,
 				b.vector.data, b.vector.stride, a.vector.size);
 			gsl_matrix_set (cov, i, j, v);
		}
	}
}


/*gsl_vector_includes() - returns 1 or 0 depending on result */
int gsl_vector_includes(gsl_vector *aVec, double aValue){
  int i;
  
  for(i = 0; i < aVec->size; i++){
    if (gsl_vector_get(aVec, i) == aValue){
      return(1);
    }
  }
  return(0);
}


void gsl_matrix_prettyPrint(gsl_matrix *m){
	int i, j;
	for(i=0;i<m->size1;i++){
		for(j=0;j<m->size2;j++){
			printf("%.6f ",gsl_matrix_get(m,i,j));
		}
		printf("\n");
	}
	printf("///////////\n");
}

double gsl_vector_sum(gsl_vector *aVec, int aVecSize){
	double sum;
	int i;
	
	sum = 0;
	for (i = 0; i < aVecSize; i++){
		sum += gsl_vector_get(aVec, i);
		}
	return sum;
}

double gsl_vector_dot_product(gsl_vector *vec1, gsl_vector *vec2){
  int i;
  double dot = 0.0;

  assert(vec1->size == vec2->size);
  for(i = 0; i < vec1->size; i++){
    dot += gsl_vector_get(vec1, i) *  gsl_vector_get(vec2, i);
  }
  return(dot);
}

void gsl_vector_outer_product(gsl_vector *vec1, gsl_vector *vec2, gsl_matrix *result){
	int i,j;
	
	assert(vec1->size == result->size1 || vec2->size == result->size2);
	for(i = 0; i < vec1->size; i++){
		for(j = 0; j < vec2->size; j++){
			gsl_matrix_set(result,i,j,gsl_vector_get(vec1,i) * gsl_vector_get(vec2,j));
		}
	}
}


double gsl_matrix_row_sum(gsl_matrix *mat, int iRow, int rowSize){
	gsl_vector_view aRow;
	double sum = 0;
	int i;
	
	aRow = gsl_matrix_row(mat, iRow);
	for (i = 0; i< rowSize; i++){
		sum += gsl_vector_get(&aRow.vector, i);
		}
	return sum;
}

//gsl_matrix_row_max -- returns the max value of a row
double gsl_matrix_row_max(gsl_matrix *mat, int iRow){
	double max = 0;
	double tmp;
	int i;
	
	for (i = 0; i< mat->size2; i++){
		tmp = gsl_matrix_get(mat,iRow,i);
		if(max < tmp) max = tmp ;
	}
	return max;
}

//gsl_matrix_sum -- returns sum of all values in matrix
double gsl_matrix_sum(gsl_matrix *mat){
	double sum = 0;
	int i, j;
	
	for (i = 0; i< mat->size1; i++){
		for (j = 0; j< mat->size2; j++){
		sum += gsl_matrix_get(mat, i,j);
		}
	}
	return sum;
}

//gsl_matrix_bootstrap-- fill matrix boot with samples with replacement from orignal
// needs an initialize random number generator
void gsl_matrix_bootstrap(gsl_matrix *orig, gsl_matrix *boot, gsl_rng *rng){
	
	size_t i, j, k, K;
	K = (size_t) orig->size1 * orig->size2;
	double probs[K], sum;
	int randM[K];
	int M;
	
	sum = gsl_matrix_sum(orig);
	M = (int) sum;
	k = 0;
	for (i = 0; i < orig->size1; i++) {
 		for (j = 0; j < orig->size2; j++) {
			probs[k++] = gsl_matrix_get(orig,i,j) / M;
		}
	}
	gsl_ran_multinomial (rng,K,M, probs, randM);
	k = 0;
	for (i = 0; i < orig->size1; i++) {
 		for (j = 0; j < orig->size2; j++) {
			gsl_matrix_set(boot,i,j,randM[k++]);
		}
	}
}
/* efficiently compute log of sum of values, which themselves are
   stored as logs: that is, return log(sum_i exp(l_i)).  The largest
   of the elements of l (call it maxval) is factored out, so that
   log(sum_i(exp(l_i))) = maxval + log(1 + sum_i(exp(l_i-maxval))),
   where the new sum is taken over 2 <= i < n.  All of the quantities
   in the exp must be negative, and those smaller than some reasonable
   threshold can be ignored. [Thanks to David Haussler for showing me
   this trick].  This code was adapted from code kindly given to me by Adam Siepel */

double log_sum(gsl_vector *vec) {
  double maxval, expsum;
  int k;

  if (vec->size > 1){
  	gsl_sort_vector(vec);
	gsl_vector_reverse(vec);
	}
	
  maxval = gsl_vector_get(vec, 0);
  expsum = 1;
  k = 1;

  while (k < vec->size && gsl_vector_get(vec, k) - maxval > SUM_LOG_THRESHOLD){
    expsum += exp(gsl_vector_get(vec, k++) - maxval);
	}
  return maxval + log(expsum);
}


/*the following methods are wrappers for use with LAPACK */

void set_lapack_entry(double *A, int i, int j, int nrows, double val){
  A[j*nrows+i] = val;
}

double get_lapack_entry(double *a, int i, int j,int nrows){
  printf("ha");
  return a[j*nrows+i];
}

/*this is meant to be a conversion utility */
double *gsl_matrix_2_lapack(gsl_matrix *m){
  double *lapack_mat;
  int i, j;

  lapack_mat = malloc(m->size1 * m->size2 * sizeof(double));
  for(i = 0; i < m->size1; i++){
    for(j = 0; j < m->size2; j++){
      set_lapack_entry(lapack_mat, i, j,m->size1, gsl_matrix_get(m, i, j));
    }
  }
  return lapack_mat;
}

gsl_matrix *lapack_2_gsl_matrix(double *A, int nrows, int ncols){
  int i, j;
  gsl_matrix *new;
  
  new = gsl_matrix_alloc(nrows,ncols);
  for(i=0;i<nrows;i++){
    for(j=0;j<ncols;j++){
      gsl_matrix_set(new,i,j,get_lapack_entry(A,i,j,nrows));
    }
  }
  return new;
}

// gsl_matrix *gsl_matrix_power_logs(gsl_matrix *aMatrix, int power){
//   gsl_matrix *powerMat, *logMat;
//   int i, j;
//   double x;
// 
//   
//   powerMat = gsl_matrix_power(aMatrix, power);
//   logMat = gsl_matrix_alloc(aMatrix->size1, aMatrix->size2);
//   for(i = 0; i < aMatrix->size1; i++){
//     for(j = 0; j < aMatrix->size2; j++){
//       x = gsl_matrix_get(powerMat, i, j);
//       gsl_matrix_set(logMat, i, j, log(x));
//     }
//   }
//   gsl_matrix_free(powerMat);
//   return(logMat);
// }
  
// gsl_matrix *gsl_matrix_power(gsl_matrix *aMatrix, int power){
//   gsl_matrix *revect, *levect;
//   double *dataMat, *wr, *wi, *vl, *vr, *aMat, *bMat, *cMat, *dMat, *work;
//   int lwork, i, j, n ;
//   assert(aMatrix->size1 == aMatrix->size2);
//   n = aMatrix->size1;
// 
//   if (power == 1){
//     gsl_matrix * solMat = gsl_matrix_alloc(n, n);
//     for(i=0; i < n; i++){
//       for(j = 0; j< n; j++){
// 	gsl_matrix_set(solMat, i, j, gsl_matrix_get(aMatrix,i,j));
//       }
//     }
//     return(solMat);
//   }
// 
//   
//   //allocate space for eigenvector matrices, and other result holders
//   levect = gsl_matrix_alloc(n, n);
//   revect = gsl_matrix_alloc(n, n);
//   wr = malloc(n * sizeof(double));
//   wi =  malloc(n * sizeof(double));
//   vl = malloc(n * n * sizeof(double));
//   vr = malloc(n * n * sizeof(double));
//   lwork = 20 * n;
//   work = malloc(lwork * sizeof(double));
//  
//   //move data to LAPACK (col major) format
//   dataMat =  gsl_matrix_2_lapack(aMatrix);
//  
//   //run LAPACK routine dgeev
// //  #ifdef OSX
//   dgeev('V','V',n, dataMat, n, wr, wi, vl, n, vr, n, work, lwork);
//   //#else
//   //dgeev_('V','V',n, dataMat, n, wr, wi, vl, n, vr, n, work, lwork, info);
//   //#endif
//   //collect eigenvalues into matrices
//   for(i=0; i < n; i++){
//     for(j = 0; j< n; j++){
//       gsl_matrix_set(revect, i, j, vr[j*2+i]);
//       gsl_matrix_set(levect, j, i, vl[j*2+i]);
//     }
//   }
//   
//   /* rescale such that revect and levect are inverses */
//   /* see Press et al. (Numerical Recipes) for a very clear explanation
//      of the relationship between the left and right eigenvectors  */
//    for(i = 0; i < n; i++){
//     double dotProd = 0.0;
//     /* compute dot product of row i in levect and col i in revect */
//     for(j = 0; j< n; j++){
//       dotProd += gsl_matrix_get(levect, i, j) * gsl_matrix_get(revect,j,i); 
//       }
//     
//       /* now rescale levect */
//     for(j = 0; j< n; j++){
//       double old = gsl_matrix_get(levect, i, j);
//       double scaled = old / dotProd;
//       gsl_matrix_set(levect, i, j, scaled);
//       } 
//   }
// 
//    //prepare for matrix multiplication
//    gsl_matrix * solMat = gsl_matrix_alloc(n, n);
//    gsl_matrix * diag = gsl_matrix_alloc(n, n);
//    gsl_matrix_set_all(diag, 0);
//    for(j = 0; j< n; j++){
//      gsl_matrix_set(diag, j,j, pow(wr[j],power));
//    }
//    
//    aMat = gsl_matrix_2_lapack(diag);
//    bMat = gsl_matrix_2_lapack(revect);
//    cMat = malloc(n * n * sizeof(double));
// 
//    //compute a^t = s^-1 * diag^t * s
//    //cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, n, n, n, 1.0, aMat, n, bMat, n, 0.0, cMat, n);
//    free(bMat);
// 
//    bMat = gsl_matrix_2_lapack(levect);
//    dMat =  malloc(n * n * sizeof(double));
//    //cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, n, n, n, 1, cMat, n, bMat, n, 0.0, dMat, n); 
// 
//    //put results into gsl_matrix
//    for(i = 0; i < n; i++){
//      for(j = 0; j < n; j++){
//        gsl_matrix_set(solMat, i, j, dMat[j*n+i]);
//      }
//    }
//   
//    //cleanup
//    gsl_matrix_free(levect);
//    gsl_matrix_free(revect);
//    gsl_matrix_free(diag);
//    free(work);
//    free(vl);
//    free(vr);
//    free(wi);
//    free(wr);
//    free(aMat);
//    free(bMat);
//    free(cMat);
//    free(dMat);
//    free(dataMat);
//    return(solMat);
// }

//returns the trace of a matrix-- sum of diag
double gsl_matrix_trace(gsl_matrix *x){
	int i;
	double sum = 0.0;
	for(i=0;i<x->size1;i++){
		sum += gsl_matrix_get(x,i,i);
	}
	return(sum);
}


//copies the lower triangle of src in to dest, sets others to zero
void gsl_matrix_lower_tri_copy(gsl_matrix *src, gsl_matrix *dest){
	int i, j;
	assert(src->size1 == dest->size1);
	assert(src->size2 == dest->size2);
	for(i = 0; i < src->size1; i++){
		for(j = 0; j < src->size2; j++){
			if(i<=j){
				gsl_matrix_set(dest,i,j,gsl_matrix_get(src,i,j));
			}
			else{
				gsl_matrix_set(dest,i,j,0.0);
			}
		}
	}
}

//copies the upper triangle of src in to dest, sets others to zero
void gsl_matrix_upper_tri_copy(gsl_matrix *src, gsl_matrix *dest){
	int i, j;
	assert(src->size1 == dest->size1);
	assert(src->size2 == dest->size2);
	for(i = 0; i < src->size1 ;i++){
		for(j = 0; j < src->size2; j++){
			if(i>=j){
				gsl_matrix_set(dest,i,j,gsl_matrix_get(src,i,j));
			}
			else{
				gsl_matrix_set(dest,i,j,0.0);
			}
		}
	}
}

//zeros out everything but lower tri-- in place
void gsl_matrix_lower_tri(gsl_matrix *x){
	int i, j;
	for(i = 0; i < x->size1; i++){
		for(j = 0; j < x->size2; j++){
			if(i>j)	gsl_matrix_set(x,i,j,0.0);
		}
	}
}

//zeros out everything but upper tri-- in place
void gsl_matrix_upper_tri(gsl_matrix *x){
		int i, j;
		for(i = 0; i < x->size1; i++){
			for(j = 0; j < x->size2; j++){
				if(i<j)	gsl_matrix_set(x,i,j,0.0);
			}
		}
}

//fillMatrixFromArray-- fills up a matrix based on nrow and ncol given *double
void fillMatrixFromArray(double *numbers, gsl_matrix *dest, int nrow, int ncol){
	int i, j, count = 0;
	assert(nrow == dest->size1);
	for(i = 0; i < nrow; i++){
		for(j=0;j<ncol;j++){
			gsl_matrix_set(dest,i,j,numbers[count++]);
		}
	}
}

//fillMatrixFromArray-- fills up a matrix based on nrow and ncol given gsl_vector
void fillMatrixFromVector(gsl_vector *numbers, gsl_matrix *dest, int nrow, int ncol){
	int i, j, count = 0;
	assert(nrow == dest->size1);
	for(i = 0; i < nrow; i++){
		for(j=0;j<ncol;j++){
			gsl_matrix_set(dest,i,j,gsl_vector_get(numbers,count++));
		}
	}
}

//fillArrayFromMatrix-- fills up a vector  based gsl_matrix
void fillArrayFromMatrix(gsl_matrix *src, gsl_vector *dest){
	int i, j, count = 0;
	
	for(i = 0; i < src->size1; i++){
		for(j=0;j<src->size2;j++){
			gsl_vector_set(dest,count++, gsl_matrix_get(src,i,j));
		}
	}
}



//fillMatrixFromCholArray-- fills up a matrix based on nrow and ncol given arraay
// representing the Cholesky decomp
void fillMatrixFromCholArray(double *numbers, gsl_matrix *dest, int nrow, int ncol){
	int i, j, count = 0;

	for(i = 0; i < nrow; i++){
		for(j=0;j<ncol;j++){
			if(i<=j)gsl_matrix_set(dest,i,j,numbers[count++]);
		}
	}
}

void fillCholArrayFromMatrix(gsl_matrix *src, gsl_vector *dest){
	int i, j, count = 0;
	
	for(i = 0; i < src->size1; i++){
		for(j=0;j<src->size2;j++){
			if(i<=j){
				 gsl_vector_set(dest,count++, gsl_matrix_get(src,i,j));
				//printf("i:%d j: %d val: %f\n",i,j,gsl_matrix_get(src,i,j) );
			}
			
		}
	}
	
}








//////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef OSX

/*some wrappers for use with the OS X vecLib LAPACK implementation */
/*wrapper for the lapack eigensystem routine dgeev_*/
/*
int dgeev(char jobvl, char jobvr, int n,  double *a, int lda, double *wr, double *wi, double *vl, 
	  int ldvl,double *vr, int ldvr, double *work, int lwork){
 
  __CLPK_integer _n, _lda, _ldvl, _ldvr, _lwork, info;
  
  _n = (__CLPK_integer) n;
  _lda =  (__CLPK_integer) lda;
  _ldvl =(__CLPK_integer) ldvl;
  _ldvr =(__CLPK_integer) ldvr;
  _lwork =(__CLPK_integer) lwork;
  
  dgeev_(&jobvl, &jobvr, &_n, a, &_lda, wr, wi, vl, &_ldvl, vr, &_ldvr, work, &_lwork, &info);
  return info;
} */

/*wrapper for the lapack eigensystem routine dsyev_ 
int dsyev(char jobvl, char uplo, int n,  double *a, int lda, double *w, double *work, int lwork){
 
  __CLPK_integer _n, _lda, _lwork, info;
  
  _n = (__CLPK_integer) n;
  _lda =  (__CLPK_integer) lda;
  _lwork =(__CLPK_integer) lwork;
  
  dsyev_(&jobvl, &uplo, &_n, a, &_lda, w, work, &_lwork, &info);
  return info;
}

int dgehrd(int n, int ilo, int ihi, double *a, int lda, double *tau, double *work, int lwork){
  __CLPK_integer _n, _ilo, _ihi, _lda,  _lwork, info;

   _n = (__CLPK_integer) n;
  _lda =  (__CLPK_integer) lda;
  _lwork =(__CLPK_integer) lwork;
  _ilo =(__CLPK_integer) ilo;
  _ihi =(__CLPK_integer) ihi;
  
  dgehrd_(&_n, &_ilo, &_ihi, a, &_lda, tau, work, &_lwork, &info);
  return info;
}

int dgebal(char job, int n,  double *a, int lda, int ilo, int ihi, double *scale){
  __CLPK_integer _n, _ilo, _ihi, _lda, info;

   _n = (__CLPK_integer) n;
  _lda =  (__CLPK_integer) lda;
  _ilo =(__CLPK_integer) ilo;
  _ihi =(__CLPK_integer) ihi;
  
  dgebal_(&job, &_n, a, &_lda, &_ilo, &_ihi, scale, &info);
  return info;
}

int dhseqr(char job, char compz, int n, int ilo, int ihi,double *h, int ldh, double *wr, double *wi, double *z, int ldz, double *work, int lwork){
  __CLPK_integer _n, _ilo, _ihi, _ldh, _ldz, _lwork, info;

   _n = (__CLPK_integer) n;
  _ldh =  (__CLPK_integer) ldh;
  _ilo =(__CLPK_integer) ilo;
  _ihi =(__CLPK_integer) ihi;
  _ldz =  (__CLPK_integer) ldz;
  _lwork = (__CLPK_integer) lwork;
  dhseqr_(&job, &compz, &_n, &_ilo, &_ihi, h, &_ldh, wr, wi, z,&_ldz,  work, &_lwork, &info);
  return info;
}
#endif

int dgeev(char jobvl, char jobvr, int n,  double *a, int lda, double *wr, double *wi, double *vl,
          int ldvl,double *vr, int ldvr, double *work, int lwork){

  long int _n, _lda, _ldvl, _ldvr, _lwork, info;

  _n =  n;
  _lda =   lda;
  _ldvl = ldvl;
  _ldvr = ldvr;
  _lwork = lwork;

  dgeev_(&jobvl, &jobvr, &_n, a, &_lda, wr, wi, vl, &_ldvl, vr, &_ldvr, work, &_lwork, &info);
  return info;
}

int dsyev(char jobvl, char uplo, int n,  double *a, int lda, double *w, double *work, int lwork){

  long int _n, _lda, _lwork, info;

  _n =  n;
  _lda =   lda;
  _lwork = lwork;

  dsyev_(&jobvl, &uplo, &_n, a, &_lda, w, work, &_lwork, &info);
  return info;
}

int dgehrd(int n, int ilo, int ihi, double *a, int lda, double *tau, double *work, int lwork){
  long int _n, _ilo, _ihi, _lda,  _lwork, info;

   _n =  n;
  _lda =   lda;
  _lwork = lwork;
  _ilo = ilo;
  _ihi = ihi;

  dgehrd_(&_n, &_ilo, &_ihi, a, &_lda, tau, work, &_lwork, &info);
  return info;
}

int dgebal(char job, int n,  double *a, int lda, int ilo, int ihi, double *scale){
  long int _n, _ilo, _ihi, _lda, info;

   _n =  n;
  _lda =   lda;
  _ilo = ilo;
  _ihi = ihi;

  dgebal_(&job, &_n, a, &_lda, &_ilo, &_ihi, scale, &info);
  return info;
}

int dhseqr(char job, char compz, int n, int ilo, int ihi,double *h, int ldh, double *wr, double *wi, double *z, int ldz, double *work, int lwork){
  long int _n, _ilo, _ihi, _ldh, _ldz, _lwork, info;

   _n =  n;
  _ldh =   ldh;
  _ilo = ilo;
  _ihi = ihi;
  _ldz =   ldz;
  _lwork =  lwork;
  dhseqr_(&job, &compz, &_n, &_ilo, &_ihi, h, &_ldh, wr, wi, z,&_ldz,  work, &_lwork, &info);
  return info;
}
*/
#endif
