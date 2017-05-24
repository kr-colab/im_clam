/*
/
/ adkGSL.h
/
/ homebrewed addons to gsl */

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <math.h>
#include <assert.h>

int gsl_vector_includes(gsl_vector *aVec, double aValue);
double gsl_vector_sum(gsl_vector *aVec, int aVecSize);
double gsl_matrix_row_sum(gsl_matrix *mat, int row, int rowSize);
double gsl_matrix_row_max(gsl_matrix *mat, int iRow);

double log_sum(gsl_vector *vec);
double  *gsl_matrix_2_lapack(gsl_matrix *m);
double gsl_vector_dot_product(gsl_vector *vec1, gsl_vector *vec2);
void gsl_vector_outer_product(gsl_vector *vec1, gsl_vector *vec2, gsl_matrix *result);

gsl_matrix *gsl_matrix_power(gsl_matrix *A, int k);
gsl_matrix *gsl_matrix_power_logs(gsl_matrix *A, int k);
void gsl_matrix_lower_tri_copy(gsl_matrix *src, gsl_matrix *dest);
void gsl_matrix_upper_tri_copy(gsl_matrix *src, gsl_matrix *dest);
void gsl_matrix_lower_tri(gsl_matrix *x);
void gsl_matrix_upper_tri(gsl_matrix *x);
double gsl_matrix_trace(gsl_matrix *x);
void gsl_matrix_covariance(gsl_matrix *data, gsl_matrix *cov);
void gsl_matrix_bootstrap(gsl_matrix *orig, gsl_matrix *boot, gsl_rng *rng);

void gsl_matrix_prettyPrint(gsl_matrix *m);

void fillMatrixFromArray(double *numbers, gsl_matrix *dest, int nrow, int ncol);
void fillMatrixFromVector(gsl_vector *numbers, gsl_matrix *dest, int nrow, int ncol);
void fillArrayFromMatrix(gsl_matrix *src, gsl_vector *dest);

void fillCholArrayFromMatrix(gsl_matrix *src, gsl_vector *dest);
void fillMatrixFromCholArray(double *numbers, gsl_matrix *dest, int nrow, int ncol);

int gsl_matrix_invert_lapack(gsl_matrix *m);

void gsl_matrix_prettyPrint(gsl_matrix *m);
double gsl_matrix_sum(gsl_matrix *mat);

//lapack wrappers
//int dgeev(char jobvl, char jobvr, int n,  double *a, int lda, double *wr, double *wi, double *vl, 
//	  int ldvl,double *vr, int ldvr, double *work, int lwork);
//int dhseqr(char job, char compz, int n, int ilo, int ihi,double *h, int ldh, double *wr, double *wi, double *z, int ldz, double *work, int lwork);

