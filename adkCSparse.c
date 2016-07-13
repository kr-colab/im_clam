#include <stdio.h>
#include <stdlib.h>
#include "cs.h"
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include "adkCSparse.h"
#include "adkGSL.h"



//gslMatrixToCSsparse -- returns new cs_sparse object. compressed
struct cs_di_sparse *gslMatrixToCSsparse(gsl_matrix *gm, int numNonZero){
	int i, j;
	struct cs_di_sparse *triplet, *compmat;
	double tmp;
	

	//now allocate new cs_sparse obj
	triplet = cs_spalloc(gm->size1, gm->size2, numNonZero, 1, 1);
	for(i=0;i<gm->size1;i++){
		for(j=0;j<gm->size2;j++){
			tmp = gsl_matrix_get(gm,i,j);
			if(tmp != 0.0) cs_entry(triplet,i,j,tmp);
		}
	}
	compmat = cs_compress(triplet);
	cs_spfree(triplet);
	return(compmat);
}


//gslMatrixToCSsparse -- returns new cs_sparse object. compressed
struct cs_di_sparse *arraysToCSsparse(int rank, int numNonZero, int *dim1, int *dim2,double *x){
	int i;
	struct cs_di_sparse *triplet, *compmat;
	
	//now allocate new cs_sparse obj
	triplet = cs_spalloc(rank, rank, numNonZero, 1, 1);
	for(i=0;i<numNonZero;i++) cs_entry(triplet,dim1[i],dim2[i],x[i]);
	compmat = cs_compress(triplet);
	cs_spfree(triplet);
	return(compmat);
}


//cs_max -- returns the maximum value of a cs_sparse object
double cs_max(struct cs_di_sparse *mat){
	int j,p;
	double max=0.0;
	double val;
	
	if (mat->nz < 0){
		for (j = 0 ; j < mat->n ; j++){
			for(p = mat->p[j] ; p < mat->p[j+1] ; p++){
				val = mat->x[p];
				if(val >max){
					max=val;
				}
			}
		}
	}
	else{ //below is untested
		for(j=0;j<mat->nzmax;j++){
			if(mat->x[j] > max){
				max = mat->x[j];
			}
		}
	}
	return(max);
}

//cs_scale -- scales each value by some constant
void cs_scale(struct cs_di_sparse *mat, double c){
	int i;
	for(i=0;i<mat->nzmax;i++){
		mat->x[i] *= c;
	}
}

//cs_scale -- scales each value by some constant
void cs_print_adk(struct cs_di_sparse *mat){
	int i,j,col,row,p;
	double val;
	gsl_matrix *gm;
	gm = gsl_matrix_alloc(mat->m,mat->n);
	gsl_matrix_set_zero(gm);


	if (mat->nz < 0){
		for (j = 0 ; j < mat->n ; j++){
			col=j;
			for(p = mat->p[j] ; p < mat->p[j+1] ; p++){
				row = mat->i [p];
				val = mat->x[p];
				gsl_matrix_set(gm,row,col,val);
			}
		}
	}
	else{

		for(i=0;i<mat->nz;i++){
			gsl_matrix_set(gm,mat->p[i],mat->i[i],mat->x[i]);
		}
	}
	gsl_matrix_prettyPrint(gm);
	gsl_matrix_free(gm);

}









