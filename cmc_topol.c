//cmc_topol.c --
// this uses stateSpace generated previously to
// create the topology matrix and the move type matrix
//defined by the markov chain

//this was only for testing purposes

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_matrix.h>
 #include <gsl/gsl_blas.h>

#include "AFS.h"
#include "adkGSL.h"
#include "cs.h"
#include "time.h"
#include "adkCSparse.h"
#include "AFS_pthreads.c"

#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))
#define NTHREADS 4

int n1, n2;
const char *filename;

afsStateSpace *stateSpace;
void getParameters(int argc,const char **argv);
void usage();

int main(int argc, const char * argv[]){
	int i,N,tmpS,win,t;
	int rc;
	void *v;
	int nnz;
	matrixThreadObject_sparse *tObj[NTHREADS];
 	pthread_t threads[NTHREADS];

	getParameters(argc,argv);
	stateSpace = afsStateSpaceImportFromFile(filename);
	//afsStateSpaceRemoveAbsorbing(stateSpace);
	N = stateSpace->nstates;
	
	
	//Multithreaded construction of original transition matrix
//	printf("initializing original transition matrix using %d threads....\n",NTHREADS);
	//setup threads
	tmpS=0;
	win =(int) ceil( (((float)stateSpace->nstates) / NTHREADS) );
//	printf("nstates: %d winSize: %d\n",stateSpace->nstates,win);
	for(t=0;t<NTHREADS;t++){
		tObj[t] = malloc(sizeof( matrixThreadObject_sparse));
		tObj[t]->nnz = 0;
		tObj[t]->dim1 = malloc(sizeof(int) * stateSpace->nstates * 3);
		tObj[t]->dim2 = malloc(sizeof(int) * stateSpace->nstates * 3);
		tObj[t]->moveType = malloc(sizeof(int) * stateSpace->nstates *  3);
		tObj[t]->topol = malloc(sizeof(double) * stateSpace->nstates *  3);
		tObj[t]->start = tmpS;
		tObj[t]->stop = MIN(tmpS+win,stateSpace->nstates);
		tObj[t]->im_params = malloc(sizeof(im_lik_params));
		tObj[t]->im_params->stateSpace = stateSpace;
//		printf("start: %d stop: %d\n",tmpS,tmpS+win);
		rc = pthread_create(&threads[t], NULL, threadDoRows_sparse, (void *)tObj[t]);
		tmpS += win;
	}
//	printf("threads alloced....\n");

	//harvest jobs
	nnz=0;
	for(t=0;t<NTHREADS;t++){
		rc = pthread_join(threads[t], &v);
		tObj[t] = ((matrixThreadObject_sparse *)v);
		nnz += tObj[t]->nnz;
	}
	//now combine outputs for nnz, dim1, and dim2
	printf("nnz: %d\n",nnz);
	for(t=0;t<NTHREADS;t++){
		for(i=0;i<tObj[t]->nnz;i++){
			printf("%lf\t%d\t%d\t%d\n",tObj[t]->topol[i],tObj[t]->moveType[i],tObj[t]->dim1[i],tObj[t]->dim2[i]);
		}
	}
	//currentParams->nnz -= 1;
	////////// Multithreading done!

	// dim1 = malloc(sizeof(int) * N * 10);
	// dim2 = malloc(sizeof(int) * N * 10);
	// moveA = malloc(sizeof(int) * N * 10);
	// topA = malloc(sizeof(double) * N * 10);
	// for(i=0;i< (N*10);i++){
	// 	moveA[i]=0;
	// 	topA[i]=0;
	// }
	// nnz = coalMarkovChainTopologyMatrix_sparse(stateSpace,topA, moveA, dim1, dim2);
	// printf("nnz: %d\n",nnz);
	// for(i=0;i<nnz;i++){
	// 	printf("%f\t%d\t%d\t%d\n",topA[i],moveA[i],dim1[i],dim2[i]);
	// }
	return(0);
}

void getParameters(int argc, const char **argv){
//	int args;
//	int i;

	if( argc <= 1){
		usage();
	}
	filename = argv[1];
	
}

void usage(){
	fprintf(stderr,"usage: cmc_topol stateSpaceFile\n");
//	fprintf(stderr,"parameters: \n");
	exit(1);
}
