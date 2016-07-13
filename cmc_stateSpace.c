//cmc_stateSpace.c --
// this program creates and outputs a coalescent markov chain's complete state space
// given a particular model

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_matrix.h> 
#include "AFS.h"

int n1, n2;

void getParameters(int argc,const char **argv);
void usage();

int main(int argc, const char * argv[]){
	int i;
	afsObject *A;
	afsStateSpace *stateSpace;
	
	
	getParameters(argc,argv);
	
	stateSpace = afsStateSpaceNew();

	
	A = afsObjectNew(n1,n2);
 	afsObjectInit(A,n1,n2);
// 	
// 	printf("%s\n",A->matString );
// //	exit(1);
	stateSpace->states[0] = A;
	stateSpace->nstates=1;
//	stateSpace->tmpStorage[A->nalleles][n1][0] = A;
//	stateSpace->nstateArray[A->nalleles][n1] =1;

	g_hash_table_insert(stateSpace->stateHash, A->matString, "TRUE");


	
	preorderEvolveState(stateSpace, A, n1, n2);

	//qsort(stateSpace->states,stateSpace->nstates, sizeof(stateSpace->states[0]), afsObjectQsortCompare);
	
	printf("n1: %d n2: %d nstates: %d\n",n1,n2,stateSpace->nstates);
	for(i=0;i<stateSpace->nstates;i++){
		afsObjectPrint(stateSpace->states[i]);
		printf("-----------\n");
		
	}


	afsStateSpaceFree(stateSpace);
	return(0);
}


void getParameters(int argc,const char **argv){
//	int args;
//	int i;

	
	if( argc < 2){
		usage();
	}
	n1 = atoi(argv[1]);
	n2 = atoi(argv[2]);
}

void usage(){
	fprintf(stderr,"usage: cmc_stateSpace n1 n2\n");
//	fprintf(stderr,"parameters: \n");
	exit(1);
}