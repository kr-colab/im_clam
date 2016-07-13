/* AFS.c
These are definitions / methods of the AFS object for a markov chain representation
of the AFS statespace during the coalecsent process

A. Kern 2012
*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_matrix.h> 
#include <gsl/gsl_matrix_int.h> 
#include <gsl/gsl_cblas.h> 
#include <gsl/gsl_sf.h> 
#include "AFS.h"
#include "adkGSL.h"
//#include "cmc_island_mle.h"
#include "multimin.h"
#include "sites.h"
#include "cs.h"
#include "adkCSparse.h"
#include <glib.h>


afsObject *afsObjectNew(int n1, int n2){
	int i;
	
	afsObject *tmp;
	tmp = (afsObject *)malloc(sizeof(afsObject));
	tmp->popMats =(gsl_matrix_int **) malloc(2* sizeof(gsl_matrix_int));
	tmp->popMats[0] = (gsl_matrix_int *) gsl_matrix_int_alloc(n1+1,n2+1);
	tmp->popMats[1] = (gsl_matrix_int *) gsl_matrix_int_alloc(n1+1,n2+1);
	tmp->npops=2; //unneeded right now but maybe will generalize in the future
	tmp->matString = malloc(sizeof(char) * ((n1+1)*(n2+1)*2) + 1);
	for(i=0;i<(n1*n2*2);i++){tmp->matString[i] ='\0';}
	
	return(tmp);
}

//afsObjectNewFrom-- creates a deep copy of an afsObject
afsObject *afsObjectNewFrom(afsObject *v){
	afsObject *tmp;
	tmp = (afsObject *) malloc(sizeof(afsObject));
	tmp->popMats = (gsl_matrix_int **) malloc(2* sizeof(gsl_matrix_int));
	tmp->popMats[0] = (gsl_matrix_int *)gsl_matrix_int_alloc(v->popMats[0]->size1,v->popMats[0]->size2);
	tmp->popMats[1] = (gsl_matrix_int *)gsl_matrix_int_alloc(v->popMats[1]->size1,v->popMats[1]->size2);
	tmp->npops=v->npops;
	gsl_matrix_int_memcpy(tmp->popMats[0],v->popMats[0]);
	gsl_matrix_int_memcpy(tmp->popMats[1],v->popMats[1]);
	tmp->nalleles=afsObjectCount(tmp);
	tmp->matString = malloc(sizeof(char) * ((v->popMats[0]->size1 +1) * (v->popMats[0]->size2 +1) * 2) + 1);
//	for(i=0;i<(v->popMats[0]->size1 * v->popMats[0]->size2 * 2);i++){tmp->matString[i] ='\0';}
	tmp->matString[0] ='\0';
	return(tmp);
}

void afsObjectFree(afsObject *v){
	gsl_matrix_int_free(v->popMats[0]);
	gsl_matrix_int_free(v->popMats[1]);
	free(v->popMats);
	free(v->matString);
	free(v);
}

//sets up the initial sample of size n1 from population 1 and size n2 from
// population 2
void afsObjectInit(afsObject *v, int n1, int n2){
	int i, j; 
	char c[1000];
	gsl_matrix_int_set_zero(v->popMats[0]);
	gsl_matrix_int_set_zero(v->popMats[1]);
	gsl_matrix_int_set(v->popMats[0],1,0,n1);
	gsl_matrix_int_set(v->popMats[1],0,1,n2);
	v->aCounts[0]=n1;
	v->aCounts[1]=n2;
	v->nalleles = n1 + n2;
	for(i=0;i<= n1;i++){
		for(j=0;j<=n2;j++){
			itoa(gsl_matrix_int_get(v->popMats[0],i,j),c);
			strcat(v->matString, c);
			itoa(gsl_matrix_int_get(v->popMats[1],i,j),c);
			strcat(v->matString, c);
		}
	}
	//printf("from within afsObjInit: %s\n",v->matString);
}

void afsObjectPrint(afsObject *v){
	int i, j;
	int m1, m2;
	m1 = v->popMats[0]->size1;
	m2 = v->popMats[0]->size2;
	for(i=0;i<m1;i++){
		for(j=0;j<m2;j++) printf("%d ",gsl_matrix_int_get(v->popMats[0],i,j));
		printf("|| ");
		for(j=0;j<m2;j++) printf("%d ",gsl_matrix_int_get(v->popMats[1],i,j));
		printf("\n");
	}
}

int nDescPop1(afsObject *v){
	int sum, i, j;
	int m1, m2;
	m1 = v->popMats[0]->size1;
	m2 = v->popMats[0]->size2;
	sum=0;
	for(i=0;i<m1;i++){
		for(j=0;j<m2;j++){
			sum+= gsl_matrix_int_get(v->popMats[0],i,j) * i;
			sum+= gsl_matrix_int_get(v->popMats[1],i,j) * i;
		}
	}
	return(sum);
}

int nDescPop2(afsObject *v){
	int sum, i, j;
	int m1, m2;
	m1 = v->popMats[0]->size1;
	m2 = v->popMats[0]->size2;
	sum=0;
	for(i=0;i<m1;i++){
		for(j=0;j<m2;j++){
			sum+= gsl_matrix_int_get(v->popMats[0],i,j) * j;
			sum+= gsl_matrix_int_get(v->popMats[1],i,j) * j;
		}
	}
	return(sum);
}
/*
void gsl_matrix_prettyPrint(gsl_matrix *m){
	int i, j;
	for(i=0;i<m->size1;i++){
		for(j=0;j<m->size2;j++){
			printf("%f ",gsl_matrix_get(m,i,j));
		}
		printf("\n");
	}
	printf("///////////\n");
}
*/
//compatibleAfsGivenInitial tests if correct number of ancestral lineages
// are present and if each entry is >= 0, i.e. no bogus negative afs cells
int compatibleAfsGivenInitial(afsObject *v, int m1, int m2){


	if (nDescPop1(v) == m1 && nDescPop2(v) == m2 && gsl_matrix_int_min(v->popMats[0]) >= 0 
		&& gsl_matrix_int_min(v->popMats[1]) >= 0)
		return(1);
	else
		return(0);
}

//nonZeroEntries -- fills a 2D list of coords for nonZero entries in matrix
void nonZeroEntries(gsl_matrix_int *m, int **coords, int *count){
	int i, j;

	*count = 0;
	for(i=0;i<m->size1;i++){
		for(j=0;j<m->size2;j++){
			if(gsl_matrix_int_get(m,i,j) != 0){
				coords[*count][0]=i;
				coords[*count][1]=j;
				*count+=1;
			}
		}
	}
} 

//countNegativeEntries-- returns the count of matrix entries < 0
int countNegativeEntries(gsl_matrix_int *m){
	int i, j, sum=0;
	for(i=0;i<m->size1;i++){
		for(j=0;j<m->size2;j++){
			if(gsl_matrix_int_get(m,i,j) < 0)
				sum+=1;
			}
		}
		return(sum);
}

//entriesGreaterThan -- fills a 2D list of coords for nonZero entries in matrix
void entriesGreaterThan(gsl_matrix_int *m, int **coords, int *count, int cut){
	int i, j;

	*count = 0;
	for(i=0;i<m->size1;i++){
		for(j=0;j<m->size2;j++){
			if(gsl_matrix_int_get(m,i,j) > cut){
				coords[*count][0]=i;
				coords[*count][1]=j;
				*count+=1;
			}
		}
	}
}

//entriesLessThan -- fills a 2D list of coords for nonZero entries in matrix
void entriesLessThan(gsl_matrix_int *m, int **coords, int *count, int cut){
	int i, j;

	*count = 0;
	for(i=0;i<m->size1;i++){
		for(j=0;j<m->size2;j++){
			if(gsl_matrix_int_get(m,i,j) < cut){
				coords[*count][0]=i;
				coords[*count][1]=j;
				*count+=1;
			}
		}
	}
}

/* factorials */

double xchoosey(int n, int k){
  if(n<k) return(0);
  if(n == k || k == 0) return(1);
  if(k == 1) return(n);
if(k==2) return((float)n*(n-1.0)/2.0);
  return(gsl_sf_choose(n,k));
}


//matrixSum -- just convenience function
int matrixSum(gsl_matrix_int *m){
	int i,j, sum=0;
	for(i=0;i<m->size1;i++)
		for(j=0;j<m->size2;j++)
		sum+= gsl_matrix_int_get(m,i,j);
	return(sum);
}

//matrixSumDouble -- just convenience function
double matrixSumDouble(gsl_matrix *m){
	int i,j; 
	double sum=0;
	for(i=0;i<m->size1;i++)
		for(j=0;j<m->size2;j++)
		sum+= gsl_matrix_get(m,i,j);
	return(sum);
}

//matrixAllGTZ-- again a silly convenience function
int matrixAllGTZ(gsl_matrix_int *m){
	int i,j;
	for(i=0;i<m->size1;i++)
		for(j=0;j<m->size2;j++)
		if(gsl_matrix_int_get(m,i,j) <=0)
		return(0);
	return(1);
}

//afsObjectsEqual -- compares objects based on matrices and nothing more
int afsObjectsEqual(afsObject *a, afsObject *b){
	if(gsl_matrix_int_equal(a->popMats[0],b->popMats[0]) 
		&& gsl_matrix_int_equal(a->popMats[1],b->popMats[1]))
		return(1);
	else
		return(0);
}

int afsObjectQsortCompare(const void *p1, const void *p2){
	  afsObject **x = (afsObject **)p1;
	  afsObject **y = (afsObject **)p2;
	
	if((*x)->nalleles < (*y)->nalleles)
		return 1;
   	else if ((*x)->nalleles > (*y)->nalleles)
      		return -1;
   	else
      		return 0;
}

//subtracts matrices from one another to aid in ID of operations
afsObject *afsObjectDelta(afsObject *a, afsObject *b){
	afsObject *tmp;

	tmp = afsObjectNewFrom(a);

	gsl_matrix_int_memcpy(tmp->popMats[0],a->popMats[0]);
	gsl_matrix_int_memcpy(tmp->popMats[1],a->popMats[1]);
	gsl_matrix_int_sub(tmp->popMats[0],b->popMats[0]);
	gsl_matrix_int_sub(tmp->popMats[1],b->popMats[1]);
	return(tmp);
}

//subtracts matrices from one another to aid in ID of operations
void afsObjectDeltaPre(afsObject *a, afsObject *b,afsObject *tmp){
	gsl_matrix_int_memcpy(tmp->popMats[0],a->popMats[0]);
	gsl_matrix_int_memcpy(tmp->popMats[1],a->popMats[1]);
	gsl_matrix_int_sub(tmp->popMats[0],b->popMats[0]);
	gsl_matrix_int_sub(tmp->popMats[1],b->popMats[1]);
}

//afsObjectCount-- returns the total number of lineages in the afs state
int afsObjectCount(afsObject *a){
	return(matrixSum(a->popMats[0])+matrixSum(a->popMats[1]));	
}

afsStateSpace *afsStateSpaceNew(){

	afsStateSpace *tmp;
	tmp =(afsStateSpace *) malloc(sizeof(afsStateSpace));
	tmp->nstates=0;
	tmp->states = (afsObject **) malloc(MAXSTATES * sizeof(afsObject));
	tmp->stateHash = g_hash_table_new(g_str_hash, g_str_equal);

	return(tmp);
}

afsStateSpace *afsStateSpaceCopy(afsStateSpace *S){
	int i;
	afsStateSpace *tmp;
	tmp =(afsStateSpace *) malloc(sizeof(afsStateSpace));
	tmp->nstates=S->nstates;
	tmp->states = (afsObject **) malloc(MAXSTATES * sizeof(afsObject));
	for(i=0;i<S->nstates;i++) tmp->states[i] = S->states[i];
//	tmp->stateHash = g_hash_table_new(g_str_hash, g_str_equal);


	return(tmp);
}
void afsStateSpaceFree(afsStateSpace *a){
	int i;
	for(i=0;i<a->nstates;i++)afsObjectFree(a->states[i]);
	free(a->states);
	// for(i=0;i<MAXSIZE;i++){
	// 	free(a->nstateArray[i]);
	// 	for(j=0;j<MAXSIZE;j++){
	// 		free(a->tmpStorage[i][j]);
	// 	}
	// 	free(a->tmpStorage[i]);
	// }
	// free(a->tmpStorage);
	// free(a->nstateArray);
	g_hash_table_destroy(a->stateHash);
	
	free(a);
}

//afsStateSpaceContainsState-- linear search for state
int afsStateSpaceContainsState(afsStateSpace *S, afsObject *aState){
	int i;
	for(i = S->nstates -1 ; i > - 1;i--){
		if(S->states[i]->nalleles == aState->nalleles)
			if(S->states[i]->aCounts[0] == aState->aCounts[0])
				if(afsObjectsEqual(S->states[i],aState))
					return(1);		
	}
	return(0);
}

/////// Below is not needed with hash implementation
//afsStateSpaceContainsStateTmpStorage-- same as above but works off of stateSpace->tmpStorage
// int afsStateSpaceContainsStateTmpStorage(afsStateSpace *S, afsObject *aState){
// 	int i;
// 	for(i = S->nstateArray[aState->nalleles][aState->aCounts[0]] -1 ; i > - 1;i--){
// 			if(afsObjectsEqual(S->tmpStorage[aState->nalleles][aState->aCounts[0]][i],aState))
// 				return(1);		
// 	}
// 	return(0);
// }

//afsStateSpaceRemoveAbsorbing-- removes all absorbing states
void afsStateSpaceRemoveAbsorbing(afsStateSpace *S){
	int i,ii;
	afsObject *tmp;
	i = 0;
	while(i < S->nstates){
		if(S->states[i]->nalleles == 1){
			tmp = S->states[i];
			for(ii=i;ii<S->nstates - 1;ii++) S->states[ii] = S->states[ii+1];
			S->nstates--;
			i-=1;
			afsObjectFree(tmp);
		}
	 	i+=1; 
	}
}

//afsStateSpaceMapPopn-- returns a vector mapping each state to its state if popn1 collapsed into popn0
void afsStateSpaceMapPopn(afsStateSpace *S, int *map){
	int i, j,valCount,valInd[S->nstates / 2];
	afsObject *tmp, *test;
	
	tmp = afsObjectNewFrom(S->states[0]);
	gsl_matrix_int_set_zero(tmp->popMats[1]);
	//find and count valid states
	valCount = 0;
	for(i=0;i<S->nstates;i++){
		if(S->states[i]->aCounts[1]==0){
			valInd[valCount]=i;
			valCount++;
		}
	}
	for(i=0;i<S->nstates;i++){
		if(S->states[i]->aCounts[1]==0){
			map[i]=i;
		}
		else{
			//create folded state in tmp
			gsl_matrix_int_memcpy(tmp->popMats[0],S->states[i]->popMats[0]);
			gsl_matrix_int_add(tmp->popMats[0],S->states[i]->popMats[1]);
			//go through valid ones to find it
			for(j=0;j<valCount;j++){
				test= S->states[valInd[j]];
				if(S->states[i]->nalleles == test->nalleles)
					if(gsl_matrix_int_equal(tmp->popMats[0],test->popMats[0]))
					map[i]=valInd[j];		
			}
		}
	}
}

//afsStateSpaceMapAndReducePopn-- returns a vector mapping each state to its state if popn1 collapsed into popn0 (map)
//adds valid states to new reducedSpace, also fills reverseMap to larger space indexing
//reducedSpace is passed as ptr, assumed to be empty, initial stateSpace
void afsStateSpaceMapAndReducePopn(afsStateSpace *S, int *map, afsStateSpace *reducedSpace, int *reverseMap){
	int i, j,valCount,valInd[S->nstates / 2];
	afsObject *tmp, *test;
	

	
	tmp = afsObjectNewFrom(S->states[0]);
	gsl_matrix_int_set_zero(tmp->popMats[1]);
	//find and count valid states
	valCount = 0;
	for(i=0;i<S->nstates;i++){
		if(S->states[i]->aCounts[1]==0){
			valInd[valCount]=i;
			reverseMap[valCount]=i;
			reducedSpace->states[valCount]=S->states[i];
			reducedSpace->nstates++;
			valCount++;
		}
	}
	for(i=0;i<S->nstates;i++){
		if(S->states[i]->aCounts[1]==0){
			map[i]=i;
		}
		else{
			//create folded state in tmp
			gsl_matrix_int_memcpy(tmp->popMats[0],S->states[i]->popMats[0]);
			gsl_matrix_int_add(tmp->popMats[0],S->states[i]->popMats[1]);
			//go through valid ones to find it
			for(j=0;j<valCount;j++){
				test= S->states[valInd[j]];
				if(S->states[i]->nalleles == test->nalleles)
					if(gsl_matrix_int_equal(tmp->popMats[0],test->popMats[0]))
					map[i]=valInd[j];		
			}
		}
	}
}


//afsStateSpaceRemovePopulation-- removes all states with members in specified popn
void afsStateSpaceRemovePopulation(afsStateSpace *S, int popn){
	int i,ii;
	afsObject *tmp;
	i = 0;
	while(i < S->nstates){
		if(S->states[i]->aCounts[1] > 0){
			tmp = S->states[i];
			for(ii=i;ii<S->nstates - 1;ii++) S->states[ii] = S->states[ii+1];
			S->nstates--;
			i-=1;
			afsObjectFree(tmp);
		}
	 	i+=1; 
	}
}

afsStateSpace *afsStateSpaceImportFromFile(const char *fileName){
	int m1,m2,ns;
	int i,j,l;
	FILE *infile;
	afsStateSpace *newStateSpace;
//	afsObject *newObj;
	int linebuf = 10000;
	char line[linebuf+1],*token, search[]=" ";

	//open file
	infile = fopen(fileName, "r");
	if (infile == NULL){
		fprintf(stderr,"Error opening infile! ARRRRR!!!!\n");
		exit(1);
	}
	//initialization
	//define state space size and shape
	fscanf(infile,"n1: %d n2: %d nstates: %d", &m1, &m2, &ns);
	newStateSpace = afsStateSpaceNew();
	newStateSpace->nstates = ns;
	fgets(line,linebuf,infile);

	for(l=0;l<ns;l++){
		newStateSpace->states[l] = afsObjectNew(m1,m2);
		for(i=0;i<m1+1;i++){
			fgets(line,linebuf,infile);
	//		printf("%s\n",line);
	//		printf("----\n");
			token = strtok(line,search);


			gsl_matrix_int_set(newStateSpace->states[l]->popMats[0],i,0,atoi(token));
			for(j=1;j<m2+1;j++){
				token = strtok(NULL,search);
		//		printf("%s\n",token);
				gsl_matrix_int_set(newStateSpace->states[l]->popMats[0],i,j,atoi(token));
			}
			token = strtok(NULL,search);
			for(j=0;j<m2+1;j++){
				token = strtok(NULL,search);
				gsl_matrix_int_set(newStateSpace->states[l]->popMats[1],i,j,atoi(token));

			}

		}
		fgets(line,linebuf,infile);
		newStateSpace->states[l]->aCounts[0] = matrixSum(newStateSpace->states[l]->popMats[0]);
		newStateSpace->states[l]->aCounts[1] = matrixSum(newStateSpace->states[l]->popMats[1]);
		newStateSpace->states[l]->nalleles = afsObjectCount(newStateSpace->states[l]);


	//	afsObjectPrint(newStateSpace->states[l]);
	//	printf("+++++\n");
	}
	fclose(infile);
	return(newStateSpace);

}

//mcMatsImportFromFile -- reads in the Topology and MoveType matrices output by cmc_printMCMats program
// this function only fills preallocated arrays and sets nnz
void mcMatsImportFromFile(const char *fileName,int *nnz,double *topol, int *moveA, int *dim1, int *dim2){
	int i;
	FILE *infile;

//	afsObject *newObj;
	int linebuf = 10000;
	char line[linebuf+1],*token, search[]=" ";

	//open file
	infile = fopen(fileName, "r");
	if (infile == NULL){
		fprintf(stderr,"Error opening infile! ARRRRR!!!!\n");
		exit(1);
	}
	//initialization
	fscanf(infile,"nnz: %d", nnz);
	for(i=0;i<*nnz;i++){
		fscanf(infile,"%lf\t%d\t%d\t%d\n",&topol[i],&moveA[i],&dim1[i],&dim2[i]);
	}
	fclose(infile);

}



//acceptOrRejectProposedState-- tests if candidate state in valid and not yet represented
// in the afsStateSpace object. Adds to array if true and iterates nstates. Free's candidate if
// false

void acceptOrRejectProposedState(afsStateSpace *S, afsObject *child, int n1, int n2){
	int i,j;
	char c[10];
	
	child->nalleles=afsObjectCount(child);
	child->aCounts[0]=matrixSum(child->popMats[0]);
	child->aCounts[1]=matrixSum(child->popMats[1]);

	for(i=0;i <= n1;i++){
		for(j=0;j <= n2;j++){
			itoa(gsl_matrix_int_get(child->popMats[0],i,j),c);
			strcat(child->matString, c);
			itoa(gsl_matrix_int_get(child->popMats[1],i,j),c);
			strcat(child->matString, c);
		}
	}
//	printf("%s\n",child->matString);
	//Below block is the old way. new way is hash implementation
	// if(compatibleAfsGivenInitial(child, n1,n2) == 1
	// && afsStateSpaceContainsStateTmpStorage(S,child) == 0){
	// 		//add in new state
	// 
	// 	S->states[S->nstates]=child;
	// 	S->nstates++;
	// 	S->tmpStorage[child->nalleles][child->aCounts[0]][S->nstateArray[child->nalleles][child->aCounts[0]]] = child;
	// 	S->nstateArray[child->nalleles][child->aCounts[0]]++;
			//check to make sure state array can keep up!
//			if(S->nstates==MAXSTATES){
//				MAXSTATES = 2*MAXSTATES;
//				S->states = realloc(S->states, MAXSTATES * sizeof(child));
//			}
	if(g_hash_table_lookup(S->stateHash, child->matString) == NULL 
		&& compatibleAfsGivenInitial(child, n1,n2) == 1){
		//add in new state
		S->states[S->nstates]=child;
		S->nstates++;
		g_hash_table_insert(S->stateHash, child->matString, "TRUE");
	
		//recurse!
		preorderEvolveState(S, child, n1,n2);
	}
	else{
		afsObjectFree(child);
	}
}



//preorderEvolveState -- this is the workhorse that fills the entire state space
// currently this only handles coalescent and migration events for 2 populations
void preorderEvolveState(afsStateSpace *S, afsObject *aState, int n1, int n2){
	int i,j, nlin1,nlin2, tmp1, tmp2;
	int **activePos1, **activePos2; //first dimension here relates to 2 popns
	afsObject *child;



	activePos1 = malloc(100*sizeof(int *));
	activePos2 = malloc(100*sizeof(int *));
	for(i=0;i<100;i++){
		activePos1[i] = malloc(2*sizeof(int));
		activePos1[i][0] = activePos1[i][1] = 0;
		activePos2[i] = malloc(2*sizeof(int));
		activePos2[i][0] = activePos2[i][1] = 0;

	}

	nonZeroEntries(aState->popMats[0], activePos1, &nlin1);
	nonZeroEntries(aState->popMats[1], activePos2, &nlin2);

	//coalescent in Popn 1?
	if(matrixSum(aState->popMats[0]) > 1){
		//coalescent is possible; permutations of lineages
		for(i=0;i<nlin1;i++){
			for(j=0;j<nlin1;j++){

				if( i != j || gsl_matrix_int_get(aState->popMats[0],activePos1[i][0],activePos1[i][1]) > 1){	
					child = afsObjectNewFrom(aState);


					//get and set below is ugly... basically subtract one from each of coal lineages add one 
					//to new lineage
					gsl_matrix_int_set(child->popMats[0],activePos1[i][0],
						activePos1[i][1], gsl_matrix_int_get(child->popMats[0],
						activePos1[i][0],activePos1[i][1]) - 1);
					gsl_matrix_int_set(child->popMats[0],activePos1[j][0],
						activePos1[j][1], gsl_matrix_int_get(child->popMats[0],
						activePos1[j][0],activePos1[j][1]) - 1);
					tmp1 = activePos1[i][0] + activePos1[j][0];
					tmp2 = activePos1[i][1] + activePos1[j][1];
					gsl_matrix_int_set(child->popMats[0],tmp1,
						tmp2, gsl_matrix_int_get(child->popMats[0],
						tmp1,tmp2) + 1);

					//now check state for compatibilty
					acceptOrRejectProposedState(S,child, n1, n2);
				}
			}
		}
	}

	//coalescent in Popn 2?
	if(matrixSum(aState->popMats[1]) > 1){
		//coalescent is possible; permutations of lineages
		for(i=0;i<nlin1;i++){
			for(j=0;j<nlin1;j++){

				if( i != j || gsl_matrix_int_get(aState->popMats[1],activePos2[i][0],activePos2[i][1]) > 1){	
					child = afsObjectNewFrom(aState);


					//get and set below is ugly... basically subtract one from each of coal lineages add one 
					//to new lineage
					gsl_matrix_int_set(child->popMats[1],activePos2[i][0],
						activePos2[i][1], gsl_matrix_int_get(child->popMats[1],
						activePos2[i][0],activePos2[i][1]) - 1);
					gsl_matrix_int_set(child->popMats[1],activePos2[j][0],
						activePos2[j][1], gsl_matrix_int_get(child->popMats[1],
						activePos2[j][0],activePos2[j][1]) - 1);
					tmp1 = activePos2[i][0] + activePos2[j][0];
					tmp2 = activePos2[i][1] + activePos2[j][1];
					gsl_matrix_int_set(child->popMats[0],tmp1,
						tmp2, gsl_matrix_int_get(child->popMats[0],
						tmp1,tmp2) + 1);

					//now check state for compatibilty
					acceptOrRejectProposedState(S,child, n1, n2);
				}
			}
		}
	}

	//migration from Pop1 to Popn 2
	for(i=0;i<nlin1;i++){
		child = afsObjectNewFrom(aState);
		gsl_matrix_int_set(child->popMats[0],activePos1[i][0],
			activePos1[i][1], gsl_matrix_int_get(child->popMats[0],
			activePos1[i][0],activePos1[i][1]) - 1);
		gsl_matrix_int_set(child->popMats[1],activePos1[i][0],
			activePos1[i][1], gsl_matrix_int_get(child->popMats[1],
			activePos1[i][0],activePos1[i][1]) + 1);
		//now check state for compatibilty
		//now check state for compatibilty
		acceptOrRejectProposedState(S,child, n1, n2);
	}	

	//migration from Pop2 to Popn 1
	for(i=0;i<nlin2;i++){
		child = afsObjectNewFrom(aState);
		gsl_matrix_int_set(child->popMats[1],activePos2[i][0],
			activePos2[i][1], gsl_matrix_int_get(child->popMats[1],
			activePos2[i][0],activePos2[i][1]) - 1);
		gsl_matrix_int_set(child->popMats[0],activePos2[i][0],
			activePos2[i][1], gsl_matrix_int_get(child->popMats[0],
			activePos2[i][0],activePos2[i][1]) + 1);
		//now check state for compatibilty
		//now check state for compatibilty
		acceptOrRejectProposedState(S,child, n1, n2);
	}

	//cleanup and free
	for(i=0;i<100;i++){
		free(activePos1[i]); 
		free(activePos2[i]); 
	}
	free(activePos1);
	free(activePos2);
}

/////coalMarkovChainTopologyMatrix -- using a previously calculated state space
//this fills a matrix of topological moves possible in the markov chain, 
//along with a lookup table of moves that are parameter dependant
//
// codes for lookup :: 
// -1 : absorb
// 0 : coal popn0
// 1 : coal popn1
// 2 : mig popn0
// 3 : mig popn2
// 666: incompatible move

void coalMarkovChainTopologyMatrix(afsStateSpace *S,gsl_matrix *topol, gsl_matrix_int *moveType, int n1, int n2){
	int i,j,k,steps,x,y;
	double top,bottom;
	afsObject *delta;
	int **activePos1, **activePos2; //first dimension here relates to 2 popns
	int nlin1,nlin2,compat;

	if(S->nstates != topol->size1 || S->nstates != moveType->size1){
		fprintf(stderr,"Error: StateSpace and matrices are not of equal size\n");
		exit(1);
	}

	//arrays for finding positions of entries
	activePos1 = malloc(100*sizeof(int *));
	activePos2 = malloc(100*sizeof(int *));
	for(i=0;i<100;i++){
		activePos1[i] = malloc(2*sizeof(int));
		activePos1[i][0] = activePos1[i][1] = 0;
		activePos2[i] = malloc(2*sizeof(int));
		activePos2[i][0] = activePos2[i][1] = 0;

	}

	for(i=0;i<S->nstates;i++){
		for(j=0;j<S->nstates;j++){
			gsl_matrix_int_set(moveType,i,j,666);

			//is ith state root state?
			if(S->states[i]->nalleles == 1){
				//set correct absorbing conditions
				if(i==j){
					gsl_matrix_set(topol,i,j,1.0);
					gsl_matrix_int_set(moveType,i,j,-1);
				}
				else{
					gsl_matrix_set(topol,i,j,0.0);
					gsl_matrix_int_set(moveType,i,j,-1);
				}		
			}
			else{ //ith not root; what is the move?
				delta = afsObjectDelta(S->states[i],S->states[j]);
				steps = S->states[i]->nalleles - S->states[j]->nalleles ;
			//	printf("i:%d j:%d steps:%d\n",i,j,steps);
				switch(steps){
					case 0:
					// 0 "steps", so no changes in nalleles between two
					//check counts and positions	
					if(abs(matrixSum(delta->popMats[0])) == 1){
						nonZeroEntries(delta->popMats[0], activePos1, &nlin1);
						nonZeroEntries(delta->popMats[1], activePos2, &nlin2);


						//exit(1);
						//migration?
						if(nlin1 == nlin2 && activePos1[0][0] == activePos2[0][0] &&
						(gsl_matrix_int_min(delta->popMats[0]) >= 0 || gsl_matrix_int_min(delta->popMats[1]) >= 0))
						{

							//migration popn1?
							if(matrixSum(delta->popMats[0]) == 1){
								gsl_matrix_set(topol,i,j,
									(double) gsl_matrix_int_get(S->states[i]->popMats[0],
									activePos1[0][0], activePos1[0][1]) / S->states[i]->aCounts[0]);
								gsl_matrix_int_set(moveType,i,j,2);
							}
							else{ //migration popn2
						//		printf("here %d %d %d\n",activePos2[0][0],activePos2[0][1], S->states[i]->aCounts[1]);
								gsl_matrix_set(topol,i,j,
									(double) gsl_matrix_int_get(S->states[i]->popMats[1],
									activePos2[0][0],activePos2[0][1]) / S->states[i]->aCounts[1]);
								gsl_matrix_int_set(moveType,i,j,3);
							}
						}
					}
					break;

					case 1://coalescent?
					//coal popn1?
					if(matrixSum(delta->popMats[0]) == 1 && countNegativeEntries(delta->popMats[1]) == 0 && S->states[j]->aCounts[0] >= 1 
						/*((matrixAllGTZ(delta->popMats[0]) || matrixAllGTZ(delta->popMats[1]))) */ ){

					//	afsObjectPrint(S->states[i]);
					//	printf("1///////////\n");
					//	afsObjectPrint(S->states[j]);
					//	printf("2///////////\n");
						entriesGreaterThan(delta->popMats[0], activePos1, &nlin1,0);//old lineages
						entriesLessThan(delta->popMats[0], activePos2, &nlin2,0);//new lineages
						if(nlin2 != 0 && nlin2 < 2){
							
							compat=0;
							if(nlin1 == 1){
								if(2*activePos1[0][0] == activePos2[0][0] && 2*activePos1[0][1] == activePos2[0][1])
									compat = 1;
							}
							if(nlin1==2){
								if(activePos1[0][0] + activePos1[1][0] == activePos2[0][0] &&
									activePos1[0][1] + activePos1[1][1] == activePos2[0][1])
									compat = 1;
							}
							if(compat != 1){
								gsl_matrix_set(topol,i,j,0.0);
								gsl_matrix_int_set(moveType,i,j,666);
							}
							else{
							
								top = 1.0;
								bottom = S->states[i]->aCounts[0] * (S->states[i]->aCounts[0]-1) /2.0 ;
								for(k=0;k<nlin1;k++){
									x = activePos1[k][0];
									y = activePos1[k][1];
									top *= (float) xchoosey(gsl_matrix_int_get(S->states[i]->popMats[0],x,y),
										gsl_matrix_int_get(delta->popMats[0],x,y));
								
								//	printf("x: %d y: %d bottom: %0.2f top:%0.2f aCounts[0]: %d delta: %d foo:%d choose: %lf\n",x,y,bottom,top,S->states[i]->aCounts[0],gsl_matrix_int_get(delta->popMats[0],x,y),gsl_matrix_int_get(S->states[i]->popMats[0],x,y),xchoosey(3,2));
								
								}
								gsl_matrix_set(topol,i,j,top/bottom);
								gsl_matrix_int_set(moveType,i,j,0);
							}
						}
					} 
					else{
						if(matrixSum(delta->popMats[1]) == 1 && countNegativeEntries(delta->popMats[0]) == 0 && S->states[j]->aCounts[1] >= 1 
							/*((matrixAllGTZ(delta->popMats[0]) || matrixAllGTZ(delta->popMats[1]))) */ ){

//							afsObjectPrint(S->states[i]);
//							printf("1///////////\n");
//							afsObjectPrint(S->states[j]);
//							printf("2///////////\n");
							entriesGreaterThan(delta->popMats[1], activePos1, &nlin1,0);//old lineages
							entriesLessThan(delta->popMats[1], activePos2, &nlin2,0);//new lineages
							if(nlin2 != 0 && nlin2 < 2){								
								compat=0;
								if(nlin1 == 1){
									if(2*activePos1[0][0] == activePos2[0][0] && 2*activePos1[0][1] == activePos2[0][1])
										compat = 1;
								}
								if(nlin1==2){
									if(activePos1[0][0] + activePos1[1][0] == activePos2[0][0] &&
										activePos1[0][1] + activePos1[1][1] == activePos2[0][1])
										compat = 1;
								}
								if(compat != 1){
									
									gsl_matrix_set(topol,i,j,0.0);
									gsl_matrix_int_set(moveType,i,j,666);
								}
								else{
									top = 1.0;
									
									bottom = S->states[i]->aCounts[1] * (S->states[i]->aCounts[1]-1) /2.0 ;
									for(k=0;k<nlin1;k++){
										x = activePos1[k][0];
										y = activePos1[k][1];

										top *= xchoosey(gsl_matrix_int_get(S->states[i]->popMats[1],x,y),
											 	gsl_matrix_int_get(delta->popMats[1],x,y));
										//	top *= gsl_matrix_int_get(S->states[i]->popMats[1],x,y);
									
									//	printf("x: %d y: %d bottom: %0.2f top:%0.2f aCounts[1]: %d delta: %d\n",x,y,bottom,top,S->states[i]->aCounts[1],gsl_matrix_int_get(delta->popMats[1],x,y));
									
									}
									
									gsl_matrix_set(topol,i,j,top/bottom);
									gsl_matrix_int_set(moveType,i,j,1);
								}
							}
						}
					}
					break;
					

				}
				afsObjectFree(delta);

			}
		}

	}
	//cleanup and free
	for(i=0;i<100;i++){
		free(activePos1[i]); 
		free(activePos2[i]); 
	}
	free(activePos1);
	free(activePos2);
}

//fills and transition matrix along with a vector of rates for use with calculating expected times
void fillTransitionMatrix(afsStateSpace *S, gsl_matrix *topol, gsl_matrix_int *moveType, gsl_matrix *transMat, gsl_vector *rates,
	double theta1, double theta2, double mig1, double mig2){
	int i, j;
	double  c0r,c1r,m0r,m1r,totR;
	
	gsl_matrix_set_zero(transMat);
	for(i=0;i<topol->size1;i++){
		c0r = S->states[i]->aCounts[0] * (S->states[i]->aCounts[0]-1) / theta1 ;
		c1r = S->states[i]->aCounts[1] * (S->states[i]->aCounts[1]-1) / theta2;
		m0r = S->states[i]->aCounts[0] * mig1 ;
		m1r = S->states[i]->aCounts[1] * mig2 ;
		totR = c0r + c1r + m0r + m1r;
		if(S->states[i]->nalleles > 1  && c0r >= 0 && c1r >= 0)
			gsl_vector_set(rates,i,1.0/totR);
		else
			gsl_vector_set(rates,i,0);
		for(j=0;j<topol->size2;j++){
		
			switch(gsl_matrix_int_get(moveType,i,j)){
			//	case 666:
			//	gsl_matrix_set(transMat,i,j,0);
				break;
				case 0: //coal pop0
				gsl_matrix_set(transMat,i,j,c0r / totR * gsl_matrix_get(topol,i,j));
				break;
				case 1: //coal pop1
				gsl_matrix_set(transMat,i,j,c1r / totR * gsl_matrix_get(topol,i,j));
				break;
				case 2: //mig pop0
				gsl_matrix_set(transMat,i,j,m0r / totR * gsl_matrix_get(topol,i,j));
				break;
				case 3: //mig pop1
				gsl_matrix_set(transMat,i,j,m1r / totR * gsl_matrix_get(topol,i,j));
			//	printf("case3 -- tot:%f c0r: %f c1r: %f m0r: %f m1r: %f res:%f\n",
			//	totR,c0r,c1r,m0r,m1r,gsl_matrix_get(transMat,i,j));
				
				break;
				case -1: //absorb state
				gsl_matrix_set(transMat,i,j, gsl_matrix_get(topol,i,j));
				break;
			}
		}
	}
}



//fills up an AFS matrix using expected times

void fillExpectedAFS(afsStateSpace *S, gsl_matrix *transMat,gsl_vector *rates, gsl_matrix *expAFS, int m1, int m2){
	int i, j, k, l, its, count;
	gsl_matrix *tmIt, *tmp, *times;
	gsl_vector *means;
	double tmpRes,max = 0.0;

	its = 30;
	//set up matrices make copy of transMat
	tmIt = gsl_matrix_alloc(transMat->size1,transMat->size2);
	tmp = gsl_matrix_alloc(transMat->size1,transMat->size2);
	
	times = gsl_matrix_alloc(transMat->size1, its);
	means = gsl_vector_alloc(its);
	gsl_matrix_set_zero(times);
	gsl_matrix_set_zero(expAFS);
	//initial state contribution to AFS; ie transMat^0
	gsl_matrix_set_identity(tmIt);
	
	for(j=0;j<S->nstates;j++){
		for(k=0;k<m1+1;k++){
			for(l=0;l<m2+1;l++){
				count = gsl_matrix_int_get(S->states[j]->popMats[0],k,l) +  gsl_matrix_int_get(S->states[j]->popMats[1],k,l);
				if(gsl_vector_get(rates,j) > 0){
					tmpRes = gsl_matrix_get(expAFS,k,l) + (count * gsl_vector_get(rates,j) * gsl_matrix_get(tmIt,0,j));										
					gsl_matrix_set(expAFS,k,l, tmpRes);
				}
			}
		}
		tmpRes = 0;
		for(k=0;k<S->nstates;k++)
			tmpRes += (gsl_vector_get(rates,k) * gsl_matrix_get(tmIt,k,j));
		gsl_matrix_set(times,j,0,tmpRes);
	}


	//now add transMat ^1
	gsl_matrix_memcpy(tmIt,transMat);
	for(j=0;j<S->nstates;j++){
		for(k=0;k<m1+1;k++){
			for(l=0;l<m2+1;l++){
				count = gsl_matrix_int_get(S->states[j]->popMats[0],k,l) +  gsl_matrix_int_get(S->states[j]->popMats[1],k,l);
				if(gsl_vector_get(rates,j) > 0){
					tmpRes = gsl_matrix_get(expAFS,k,l) + (count * gsl_vector_get(rates,j) * gsl_matrix_get(tmIt,0,j));
					gsl_matrix_set(expAFS,k,l, tmpRes);
				}
			}
		}
		tmpRes = 0;
		for(k=0;k<S->nstates;k++)
			tmpRes += ((gsl_vector_get(rates,k) + gsl_matrix_get(times,k,0)) * gsl_matrix_get(tmIt,k,j));
		gsl_matrix_set(times,j,1,tmpRes );
	}

	//now do matrix powers
	i = 2;
	while(i<its && max < 0.99){
	//	printf("power: %d\n",i);
	//	gsl_matrix_prettyPrint(tmIt);
		//compute matrix product
		gsl_matrix_set_zero(tmp);
		gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0,tmIt,transMat,0.0,tmp);
		gsl_matrix_memcpy(tmIt,tmp);
		for(j=0;j<S->nstates;j++){
			for(k=0;k<m1+1;k++){
				for(l=0;l<m2+1;l++){
					count = gsl_matrix_int_get(S->states[j]->popMats[0],k,l) +  gsl_matrix_int_get(S->states[j]->popMats[1],k,l);
					if(gsl_vector_get(rates,j) > 0){
						tmpRes = gsl_matrix_get(expAFS,k,l) + (count * gsl_vector_get(rates,j) * gsl_matrix_get(tmIt,0,j));
						gsl_matrix_set(expAFS,k,l, tmpRes);
					}
				}
			}
		//	for(k=0;k<S->nstates;k++)
		//		tmpRes += ((gsl_vector_get(rates,k) + gsl_matrix_get(times,k,i-1)) * gsl_matrix_get(tmIt,k,j)) ;
		//	gsl_matrix_set(times,j,i,tmpRes);
		}
		//printf("computed transition matrix power %d....\n",i);
		max = gsl_matrix_row_max(tmIt,0);
		i++;	
	}
//	gsl_matrix_prettyPrint(times);
	for(i=0;i<its;i++){
		tmpRes = 0.0;
		for(j=0;j<S->nstates;j++){
			tmpRes += gsl_matrix_get(times,j,i);
		}
		gsl_vector_set(means,i,tmpRes / (double) S->nstates);
	}
//	gsl_vector_fprintf(stdout,means,"%f");

//	gsl_matrix_prettyPrint(expAFS);

	tmpRes = matrixSumDouble(expAFS);
	gsl_matrix_scale(expAFS, 1.0 / tmpRes);
	gsl_matrix_free(tmp);
	gsl_matrix_free(tmIt);
	gsl_matrix_free(times);
	gsl_vector_free(means);
	
		
}


//fills up an AFS matrix using expected times and a precalculated visitMat
void fillExpectedAFS_visits(afsStateSpace *S, gsl_matrix *visitMat,gsl_vector *rates, gsl_matrix *expAFS){
	int  j, k, l, m1, m2, count;

	double tmpRes;

	m1=expAFS->size1;
	m2=expAFS->size2;
	
	gsl_matrix_set_zero(expAFS);
	
	for(j=0;j<S->nstates;j++){
		for(k=0;k<m1;k++){
			for(l=0;l<m2;l++){
				count = gsl_matrix_int_get(S->states[j]->popMats[0],k,l) +  gsl_matrix_int_get(S->states[j]->popMats[1],k,l);
				if(gsl_vector_get(rates,j) > 0){
					tmpRes = gsl_matrix_get(expAFS,k,l) + (count * gsl_vector_get(rates,j) * gsl_matrix_get(visitMat,0,j));										
					gsl_matrix_set(expAFS,k,l, tmpRes);
				}
			}
		}
	}
	tmpRes = matrixSumDouble(expAFS);
	gsl_matrix_scale(expAFS, 1.0 / tmpRes);
}


/////coalMarkovChainTopologyMatrix_2 -- using a previously calculated state space
//this fills a matrix of topological moves possible in the markov chain, 
//along with a lookup table of moves that are parameter dependant
// also returns list of nonzero entries and the numbers of such entries
//
// codes for lookup :: 
// -1 : absorb
// 0 : coal popn0
// 1 : coal popn1
// 2 : mig popn0
// 3 : mig popn2
// 666: incompatible move

int coalMarkovChainTopologyMatrix_2(afsStateSpace *S,gsl_matrix *topol, gsl_matrix_int *moveType, int *dim1, int *dim2){
	int i,j,k,steps,x,y;
	double top,bottom;
	afsObject *delta;
	int **activePos1, **activePos2; //first dimension here relates to 2 popns
	int nlin1,nlin2,compat;
	int nonZeroCount = 0;
	nlin1 = nlin2 = 0;
	if(S->nstates != topol->size1 || S->nstates != moveType->size1){
		fprintf(stderr,"Error: StateSpace and matrices are not of equal size\n");
		exit(1);
	}

	//arrays for finding positions of entries
	activePos1 = malloc(200*sizeof(int *));
	activePos2 = malloc(200*sizeof(int *));
	for(i=0;i<200;i++){
		activePos1[i] = malloc(2*sizeof(int));
		activePos1[i][0] = activePos1[i][1] = 0;
		activePos2[i] = malloc(2*sizeof(int));
		activePos2[i][0] = activePos2[i][1] = 0;

	}

//	printf("\033[2J");        /*  clear the screen  */
//	printf("\033[H");         /*  position cursor at top-left corner */
	for(i=0;i<S->nstates;i++){
		
		
		for(j=0;j<S->nstates;j++){
//			printf("%.4f percent done with constructing initial transition matrix...\n",(float)tCount++ / (float)tot);
//			sleep(0.001);
//			printf("\033[H");  
			gsl_matrix_int_set(moveType,i,j,666);

			//is ith state root state?
			if(S->states[i]->nalleles == 1){
				//set correct absorbing conditions
				if(i==j){
					gsl_matrix_set(topol,i,j,1.0);
					gsl_matrix_int_set(moveType,i,j,-1);
					dim1[nonZeroCount] = i;
					dim2[nonZeroCount] = j;
					nonZeroCount += 1;
				}
				else{
					gsl_matrix_set(topol,i,j,0.0);
					gsl_matrix_int_set(moveType,i,j,-1);
				}		
			}
			else{ //ith not root; what is the move?
				steps = S->states[i]->nalleles - S->states[j]->nalleles ;
				if(steps < 2){
					delta = afsObjectDelta(S->states[i],S->states[j]);
					switch(steps){
						case 0:
						// 0 "steps", so no changes in nalleles between two
						//check counts and positions	
						if(abs(matrixSum(delta->popMats[0])) == 1){
							nonZeroEntries(delta->popMats[0], activePos1, &nlin1);
							nonZeroEntries(delta->popMats[1], activePos2, &nlin2);

							//migration?
							if(nlin1 == nlin2 && activePos1[0][0] == activePos2[0][0] &&
							(gsl_matrix_int_min(delta->popMats[0]) >= 0 || gsl_matrix_int_min(delta->popMats[1]) >= 0))
							{

								//migration popn1?
								if(matrixSum(delta->popMats[0]) == 1){
									gsl_matrix_set(topol,i,j,
										(double) gsl_matrix_int_get(S->states[i]->popMats[0],
										activePos1[0][0], activePos1[0][1]) / S->states[i]->aCounts[0]);
									gsl_matrix_int_set(moveType,i,j,2);
									dim1[nonZeroCount] = i;
									dim2[nonZeroCount] = j;
									nonZeroCount += 1;
								}
								else{ //migration popn2
									gsl_matrix_set(topol,i,j,
										(double) gsl_matrix_int_get(S->states[i]->popMats[1],
										activePos2[0][0],activePos2[0][1]) / S->states[i]->aCounts[1]);
									gsl_matrix_int_set(moveType,i,j,3);
									dim1[nonZeroCount] = i;
									dim2[nonZeroCount] = j;
									nonZeroCount += 1;
								}
							}
						}
						break;

						case 1://coalescent?
						//coal popn1?
						if(matrixSum(delta->popMats[0]) == 1 && countNegativeEntries(delta->popMats[1]) == 0 && S->states[j]->aCounts[0] >= 1 ){
							entriesGreaterThan(delta->popMats[0], activePos1, &nlin1,0);//old lineages
							entriesLessThan(delta->popMats[0], activePos2, &nlin2,0);//new lineages
							if(nlin2 != 0 && nlin2 < 2){
							
								compat=0;
								if(nlin1 == 1){
									if(2*activePos1[0][0] == activePos2[0][0] && 2*activePos1[0][1] == activePos2[0][1])
										compat = 1;
								}
								if(nlin1==2){
									if(activePos1[0][0] + activePos1[1][0] == activePos2[0][0] &&
										activePos1[0][1] + activePos1[1][1] == activePos2[0][1])
										compat = 1;
								}
								if(compat != 1){
									gsl_matrix_set(topol,i,j,0.0);
									gsl_matrix_int_set(moveType,i,j,666);
								}
								else{
							
									top = 1.0;
									bottom = S->states[i]->aCounts[0] * (S->states[i]->aCounts[0]-1) /2.0 ;
									for(k=0;k<nlin1;k++){
										x = activePos1[k][0];
										y = activePos1[k][1];
										top *= (float) xchoosey(gsl_matrix_int_get(S->states[i]->popMats[0],x,y),
											gsl_matrix_int_get(delta->popMats[0],x,y));
									}
									gsl_matrix_set(topol,i,j,top/bottom);
									gsl_matrix_int_set(moveType,i,j,0);
									dim1[nonZeroCount] = i;
									dim2[nonZeroCount] = j;
									nonZeroCount += 1;
								}
							}
						} 
						else{
							if(matrixSum(delta->popMats[1]) == 1 && countNegativeEntries(delta->popMats[0]) == 0 && S->states[j]->aCounts[1] >= 1 ){
								entriesGreaterThan(delta->popMats[1], activePos1, &nlin1,0);//old lineages
								entriesLessThan(delta->popMats[1], activePos2, &nlin2,0);//new lineages
								if(nlin2 != 0 && nlin2 < 2){								
									compat=0;
									if(nlin1 == 1){
										if(2*activePos1[0][0] == activePos2[0][0] && 2*activePos1[0][1] == activePos2[0][1])
											compat = 1;
									}
									if(nlin1==2){
										if(activePos1[0][0] + activePos1[1][0] == activePos2[0][0] &&
											activePos1[0][1] + activePos1[1][1] == activePos2[0][1])
											compat = 1;
									}
									if(compat != 1){
									
										gsl_matrix_set(topol,i,j,0.0);
										gsl_matrix_int_set(moveType,i,j,666);
									}
									else{
										top = 1.0;
									
										bottom = S->states[i]->aCounts[1] * (S->states[i]->aCounts[1]-1) /2.0 ;
										for(k=0;k<nlin1;k++){
											x = activePos1[k][0];
											y = activePos1[k][1];

											top *= xchoosey(gsl_matrix_int_get(S->states[i]->popMats[1],x,y),
												 	gsl_matrix_int_get(delta->popMats[1],x,y));
										}
									
										gsl_matrix_set(topol,i,j,top/bottom);
										gsl_matrix_int_set(moveType,i,j,1);
										dim1[nonZeroCount] = i;
										dim2[nonZeroCount] = j;
										nonZeroCount += 1;
									}
								}
							}
						}
						break;
					}
					afsObjectFree(delta);
				}
			}
		}

	}
	//cleanup and free
	for(i=0;i<200;i++){
		free(activePos1[i]); 
		free(activePos2[i]); 
	}
	free(activePos1);
	free(activePos2);
	return(nonZeroCount);
}

/////coalMarkovChainTopologyMatrix_sparse -- using a previously calculated state space
//fills requirements for triplet style sparse matrix for the instantaneous rate matrix
// in the markov chain along with a lookup table of moves that are parameter dependant
// also returns list of nonzero entries and the numbers of such entries
//this fills a matrix of topological moves possible in the markov chain, 
// codes for lookup :: 
// -1 : absorb
// 0 : coal popn0
// 1 : coal popn1
// 2 : mig popn0
// 3 : mig popn2
// 666: incompatible move

int coalMarkovChainTopologyMatrix_sparse(afsStateSpace *S,double *topol, int *moveType, int *dim1, int *dim2){
	int i,j,k,steps,x,y;
	double top,bottom;
	afsObject *delta;
	int **activePos1, **activePos2; //first dimension here relates to 2 popns
	int nlin1,nlin2,compat;
	int nonZeroCount;
	nlin1 = nlin2 = 0;

	//arrays for finding positions of entries
	activePos1 = malloc(200*sizeof(int *));
	activePos2 = malloc(200*sizeof(int *));
	for(i=0;i<200;i++){
		activePos1[i] = malloc(2*sizeof(int));
		activePos1[i][0] = activePos1[i][1] = 0;
		activePos2[i] = malloc(2*sizeof(int));
		activePos2[i][0] = activePos2[i][1] = 0;

	}


	nonZeroCount = 0;
//	printf("\033[2J");        /*  clear the screen  */
//	printf("\033[H");         /*  position cursor at top-left corner */
	for(i=0;i<S->nstates;i++){
		for(j=0;j<S->nstates;j++){
//			printf("%.4f percent done with constructing initial transition matrix...\n",(float)tCount++ / (float)tot);
//			sleep(0.001);
//			printf("\033[H");  
			
			//is ith state root state?
			if(S->states[i]->nalleles == 1){
				//set correct absorbing conditions
				if(i==j){
					topol[nonZeroCount] = 1.0;
					moveType[nonZeroCount] = -1;
					dim1[nonZeroCount] = i;
					dim2[nonZeroCount] = j;
					nonZeroCount += 1;
				}
		
			}

			else{ //ith not root; what is the move?
				steps = S->states[i]->nalleles - S->states[j]->nalleles ;
				if(steps < 2){
					delta = afsObjectDelta(S->states[i],S->states[j]);					
					switch(steps){
						case 0:
						// 0 "steps", so no changes in nalleles between two
						//check counts and positions	
						if(abs(matrixSum(delta->popMats[0])) == 1){
							nonZeroEntries(delta->popMats[0], activePos1, &nlin1);
							nonZeroEntries(delta->popMats[1], activePos2, &nlin2);

							//migration?
							if(nlin1 == nlin2 && activePos1[0][0] == activePos2[0][0] &&
							(gsl_matrix_int_min(delta->popMats[0]) >= 0 || gsl_matrix_int_min(delta->popMats[1]) >= 0))
							{

								//migration popn1?
								if(matrixSum(delta->popMats[0]) == 1){
									topol[nonZeroCount] = (double) gsl_matrix_int_get(S->states[i]->popMats[0],
										activePos1[0][0], activePos1[0][1]) / S->states[i]->aCounts[0];
									moveType[nonZeroCount] = 2;
									dim1[nonZeroCount] = i;
									dim2[nonZeroCount] = j;
									nonZeroCount += 1;
								}
								else{ //migration popn2
									topol[nonZeroCount] = (double) gsl_matrix_int_get(S->states[i]->popMats[1],
										activePos2[0][0],activePos2[0][1]) / S->states[i]->aCounts[1];
									moveType[nonZeroCount] = 3;
									
									dim1[nonZeroCount] = i;
									dim2[nonZeroCount] = j;
									nonZeroCount += 1;
								}
							}
						}
						break;

						case 1://coalescent?
						//coal popn1?
						if(matrixSum(delta->popMats[0]) == 1 && countNegativeEntries(delta->popMats[1]) == 0 && S->states[j]->aCounts[0] >= 1 ){
							entriesGreaterThan(delta->popMats[0], activePos1, &nlin1,0);//old lineages
							entriesLessThan(delta->popMats[0], activePos2, &nlin2,0);//new lineages
							if(nlin2 != 0 && nlin2 < 2){
							
								compat=0;
								if(nlin1 == 1){
									if(2*activePos1[0][0] == activePos2[0][0] && 2*activePos1[0][1] == activePos2[0][1])
										compat = 1;
								}
								if(nlin1==2){
									if(activePos1[0][0] + activePos1[1][0] == activePos2[0][0] &&
										activePos1[0][1] + activePos1[1][1] == activePos2[0][1])
										compat = 1;
								}
								if(compat == 1){
									top = 1.0;
									bottom = S->states[i]->aCounts[0] * (S->states[i]->aCounts[0]-1) /2.0 ;
									for(k=0;k<nlin1;k++){
										x = activePos1[k][0];
										y = activePos1[k][1];
										top *= (float) xchoosey(gsl_matrix_int_get(S->states[i]->popMats[0],x,y),
											gsl_matrix_int_get(delta->popMats[0],x,y));
									}

									topol[nonZeroCount] = top/bottom;
									moveType[nonZeroCount] = 0;
									dim1[nonZeroCount] = i;
									dim2[nonZeroCount] = j;
									nonZeroCount += 1;
								}
							}
						} 
						else{
							if(matrixSum(delta->popMats[1]) == 1 && countNegativeEntries(delta->popMats[0]) == 0 && S->states[j]->aCounts[1] >= 1 ){
								entriesGreaterThan(delta->popMats[1], activePos1, &nlin1,0);//old lineages
								entriesLessThan(delta->popMats[1], activePos2, &nlin2,0);//new lineages
								if(nlin2 != 0 && nlin2 < 2){								
									compat=0;
									if(nlin1 == 1){
										if(2*activePos1[0][0] == activePos2[0][0] && 2*activePos1[0][1] == activePos2[0][1])
											compat = 1;
									}
									if(nlin1==2){
										if(activePos1[0][0] + activePos1[1][0] == activePos2[0][0] &&
											activePos1[0][1] + activePos1[1][1] == activePos2[0][1])
											compat = 1;
									}
									if(compat == 1){
										top = 1.0;
										bottom = S->states[i]->aCounts[1] * (S->states[i]->aCounts[1]-1) /2.0 ;
										for(k=0;k<nlin1;k++){
											x = activePos1[k][0];
											y = activePos1[k][1];

											top *= xchoosey(gsl_matrix_int_get(S->states[i]->popMats[1],x,y),
												 	gsl_matrix_int_get(delta->popMats[1],x,y));
										}
									
										topol[nonZeroCount] = top/bottom;
										moveType[nonZeroCount] = 1;
										dim1[nonZeroCount] = i;
										dim2[nonZeroCount] = j;
										nonZeroCount += 1;
									}
								}
							}
						}
						break;
					}
					afsObjectFree(delta);
				}
			}
		}

	}
	//cleanup and free
	for(i=0;i<200;i++){
		free(activePos1[i]); 
		free(activePos2[i]); 
	}
	free(activePos1);
	free(activePos2);
	return(nonZeroCount);
}

//readTopolFile-- this function allocs sparse matrix representation of topology and movetype matrices
// and then fills them from reading from a file.
void readTopolFile(const char *fileName, double *topA, int *moveA, int *dim1, int *dim2){
	int i, nnz;
	double tmp1;
	int tmp2,tmp3,tmp4;
	FILE *infile;
	int linebuf = 10000;
	char line[linebuf+1];

	//open file
	infile = fopen(fileName, "r");
	if (infile == NULL){
		fprintf(stderr,"Error opening infile! ARRRRR!!!!\n");
		exit(1);
	}
	//initialization
	fscanf(infile,"nnz: %d", &nnz);
	topA = malloc(sizeof(double)*nnz);
	moveA = malloc(sizeof(int)*nnz);
	dim1 = malloc(sizeof(int)*nnz);
	dim2 = malloc(sizeof(int)*nnz);
	
	
	
	for(i=0;i<nnz;i++){
		fscanf(infile,"%lf\t%d\t%d\t%d", &tmp1, &tmp2, &tmp3,&tmp4);
		topA[i] = tmp1;
		moveA[i]=tmp2;
		dim1[i]=tmp3;
		dim2[i]=tmp4;
	}
}

//fills a sparse transition matrix along with a vector of rates for use with calculating expected times
struct cs_di_sparse *fillTransitionMatrix_2(afsStateSpace *S, gsl_matrix *topol, gsl_matrix_int *moveType, gsl_vector *rates, int nzCount,
	int *dim1, int *dim2, double theta1, double theta2, double mig1, double mig2){
	int i, j,k;
	int N = topol->size1;
	double  c0r[N],c1r[N],m0r[N],m1r[N],totR[N];
	
	
	struct cs_di_sparse *triplet, *tmat;

	

	//now allocate new cs_sparse obj
	triplet = cs_spalloc(topol->size1, topol->size2, nzCount + topol->size2, 1, 1); //alloc sparse mat with extra space for nonzero identity mats
	for(i=0;i<N;i++){
		c0r[i] = S->states[i]->aCounts[0] * (S->states[i]->aCounts[0]-1) / theta1 ;
		c1r[i] = S->states[i]->aCounts[1] * (S->states[i]->aCounts[1]-1) / theta2;
		m0r[i] = S->states[i]->aCounts[0] * mig1 ;
		m1r[i] = S->states[i]->aCounts[1] * mig2 ;
		totR[i] = c0r[i] + c1r[i] + m0r[i] + m1r[i];
		if(S->states[i]->nalleles > 1  && c0r[i] >= 0 && c1r[i] >= 0)
			gsl_vector_set(rates,i,1.0/totR[i]);
		else
			gsl_vector_set(rates,i,0);
	}
	
	for(k=0;k<nzCount;k++){
		i = dim1[k];
		j= dim2[k];
		switch(gsl_matrix_int_get(moveType,i,j)){
			case 0: //coal pop0
			cs_entry(triplet,i,j,c0r[i] / totR[i] * gsl_matrix_get(topol,i,j));
			break;
			case 1: //coal pop1
			cs_entry(triplet,i,j,c1r[i] / totR[i] * gsl_matrix_get(topol,i,j));
			break;
			case 2: //mig pop0
			cs_entry(triplet,i,j,m0r[i] / totR[i] * gsl_matrix_get(topol,i,j));
			break;
			case 3: //mig pop1
			cs_entry(triplet,i,j,m1r[i] / totR[i] * gsl_matrix_get(topol,i,j));
			break;
		}
	}
	tmat = cs_compress(triplet);
	cs_spfree(triplet);
	return(tmat);

}


//fills up an AFS matrix using expected times and a precalculated visitMat
void fillExpectedAFS_visits_2(afsStateSpace *S, double *visitMat ,gsl_vector *rates, gsl_matrix *expAFS){
	int  j, k, l, m1, m2, count;

	double tmpRes;

	m1=expAFS->size1;
	m2=expAFS->size2;
	
	gsl_matrix_set_zero(expAFS);
	
	for(j=0;j<S->nstates;j++){
		for(k=0;k<m1;k++){
			for(l=0;l<m2;l++){
				count = gsl_matrix_int_get(S->states[j]->popMats[0],k,l) +  gsl_matrix_int_get(S->states[j]->popMats[1],k,l);
				tmpRes = gsl_matrix_get(expAFS,k,l) + (count * gsl_vector_get(rates,j) * visitMat[j]);	
			//	printf("tmpRes[%d][%d]: %f rates[j]: %f visits[j]: %f count: %d\n",tmpRes,k,l,gsl_vector_get(rates,j),visitMat[j],count);
				gsl_matrix_set(expAFS,k,l, tmpRes);
			}
		}
	}
	tmpRes = matrixSumDouble(expAFS);
	gsl_matrix_scale(expAFS, 1.0 / tmpRes);
}
/////////////////// Likelihood stuff 

//fillLogAFS_ptr -- wrapper around above that uses the island_like_param struct to hold all revelant info
void fillLogAFS_ptr(void * p){
// 	struct island_lik_params * params = (struct island_lik_params *) p;
// 	int i,j;
// 	
// 	//initialize some vectors
// 	gsl_vector_set_zero(params->rates);
// 	//fill transMat
// 	fillTransitionMatrix(params->stateSpace, params->topol,
// 		params->moveType,params->transMat,params->rates,
// 		gsl_vector_get(params->paramVector,0), gsl_vector_get(params->paramVector,1),
// 		gsl_vector_get(params->paramVector,2), gsl_vector_get(params->paramVector,3));
// 	gsl_matrix_set_identity(params->visitMat);
// 	gsl_matrix_sub(params->visitMat, params->transMat);
// 	gsl_matrix_invert_lapack(params->visitMat);
// 	fillExpectedAFS_visits(params->stateSpace, params->visitMat,params->rates, params->expAFS);
// //	fillExpectedAFS(params->stateSpace, params->transMat, params->rates, 
// //		params->expAFS, params->n1, params->n2);
// //	gsl_matrix_prettyPrint(params->expAFS);
// 
// 	for(i=0;i<params->expAFS->size1;i++){
// 		for(j=0;j<params->expAFS->size2;j++){
// 			gsl_matrix_set(params->expAFS,i,j,log(gsl_matrix_get(params->expAFS,i,j)));
// 		}
// 	}
//	gsl_matrix_prettyPrint(params->expAFS);
}

//fillLogAFS_ptr_sp -- uses cxsparse library for inversion
void fillLogAFS_ptr_sp(void * p){
	struct island_lik_params * params = (struct island_lik_params *) p;
	int i,j,nz;
 	struct cs_di_sparse *spMat, *mt;
	//initialize some vectors
	gsl_vector_set_zero(params->rates);
	//fill transMat
	fillTransitionMatrix(params->stateSpace, params->topol,
		params->moveType,params->transMat,params->rates,
		gsl_vector_get(params->paramVector,0), gsl_vector_get(params->paramVector,1),
		gsl_vector_get(params->paramVector,2), gsl_vector_get(params->paramVector,3));
	gsl_matrix_set_identity(params->visitMat);
	gsl_matrix_sub(params->visitMat, params->transMat);
	
	nz=0;
	for(i=0;i<params->visitMat->size1;i++){
		for(j=0;j<params->visitMat->size2;j++){
			nz += (gsl_matrix_get(params->visitMat,i,j) == 0.0) ? 0 : 1;
		}
	}
	//create unit array for solve
	(params->cs_work)[0] = 1;
	for(i=1; i<params->visitMat->size1; i++) (params->cs_work)[i] = 0;
	
	spMat = gslMatrixToCSsparse(params->visitMat,nz);
	mt = cs_transpose(spMat,1);
	cs_lusol(0, mt, params->cs_work, 1e-12);

	fillExpectedAFS_visits_2(params->stateSpace, params->cs_work,params->rates, params->expAFS);
	for(i=0;i<params->expAFS->size1;i++){
		for(j=0;j<params->expAFS->size2;j++){
			gsl_matrix_set(params->expAFS,i,j,log(gsl_matrix_get(params->expAFS,i,j)));
		}
	}
	cs_spfree(spMat);
	cs_spfree(mt);
}

//fillLogAFS_ptr_sp2 -- uses cxsparse library for inversion -- will replace above
void fillLogAFS_ptr_sp2(void * p){
	struct island_lik_params * params = (struct island_lik_params *) p;
	int i,j, N = params->stateSpace->nstates;
 	cs *spMat, *mt, *ident, *eye, *tmpMat;
	//initialize some vectors
	gsl_vector_set_zero(params->rates);
	//fill transMat
	tmpMat = fillTransitionMatrix_2(params->stateSpace, params->topol,
			params->moveType,params->rates,params->nzCount, params->dim1, params->dim2,
			1.0, gsl_vector_get(params->paramVector,0),
			gsl_vector_get(params->paramVector,1), gsl_vector_get(params->paramVector,2));
		
	//add negative ident
	ident = cs_spalloc(N,N,N,1,1);
	for(i=0;i<N;i++) cs_entry(ident,i,i,1);
	eye = cs_compress(ident);
	spMat = cs_add(eye,tmpMat,1.0,-1.0);

	mt = cs_transpose(spMat,1);
	//cs_print(mt,0);
	
	//create unit array for solve
	(params->cs_work)[0] = 1;
	for(i=1; i<params->visitMat->size1; i++) (params->cs_work)[i] = 0;
	
	cs_dropzeros(mt);
	cs_lusol(0, mt, params->cs_work, 1e-12);

	fillExpectedAFS_visits_2(params->stateSpace, params->cs_work,params->rates, params->expAFS);

	
	for(i=0;i<params->expAFS->size1;i++){
		for(j=0;j<params->expAFS->size2;j++){
			gsl_matrix_set(params->expAFS,i,j,log(gsl_matrix_get(params->expAFS,i,j)));
		}
	}
	cs_spfree(spMat);
	cs_spfree(mt);
	cs_spfree(ident);
	cs_spfree(eye);
	
//	cs_spfree(cMat);
	cs_spfree(tmpMat);
}

void fillLogAFS_ptr_sp_hmm(int state,void * p){
	struct island_lik_params * params = (struct island_lik_params *) p;
	int i,j, N = params->stateSpace->nstates;
	int migInd = (params->nstates + 1) + (2 * state);
	cs *spMat, *mt, *ident, *eye;
	struct cs_di_sparse *tmpMat;
	//initialize some vectors
	gsl_vector_set_zero(params->rates);
	//fill transMat
//	printf("params: %f %f %f\n", gsl_vector_get(params->paramVector,0), gsl_vector_get(params->paramVector,migInd), gsl_vector_get(params->paramVector,migInd+1));
	tmpMat = fillTransitionMatrix_2(params->stateSpace, params->topol,
			params->moveType,params->rates,params->nzCount, params->dim1, params->dim2,
			1.0, gsl_vector_get(params->paramVector,0),
			gsl_vector_get(params->paramVector,migInd), gsl_vector_get(params->paramVector,migInd+1)); //find right lookuyp
		
	//add negative ident
	ident = cs_spalloc(N,N,N,1,1);
	for(i=0;i<N;i++) cs_entry(ident,i,i,1);
	eye = cs_compress(ident);
	spMat = cs_add(eye,tmpMat,1.0,-1.0);
	mt = cs_transpose(spMat,1);
	
	//create unit array for solve
	(params->cs_work)[0] = 1;
	for(i=1; i<params->visitMat->size1; i++) (params->cs_work)[i] = 0;
	
	cs_dropzeros(mt);
	cs_lusol(0, mt, params->cs_work, 1e-12);

	fillExpectedAFS_visits_2(params->stateSpace, params->cs_work,params->rates, params->expAFS);

	
	for(i=0;i<params->expAFS->size1;i++){
		for(j=0;j<params->expAFS->size2;j++){
			gsl_matrix_set(params->expAFS,i,j,log(gsl_matrix_get(params->expAFS,i,j)));
		}
	}
	cs_spfree(spMat);
	cs_spfree(mt);
	cs_spfree(ident);
	cs_spfree(eye);
	cs_spfree(tmpMat);
}


 //likelihood function for HKA given the data; with constraints
// used with multimin minimizer
void island_lik_constraint(size_t n, const double *pGuess, void * p, double *result){
	struct island_lik_params * params = (struct island_lik_params *) p;
	double lik;
	int i,j;


	for(i=0;i<n;i++){
		gsl_vector_set(params->paramVector,i,pGuess[i]);
	//	printf("params[%d]: %f\n",i,gsl_vector_get(params->paramVector,i));
	}
	lik = 0.0;
	for(j=0;j < params->nstates;j++){
		//calc log probs
		fillLogAFS_ptr_sp2(p);
		for(i = 0; i < params->snpNumber; i++){
			//printf("here data[%d] p1: %d p2: %d\n",i,params->data[i].p1,params->data[i].p2);
			lik += gsl_matrix_get(params->expAFS,params->data[i].p1,params->data[i].p2);
			//printf("lik:%f %f\n",lik,gsl_matrix_get(params->expAFS,params->data[i].p1,params->data[i].p2));
		
		}
	}
//	printf("lik: %f\n",lik);
	*result=-lik;
}


void NumericalGradient_island_lik(size_t n, const double *v, void *params, double *df){
	double CurrentF, NewF ;
	int i,j, counter;
	int NumVars = n;
	double VarList[n], temp[n];
	double epsilon=1e-5;

//get initial values
	island_lik_constraint(n,v,params,&CurrentF);

	for ( i=0; i < NumVars ; i++){
		VarList[i] = v[i];
	}

	for ( counter=0; counter < NumVars ; counter++){
		for( j=0; j< NumVars; j++)
			temp[j] = VarList[j];
		
		temp[counter] +=  epsilon;
		island_lik_constraint(n,temp,params,&NewF);
		df[counter] = (NewF -  CurrentF)/epsilon;
	}

}

void island_lik_fdf(const size_t n, const double *x,void *fparams,double *fval,double *grad){
	island_lik_constraint(n,x,fparams,fval);
	NumericalGradient_island_lik(n,x, fparams, grad);
}


 //likelihood function for model given the data; with constraints
// used with multimin minimizer
void island_lik_constraint_hmm(size_t n, const double *pGuess, void * p, double *result){
	struct island_lik_params * params = (struct island_lik_params *) p;
	double lik;
	int i,j;

	for(i=0;i<n;i++){
		gsl_vector_set(params->paramVector,i,pGuess[i]);
	//	printf("params[%d]: %f\n",i,gsl_vector_get(params->paramVector,i));
	}
	lik = 0.0;
	for(j=0;j < params->nstates;j++){
		//calc log probs
		fillLogAFS_ptr_sp_hmm(j,p);
		for(i = 0; i < params->snpNumber; i++){
			
			if(params->data[i].p1==0 && params->data[i].p2==0){
		//		printf("%d there\n",i);
				lik += log( 1.0 - gsl_vector_get(params->paramVector,j+1)) *  gsl_matrix_get(params->posts,j,i);
			}
			else{
		//		printf("here data[%d] p1: %d p2: %d\n",i,params->data[i].p1,params->data[i].p2);
				lik += (log(gsl_vector_get(params->paramVector,j+1)) + gsl_matrix_get(params->expAFS,params->data[i].p1,params->data[i].p2))
				 	*  gsl_matrix_get(params->posts,j,i);
	//			printf("lik:%f %f\n",lik,gsl_matrix_get(params->expAFS,params->data[i].p1,params->data[i].p2));
			}
		}
	}
//	printf("lik: %f\n",lik);
	*result=-lik;
//	exit(1);
}


//maximize_island_lik_const-- maximizes the island_lik_constraint function using constraints via
// multimin.c addon
void maximize_island_lik_const(double *lik, void * p){
	struct island_lik_params * params = (struct island_lik_params *) p;
	int np = 3;
	int i;
	double x[np];
	double minimum=6666666.6;
	//use Nedler-Mead Simplex
	struct multimin_params optim_par = {.01,1e-4,1000,1e-4,1e-10,4,1};
	//define constraints
	unsigned type[np];
	double xmin[np];
	double xmax[np];
	for(i=0;i<np;i++){
			type[i] = 3;
			xmin[i] = 0.00;
			xmax[i] = 20;
			
	}

	for(i=0;i<np;i++){
		x[i] = gsl_vector_get(params->paramVector,i);
	}
	printf("Initial Parameter Values:\n");
	for(i=0;i<np;i++) printf("%f ",x[i]);
	printf("\n");
	multimin(np,x,&minimum,type,xmin,xmax,&island_lik_constraint,&NumericalGradient_island_lik,&island_lik_fdf,p,optim_par);
	for (i = 0; i < np; i++){
		gsl_vector_set(params->paramVector,i,x[i]);
		gsl_vector_set(params->mlParams,i,x[i]);
	//	printf("%f ",gsl_vector_get(params->paramVector,i));
	}
	*lik = -minimum;
	params->maxLik = -minimum;
	params->optimized=1;
//	printf("\n");

}


//maximize_island_lik_hmm_const-- maximizes the island_lik_constraint function using constraints via
// multimin.c addon
void maximize_island_lik_hmm_const(double *lik, void * p){
	struct island_lik_params * params = (struct island_lik_params *) p;
	int np = params->nParams;
	int i;
	double x[np];
	double minimum=666666666.6;
	//use Nedler-Mead Simplex
	struct multimin_params optim_par = {.01,1e-4,500,1e-3,1e-5,4,0};
	//define constraints
	unsigned type[np];
	double xmin[np];
	double xmax[np];
	//N2 constraint
	type[0] = 3; xmin[0] = 0.01; xmax[0] = 10.0;
	// Theta / site constraints
	for(i=1;i<=params->nstates;i++){
			type[i] = 3;
			xmin[i] = 0.0001;
			xmax[i] = 0.8;
	}
	//migration params
	for(i=params->nstates + 1;i<np;i++){
			type[i] = 3;
			xmin[i] = 0.00001;
			xmax[i] = 20;
			
	}

	for(i=0;i<np;i++)
		x[i] = gsl_vector_get(params->paramVector,i);
	
		
	multimin((size_t)np,x,&minimum,type,xmin,xmax,&island_lik_constraint_hmm,NULL,NULL,p,optim_par);
	for (i = 0; i < np; i++){
		gsl_vector_set(params->paramVector,i,x[i]);
//		printf("%f ",x[i]);
	}
	*lik = -minimum;


}


 //likelihood function for island model given the data; with constraints
// used with multimin minimizer; gets the confidence interval point
void island_lik_constraint_CI(size_t n, const double *pGuess, void * p, double *result){
	struct island_lik_params * params = (struct island_lik_params *) p;
	double lik;
	int i,j;


	for(i=0;i<n;i++){
		gsl_vector_set(params->paramVector,i,pGuess[i]);
	//	printf("params[%d]: %f\n",i,gsl_vector_get(params->paramVector,i));
	}
	lik = 0.0;
	for(j=0;j < params->nstates;j++){
		//calc log probs
		fillLogAFS_ptr_sp2(p);
		for(i = 0; i < params->snpNumber; i++){
			//printf("here data[%d] p1: %d p2: %d\n",i,params->data[i].p1,params->data[i].p2);
			lik += gsl_matrix_get(params->expAFS,params->data[i].p1,params->data[i].p2);
			//printf("lik:%f %f\n",lik,gsl_matrix_get(params->expAFS,params->data[i].p1,params->data[i].p2));

		}
	}
//	printf("lik: %f\n",lik);

	*result= lik - params->maxLik - 1.92; 
}



//maximize_island_lik_CIUpper-- maximizes the island_lik_constraint function using constraints via
// multimin.c addon
void maximize_island_lik_CIUpper(double *lik, void * p){
	struct island_lik_params * params = (struct island_lik_params *) p;
	int np = 3;
	int i;
	double x[np];
	double minimum=6666666.6;
	//use Nedler-Mead Simplex
	struct multimin_params optim_par = {.01,1e-4,1000,1e-4,1e-10,4,0};
	//define constraints
	unsigned type[np];
	double xmin[np];
	double xmax[np];
	if(params->optimized == 1){
		for(i=0;i<np;i++){
				type[i] = 3;
				xmin[i] = gsl_vector_get(params->mlParams,i);
				xmax[i] = 40;
			
		}

		for(i=0;i<np;i++){
			x[i] = gsl_vector_get(params->mlParams,i) + 0.01;
			printf("%f %f\n",x[i],xmin[i]);
		}
		printf("Upper Initial Parameter Values:\n");
		for(i=0;i<np;i++) printf("%f ",x[i]);
		printf("\n");
		multimin(np,x,&minimum,type,xmin,xmax,&island_lik_constraint_CI,&NumericalGradient_island_lik,&island_lik_fdf,p,optim_par);
		for (i = 0; i < np; i++){
			gsl_vector_set(params->paramCIUpper,i,x[i]);
		//	printf("%f ",gsl_vector_get(params->paramVector,i));
		}
		*lik = -minimum;
	}
	else{
		printf("can't compute CIs on non-optimized params\n");
	}
//	printf("\n");

}

// multimin.c addon
void maximize_island_lik_CILower(double *lik, void * p){
	struct island_lik_params * params = (struct island_lik_params *) p;
	int np = 3;
	int i;
	double x[np];
	double minimum=6666666.6;
	//use Nedler-Mead Simplex
	struct multimin_params optim_par = {.1,1e-4,1000,1e-4,1e-10,4,0};
	//define constraints
	unsigned type[np];
	double xmin[np];
	double xmax[np];
	if(params->optimized == 1){
		for(i=0;i<np;i++){
				type[i] = 3;
				xmin[i] = 0.00001;
				xmax[i] = gsl_vector_get(params->mlParams,i);
			
		}

		for(i=0;i<np;i++){
			x[i] = gsl_vector_get(params->mlParams,i) - 0.5;
		}
		printf("Lower Initial Parameter Values:\n");
		for(i=0;i<np;i++) printf("%f ",x[i]);
		printf("\n");
		multimin(np,x,&minimum,type,xmin,xmax,&island_lik_constraint_CI,&NumericalGradient_island_lik,&island_lik_fdf,p,optim_par);
		for (i = 0; i < np; i++){
			gsl_vector_set(params->paramCILower,i,x[i]);
		//	printf("%f ",gsl_vector_get(params->paramVector,i));
		}
		*lik = -minimum;
	}
	else{
		printf("can't compute CIs on non-optimized params\n");
	}
//	printf("\n");

}

void itoa(int n, char s[])
{
    int i, sign;

    if ((sign = n) < 0)  /* record sign */
        n = -n;          /* make n positive */
    i = 0;
    do {       /* generate digits in reverse order */
        s[i++] = n % 10 + '0';   /* get next digit */
    } while ((n /= 10) > 0);     /* delete it */
    if (sign < 0)
        s[i++] = '-';
    s[i] = '\0';
    reverse(s);
}

 /* reverse:  reverse string s in place */
 void reverse(char s[])
 {
     int i, j;
     char c;
 
     for (i = 0, j = strlen(s)-1; i<j; i++, j--) {
         c = s[i];
         s[i] = s[j];
         s[j] = c;
     }
 }

/*
////////////////////////////////////////////////////////
//Population splitting stuff
//fills up an AFS matrix using expected times
// Currently Does not work....

void fillExpectedAFS_IM(afsStateSpace *S, gsl_matrix *transMat,gsl_vector *rates, gsl_matrix *expAFS, int m1, int m2){
	int i, j, k, l, its, count;
	gsl_matrix *tmIt, *tmp, *times;
	
	double tmpRes, meanTime;

	its = 20;
	
	//set up matrices make copy of transMat
	tmIt = gsl_matrix_alloc(transMat->size1,transMat->size2);
	tmp = gsl_matrix_alloc(transMat->size1,transMat->size2);
	
	times = gsl_matrix_alloc(transMat->size1, its);
	means = gsl_vector_alloc(its);
	gsl_matrix_set_zero(times);
	
	//initial state contribution to AFS; ie transMat^0
	gsl_matrix_set_identity(tmIt);
	meanTime = 0;
	for(j=0;j<S->nstates;j++){
		for(k=0;k<m1+1;k++){
			for(l=0;l<m2+1;l++){
				count = gsl_matrix_int_get(S->states[j]->popMats[0],k,l) +  gsl_matrix_int_get(S->states[j]->popMats[1],k,l);
				if(gsl_vector_get(rates,j) > 0){
					tmpRes = gsl_matrix_get(expAFS,k,l) + (count * gsl_vector_get(rates,j) * gsl_matrix_get(tmIt,0,j));
					gsl_matrix_set(expAFS,k,l, tmpRes);
				}
			}
		}
		tmpRes = 0;
		for(k=0;k<S->nstates;k++)
			tmpRes += (gsl_vector_get(rates,k) * gsl_matrix_get(tmIt,k,j));
		gsl_matrix_set(times,j,0,tmpRes);
		meanTime+=tmpRes;
	}
	meanTime /= S->nstates;

	//now add transMat ^1
	gsl_matrix_memcpy(tmIt,transMat);
	for(j=0;j<S->nstates;j++){
		for(k=0;k<m1+1;k++){
			for(l=0;l<m2+1;l++){
				count = gsl_matrix_int_get(S->states[j]->popMats[0],k,l) +  gsl_matrix_int_get(S->states[j]->popMats[1],k,l);
				if(gsl_vector_get(rates,j) > 0){
					tmpRes = gsl_matrix_get(expAFS,k,l) + (count * gsl_vector_get(rates,j) * gsl_matrix_get(tmIt,0,j));
					gsl_matrix_set(expAFS,k,l, tmpRes);
				}
			}
		}
		tmpRes = 0;
		for(k=0;k<S->nstates;k++)
			tmpRes += ((gsl_vector_get(rates,k) + gsl_matrix_get(times,k,0)) * gsl_matrix_get(tmIt,k,j));
		gsl_matrix_set(times,j,1,tmpRes );
	}

	//now do matrix powers
	for(i=2;i<its;i++){
	//	printf("power: %d\n",i);
	//	gsl_matrix_prettyPrint(tmIt);
		//compute matrix product
		gsl_matrix_set_zero(tmp);
		gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0,tmIt,transMat,0.0,tmp);
		gsl_matrix_memcpy(tmIt,tmp);
		for(j=0;j<S->nstates;j++){
			for(k=0;k<m1+1;k++){
				for(l=0;l<m2+1;l++){
					count = gsl_matrix_int_get(S->states[j]->popMats[0],k,l) +  gsl_matrix_int_get(S->states[j]->popMats[1],k,l);
					if(gsl_vector_get(rates,j) > 0){
						tmpRes = gsl_matrix_get(expAFS,k,l) + (count * gsl_vector_get(rates,j) * gsl_matrix_get(tmIt,0,j));
						gsl_matrix_set(expAFS,k,l, tmpRes);
					}
				}
			}
			for(k=0;k<S->nstates;k++)
				tmpRes += ((gsl_vector_get(rates,k) + gsl_matrix_get(times,k,i-1)) * gsl_matrix_get(tmIt,k,j)) ;
			gsl_matrix_set(times,j,i,tmpRes);
		}
		
	}

	for(i=0;i<its;i++){
		tmpRes = 0.0;
		for(j=0;j<S->nstates;j++){
			tmpRes += gsl_matrix_get(times,j,i);
		}
		gsl_vector_set(means,i,tmpRes / (double) S->nstates);
	}
//	gsl_vector_fprintf(stdout,means,"%f");
	
//	exit(1);
//	gsl_matrix_prettyPrint(expAFS);

	tmpRes = matrixSumDouble(expAFS);
	gsl_matrix_scale(expAFS, 1.0 / tmpRes);
	gsl_matrix_free(tmp);
	gsl_matrix_free(tmIt);
	
		
}


void overlayPopnSplitting(afsStateSpace *S,gsl_matrix *topol, gsl_matrix_int *moveType, int n1, int n2){
	int i,j,k,steps,x,y;
	double tmp,top,bottom;
	afsObject *delta;

	int nlin1,nlin2, position,compat;



	gsl_matrix_int_fprintf(stdout,moveType,"%d");	
	
	for(i=0;i<S->nstates;i++){
		for(j=0;j<S->nstates;j++){
			//is ith state root state?
			if(S->states[i]->nalleles != 1){
				
				delta = afsObjectDelta(S->states[i],S->states[j]);
				steps = S->states[i]->nalleles - S->states[j]->nalleles ;
			//	printf("i:%d j:%d steps:%d\n",i,j,steps);
				if(steps==0 && abs(matrixSum(delta->popMats[0])) >= 1){
					// 0 "steps", so no changes in nalleles between two
					//check counts and positions
					
					if((S->states[j]->aCounts[1] == 0 || S->states[j]->aCounts[0] == 0)
					 	&& (S->states[i]->aCounts[0] > 0 && S->states[i]->aCounts[1] > 0 ) ){
						//population split or potentially migration compatible
						printf("here\n");
						switch(gsl_matrix_int_get(moveType,i,j)){
							case 666:
							gsl_matrix_int_set(moveType,i,j,4); // case 4 indicates only split possible
							break;
							case 2: //mig popn0
							gsl_matrix_int_set(moveType,i,j,5); // case 5 indicates mig popn0 or split possible
							break;
							case 3: //mig popn1
							gsl_matrix_int_set(moveType,i,j,6); // case 6 indicates mig popn0 or split possible
							break;
						}
					}
				}
			}
		}
					
	}

	gsl_matrix_int_fprintf(stdout,moveType,"%d");	
}

*/
