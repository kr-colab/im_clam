
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_matrix.h> 
#include <gsl/gsl_cblas.h> 
#include <gsl/gsl_sf.h> 
#include "adkGSL.h"
#include "AFS.h"
#include "cs.h"
#include "adkCSparse.h"
#include "AFS_ctmc.h"
#include <pthread.h>


void coalMarkovChainTopologyMatrix_pthread(afsStateSpace *S,gsl_matrix *topol, gsl_matrix_int *moveType, int *dim1, int *dim2, int start, int stop, int *nnz);
void coalMarkovChainTopologyMatrix_sparse_pthread(afsStateSpace *S,double *topol, int *moveType, int *dim1, int *dim2, int start, int stop, int *nnz);


typedef struct matrixThreadObject
{
	int nnz,start, stop;
	int *dim1,*dim2; 
	struct island_lik_params *island_params;
}
matrixThreadObject;


typedef struct matrixThreadObject_sparse
{
	int nnz,start, stop, nnzA;
	int *dim1,*dim2, *moveType;
	double *topol;
	struct im_lik_params *im_params;
}
matrixThreadObject_sparse;


//gonna try to use a global declared in island.c
void *threadDoRows(void* arg){
	matrixThreadObject *obj = (matrixThreadObject*)(arg);
	
	coalMarkovChainTopologyMatrix_pthread(obj->island_params->stateSpace,obj->island_params->topol, obj->island_params->moveType, 
	obj->dim1,obj->dim2,obj->start,obj->stop,&(obj->nnz));
	//obj->nnz -= 1;
	
	pthread_exit((void* )obj);
}

//sparse matrix version
void *threadDoRows_sparse(void* arg){
	matrixThreadObject_sparse *obj = (matrixThreadObject_sparse*)(arg);
	
	coalMarkovChainTopologyMatrix_sparse_pthread(obj->im_params->stateSpace,obj->topol, obj->moveType, 
	obj->dim1,obj->dim2,obj->start,obj->stop,&(obj->nnz));
	//obj->nnz -= 1;
	
	pthread_exit((void* )obj);
}

//sparse matrix version on reduced state space
void *threadDoRows_sparse_reduced(void* arg){
	matrixThreadObject_sparse *obj = (matrixThreadObject_sparse*)(arg);
	
	coalMarkovChainTopologyMatrix_sparse_pthread(obj->im_params->reducedStateSpace,obj->topol, obj->moveType, 
	obj->dim1,obj->dim2,obj->start,obj->stop,&(obj->nnzA));
	//obj->nnz -= 1;
	
	pthread_exit((void* )obj);
}

void coalMarkovChainTopologyMatrix_pthread(afsStateSpace *S,gsl_matrix *topol, gsl_matrix_int *moveType, int *dim1, int *dim2, int start, int stop, int *nnz){
	int i,j,k,steps,x,y;
	double top,bottom;
	afsObject *delta;
	int **activePos1, **activePos2; //first dimension here relates to 2 popns
	int nlin1,nlin2,compat;
	int nonZeroCount, tot, tCount;
	int size = 200;
	nlin1 = nlin2 = 0;
	nonZeroCount = 0;
	if(S->nstates != topol->size1 || S->nstates != moveType->size1){
		fprintf(stderr,"Error: StateSpace and matrices are not of equal size\n");
		exit(1);
	}
	*nnz = 0;
	//arrays for finding positions of entries
	activePos1 = malloc(size*sizeof(int *));
	activePos2 = malloc(size*sizeof(int *));
	for(i=0;i<size;i++){
		activePos1[i] = malloc(2*sizeof(int));
		activePos1[i][0] = activePos1[i][1] = 0;
		activePos2[i] = malloc(2*sizeof(int));
		activePos2[i][0] = activePos2[i][1] = 0;

	}

	tot = S->nstates * S->nstates;
	tCount = 0;
	delta = afsObjectNewFrom(S->states[0]);
	for(i=start;i<stop;i++){
		for(j=0;j<S->nstates;j++){
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
				//	delta = afsObjectDelta(S->states[i],S->states[j]);
					afsObjectDeltaPre(S->states[i],S->states[j],delta);
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
					//afsObjectFree(delta);
				}
			}
		}

	}
	//cleanup and free
	for(i=0;i<size;i++){
		free(activePos1[i]); 
		free(activePos2[i]); 
	}
	free(activePos1);
	free(activePos2);
	*nnz = nonZeroCount;
}

void coalMarkovChainTopologyMatrix_sparse_pthread(afsStateSpace *S,double *topol, int *moveType, int *dim1, int *dim2, int start, int stop, int *nnz){
	int i,j,k,steps,x,y;
	double top,bottom;
	afsObject *delta;
	int **activePos1, **activePos2; //first dimension here relates to 2 popns
	int nlin1,nlin2,compat;
	int nonZeroCount, tot, tCount;
	int size = 400;
	nlin1 = nlin2 = 0;

	nonZeroCount = 0;
	*nnz = 0;
	//arrays for finding positions of entries
	activePos1 = malloc(size*sizeof(int *));
	activePos2 = malloc(size*sizeof(int *));
	for(i=0;i<size;i++){
		activePos1[i] = malloc(2*sizeof(int));
		activePos1[i][0] = activePos1[i][1] = 0;
		activePos2[i] = malloc(2*sizeof(int));
		activePos2[i][0] = activePos2[i][1] = 0;
	}

	tot = S->nstates * S->nstates;
	tCount = 0;
	delta = afsObjectNewFrom(S->states[0]);
	for(i=start;i<stop;i++){
		for(j=0;j<S->nstates;j++){
		//	printf("nonzerocount=%d\n",nonZeroCount);
			
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
								(gsl_matrix_int_min(delta->popMats[0]) >= 0 || gsl_matrix_int_min(delta->popMats[1]) >= 0)){

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
					//afsObjectFree(delta);
				}
			}
		}
	}
	//cleanup and free
	for(i=0;i<size;i++){
		free(activePos1[i]); 
		free(activePos2[i]); 
	}
	free(activePos1);
	free(activePos2);
	*nnz = nonZeroCount;
}