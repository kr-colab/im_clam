#include <petscksp.h>
#include "AFS.h"
#include "AFS_ctmc_petsc.h"
#include "im_clam.h"
#include "sparseExp.h"
#include "cs.h"
#include "AFS_ctmc.h"
#include "nlopt.h"
#include "adkGSL.h"
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_blas.h>
//Functions for jointly computing CTMC transitions and embedded chain

//this function fills two petsc matrices: 1) embedded Discrete MC trans Mat, 2) CTMC trans mat
// also fills the rates vector

void fillPetscTransMats(afsStateSpace *S, double *topol, int *moveType, int *nzCount,
	int *dim1, int *dim2, double theta1, double theta2, double mig1, double mig2, Mat *DTMC, Mat *CTMC, gsl_vector *rates){
	int i, j,k, newnz, abcount=0;
	int N = S->nstates,abFlag;
	double  c0r[N],c1r[N],m0r[N],m1r[N],totR[N],rowSums[N], tmp;;
	
	gsl_vector_set_zero(rates);
	//set rates vector
	for(i=0;i<N;i++){
		c0r[i] = S->states[i]->aCounts[0] * (S->states[i]->aCounts[0]-1)  / theta1 ;
		if(theta2 == 0)
			c1r[i] = 0.0;
		else
			c1r[i] = S->states[i]->aCounts[1] * (S->states[i]->aCounts[1]-1)  / theta2;
	
		m0r[i] = S->states[i]->aCounts[0] * mig1 ;
		m1r[i] = S->states[i]->aCounts[1] * mig2 ;
		totR[i] = c0r[i] + c1r[i] + m0r[i] + m1r[i];
	
		
		rowSums[i] = 0.0;
		if(S->states[i]->nalleles > 1  && c0r[i] >= 0 && c1r[i] >= 0 && totR[i] != 0.0)
			gsl_vector_set(rates,i,1.0/totR[i]);
		else
			gsl_vector_set(rates,i,0);
	}

	newnz=0;
	for(k=0;k<*nzCount;k++){
		i = dim1[k];
		j= dim2[k];
//		printf("i: %d, j:%d, movetype: %d\n",i,j,moveType[k]);
		if(S->states[i]->nalleles ==1 || S->states[j]->nalleles ==1 )abFlag=1;
	//	if(S->states[i]->nalleles ==1 )abFlag=1;
		else abFlag = 0;
		switch(moveType[k]){
			case 0: //coal pop0
			if(totR[i]==0)tmp=0;
			else tmp = c0r[i] / totR[i] * topol[k];

			MatSetValue(*CTMC,i,j,tmp*totR[i],INSERT_VALUES);	
			rowSums[i] += tmp* totR[i];
			if(abFlag==0)MatSetValue(*DTMC,i,j,tmp,INSERT_VALUES);				
			break;
			case 1: //coal pop1
			if(totR[i]==0)tmp=0;
			else tmp = c1r[i] / totR[i] * topol[k];
			
			MatSetValue(*CTMC,i,j,tmp*totR[i],INSERT_VALUES);	
			rowSums[i] += tmp* totR[i];
			if(abFlag==0)MatSetValue(*DTMC,i,j,tmp,INSERT_VALUES);			
			break;
			case 2: //mig pop0
			if(totR[i]==0)tmp=0;
			else tmp = m0r[i] / totR[i] * topol[k];
		
			MatSetValue(*CTMC,i,j,tmp*totR[i],INSERT_VALUES);	
			rowSums[i] += tmp* totR[i];
			if(abFlag==0)MatSetValue(*DTMC,i,j,tmp,INSERT_VALUES);			
			break;
			case 3: //mig pop1
			if(totR[i]==0)tmp=0;
			else tmp = m1r[i] / totR[i] * topol[k];			
		
			MatSetValue(*CTMC,i,j,tmp*totR[i],INSERT_VALUES);	
			rowSums[i] += tmp* totR[i];
			if(abFlag==0)MatSetValue(*DTMC,i,j,tmp,INSERT_VALUES);		
			break;
			case -1:
			abcount+=1;
			break;			
		}
	}
	newnz=*nzCount ;
	//now add diagonal elements; topol pre-set to -1
	for(k=0;k<N;k++){
		if((rowSums[k]) != 0.0){
			MatSetValue(*CTMC,k,k,0.0-rowSums[k],INSERT_VALUES);				
			newnz++;
		}
		MatSetValue(*DTMC,k,k,0.0,INSERT_VALUES);
	}
	*nzCount = newnz;
	/* Assemble the matrix */
  	MatAssemblyBegin(*CTMC,MAT_FINAL_ASSEMBLY);
  	MatAssemblyEnd(*CTMC,MAT_FINAL_ASSEMBLY);
  	MatAssemblyBegin(*DTMC,MAT_FINAL_ASSEMBLY);
  	MatAssemblyEnd(*DTMC,MAT_FINAL_ASSEMBLY);
}

//fillPetscCsparseTransMats--this function fills one petsc matrix with CTMC transmat and returns pointer to CSparse matrix for DTMC trans mat
// also fills the rates vector
//note that the DTMC matrix needs to be cleared before entry to this function and compressed afterward.

struct cs_di_sparse * fillPetscCsparseTransMats(afsStateSpace *S, double *topol, int *moveType, int *nzCount,
	int *dim1, int *dim2, double theta1, double theta2, double mig1, double mig2, Mat *CTMC, gsl_vector *rates){
	int i, j,k, newnz, abcount=0;
	int N = S->nstates,abFlag;
	double  c0r[N],c1r[N],m0r[N],m1r[N],totR[N],rowSums[N], tmp;;
	struct cs_di_sparse *triplet, *tmat;
	
	//allocate new cs_sparse obj
	triplet = cs_spalloc(N, N, *nzCount + N , 1, 1); //alloc sparse mat with extra space for nonzero identity mats
	
	gsl_vector_set_zero(rates);
	//set rates vector
	for(i=0;i<N;i++){
		c0r[i] = S->states[i]->aCounts[0] * (S->states[i]->aCounts[0]-1)  / theta1 ;
		if(theta2 == 0)
			c1r[i] = 0.0;
		else
			c1r[i] = S->states[i]->aCounts[1] * (S->states[i]->aCounts[1]-1)  / theta2;
		m0r[i] = S->states[i]->aCounts[0] * mig1 ;
		m1r[i] = S->states[i]->aCounts[1] * mig2 ;
		totR[i] = c0r[i] + c1r[i] + m0r[i] + m1r[i];
	
		
		rowSums[i] = 0.0;
		if(S->states[i]->nalleles > 1  && c0r[i] >= 0 && c1r[i] >= 0 && totR[i] != 0.0)
			gsl_vector_set(rates,i,1.0/totR[i]);
		else
			gsl_vector_set(rates,i,0);
	}

	newnz=0;
	for(k=0;k<*nzCount;k++){
		i = dim1[k];
		j= dim2[k];
//		printf("i: %d, j:%d, movetype: %d\n",i,j,moveType[k]);
		if(S->states[i]->nalleles ==1 || S->states[j]->nalleles ==1 )abFlag=1;
	//	if(S->states[i]->nalleles ==1 )abFlag=1;
		else abFlag = 0;
		switch(moveType[k]){
			case 0: //coal pop0
			if(totR[i]==0)tmp=0;
			else tmp = c0r[i] / totR[i] * topol[k];
		//	printf("here--- i:%d j:%d tmp*totR[i]:%f",i,j,tmp*totR[i]);
			MatSetValue(*CTMC,i,j,tmp*totR[i],INSERT_VALUES);	
			rowSums[i] += tmp* totR[i];
			if(abFlag==0)cs_entry(triplet,i,j,tmp);				
			break;
			case 1: //coal pop1
			if(totR[i]==0)tmp=0;
			else tmp = c1r[i] / totR[i] * topol[k];
			
			MatSetValue(*CTMC,i,j,tmp*totR[i],INSERT_VALUES);	
			rowSums[i] += tmp* totR[i];
			if(abFlag==0)cs_entry(triplet,i,j,tmp);			
			break;
			case 2: //mig pop0
			if(totR[i]==0)tmp=0;
			else tmp = m0r[i] / totR[i] * topol[k];
		
			MatSetValue(*CTMC,i,j,tmp*totR[i],INSERT_VALUES);	
			rowSums[i] += tmp* totR[i];
			if(abFlag==0)cs_entry(triplet,i,j,tmp);			
			break;
			case 3: //mig pop1
			if(totR[i]==0)tmp=0;
			else tmp = m1r[i] / totR[i] * topol[k];			
		
			MatSetValue(*CTMC,i,j,tmp*totR[i],INSERT_VALUES);	
			rowSums[i] += tmp* totR[i];
			if(abFlag==0)cs_entry(triplet,i,j,tmp);		
			break;
			case -1:
			abcount+=1;
			break;			
		}
	}
	newnz=*nzCount ;
	//now add diagonal elements; topol pre-set to -1
	for(k=0;k<N;k++){
		if((rowSums[k]) != 0.0){
			MatSetValue(*CTMC,k,k,0.0-rowSums[k],INSERT_VALUES);				
			newnz++;
		}
	}
	*nzCount = newnz;
	/* Assemble the matrix */
  	MatAssemblyBegin(*CTMC,MAT_FINAL_ASSEMBLY);
  	MatAssemblyEnd(*CTMC,MAT_FINAL_ASSEMBLY);
	tmat = cs_compress(triplet);
	cs_spfree(triplet);
	return(tmat);
}

//fillPetscCsparseTransMats_prealloc--this function fills one petsc matrix with CTMC transmat and returns pointer to CSparse matrix for DTMC trans mat
// also fills the rates vector
struct cs_di_sparse * fillPetscCsparseTransMats_prealloc(afsStateSpace *S, double *topol, int *moveType, int *nzCount,
	int *dim1, int *dim2, double theta1, double theta2, double mig1, double mig2,
	Mat *CTMC, gsl_vector *rates, struct cs_di_sparse *triplet){
		
	int i, j,k, newnz, abcount=0;
	int N = S->nstates,abFlag;
	double  *c0r,*c1r,*m0r,*m1r,*totR,*rowSums, tmp;
	struct cs_di_sparse *tmat;
	
	c0r = malloc(sizeof(double)*N);
	c1r = malloc(sizeof(double)*N);
	m0r = malloc(sizeof(double)*N);
	m1r = malloc(sizeof(double)*N);
	totR = malloc(sizeof(double)*N);
	rowSums = malloc(sizeof(double)*N);
	
	gsl_vector_set_zero(rates);
	//set rates vector
	for(i=0;i<N;i++){
		c0r[i] = S->states[i]->aCounts[0] * (S->states[i]->aCounts[0]-1)  / theta1 ;
		if(theta2 == 0)
			c1r[i] = 0.0;
		else
			c1r[i] = S->states[i]->aCounts[1] * (S->states[i]->aCounts[1]-1)  / theta2;
		m0r[i] = S->states[i]->aCounts[0] * mig1 ;
		m1r[i] = S->states[i]->aCounts[1] * mig2 ;
		totR[i] = c0r[i] + c1r[i] + m0r[i] + m1r[i];
	
		
		rowSums[i] = 0.0;
		if(S->states[i]->nalleles > 1  && c0r[i] >= 0 && c1r[i] >= 0 && totR[i] != 0.0)
			gsl_vector_set(rates,i,1.0/totR[i]);
		else
			gsl_vector_set(rates,i,0);
	}

	newnz=0;
	for(k=0;k<*nzCount;k++){
		i = dim1[k];
		j= dim2[k];
//		printf("i: %d, j:%d, movetype: %d\n",i,j,moveType[k]);
		if(S->states[i]->nalleles ==1 || S->states[j]->nalleles ==1 )abFlag=1;
	//	if(S->states[i]->nalleles ==1 )abFlag=1;
		else abFlag = 0;
		switch(moveType[k]){
			case 0: //coal pop0
			if(totR[i]==0)tmp=0;
			else tmp = c0r[i] / totR[i] * topol[k];
		//	printf("here--- i:%d j:%d tmp*totR[i]:%f",i,j,tmp*totR[i]);
			MatSetValue(*CTMC,i,j,tmp*totR[i],INSERT_VALUES);	
			rowSums[i] += tmp* totR[i];
			if(abFlag==0)cs_entry(triplet,i,j,tmp);				
			break;
			case 1: //coal pop1
			if(totR[i]==0)tmp=0;
			else tmp = c1r[i] / totR[i] * topol[k];
			
			MatSetValue(*CTMC,i,j,tmp*totR[i],INSERT_VALUES);	
			rowSums[i] += tmp* totR[i];
			if(abFlag==0)cs_entry(triplet,i,j,tmp);			
			break;
			case 2: //mig pop0
			if(totR[i]==0)tmp=0;
			else tmp = m0r[i] / totR[i] * topol[k];
		
			MatSetValue(*CTMC,i,j,tmp*totR[i],INSERT_VALUES);	
			rowSums[i] += tmp* totR[i];
			if(abFlag==0)cs_entry(triplet,i,j,tmp);			
			break;
			case 3: //mig pop1
			if(totR[i]==0)tmp=0;
			else tmp = m1r[i] / totR[i] * topol[k];			
		
			MatSetValue(*CTMC,i,j,tmp*totR[i],INSERT_VALUES);	
			rowSums[i] += tmp* totR[i];
			if(abFlag==0)cs_entry(triplet,i,j,tmp);		
			break;
			case -1:
			abcount+=1;
			break;			
		}
	}
	newnz=*nzCount ;
	//now add diagonal elements; topol pre-set to -1
	for(k=0;k<N;k++){
		if((rowSums[k]) != 0.0){
			MatSetValue(*CTMC,k,k,0.0-rowSums[k],INSERT_VALUES);				
			newnz++;
		}
	}
	*nzCount = newnz;
	/* Assemble the matrix */
  	MatAssemblyBegin(*CTMC,MAT_FINAL_ASSEMBLY);
  	MatAssemblyEnd(*CTMC,MAT_FINAL_ASSEMBLY);
	tmat = cs_compress(triplet);
	
	free(c0r);
	free(c1r);
	free(m0r);
	free(m1r);
	free(totR);
	free(rowSums);	
	return(tmat);
}

//calcLogAFS_IM -- returns the expects log AFS from an IM model
void calcLogAFS_IM(void * p){
	struct clam_lik_params * params = (struct clam_lik_params *) p;
 	PetscInt i,j, N = params->stateSpace->nstates;
  	cs *spMat, *mt, *tmpMat ;
  	double timeV, sum, thetaA, theta2, mig1, mig2;
  	PetscInt Na;
  	gsl_vector *tmpStates;
 	css *S ;
	csn *NN ;
	PetscInt n,Ntmp ;
	double *xx;
 	PetscInt        iStart,iEnd, *idx;
	const PetscInt *idx2;
	PetscScalar    *tmpArray;
	const PetscScalar *tmpArrayC;
	PetscErrorCode ierr;
	MFNConvergedReason reason;
	
 	//initialize some vectors
 	gsl_vector_set_zero(params->rates);
	MatZeroEntries(params->C);
	MatZeroEntries(params->C2);
 	tmpStates = gsl_vector_alloc(N);
	Ntmp=N;
	Na=params->Na;
 	//For straight MLE the paramVector takes the form [N2,NA,m1,m2,t]
 	theta2 = gsl_vector_get(params->paramVector,0);
 	thetaA = gsl_vector_get(params->paramVector,1);
 	mig1 = gsl_vector_get(params->paramVector,2);
 	mig2 = gsl_vector_get(params->paramVector,3);
 	timeV = gsl_vector_get(params->paramVector,4);
// 	printf("params-> %f %f %f %f %f\n",theta2,thetaA,mig1,mig2,timeV);
// 	printf("nnz %d nnzA%d \n",params->nnz,params->nnzA);
 
	//fill transMat
 //	tmpMat = fillPetscCsparseTransMats_prealloc(params->stateSpace, params->top, params->move, &params->nnz,
 //		params->dim1, params->dim2, 1, theta2, mig1, mig2,  &params->C, params->rates,params->triplet);
 
	tmpMat = fillPetscCsparseTransMats(params->stateSpace, params->top, params->move, &params->nnz,
	 		params->dim1, params->dim2, 1, theta2, mig1, mig2,  &params->C, params->rates);	
	//using CSparse
	// S = (I-P)^-1
	//add negative ident
//	ident = cs_spalloc(N,N,N,1,1);
//	for(i=0;i<N;i++) cs_entry(ident,i,i,1);
//	eye = cs_compress(ident);
	spMat = cs_add(params->eye,tmpMat,1.0,-1.0);
	//cs_print_adk(spMat);
	mt = cs_transpose(spMat,1);
	//cs_print_adk(mt);
	

	n = mt->n ;
	S = cs_sqr (0, mt, 0) ;              /* ordering and symbolic analysis */
	NN = cs_lu (mt, S, 1e-12) ;                 /* numeric LU factorization */
	xx = cs_malloc (n, sizeof (double)) ;    /* get workspace */
	
	MatGetOwnershipRange(params->denseMat1,&iStart,&iEnd);
	PetscMalloc1(N,&idx);
	for(j=0;j<N;j++) idx[j]=j;

	////Compute Entire Inverse Mat
	for(j=iStart;j<iEnd;j++){
	//create unit array for solve
		for(i=0; i<N; i++)params->b[i] = 0.0;
		params->b[j]=1.0;
		//factor outside loop leads to ~40x speedup
		cs_ipvec (NN->pinv, params->b, xx, n) ;       /* x = b(p) */
		cs_lsolve (NN->L, xx) ;               /* x = L\x */
		cs_usolve (NN->U, xx) ;               /* x = U\x */
		cs_ipvec (S->q, xx, params->b, n) ;          /* b(q) = x */
	//	if(!rank)printf("row %d\n",j);
		MatSetValues(params->denseMat1,1,&j,N,idx,params->b,INSERT_VALUES);
	}
	MatAssemblyBegin(params->denseMat1,MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(params->denseMat1,MAT_FINAL_ASSEMBLY);
	//temporary clean up
	free(xx);
	cs_spfree(spMat);
	cs_spfree(mt);
//	cs_spfree(ident);
//	cs_spfree(eye);
	cs_spfree(tmpMat);
	cs_nfree(NN);
	cs_sfree(S);
	PetscFree(idx);

//	MatView(params->denseMat1,PETSC_VIEWER_STDOUT_WORLD);

	//broadcast the invMat[0,] as vector to each processor
	if(params->rank==0){
		MatGetRow(params->denseMat1,0,&N,&idx2,&tmpArrayC);
		for (j = 0; j < N; j++) params->b[j] = tmpArrayC[j];
		MatRestoreRow(params->denseMat1,0,&Ntmp,&idx2,&tmpArrayC);	//have to use Ntmp here as MatRestoreRow does something funky
	}
	MPI_Bcast(params->b,N,MPI_DOUBLE,0,PETSC_COMM_WORLD);
	

	//Get Island Time 0-INF Unnormal
	gsl_matrix_set_zero(params->expAFS);
	fillExpectedAFS_unnorm(params->stateSpace, params->b,params->rates,params->expAFS);
	// printf("////////////////////Island Time 0 - INF unnormalized\n");
	// for(i=0;i< (params->n1+1) ;i++){
	// 	for(j=0;j< (params->n2+1);j++){
	// 		printf("%.5f ",gsl_matrix_get(params->expAFS,i,j));
	// 	}
	// 	printf("\n");
	// }
	
	//Matrix Exponentiation to get state vector at time t
	//SLEPC Stuff ahead!
	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
	Create the solver and set various options
	- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

	/* 
	Set operator matrix, the function to compute, and other options
	*/
	//transpose
	MatTranspose(params->C,MAT_REUSE_MATRIX,&params->C);
	ierr = MFNSetOperator(params->mfn,params->C);CHKERRV(ierr);
	/*
	Set solver parameters at runtime
	*/
	
	ierr = MFNSetScaleFactor(params->mfn,timeV); CHKERRV(ierr);
	ierr = MFNSetFromOptions(params->mfn);CHKERRV(ierr);

	/* set v = e_1 */
	ierr = MatGetVecs(params->C,PETSC_NULL,&params->y);CHKERRV(ierr);
	ierr = MatGetVecs(params->C,&params->v, PETSC_NULL);CHKERRV(ierr);
	ierr = VecSetValue(params->v,0,1.0,INSERT_VALUES);CHKERRV(ierr);

	ierr = VecAssemblyBegin(params->v);CHKERRV(ierr);
	ierr = VecAssemblyEnd(params->v);CHKERRV(ierr);
	ierr = MFNSolve(params->mfn,params->v,params->y);CHKERRV(ierr);
	ierr = MFNGetConvergedReason(params->mfn,&reason);CHKERRV(ierr);
	if (reason!=MFN_CONVERGED_TOL) ierr = 1; CHKERRV(ierr);

	//retranspose for later use
	MatTranspose(params->C,MAT_REUSE_MATRIX,&params->C);

	//state vector at time t stored in y
	//scatter y and store in expoArray
	VecScatterCreateToAll(params->y,&params->ctx,&params->v_seq);
	VecScatterBegin(params->ctx,params->y,params->v_seq,INSERT_VALUES,SCATTER_FORWARD);
	VecScatterEnd(params->ctx,params->y,params->v_seq,INSERT_VALUES,SCATTER_FORWARD);
	VecGetArray(params->v_seq,&tmpArray);
	for (j = 0; j < N; j++){
		params->expoArray[j] = tmpArray[j];
		//printf("%.5f\n",expoArray[j]);	
	} 

	VecRestoreArray(params->v_seq,&tmpArray);
	VecScatterDestroy(&params->ctx);
	VecDestroy(&params->v_seq);
	VecDuplicate(params->y,&params->x);
	
	//VecView(y, PETSC_VIEWER_STDOUT_WORLD);
	MatMultTranspose(params->denseMat1,params->y,params->x);
	VecScatterCreateToAll(params->x,&params->ctx,&params->v_seq);
	VecScatterBegin(params->ctx,params->x,params->v_seq,INSERT_VALUES,SCATTER_FORWARD);
	VecScatterEnd(params->ctx,params->x,params->v_seq,INSERT_VALUES,SCATTER_FORWARD);
	VecGetArray(params->v_seq,&tmpArray);
	for (j = 0; j < N; j++){
		params->b[j] = tmpArray[j];
		//printf("rank %d expoArray[%d]:%.5f\n",rank,j,expoArray[j]);	
		//printf("rank %d b[%d]:%.5f\n",rank,j,b[j]);	
	} 

	VecRestoreArray(params->v_seq,&tmpArray);
	VecScatterDestroy(&params->ctx);
	VecDestroy(&params->v_seq);


	fillExpectedAFS_unnorm(params->stateSpace, params->b,params->rates,params->expAFS2);
	//now subtract older AFS (t_t - t_INF) from total AFS (t_0 - t_INF)
	/////////////////
	///////

	gsl_matrix_sub(params->expAFS,params->expAFS2);
	
	//////////////////////
	//Final phase, collapse popns, reset rates
	//now collapse populations
	//state vector at time t stored in st[] for largerStateSpace; use reverseMap
	for(i=0;i<N;i++)params->st[i]=0.0;
	for(i=0;i<N;i++){params->st[params->map[i]]+=params->expoArray[i];
	//	printf("expoArray[%d]:%f\n",i,expoArray[i]);
	}
	for(i=0; i<params->reducedStateSpace->nstates; i++){
		VecSetValue(params->ancStateVec,i,params->st[params->reverseMap[i]],INSERT_VALUES );
	}

	//Re-Fill trans mats for Ancestral Params

	params->nnzA = coalMarkovChainTopologyMatrix_sparse(params->reducedStateSpace,params->topA,params->moveA, params->dim1A, params->dim2A);
	
	tmpMat=fillPetscCsparseTransMats(params->reducedStateSpace, params->topA, params->moveA, &params->nnzA,
		params->dim1A, params->dim2A, thetaA, 0, 0, 0, &params->C2, params->rates);
	//////////////////////
	//
	// S = (I-P)^-1
	//
	/////////
	//////
	//add negative ident
//	ident = cs_spalloc(Na,Na,Na,1,1);
//	for(i=0;i<Na;i++) cs_entry(ident,i,i,1);	
//	eye = cs_compress(ident);	
	spMat = cs_add(params->eyeAnc,tmpMat,1.0,-1.0);

	//cs_print_adk(spMat);
	mt = cs_transpose(spMat,1);
	//VecView(y,PETSC_VIEWER_STDOUT_WORLD);
	n = mt->n ;
	S = cs_sqr (0, mt, 0) ;              /* ordering and symbolic analysis */
	NN = cs_lu (mt, S, 1e-12) ;                 /* numeric LU factorization */
	xx = cs_malloc (n, sizeof (double)) ;    /* get workspace */
	MatGetOwnershipRange(params->denseMat2,&iStart,&iEnd);
//	printf("here\n");
	PetscFree(idx);
	PetscMalloc1(Na,&idx);
	for(j=0;j<Na;j++) idx[j]=j;
	////Compute Entire Inverse Mat
	for(j=iStart;j<iEnd;j++){
	//create unit array for solve
		for(i=0; i<Na; i++)params->b[i] = 0.0;
		params->b[j]=1.0;
		//factor outside loop leads to ~40x speedup
		cs_ipvec (NN->pinv, params->b, xx, n) ;       /* x = b(p) */
		cs_lsolve (NN->L, xx) ;               /* x = L\x */
		cs_usolve (NN->U, xx) ;               /* x = U\x */
		cs_ipvec (S->q, xx, params->b, n) ;          /* b(q) = x */
	//	if(!rank)printf("row %d\n",j);
		MatSetValues(params->denseMat2,1,&j,Na,idx,params->b,INSERT_VALUES);
	}
	MatAssemblyBegin(params->denseMat2,MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(params->denseMat2,MAT_FINAL_ASSEMBLY);
	free(xx);
	
	MatMultTranspose(params->denseMat2,params->ancStateVec,params->ancResVec);
	VecScatterCreateToAll(params->ancResVec,&params->ctx,&params->v_seq);
	VecScatterBegin(params->ctx,params->ancResVec,params->v_seq,INSERT_VALUES,SCATTER_FORWARD);
	VecScatterEnd(params->ctx,params->ancResVec,params->v_seq,INSERT_VALUES,SCATTER_FORWARD);


	VecGetArray(params->v_seq,&tmpArray);
	for (j = 0; j < Na; j++){
		params->expoArray[j] = tmpArray[j];
		//printf("%.5f\n",expoArray[j]);	
	}

	VecRestoreArray(params->v_seq,&tmpArray);
	VecScatterDestroy(&params->ctx);
	VecDestroy(&params->v_seq);

	//get contribution starting at st
	gsl_matrix_set_zero(params->expAFS2);
	fillExpectedAFS_unnorm(params->reducedStateSpace, params->expoArray,params->rates,params->expAFS2);


	// printf("////////////////////--Ancestral bit\n");
	// for(i=0;i< (n1+1) ;i++){
	// 	for(j=0;j< (n2+1);j++){
	// 		printf("%.5f ",gsl_matrix_get(expAFS2,i,j));
	// 	}
	// 	printf("\n");
	// }

	//if(params->rank==0)printf("////////////////////normalized IM:\n");
			
	gsl_matrix_add(params->expAFS,params->expAFS2);
	sum = matrixSumDouble(params->expAFS);
	gsl_matrix_scale(params->expAFS, 1.0 / sum);
	params->meanTreeLength = sum;
	// for(i=0;i< (params->n1+1) ;i++){
	// 		for(j=0;j< (params->n2+1);j++){
	// 			if(params->rank==0)printf("%.5f ",gsl_matrix_get(params->expAFS,i,j));
	// 		}
	// 		if(params->rank==0)printf("\n");
	// 	}
	// 	
	
	//clean up
	PetscFree(idx);
	VecDestroy(&params->y);
	VecDestroy(&params->x);
	VecDestroy(&params->v);
	PetscFree(tmpArray);
	cs_spfree(spMat);
	cs_spfree(mt);
//	cs_spfree(ident);
//	cs_spfree(eye);
	cs_spfree(tmpMat);
	gsl_vector_free(tmpStates);
	cs_nfree(NN);
	cs_sfree(S);
	
}


//calcLogAFS_IM_allPETSC -- returns the expects log AFS from an IM model; only uses PETSC
void calcLogAFS_IM_allPETSC(void * p){
	struct clam_lik_params * params = (struct clam_lik_params *) p;
	PetscInt i,j, N = params->stateSpace->nstates;
	cs *spMat, *mt, *ident, *eye, *tmpMat ;
	double timeV, sum, thetaA, theta2, mig1, mig2;
	PetscInt Na;
	gsl_vector *tmpStates;
	css *S ;
	csn *NN ;
	PetscInt n,Ntmp ;
	double *xx;
	PetscInt        iStart,iEnd, *idx;
	const PetscInt *idx2;
	PetscScalar    *tmpArray;
	PetscScalar negOne = -1.0;
//	PetscScalar one = 1.0;

	const PetscScalar *tmpArrayC;
	PetscErrorCode ierr;
	IS perm,iperm;
	MatFactorInfo  info;
	MFNConvergedReason reason;


	//initialize some vectors
	gsl_vector_set_zero(params->rates);
	MatZeroEntries(params->C);
	MatZeroEntries(params->C2);
	MatZeroEntries(params->D);
	tmpStates = gsl_vector_alloc(N);
	Ntmp=N;
	Na=params->Na;

	//For straight MLE the paramVector takes the form [N2,NA,m1,m2,t]
	theta2 = gsl_vector_get(params->paramVector,0);
	thetaA = gsl_vector_get(params->paramVector,1);
	mig1 = gsl_vector_get(params->paramVector,2);
	mig2 = gsl_vector_get(params->paramVector,3);
	timeV = gsl_vector_get(params->paramVector,4);
// 	printf("params-> %f %f %f %f %f\n",theta2,thetaA,mig1,mig2,timeV);
// 	printf("nnz %d nnzA%d \n",params->nnz,params->nnzA);

	//fill transMat
	fillPetscTransMats(params->stateSpace, params->top, params->move, &params->nnz,
		params->dim1, params->dim2, 1, theta2, mig1, mig2, &params->D, &params->C, params->rates);

	//using CSparse
	// S = (I-P)^-1
	//subtract DTMC mat from identity
	
	MatZeroEntries(params->D_copy);
	MatCopy(params->ident,params->D_copy,DIFFERENT_NONZERO_PATTERN);
	MatAXPY(params->D_copy,negOne,params->D,DIFFERENT_NONZERO_PATTERN);
//	MatTranspose(params->D_copy,MAT_REUSE_MATRIX,&params->D);

	ierr = MatGetOrdering(params->D_copy,  MATORDERINGNATURAL,  &perm,  &iperm);
	     
	ierr = MatFactorInfoInitialize(&info); 
	ierr = MatGetFactor(params->D_copy,MATSOLVERSUPERLU_DIST,MAT_FACTOR_LU,&params->F); 
	PetscInt icntl_7 = 5;
	//ierr = MatMumpsSetIcntl(params->F,7,icntl_7);
	info.fill = 5.0; 
	ierr = MatLUFactorSymbolic(params->F,params->D_copy,perm,iperm,&info);
	ierr = MatLUFactorNumeric(params->F,params->D_copy,&info);


	//ierr = MatLUFactor(params->D_copy, perm, iperm, &info); 

	////Compute Entire Inverse Mat
	ierr = MatMatSolve(params->F,params->denseIdent,params->denseMat1); 


	//MatView(params->denseMat1,PETSC_VIEWER_STDOUT_SELF);	

//	MatGetOwnershipRange(params->denseMat1,&iStart,&iEnd);
//	PetscMalloc1(N,&idx);
//	for(j=0;j<N;j++) idx[j]=j;


// 	for(j=0;j<N;j++){
// 	//create unit array for solve
// 		VecZeroEntries(params->bInv);
// 		VecZeroEntries(params->xInv);
// 		VecSetValue(params->bInv,j,one,INSERT_VALUES);
// 		VecAssemblyBegin(params->bInv);
// 		VecAssemblyEnd(params->bInv);
// 
// 		KSPSolve(params->ksp,params->bInv,params->xInv);
// 		VecGetValues(params->x,N,idx,hold);
// 		MatSetValues(params->denseMat1,1,&j,N,idx,hold,INSERT_VALUES);
// 	}
// 	MatAssemblyBegin(params->denseMat1,MAT_FINAL_ASSEMBLY);
// 	MatAssemblyEnd(params->denseMat1,MAT_FINAL_ASSEMBLY);
// 	//temporary clean up
// //	free(xx);
// //	cs_spfree(spMat);
// //	cs_spfree(mt);
// //	cs_spfree(ident);
// //	cs_spfree(eye);
// //	cs_spfree(tmpMat);
// //	cs_nfree(NN);
// //	cs_sfree(S);
// 	PetscFree(idx);
// 	KSPDestroy(&params->ksp);
//	MatView(params->denseMat1,PETSC_VIEWER_STDOUT_WORLD);


	MatDestroy(&params->F);
	//broadcast the invMat[0,] as vector to each processor
	if(params->rank==0){
		MatGetRow(params->denseMat1,0,&N,&idx2,&tmpArrayC);
		for (j = 0; j < N; j++) params->b[j] = tmpArrayC[j];
		MatRestoreRow(params->denseMat1,0,&Ntmp,&idx2,&tmpArrayC);	//have to use Ntmp here as MatRestoreRow does something funky
	}
	MPI_Bcast(params->b,N,MPI_DOUBLE,0,PETSC_COMM_WORLD);


	//Get Island Time 0-INF Unnormal
	gsl_matrix_set_zero(params->expAFS);
	fillExpectedAFS_unnorm(params->stateSpace, params->b,params->rates,params->expAFS);
	// printf("////////////////////Island Time 0 - INF unnormalized\n");
	// for(i=0;i< (params->n1+1) ;i++){
	// 	for(j=0;j< (params->n2+1);j++){
	// 		printf("%.5f ",gsl_matrix_get(params->expAFS,i,j));
	// 	}
	// 	printf("\n");
	// }

	//Matrix Exponentiation to get state vector at time t
	//SLEPC Stuff ahead!
	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
	Create the solver and set various options
	- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

	/* 
	Set operator matrix, the function to compute, and other options
	*/
	//transpose
	MatTranspose(params->C,MAT_REUSE_MATRIX,&params->C);
	ierr = MFNSetOperator(params->mfn,params->C);CHKERRV(ierr);
	/*
	Set solver parameters at runtime
	*/

	ierr = MFNSetScaleFactor(params->mfn,timeV); CHKERRV(ierr);
	ierr = MFNSetFromOptions(params->mfn);CHKERRV(ierr);

	/* set v = e_1 */
	ierr = MatGetVecs(params->C,PETSC_NULL,&params->y);CHKERRV(ierr);
	ierr = MatGetVecs(params->C,&params->v, PETSC_NULL);CHKERRV(ierr);
	ierr = VecSetValue(params->v,0,1.0,INSERT_VALUES);CHKERRV(ierr);

	ierr = VecAssemblyBegin(params->v);CHKERRV(ierr);
	ierr = VecAssemblyEnd(params->v);CHKERRV(ierr);
	ierr = MFNSolve(params->mfn,params->v,params->y);CHKERRV(ierr);
	ierr = MFNGetConvergedReason(params->mfn,&reason);CHKERRV(ierr);
	if (reason!=MFN_CONVERGED_TOL) ierr = 1; CHKERRV(ierr);

	//retranspose for later use
	MatTranspose(params->C,MAT_REUSE_MATRIX,&params->C);

	//state vector at time t stored in y
	//scatter y and store in expoArray
	VecScatterCreateToAll(params->y,&params->ctx,&params->v_seq);
	VecScatterBegin(params->ctx,params->y,params->v_seq,INSERT_VALUES,SCATTER_FORWARD);
	VecScatterEnd(params->ctx,params->y,params->v_seq,INSERT_VALUES,SCATTER_FORWARD);
	VecGetArray(params->v_seq,&tmpArray);
	for (j = 0; j < N; j++){
		params->expoArray[j] = tmpArray[j];
		//printf("%.5f\n",expoArray[j]);	
	} 

	VecRestoreArray(params->v_seq,&tmpArray);
	VecScatterDestroy(&params->ctx);
	VecDestroy(&params->v_seq);
	VecDuplicate(params->y,&params->x);

	//VecView(y, PETSC_VIEWER_STDOUT_WORLD);
	MatMultTranspose(params->denseMat1,params->y,params->x);
	VecScatterCreateToAll(params->x,&params->ctx,&params->v_seq);
	VecScatterBegin(params->ctx,params->x,params->v_seq,INSERT_VALUES,SCATTER_FORWARD);
	VecScatterEnd(params->ctx,params->x,params->v_seq,INSERT_VALUES,SCATTER_FORWARD);
	VecGetArray(params->v_seq,&tmpArray);
	for (j = 0; j < N; j++){
		params->b[j] = tmpArray[j];
		//printf("rank %d expoArray[%d]:%.5f\n",rank,j,expoArray[j]);	
		//printf("rank %d b[%d]:%.5f\n",rank,j,b[j]);	
	} 

	VecRestoreArray(params->v_seq,&tmpArray);
	VecScatterDestroy(&params->ctx);
	VecDestroy(&params->v_seq);


	fillExpectedAFS_unnorm(params->stateSpace, params->b,params->rates,params->expAFS2);
	//now subtract older AFS (t_t - t_INF) from total AFS (t_0 - t_INF)
	/////////////////
	///////

	gsl_matrix_sub(params->expAFS,params->expAFS2);

	//////////////////////
	//Final phase, collapse popns, reset rates
	//now collapse populations
	//state vector at time t stored in st[] for largerStateSpace; use reverseMap
	for(i=0;i<N;i++)params->st[i]=0.0;
	for(i=0;i<N;i++){params->st[params->map[i]]+=params->expoArray[i];
	//	printf("expoArray[%d]:%f\n",i,expoArray[i]);
	}
	for(i=0; i<params->reducedStateSpace->nstates; i++){
		VecSetValue(params->ancStateVec,i,params->st[params->reverseMap[i]],INSERT_VALUES );
	}

	//Re-Fill trans mats for Ancestral Params

	params->nnzA = coalMarkovChainTopologyMatrix_sparse(params->reducedStateSpace,params->topA,params->moveA, params->dim1A, params->dim2A);

	tmpMat=fillPetscCsparseTransMats(params->reducedStateSpace, params->topA, params->moveA, &params->nnzA,
		params->dim1A, params->dim2A, thetaA, 0, 0, 0, &params->C2, params->rates);
	//////////////////////
	//
	// S = (I-P)^-1
	//
	/////////
	//////
	//add negative ident
	ident = cs_spalloc(Na,Na,Na,1,1);
	for(i=0;i<Na;i++) cs_entry(ident,i,i,1);	
	eye = cs_compress(ident);	
	spMat = cs_add(eye,tmpMat,1.0,-1.0);

	//cs_print_adk(spMat);
	mt = cs_transpose(spMat,1);
	//VecView(y,PETSC_VIEWER_STDOUT_WORLD);
	n = mt->n ;
	S = cs_sqr (0, mt, 0) ;              /* ordering and symbolic analysis */
	NN = cs_lu (mt, S, 1e-12) ;                 /* numeric LU factorization */
	xx = cs_malloc (n, sizeof (double)) ;    /* get workspace */
	MatGetOwnershipRange(params->denseMat2,&iStart,&iEnd);
//	printf("here\n");
//	PetscFree(idx);
	PetscMalloc1(Na,&idx);
	for(j=0;j<Na;j++) idx[j]=j;
	////Compute Entire Inverse Mat
	for(j=iStart;j<iEnd;j++){
	//create unit array for solve
		for(i=0; i<Na; i++)params->b[i] = 0.0;
		params->b[j]=1.0;
		//factor outside loop leads to ~40x speedup
		cs_ipvec (NN->pinv, params->b, xx, n) ;       /* x = b(p) */
		cs_lsolve (NN->L, xx) ;               /* x = L\x */
		cs_usolve (NN->U, xx) ;               /* x = U\x */
		cs_ipvec (S->q, xx, params->b, n) ;          /* b(q) = x */
	//	if(!rank)printf("row %d\n",j);
		MatSetValues(params->denseMat2,1,&j,Na,idx,params->b,INSERT_VALUES);
	}
	MatAssemblyBegin(params->denseMat2,MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(params->denseMat2,MAT_FINAL_ASSEMBLY);
	free(xx);

	MatMultTranspose(params->denseMat2,params->ancStateVec,params->ancResVec);
	VecScatterCreateToAll(params->ancResVec,&params->ctx,&params->v_seq);
	VecScatterBegin(params->ctx,params->ancResVec,params->v_seq,INSERT_VALUES,SCATTER_FORWARD);
	VecScatterEnd(params->ctx,params->ancResVec,params->v_seq,INSERT_VALUES,SCATTER_FORWARD);


	VecGetArray(params->v_seq,&tmpArray);
	for (j = 0; j < Na; j++){
		params->expoArray[j] = tmpArray[j];
		//printf("%.5f\n",expoArray[j]);	
	}

	VecRestoreArray(params->v_seq,&tmpArray);
	VecScatterDestroy(&params->ctx);
	VecDestroy(&params->v_seq);

	//get contribution starting at st
	gsl_matrix_set_zero(params->expAFS2);
	fillExpectedAFS_unnorm(params->reducedStateSpace, params->expoArray,params->rates,params->expAFS2);


	// printf("////////////////////--Ancestral bit\n");
	// for(i=0;i< (n1+1) ;i++){
	// 	for(j=0;j< (n2+1);j++){
	// 		printf("%.5f ",gsl_matrix_get(expAFS2,i,j));
	// 	}
	// 	printf("\n");
	// }

	//if(params->rank==0)printf("////////////////////normalized IM:\n");

	gsl_matrix_add(params->expAFS,params->expAFS2);
	sum = matrixSumDouble(params->expAFS);
	gsl_matrix_scale(params->expAFS, 1.0 / sum);
	// for(i=0;i< (params->n1+1) ;i++){
	// 		for(j=0;j< (params->n2+1);j++){
	// 			if(params->rank==0)printf("%.5f ",gsl_matrix_get(params->expAFS,i,j));
	// 		}
	// 		if(params->rank==0)printf("\n");
	// 	}
	// 	

	//clean up
	PetscFree(idx);
	VecDestroy(&params->y);
	VecDestroy(&params->x);
	VecDestroy(&params->v);
	ierr = ISDestroy(&perm);
  	ierr = ISDestroy(&iperm);
	PetscFree(tmpArray);
	cs_spfree(spMat);
	cs_spfree(mt);
	cs_spfree(ident);
	cs_spfree(eye);
	cs_spfree(tmpMat);
	gsl_vector_free(tmpStates);
	cs_nfree(NN);
	cs_sfree(S);

}
double calcLikNLOpt(unsigned n, const double *point, double *gradients, void *p){
	int i,j;
	struct clam_lik_params * params = (struct clam_lik_params *) p;
	double localNNZ = params->nnz;
	double output = 666.0; 
	double lik=0.0;
	double x[5];
	
	for(i = 0;i<5;i++) x[i] = point[i];
	MPI_Bcast(x,5,MPI_DOUBLE,0,PETSC_COMM_WORLD);
	
	for(i = 0;i<5;i++) gsl_vector_set(params->paramVector, i, x[i]);
	
	//fill in the expAFS table
	calcLogAFS_IM(p);
	params->nnz = localNNZ;
	//compute lik
	for(i=0;i<params->obsData->size1;i++){
		for(j=0;j<params->obsData->size2;j++){
			if (gsl_matrix_get(params->expAFS,i,j) > 0.0) //corners of AFS are zero prob
				lik += gsl_matrix_get(params->obsData,i,j) * (log(gsl_matrix_get(params->expAFS,i,j)));
		}
	}
	output = -1.0* lik;
	MPI_Bcast(&output,1,MPI_DOUBLE,0,PETSC_COMM_WORLD);
	if(params->rank == 0 && vbse){
		for(i = 0;i<5;i++) printf("x[%d]: %lf\t",i,x[i]);
		printf("lik: %lf\n",lik);
	}
	params->fEvals += 1;
//	gsl_matrix_prettyPrint(params->expAFS);
	return(output);
}

double calcLikNLOpt_gradients(unsigned n, const double *point, double *gradients, void *p){
	int i,j;
	//struct clam_lik_params * params = (struct clam_lik_params *) p;
	//double localNNZ = params->nnz;
	double output = 666.0; 
	double lik=0.0;
	double newF;
	double x[5];
	double temp[5];
	double epsilon=1e-5;
	
	
	for(i = 0;i<5;i++) x[i] = point[i];
	MPI_Bcast(x,5,MPI_DOUBLE,0,PETSC_COMM_WORLD);
	
	lik = calcLikNLOpt(n, point, NULL, p);
	if (gradients != NULL){
		for(i=0; i<n ;i++){
			for( j=0; j< n; j++)
				temp[j] = x[j];
		
			temp[i] +=  epsilon;
			newF = calcLikNLOpt(n,temp,NULL,p);
			gradients[i] = (newF -  lik)/epsilon;
		}
	}
	
	output = lik;
	MPI_Bcast(&output,1,MPI_DOUBLE,0,PETSC_COMM_WORLD);
//	gsl_matrix_prettyPrint(params->expAFS);
	return(output);
}

//maximizeLik-- maximizes the logLikFunction function using constraints via
// nlOpt library

void maximizeLikNLOpt(double *lik, void *p, double *mle){
	//struct clam_lik_params * params = (struct clam_lik_params *) p;
	int np = 5;
	int err;
	//double x[5] = {1.0,1.0,1.0,1.0,1.5};
	double minimum=666.6;
	nlopt_opt opt;
	
	opt = nlopt_create(NLOPT_LD_LBFGS, np);
	//define constraints


	nlopt_set_lower_bounds(opt, lowerBounds);
	nlopt_set_upper_bounds(opt, upperBounds);
	nlopt_set_min_objective(opt, calcLikNLOpt_gradients, (void *) p);
	
	nlopt_set_xtol_rel(opt, 1e-6);
	//set the seed for deterministic sequence in COBYLA alg
	//nlopt_srand(1232333);
	err = nlopt_optimize(opt, mle, &minimum);
	if (err < 0) {
	    fprintf(stderr,"nlopt failed with code %d!\n",err);
	}
//	MPI_Barrier(MPI_COMM_WORLD);
	*lik=minimum;
	nlopt_destroy(opt);
	
}

//maximizeLikNLOpt_twoStage-- begins with coarse global, then goes local

void maximizeLikNLOpt_twoStage(double *lik, void *p, double *mle){
	//struct clam_lik_params * params = (struct clam_lik_params *) p;
	int np = 5;
	int err;
	//double x[5] = {1.0,1.0,1.0,1.0,1.5};
	double minimum=666.6;
	nlopt_opt opt;
	
	opt = nlopt_create(NLOPT_GN_DIRECT_L_RAND, np);


	nlopt_set_lower_bounds(opt, lowerBounds);
	nlopt_set_upper_bounds(opt, upperBounds);
	nlopt_set_min_objective(opt, calcLikNLOpt_gradients, (void *) p);
	


	
	//nlopt_set_xtol_rel(opt, 0.25);
	//set the seed for deterministic sequence in COBYLA alg
	//nlopt_srand(1232333);
	err = nlopt_optimize(opt, mle, &minimum);
	if (err < 0) {
	    fprintf(stderr,"nlopt failed with code %d!\n",err);
	}
//	MPI_Barrier(MPI_COMM_WORLD);
	*lik=minimum;
	nlopt_destroy(opt);
	
	
}

//maximizeLikNLOpt_MLSL-- multi-level restarts

void maximizeLikNLOpt_MLSL(double *lik, void *p, double *mle){
	int np = 5;
	int err;
	//double x[5] = {1.0,1.0,1.0,1.0,1.5};
	double minimum=666.6;
	nlopt_opt opt, l_opt;
	
	opt = nlopt_create(NLOPT_G_MLSL, np);


	nlopt_set_lower_bounds(opt, lowerBounds);
	nlopt_set_upper_bounds(opt, upperBounds);
	nlopt_set_min_objective(opt, calcLikNLOpt_gradients, (void *) p);
	
	l_opt = nlopt_create(NLOPT_LD_LBFGS, np);
	nlopt_set_lower_bounds(l_opt, lowerBounds);
	nlopt_set_upper_bounds(l_opt, upperBounds);
	nlopt_set_min_objective(l_opt, calcLikNLOpt_gradients, (void *) p);
	//define constraints


	nlopt_set_local_optimizer(opt, l_opt);
	
	//nlopt_set_xtol_rel(opt, 0.25);
	//set the seed for deterministic sequence in COBYLA alg
	//nlopt_srand(1232333);
	err = nlopt_optimize(opt, mle, &minimum);
	if (err < 0) {
	    fprintf(stderr,"nlopt failed with code %d!\n",err);
	}
//	MPI_Barrier(MPI_COMM_WORLD);
	*lik=minimum;
	nlopt_destroy(l_opt);
	nlopt_destroy(opt);
	
	
}


///////////////// Hessian stuff ///////////////////

double hessianMatrix_element(double lik, double *mle, int i, int j, double hi, double hj, void *p){
	//struct clam_lik_params * params = (struct clam_lik_params *) p;
	double Pi0, Pj0,output,fpp,fpm,fmp,fmm;
	
	Pi0 = mle[i]; Pj0 = mle[j];
	
	//printf("mle: %f\n",lik);	
	if(i==j){
		mle[i] = Pi0 + hi;
		fpp = calcLikNLOpt(5,mle,NULL,p);
		mle[i] = Pj0 - hi;
		fmm = calcLikNLOpt(5,mle,NULL,p);
		output = ((fpp - lik) + (fmm - lik))/(hi*hi);
	}
	else{
		// f(xi + hi, xj + h)
		mle[i] = Pi0 + hi;
		mle[j] = Pj0 + hj;
		fpp = calcLikNLOpt(5,mle,NULL,p);

	        // f(xi + hi, xj - hj)
		mle[i] = Pi0 + hi;
		mle[j] = Pj0 - hj;
		fpm = calcLikNLOpt(5,mle,NULL,p);

	        // f(xi - hi, xj + hj)
		mle[i] = Pi0 - hi;
		mle[j] = Pj0 + hj;
		fmp = calcLikNLOpt(5,mle,NULL,p);

	        // f(xi - hi, xj - hj)
		mle[i] = Pi0 - hi;
		mle[j] = Pj0 - hj;
		fmm = calcLikNLOpt(5,mle,NULL,p);

		output = (fpp - fpm - fmp + fmm)/(4 * hi * hj);
	}
	mle[i]=Pi0; mle[j]=Pj0;
	return output;
}

gsl_matrix *hessian(double *mle, double lik, void *p){
	int i, j;
	double eps;
	gsl_matrix *H;
	
	
	H = gsl_matrix_alloc(5,5);
	eps = 0.01;
	for(i=0;i<5;i++){
		for (j=i;j<5;j++){
			gsl_matrix_set(H,i,j,hessianMatrix_element(lik, mle, i, j, eps, eps, p));
			gsl_matrix_set(H,j,i,gsl_matrix_get(H,i,j));
		}
	}
	return(H);
}

gsl_matrix *getFisherInfoMatrix(double *mle, double lik, void *p){
	gsl_matrix *H;
	gsl_matrix *fi;
	gsl_permutation *perm = gsl_permutation_alloc(5);
	int s;
	
	H = hessian(mle,lik,p);
	//gsl_matrix_prettyPrint(H);printf("\n");
	fi = gsl_matrix_alloc(5,5);
	gsl_matrix_scale(H,-1.0);
	gsl_linalg_LU_decomp (H, perm, &s);    
	gsl_linalg_LU_invert (H, perm, fi);
	//gsl_matrix_prettyPrint(fi);printf("\n");
	return(fi);
}

gsl_matrix *getGodambeInfoMatrix(double *mle, double lik, void *p){
	gsl_matrix *H;
	gsl_matrix *J,*Jtemp, *Jinv, *origData, *boot;
	gsl_vector *cU, *grad_temp;
	struct clam_lik_params * params = (struct clam_lik_params *) p;
	gsl_permutation *perm = gsl_permutation_alloc(5);
	int s,i, nBoot;
	
	origData = params->obsData;
	nBoot = 100;
	
	H = hessian(mle,lik,p);
	gsl_matrix_scale(H,-1.0);
	//gsl_matrix_prettyPrint(H);printf("\n");
	J = gsl_matrix_alloc(5,5);
	Jinv = gsl_matrix_alloc(5,5);
	Jtemp = gsl_matrix_alloc(5,5);
	boot = gsl_matrix_alloc(origData->size1,origData->size2);
	
	gsl_matrix_set_zero(J);
	gsl_matrix_set_zero(Jinv);
	gsl_matrix_set_zero(Jtemp);
	cU = gsl_vector_alloc(5);
	gsl_vector_set_zero(cU);
	
	//do bootstraps, estimate J
	for(i=0;i<nBoot;i++){
		gsl_matrix_bootstrap(origData, boot, params->rng);
		params->obsData = boot;
		grad_temp = getGradient(mle,p);
		gsl_vector_outer_product(grad_temp,grad_temp,Jtemp);
		gsl_matrix_add(J,Jtemp);
		gsl_vector_add(cU,grad_temp);
		//gsl_matrix_prettyPrint(Jtemp);
		gsl_vector_free(grad_temp);
	}
	gsl_vector_scale(cU,1.0/(float)nBoot);
	gsl_matrix_scale(J,1.0/(float)nBoot);
	//printf("cU:\n");
	//gsl_vector_fprintf(stdout,cU,"%f");
	//printf("J:\n");
	//gsl_matrix_prettyPrint(J);printf("\n");
	//G = H*J^-1*H
	gsl_linalg_LU_decomp (J, perm, &s);    
	gsl_linalg_LU_invert (J, perm, Jinv);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, H, Jinv, 0.0, Jtemp);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, Jtemp, H, 0.0, J);
	//printf("godambe:\n");
	//gsl_matrix_prettyPrint(J);printf("\n");
	//need to invert this to get the Godambe information matrix
	gsl_linalg_LU_decomp (J, perm, &s);    
	gsl_linalg_LU_invert (J, perm, Jinv);
	return(Jinv);
}

gsl_vector *getGradient(double *mle, void *p){
	int i, j;
	double eps, epsVals[5], paramTemp[5], fp, fm;
	gsl_vector *gradient;
	
	gradient = gsl_vector_alloc(5);
	eps = 0.01;
	for(i=0;i<5;i++){
		if( mle[i] != 0){
			epsVals[i] = mle[i] * eps;
			//printf("eps[%d]=%f\n",i,epsVals[i]);
		}
		else{
			epsVals[i] = eps;
		}
		paramTemp[i] = mle[i];
	}
	for(i=0;i<5;i++){
		for(j=0;j<5;j++) paramTemp[j] = mle[j];
		if( mle[i] != 0){
			paramTemp[i] = mle[i] + epsVals[i];
			fp = calcLikNLOpt(5,paramTemp,NULL,p);
			paramTemp[i] = mle[i] - epsVals[i];
			fm = calcLikNLOpt(5,paramTemp,NULL,p);
			gsl_vector_set(gradient,i,(fp-fm)/(2*epsVals[i]));
		}
		else{
			paramTemp[i] = mle[i] + epsVals[i];
			fp = calcLikNLOpt(5,paramTemp,NULL,p);
			paramTemp[i] = mle[i];
			fm = calcLikNLOpt(5,paramTemp,NULL,p);
			gsl_vector_set(gradient,i,(fp-fm)/epsVals[i]);
		}
	}
	return(gradient);
}

