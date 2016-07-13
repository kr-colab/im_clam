#include "AFS.h"
#include "cs.h"
#include "AFS_ctmc.h"
#include "sparseExp.h"
#include "adkCSparse.h"
#include "adkGSL.h"
//#include "cuba.h"
#include <gsl/gsl_blas.h>


void fillExpectedAFS_unnorm(afsStateSpace *S, double *visitMat ,gsl_vector *rates, gsl_matrix *expAFS);



//fills a sparse transition matrix for Continuous Time model
struct cs_di_sparse *fillTransitionMatrix_ctmc(afsStateSpace *S, double *topol, int *moveType, int *nzCount,
	int *dim1, int *dim2, double theta1, double theta2, double mig1, double mig2){
	int i, j,k;
	int N = S->nstates;
	double  c0r[N],c1r[N],m0r[N],m1r[N],totR[N],rowSums[N], tmp;;
	
	
	struct cs_di_sparse *triplet;

	

	//now allocate new cs_sparse obj
	triplet = cs_spalloc(N, N, *nzCount + N , 1, 1); //alloc sparse mat with extra space for nonzero identity mats
	for(i=0;i<N;i++){
		c0r[i] = S->states[i]->aCounts[0] * (S->states[i]->aCounts[0]-1) / theta1 ;
		c1r[i] = S->states[i]->aCounts[1] * (S->states[i]->aCounts[1]-1) / theta2;
		m0r[i] = S->states[i]->aCounts[0] * mig1 ;
		m1r[i] = S->states[i]->aCounts[1] * mig2 ;
		totR[i] = c0r[i] + c1r[i] + m0r[i] + m1r[i];
		rowSums[i] = 0.0;
	}
	
	for(k=0;k<*nzCount;k++){
		i = dim1[k];
		j= dim2[k];
		switch(moveType[k]){
			case 0: //coal pop0
			tmp = c0r[i] / totR[i] * topol[k];
			cs_entry(triplet,i,j,tmp);
			rowSums[i] += tmp;
			break;
			case 1: //coal pop1
			tmp = c1r[i] / totR[i] * topol[k];
			cs_entry(triplet,i,j,tmp);
			rowSums[i] += tmp;
			break;
			case 2: //mig pop0
			tmp = m0r[i] / totR[i] * topol[k];
			cs_entry(triplet,i,j,tmp);
			rowSums[i] += tmp;
			break;
			case 3: //mig pop1
			tmp = m1r[i] / totR[i] * topol[k];
			cs_entry(triplet,i,j,tmp);
			rowSums[i] += tmp;
			break;
		}
	}
	//now add diagonal elements; topol pre-set to -1
	for(k=0;k<N;k++){
		if((rowSums[k]) != 0.0){
	 		cs_entry(triplet,k,k,0-rowSums[k]);
			*nzCount += 1;
		}
	}
//	cs_dropzeros(triplet);
//	tmat = cs_compress(triplet);
//	cs_spfree(triplet);
	return(triplet);

}

//this only fills a new set of arrays
void fillTransitionMatrixArray_ctmc(afsStateSpace *S, double *topol, int *moveType, int *nzCount,
	int *dim1, int *dim2, double theta1, double theta2, double mig1, double mig2,double *newArray,
	int *newDim1, int *newDim2){
	int i, j,k, newnz;
	int N = S->nstates;
	double  c0r[N],c1r[N],m0r[N],m1r[N],totR[N],rowSums[N], tmp;;

	for(i=0;i<N;i++){
		c0r[i] = S->states[i]->aCounts[0] * (S->states[i]->aCounts[0]-1) / theta1 ;
		c1r[i] = S->states[i]->aCounts[1] * (S->states[i]->aCounts[1]-1) / theta2;
		m0r[i] = S->states[i]->aCounts[0] * mig1 ;
		m1r[i] = S->states[i]->aCounts[1] * mig2 ;
		totR[i] = c0r[i] + c1r[i] + m0r[i] + m1r[i];
		rowSums[i] = 0.0;
	}	
	newnz=0;
	for(k=0;k<*nzCount;++k){
		i = dim1[k];
		j= dim2[k];
		switch(moveType[k]){
			case 0: //coal pop0
			tmp = c0r[i] / totR[i] * topol[k];
			newArray[k]=tmp;
			newDim1[k] = i+1;
			newDim2[k] = j+1;
			rowSums[i] += tmp;
			newnz++;
			break;
			case 1: //coal pop1
			tmp = c1r[i] / totR[i] * topol[k];
			newArray[k]=tmp;
			newDim1[k] = i+1;
			newDim2[k] = j+1;
			rowSums[i] += tmp;
			newnz++;
			break;
			case 2: //mig pop0
			tmp = m0r[i] / totR[i] * topol[k];
			newArray[k]=tmp;
			newDim1[k] = i+1;
			newDim2[k] = j+1;
			rowSums[i] += tmp;
			newnz++;
			break;
			case 3: //mig pop1
			tmp = m1r[i] / totR[i] * topol[k];
			newArray[k]=tmp;
			newDim1[k] = i+1;
			newDim2[k] = j+1;
			rowSums[i] += tmp;
			newnz++;
			break;
		}
	}

	//now add diagonal elements; topol pre-set to -1
	for(k=0;k<N;k++){
		if((rowSums[k]) != 0.0){
			newArray[newnz]=0-rowSums[k];
			newDim1[newnz] = k+1;
			newDim2[newnz] = k+1;
			newnz++;
		}
	}
	
	//sanity check
//	for(k=0;k<newnz;k++){
//		printf("row: %d col: %d val: %g\n",newDim1[k],newDim2[k],newArray[k]);
//	}
	*nzCount = newnz;
}

//Functions for jointly computing CTMC transitions and embedded chain

//this function returns pointer to embedded Discrete MC trans Mat in sparse form, while filling rates vector and filling
// arrays for CTMC representation
struct cs_di_sparse * fillTransitionMatrixArray_embed_ctmc(afsStateSpace *S, double *topol, int *moveType, int *nzCount,
	int *dim1, int *dim2, double theta1, double theta2, double mig1, double mig2,double *newArray,
	int *newDim1, int *newDim2, gsl_vector *rates){
	int i, j,k, newnz, abcount=0;
	int N = S->nstates,abFlag;
	double  c0r[N],c1r[N],m0r[N],m1r[N],totR[N],rowSums[N], tmp;;
	struct cs_di_sparse *triplet, *tmat;

	gsl_vector_set_zero(rates);
	//allocate new cs_sparse obj
	triplet = cs_spalloc(N, N, *nzCount + N , 1, 1); //alloc sparse mat with extra space for nonzero identity mats

	for(i=0;i<N;i++){
		c0r[i] = S->states[i]->aCounts[0] * (S->states[i]->aCounts[0]-1)  / theta1 ;
	//	printf("co check:%f\n",c0r[i]);
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
			newArray[k]=tmp* totR[i];
			newDim1[k] = i;
			newDim2[k] = j;
			rowSums[i] += tmp* totR[i];
			if(abFlag==0)cs_entry(triplet,i,j,tmp);			
			break;
			case 1: //coal pop1
			if(totR[i]==0)tmp=0;
			else tmp = c1r[i] / totR[i] * topol[k];
			newArray[k]=tmp* totR[i];
			newDim1[k] = i;
			newDim2[k] = j;
			rowSums[i] += tmp* totR[i];
			if(abFlag==0)cs_entry(triplet,i,j,tmp);			
			break;
			case 2: //mig pop0
			if(totR[i]==0)tmp=0;
			else tmp = m0r[i] / totR[i] * topol[k];
			newArray[k]=tmp* totR[i];
			newDim1[k] = i;
			newDim2[k] = j;
			rowSums[i] += tmp* totR[i];
			if(abFlag==0)cs_entry(triplet,i,j,tmp);			
			break;
			case 3: //mig pop1
			if(totR[i]==0)tmp=0;
			else tmp = m1r[i] / totR[i] * topol[k];			
			newArray[k]=tmp * totR[i];
			newDim1[k] = i;
			newDim2[k] = j;
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
			newArray[newnz]=0.0-rowSums[k];
			newDim1[newnz] = k;
			newDim2[newnz] = k;
			newnz++;
		}
	}

	tmat = cs_compress(triplet);
	cs_spfree(triplet);
	*nzCount = newnz;
	return(tmat);	
}










//fills up an AFS matrix using expected times and a precalculated visitMat
void fillExpectedAFS_ctmc(afsStateSpace *S, double *visitMat, gsl_matrix *expAFS){
	int  j, k, l, m1, m2, count;

	double tmpRes;

	m1=expAFS->size1;
	m2=expAFS->size2;
	
	gsl_matrix_set_zero(expAFS);
	
	for(j=0;j<S->nstates;j++){
		if(S->states[j]->nalleles != 1){
			for(k=0;k<m1;k++){
			for(l=0;l<m2;l++){
				count = gsl_matrix_int_get(S->states[j]->popMats[0],k,l) +  gsl_matrix_int_get(S->states[j]->popMats[1],k,l);
				tmpRes = gsl_matrix_get(expAFS,k,l) + (count * visitMat[j]);	
				gsl_matrix_set(expAFS,k,l, tmpRes);
				}
			}
		}
	}
//	tmpRes = matrixSumDouble(expAFS);
//	gsl_matrix_scale(expAFS, 1.0 / tmpRes);
}


//fills up an AFS matrix using expected times and a precalculated visitMat
void fillExpectedAFS_unnorm(afsStateSpace *S, double *visitMat ,gsl_vector *rates, gsl_matrix *expAFS){
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
}

//fills up an AFS matrix using expected times and a precalculated visitMat
void fillExpectedAFS_unnorm_mat(afsStateSpace *S, double **visitMat , double *initialStateProbs, gsl_vector *rates, gsl_matrix *expAFS){
	int  j, k, l,i, m1, m2, count;

	double tmpRes,sum;

	m1=expAFS->size1;
	m2=expAFS->size2;
	
	gsl_matrix_set_zero(expAFS);	
	for(j=0;j<S->nstates;j++){
		for(k=0;k<m1;k++){
			for(l=0;l<m2;l++){
				count = gsl_matrix_int_get(S->states[j]->popMats[0],k,l) +  gsl_matrix_int_get(S->states[j]->popMats[1],k,l);
				sum=0.0;
				for(i=0;i<S->nstates;i++){
					sum += visitMat[i][j] * initialStateProbs[i];
				}
				tmpRes	= gsl_matrix_get(expAFS,k,l) + (count * gsl_vector_get(rates,j) * sum);
			//	printf("tmpRes[%d][%d]: %f rates[j]: %f visits[j]: %f count: %d\n",k,l,tmpRes,gsl_vector_get(rates,j),sum,count);
				gsl_matrix_set(expAFS,k,l, tmpRes);
			}
		}
	}
}

//fillLogAFS -- uses cxsparse library for inversion -- will replace above
void fillLogAFS(void * p){
	struct im_lik_params * params = (struct im_lik_params *) p;
	struct expoMatObj anExpoMatObj;
	int i,j, N = params->stateSpace->nstates;
 	cs *spMat, *mt, *ident, *eye, *tmpMat ,*tmpMat2;
	double timeV, sum, thetaA, theta2, mig1, mig2;
	int Na;
	gsl_vector *tmpStates;
	
	
	//initialize some vectors
	gsl_vector_set_zero(params->rates);
	tmpStates = gsl_vector_alloc(N);
	
	//For straight MLE the paramVector takes the form [N2,NA,m1,m2,t]
	theta2 = gsl_vector_get(params->paramVector,0);
	thetaA = gsl_vector_get(params->paramVector,1);
	mig1 = gsl_vector_get(params->paramVector,2);
	mig2 = gsl_vector_get(params->paramVector,3);
	timeV = gsl_vector_get(params->paramVector,4);
//	printf("params-> %f %f %f %f %f\n",theta2,thetaA,mig1,mig2,timeV);
	//fill transMat
	tmpMat= fillTransitionMatrixArray_embed_ctmc(params->stateSpace, params->topA, params->moveA,
	 	&(params->nnz),params->dim1, params->dim2,1.0, theta2, mig1, mig2, 
	 	params->new, params->newdim1, params->newdim2,params->rates);
	
	//using CSparse
	// S = (I-P)^-1
	
	//add negative ident
	ident = cs_spalloc(N,N,N,1,1);
	for(i=0;i<N;i++) cs_entry(ident,i,i,1);
	eye = cs_compress(ident);
	spMat = cs_add(eye,tmpMat,1.0,-1.0);
	mt = cs_transpose(spMat,1);
	
	////Compute Entire Inverse Mat
	for(j=0;j<N;j++){
		//create unit array for solve
		for(i=0; i<N; i++)params->b[i] = 0.0;
		params->b[j]=1.0;
		cs_lusol(0, mt, params->b, 1e-12);
		for(i=0; i<N; i++){
			params->invMat[j][i]=params->b[i];
			gsl_matrix_set(params->invMatGSL,j,i,params->b[i]);
		}
		printf("inv mat row %d complete\n",j);
	}
	
	//Get Island Time 0-INF Unnormal
	gsl_matrix_set_zero(params->expAFS);
	for(i=0; i<N; i++)params->b[i] = params->invMat[0][i];
	fillExpectedAFS_unnorm(params->stateSpace, params->b,params->rates,params->expAFS);
	// printf("////////////////////Island Time 0 - INF unnormalized\n");
	// for(i=0;i< (params->expAFS->size1) ;i++){
	// 	for(j=0;j< (params->expAFS->size2);j++){
	// 		printf("%.5f ",gsl_matrix_get(params->expAFS,i,j));
	// 	}
	// 	printf("\n");
	// }
	
	//Matrix Exponentiation to get state vector at time t
	//

	//push indices
	for(i=0;i<params->nnz;i++){
		params->expodim1[i]=params->newdim1[i]+1;
		params->expodim2[i]=params->newdim2[i]+1;
		
	}
	//reset expoInit
	for(i = 0;i<N;i++){
		params->expoInit[i]=0.0;
		params->expoResult[i]=0.0;
	}
	params->expoInit[0] = 1.0;

	//sanity check ///////
//	for(i=0;i<params->nnz;i++){
//		printf("a[%d]=%f expodim1=%d expodim2=%d \n", i,params->anExpoMatObj->a[i],params->anExpoMatObj->ai[i],params->anExpoMatObj->aj[i]);
//	}
//	printf("rank=%d nnz=%d\n", params->anExpoMatObj->rank, params->anExpoMatObj->nnz);
	
//	for(i=0;i<N;i++)printf("params->anExpoMatObj->resultSpace[%d]=%f\n",i,params->anExpoMatObj->resultSpace[i]);
//	for(i=0;i<N;i++)printf("params->anExpoMatObj->initVec[%d]=%f\n",i,params->anExpoMatObj->initVec[i]);
	
	anExpoMatObj.ai = params->expodim1;
	anExpoMatObj.aj = params->expodim2;
	anExpoMatObj.a = params->new;
	anExpoMatObj.initVec = params->expoInit;
	anExpoMatObj.resultSpace = params->expoResult;
//	sparse_exponential_single_row(N, params->nnz, anExpoMatObj.ai, anExpoMatObj.aj,
//		anExpoMatObj.a, timeV,anExpoMatObj.initVec, params->expoResult, 0); 
	//state vector at time t stored in b[i] 
	for(i=0; i<N; i++){
		params->b[i] = params->expoResult[i];
//		printf("params->expoResult[%d]=%f\n",i,params->expoResult[i] );		
		gsl_vector_set(tmpStates,i,params->expoResult[i]);
	}

	gsl_blas_dgemv(CblasTrans, 1.0, params->invMatGSL,tmpStates, 0.0, params->resVec);
	

	for(i=0; i<N; i++){
		params->b[i] = gsl_vector_get(params->resVec,i);
	}
	gsl_matrix_set_zero(params->expAFS2);
	fillExpectedAFS_unnorm(params->stateSpace, params->b,params->rates,params->expAFS2);
	
	//now subtract older AFS (t_t - t_INF) from total AFS (t_0 - t_INF)
	/////////////////
	///////
//	printf("////////////////////Island time0-t=%f unnormalized\n",timeV);
	
	gsl_matrix_sub(params->expAFS,params->expAFS2);
	// for(i=0;i< (params->expAFS->size1) ;i++){
	// 	for(j=0;j< (params->expAFS->size2);j++){
	// 		printf("%.5f ",gsl_matrix_get(params->expAFS,i,j));
	// 	}
	// 	printf("\n");
	// }
	// printf("////////////////////\n");
	
	
	//mapping of states already complete when params initialized
	//state vector at time t stored in st[] for largerStateSpace; use reverseMap
	for(i=0;i<N;i++)params->st[i]=0.0;
	for(i=0;i<N;i++)params->st[params->map[i]]+=params->expoResult[i];
	for(i=0; i<params->reducedStateSpace->nstates; i++){
		gsl_vector_set(params->ancStateVec,i,params->st[params->reverseMap[i]]);
	}
	Na = params->reducedStateSpace->nstates;
	//fill up ancestral transition matrix
	gsl_vector_set_zero(params->rates);
	tmpMat2= fillTransitionMatrixArray_embed_ctmc(params->reducedStateSpace, params->topA2, params->moveA2, &(params->nnzA),
		params->dim1A, params->dim2A,thetaA, 0, 0, 0,params->new, params->newdim1, params->newdim2,params->rates);
		
	ident = cs_spalloc(Na,Na,Na,1,1);
	for(i=0;i<Na;i++) cs_entry(ident,i,i,1);
	eye = cs_compress(ident);
	spMat = cs_add(eye,tmpMat2,1.0,-1.0);
	mt = cs_transpose(spMat,1);

//	gsl_matrix_free(invMatGSL);
//	invMatGSLA = gsl_matrix_calloc(Na,Na);

	for(j=0;j<Na;j++){
		//create unit array for solve
		for(i=0; i<Na; i++)params->b[i] = 0.0;
		params->b[j]=1.0;
		cs_lusol(0, mt, params->b, 1e-12);
		for(i=0; i<Na; i++){
			gsl_matrix_set(params->invMatGSLA,j,i,params->b[i]);
		}
	}

	//get contribution starting at st
	gsl_blas_dgemv(CblasTrans, 1.0, params->invMatGSLA, params->ancStateVec, 0.0, params->ancResVec);
	for(i=0; i<Na; i++)params->b[i] = gsl_vector_get(params->ancResVec,i);
	gsl_matrix_set_zero(params->expAFS2);
	fillExpectedAFS_unnorm(params->reducedStateSpace, params->b,params->rates,params->expAFS2);
	
	////////////////////normalized IM AFS
	gsl_matrix_add(params->expAFS,params->expAFS2);
	sum = matrixSumDouble(params->expAFS);
	gsl_matrix_scale(params->expAFS, 1.0 / sum);
	
	
	///////////// Replace AFS with log(AFS)
	// for(i=0;i<params->expAFS->size1;i++){
	// 	for(j=0;j<params->expAFS->size2;j++){
	// 		gsl_matrix_set(params->expAFS,i,j,log(gsl_matrix_get(params->expAFS,i,j)));
	// 	}
	// }
	
	//clean up
	cs_spfree(spMat);
	cs_spfree(mt);
	cs_spfree(ident);
	cs_spfree(eye);
	cs_spfree(tmpMat);
	cs_spfree(tmpMat2);
	gsl_vector_free(tmpStates);
}


// gsl_matrix *afs_im_period1(void * p){
// 	struct im_lik_params * params = (struct im_lik_params *) p;
// 	int i,j, N = params->stateSpace->nstates,n1,n2;
// 	double b[N],invMat[N][N];
// 	struct cs_di_sparse *spMat, *mt, *ident, *eye, *tmpMat;
// 	gsl_matrix *expAFS2;
// 
// 	n1 = params->stateSpace->states[0]->popMats[0]->size1 - 1;
// 	n2 = params->stateSpace->states[0]->popMats[0]->size2 - 1;
// 	//Fill original trans mats	
// 	tmpMat= fillTransitionMatrixArray_embed_ctmc(params->stateSpace, params->topA, params->moveA, &(params->nnz),
// 		params->dim1, params->dim2,1.0,gsl_vector_get(params->paramVector,0),gsl_vector_get(params->paramVector,2),
// 		gsl_vector_get(params->paramVector,3),params->new, params->newdim1, params->newdim2,params->rates);
// 	//Step 1.
// 	//using CSparse
// 	// S = (I-P)^-1
// 	//add negative ident
// 	ident = cs_spalloc(N,N,N,1,1);
// 	for(i=0;i<N;i++) cs_entry(ident,i,i,1);
// 	eye = cs_compress(ident);
// 	spMat = cs_add(eye,tmpMat,1.0,-1.0);
// 	mt = cs_transpose(spMat,1);	
// 	//create unit array for solve
// 	for(i=0; i<N; i++){ 
// 		b[i] = 0.0;
// 	}
// 	b[0] = 1.0;
// 	cs_dropzeros(mt);
// 	cs_lusol(0, mt, b, 1e-12);//b now holds row of interest
// 	//get AFS contribution form time 0 to time INF
// 	gsl_matrix_set_zero(params->expAFS);
// 	fillExpectedAFS_unnorm(params->stateSpace, b,params->rates,params->expAFS);
// 
// 	printf("////////////////////time0-INF\n");
// 	for(i=0;i< (n1+1) ;i++){
// 		for(j=0;j< (n2+1);j++){
// 			printf("%.5f ",gsl_matrix_get(params->expAFS,i,j));
// 		}
// 		printf("\n");
// 	}
// 	
// 	//Step 2
// 	//Exponentiate matrix at time using CTMC representation
// 
// 	sparse_exponential_single_row(N, params->nnz, params->newdim1, params->newdim2,
// 		params->new, gsl_vector_get(params->paramVector,4),b, params->resultTmp, 0); 
// 	
// 	//recompute Inverse
// 	for(j=0;j<N;j++){
// 		//create unit array for solve
// 		b[j]=1.0;
// 		if(j>0)b[j-1]=0.0;
// 		cs_lusol(0, mt, b, 1e-12);
// 		for(i=0; i<N; i++)params->invMat[j][i]=b[i]*params->resultTmp[j];
// 	}
// 	
// 	expAFS2 = gsl_matrix_alloc(params->expAFS->size1,params->expAFS->size2);
// 	gsl_matrix_set_zero(expAFS2);
// //	fillExpectedAFS_unnorm_mat(params->stateSpace, params->invMat,params->rates,expAFS2);
// 
// 	printf("////////////////////time0-t\n");
// 	
// 	gsl_matrix_sub(params->expAFS,expAFS2);
// 	for(i=0;i< (n1+1) ;i++){
// 		for(j=0;j< (n2+1);j++){
// 			printf("%.5f ",gsl_matrix_get(params->expAFS,i,j));
// 		}
// 		printf("\n");
// 	}
// 	cs_spfree(spMat);
// 	cs_spfree(mt);
// 	cs_spfree(ident);
// 	cs_spfree(eye);
// 	cs_spfree(tmpMat);
// }

//////////////////// Uniformization Algorithm for Exponentiation
/// P. 413 from Stewart Book 
void uniformizExpSparseMat(int N, int nnz, double *vals, int *dim1, int *dim2, double tSoln, double *result){
	int i, k, ii;
	double cosi, sigma, nu, gamma, eps,y[N],tmpY[N],scale;
	struct cs_di_sparse *tmpMat, *ctmc, *ident, *eye, *new;
	
	//stuff values into more convenient matrix
	//uniformization of matrix
	tmpMat = cs_spalloc(N,N,N,1,1);
	ident = cs_spalloc(N,N,N,1,1);
	for(i=0;i<nnz;i++)cs_entry(tmpMat,dim1[i],dim2[i],vals[i]);
	ctmc=cs_compress(tmpMat);
	gamma = cs_max(ctmc);
	printf("gamma=%f\n------------------\n",gamma);
	cs_print_adk(ctmc);
	cs_scale(ctmc,1.0/gamma);
	cs_print_adk(ctmc);
	
	for(i=0;i<N;i++) cs_entry(ident,i,i,1);
	eye = cs_compress(ident);
	new = cs_add(eye,ctmc,1.0,1.0);
	cs_print_adk(new);

	//tmpMat now holds the subjugated dtmc
	
	//step one of uniformization is to compute K, number of terms in the summation
	eps=0.001;
	k=0; cosi=1.0; sigma=1.0; nu = (1-eps)/exp(-1.0*gamma*tSoln);
	while(sigma < nu){
		k+=1;
		cosi *= (gamma * tSoln)/k;
		sigma += cosi;
	}
	printf("here is k: %d\n",k);
	
	//initialize result; 
	for(i=1;i<N;i++){
		result[i]=0.0;
		y[i]=0.0;
		tmpY[i]=0.0;
	}
	result[0]=1.0;y[0]=1.0;tmpY[0]=0.0;
	printf("initialization done\n");
	for(ii=0;ii<N;ii++) printf("tmpY[%d]=%f\n",ii,tmpY[ii]);
	printf("---------------------\n");
	for(ii=0;ii<N;ii++)printf("y[%d]=%f\n",ii,y[ii]);
	printf("---------------------\n");
	
	//now approximate pi(t)
	int testVal;
	for(i=1;i<=k;i++){
		printf("k loop k=%d\n",i);
		testVal=cs_gaxpy(new, y, tmpY); //tmpY now contains the answer
		printf("test of gaxpy: %d\n",testVal);
		for(ii=0;ii<N;ii++) printf("tmpY[%d]=%f\n",ii,tmpY[ii]);
		printf("---------------------\n");
		for(ii=0;ii<N;ii++)printf("y[%d]=%f\n",ii,y[ii]);
		printf("---------------------\n");
		cs_print_adk(tmpMat);
		printf("---------------------\n");
		
		for(ii=0;ii<N;ii++){
	//		printf("tmpY[%d]=%f\n",ii,tmpY[ii]);
	//		printf("y[%d]=%f\n",ii,y[ii]);
			y[ii]=tmpY[ii] * ((gamma*tSoln)/i);
			result[ii]+=y[ii];
			tmpY[ii]=0;
		}
		
	}
	//final scaling
	scale = exp(-1.0*gamma*tSoln);
	for(i=0;i<N;i++){
		result[i] *= scale;
		printf("result[%d]=%f\n",i,result[i]);
	}	
}

