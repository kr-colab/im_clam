#ifndef CTMC_h
#define CTMC_h
//expoMatObj--basically a sparse matrix representation
typedef struct expoMatObj
{
	int rank,nnz;
	int *ai, *aj; //this is COO formated sparse matrix format
	double *a; //matrix values
	double *initVec; //initial state vector
	double *resultSpace; // N x N matrix to store possible dense output
	double intUp; //upper limit on definite integral

}
expoMatObj;

//for use in programs like cmc_island_mle
typedef struct im_lik_params{
	int snpNumber, maxSampleSize,nstates,n1,n2, hStates,nParams;
	gsl_vector *weights, *sampleSizeVector,*rates, *stateVec, *resVec, *ancStateVec, *ancResVec;
	gsl_vector *paramVector,*paramCIUpper,*paramCILower, *mlParams; 
	//For straight MLE the paramVector takes the form [N2,NA,m1,m2,t]
	gsl_matrix *log_afs, *transMat,*expAFS, *expAFS2, *invMatGSL, *invMatGSLA;
	afsStateSpace *stateSpace ,*reducedStateSpace;
	struct site *data;
	double *cs_work,*new,*topA,*resultTmp,**invMat,*b,*res, *st, *topA2;
	int *dim1, *dim2, *moveA, *newdim1, *newdim2, *expodim1, *expodim2; //these pointers ease the construction of sparse matrices
	int *dim1A, *dim2A, *moveA2;
	int nnz, nnzA;
	double maxLik;
	int *map, *reverseMap; //these provide mapping of between state spaces of two phases
	int optimized;
	double *expoResult, *expoInit;

}im_lik_params;

int integrateTimeDistribution(struct expoMatObj *foo);

struct cs_di_sparse *fillTransitionMatrix_ctmc(afsStateSpace *S, double *topol, int *moveType, int *nzCount,
	int *dim1, int *dim2, double theta1, double theta2, double mig1, double mig2);

void fillTransitionMatrixArray_ctmc(afsStateSpace *S, double *topol, int *moveType, int *nzCount,
	int *dim1, int *dim2, double theta1, double theta2, double mig1, double mig2,double *newArray,
	int *newDim1, int *newDim2);
	
	
void fillExpectedAFS_ctmc(afsStateSpace *S, double *visitMat, gsl_matrix *expAFS);
void fillExpectedAFS_unnorm(afsStateSpace *S, double *visitMat ,gsl_vector *rates, gsl_matrix *expAFS);

void fillExpectedAFS_unnorm_mat(afsStateSpace *S, double **visitMat , double *initialStateProbs, gsl_vector *rates, gsl_matrix *expAFS);

struct cs_di_sparse * fillTransitionMatrixArray_embed_ctmc(afsStateSpace *S, double *topol, int *moveType, int *nzCount,
	int *dim1, int *dim2, double theta1, double theta2, double mig1, double mig2,double *newArray,
	int *newDim1, int *newDim2, gsl_vector *rates);
	
//static int matExpIntegrand(const int *ndim, const double xx[],const int *ncomp, double ff[], void *userdata);
	
void uniformizExpSparseMat(int N, int nnz, double *vals, int *dim1, int *dim2, double tSoln, double *result);
void fillLogAFS(void * p);
	
#endif
