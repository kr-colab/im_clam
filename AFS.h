#ifndef AFS_H
#define AFS_H
#include <gsl/gsl_matrix.h> 
#include <glib.h>

//afsObject represents a two population coalescent 
typedef struct afsObject
{
	int npops,nalleles;
	int aCounts[2]; //this keeps track of popn specific allele counts
	gsl_matrix_int **popMats;
	char *matString;
	

}
afsObject;

typedef struct afsStateSpace
{
	int nstates;
	afsObject **states;
	//below are only to hasten the process of state
	//space enumeration
	//int **nstateArray;
	//afsObject ****tmpStorage;
	GHashTable *stateHash;
}
afsStateSpace;


//for use in programs like cmc_island_mle
struct island_lik_params{
	int snpNumber, maxSampleSize,nstates,n1,n2, hStates,nParams;
	gsl_vector *weights, *sampleSizeVector,*rates;
	//paramVector for HMM == [N2,thetaState1,thetaState2,...,thetaStaten, m1State1,m2State1,m1State2,m2State2,...,m1Staten,m2Staten] <n is number of hmm states>;
	gsl_vector *paramVector,*paramCIUpper,*paramCILower, *mlParams; 
	//For straight MLE the paramVector takes the form [N1,N2,m1,m2]
	gsl_matrix *log_afs, *posts, *transMat,*topol,*expAFS, *visitMat;
	gsl_matrix_int *moveType;
	afsStateSpace *stateSpace;
	struct site *data;
	double *cs_work;
	int *dim1, *dim2; //these pointers ease the construction of sparse matrices
	int nzCount;
	double maxLik;
	int optimized;

}island_lik_params;





#define MAXSTATES 100000000
#define MAXSIZE 20

afsObject *afsObjectNew(int n1, int n2);
void afsObjectFree(afsObject *v);
void afsObjectInit(afsObject *v,int n1, int n2);
void afsObjectPrint(afsObject *v);
afsObject *afsObjectNewFrom(afsObject *v);
afsObject *afsObjectDelta(afsObject *a, afsObject *b);
void afsObjectDeltaPre(afsObject *a, afsObject *b,afsObject *tmp);

int afsObjectCount(afsObject *a);

afsStateSpace *afsStateSpaceNew();
void afsStateSpaceFree(afsStateSpace *a);
int afsStateSpaceContainsState(afsStateSpace *S, afsObject *aState);
void afsStateSpaceRemoveAbsorbing(afsStateSpace *S);
int afsStateSpaceContainsStateTmpStorage(afsStateSpace *S, afsObject *aState);
void afsStateSpaceMapPopn(afsStateSpace *S, int *map);
void afsStateSpaceMapAndReducePopn(afsStateSpace *S, int *map, afsStateSpace *reducedSpace, int *reverseMap);

void afsStateSpaceRemovePopulation(afsStateSpace *S, int popn);
afsStateSpace *afsStateSpaceCopy(afsStateSpace *S);

void mcMatsImportFromFile(const char *fileName,int *nnz,double *topol, int *moveA, int *dim1, int *dim2);


int nDescPop1(afsObject *v);
int nDescPop2(afsObject *v);
int compatibleAfsGivenInitial(afsObject *v, int n1, int n2);
int afsObjectsEqual(afsObject *a, afsObject *b);
void acceptOrRejectProposedState(afsStateSpace *S, afsObject *child, int n1, int n2);
int afsObjectQsortCompare(const void *p1, const void *p2);


afsStateSpace *afsStateSpaceImportFromFile(const char *fileName);

void nonZeroEntries(gsl_matrix_int *m, int **coords, int *count);
void entriesGreaterThan(gsl_matrix_int *m, int **coords, int *count, int cut);
void entriesLessThan(gsl_matrix_int *m, int **coords, int *count, int cut);
int countNegativeEntries(gsl_matrix_int *m);

double xchoosey(int n, int k);

int matrixAllGTZ(gsl_matrix_int *m);

int matrixSum(gsl_matrix_int *m);
void gsl_matrix_prettyPrint(gsl_matrix *m);
double matrixSumDouble(gsl_matrix *m);

void preorderEvolveState(afsStateSpace *S, afsObject *aState, int n1, int n2);

void coalMarkovChainTopologyMatrix(afsStateSpace *S,gsl_matrix *topol, gsl_matrix_int *moveType, int n1, int n2);
int coalMarkovChainTopologyMatrix_2(afsStateSpace *S,gsl_matrix *topol, gsl_matrix_int *moveType, int *dim1, int *dim2);
int coalMarkovChainTopologyMatrix_sparse(afsStateSpace *S,double *topol, int *moveType, int *dim1, int *dim2);

void readTopolFile(const char *fileName, double *topA, int *moveA, int *dim1, int *dim2);


void fillTransitionMatrix(afsStateSpace *S, gsl_matrix *topol, gsl_matrix_int *moveType, gsl_matrix *transMat, gsl_vector *rates,
	double theta1, double theta2, double mig1, double mig2);
struct cs_di_sparse *fillTransitionMatrix_2(afsStateSpace *S, gsl_matrix *topol, gsl_matrix_int *moveType, gsl_vector *rates, int nzCount,
	int *dim1, int *dim2, double theta1, double theta2, double mig1, double mig2);	
	
void fillExpectedAFS(afsStateSpace *S, gsl_matrix *transMat,gsl_vector *rates, gsl_matrix *expAFS, int n1, int n2);

void fillExpectedAFS_visits(afsStateSpace *S, gsl_matrix *visitMat,gsl_vector *rates, gsl_matrix *expAFS);
void fillExpectedAFS_visits_2(afsStateSpace *S, double *visitMat,gsl_vector *rates, gsl_matrix *expAFS);

void overlayPopnSplitting(afsStateSpace *S,gsl_matrix *topol, gsl_matrix_int *moveType, int n1, int n2);


void maximize_island_lik_const(double *lik, void * p);
void maximize_island_lik_hmm_const(double *lik, void * p);

void maximize_island_lik_CILower(double *lik, void * p);
void maximize_island_lik_CIUpper(double *lik, void * p);
void island_lik_constraint_CI(size_t n, const double *pGuess, void * p, double *result);

void island_lik_constraint(size_t n, const double *pGuess, void * p, double *result);
void island_lik_constraint_hmm(size_t n, const double *pGuess, void * p, double *result);

void fillLogAFS_ptr_sp(void * p);
void fillLogAFS_ptr_sp2(void * p);
void fillLogAFS_ptr_sp_hmm(int state,void * p);

void itoa(int n, char s[]);
void reverse(char s[]);

 
void NumericalGradient_island_lik(size_t n, const double *v, void *params, double *df);
void island_lik_fdf(const size_t n, const double *x,void *fparams,double *fval,double *grad);

#endif
