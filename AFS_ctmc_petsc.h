#ifndef CTMC_PETSC
#define CTMC_PETSC

#include <petscksp.h>
extern double lowerBounds[5];
extern double upperBounds[5];

#endif

extern PetscBool vbse;



void fillPetscTransMats(afsStateSpace *S, double *topol, int *moveType, int *nzCount,
	int *dim1, int *dim2, double theta1, double theta2, double mig1, double mig2, Mat *DTMC, Mat *CTMC, gsl_vector *rates);
	
struct cs_di_sparse * fillPetscCsparseTransMats(afsStateSpace *S, double *topol, int *moveType, int *nzCount,
	int *dim1, int *dim2, double theta1, double theta2, double mig1, double mig2, Mat *CTMC, gsl_vector *rates);

struct cs_di_sparse * fillPetscCsparseTransMats_prealloc(afsStateSpace *S, double *topol, int *moveType, int *nzCount,
	int *dim1, int *dim2, double theta1, double theta2, double mig1, double mig2,
	Mat *CTMC, gsl_vector *rates, struct cs_di_sparse *triplet);
	
void calcLogAFS_IM(void * p);
void calcLogAFS_IM_allPETSC(void * p);
double calcLikNLOpt(unsigned n, const double *point, double *gradients, void *p);
double calcLikNLOpt_gradients(unsigned n, const double *point, double *gradients, void *p);

void maximizeLikNLOpt(double *lik, void *p, double *mle);
void maximizeLikNLOpt_twoStage(double *lik, void *p, double *mle);
void maximizeLikNLOpt_MLSL(double *lik, void *p, double *mle);

double hessianMatrix_element(double lik, double *mle, int i, int j, double hi, double hj, void *p);
gsl_matrix *hessian(double *mle, double lik, void *p);
gsl_matrix *getFisherInfoMatrix(double *mle, double lik, void *p);
gsl_vector *getGradient(double *mle, void *p);
gsl_matrix *getGodambeInfoMatrix(double *mle, double lik, void *p);


