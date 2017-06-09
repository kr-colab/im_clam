#include <slepcmfn.h>
#include "cs.h"
#include <gsl/gsl_rng.h>
//for use in programs like im_clam
typedef struct clam_lik_params{
	int snpNumber, maxSampleSize,nstates,n1,n2,nParams;
	gsl_vector *sampleSizeVector,*rates, *stateVec, *resVec ;
	gsl_vector *paramVector,*paramCIUpper,*paramCILower, *mlParams; 
	//For straight MLE the paramVector takes the form [N2,NA,m1,m2,t]
	gsl_matrix *log_afs, *transMat,*expAFS, *expAFS2, *invMatGSL, *invMatGSLA, *obsData;
	afsStateSpace *stateSpace ,*reducedStateSpace;
	gsl_rng *rng;
	struct site *data;
	double *cs_work,*new,*topA,*resultTmp,**invMat,*b,*res, *st, *top;
	int *dim1, *dim2, *moveA, *move, *dim1A, *dim2A; //these pointers ease the construction of sparse matrices
	int nnz, nnzA;
	double maxLik;
	int *map, *reverseMap; //these provide mapping of between state spaces of two phases
	int optimized;
	int rank;
	int Na,fEvals;
	double *expoArray;
	double meanTreeLength;
	cs_di *spMat, *eye, *eyeAnc, *triplet, *tmpMat, *tmpMat2;

	Vec  x_0, x, u,bP, u2, bP2, v, y, v_seq, *x2, ancStateVec, ancResVec;
	Vec xInv,bInv;
	Mat            ident,C, C2, C_transpose, D, D_copy, F, denseIdent, denseMat1, denseMat2, denseMat3;
	VecScatter	ctx;
	MFN                mfn;
	PetscInt        iStart,iEnd;
      /* linear solver context */
	PC             pc;           /* preconditioner context */
	
}clam_lik_params;
