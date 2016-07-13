void sparse_exponential(int n, int nnz, int *ia, int *ja, double *a, double t, double *initialArray, double *res);    
void sparse_exponential_single_row(int n, int nnz, int *ia, int *ja, double *a, double t, double *initialArray, double *res, int row);

void wrapsingledmexpv_(int * n,int* m,double * t,double* v,double * w,double* tol,
	double* anorm,double* wsp,int * lwsp,int* iwsp,int *liwsp, int * itrace,int *iflag,
			int *ia, int *ja, double *a, int *nz, double * res);

 void wrapalldmexpv_(int * n,int* m,double * t,double* v,double * w,double* tol,
                double* anorm,double* wsp,int * lwsp,int* iwsp,int *liwsp, int * itrace,int *iflag,
                int *ia, int *ja, double *a, int *nz, double * res);

void myDMEXPV_(int* n, int* m, double* t, double* v, double* w, double* tol,
	double* anorm, double* wsp, int* lwsp, int* iwsp, int* liwsp, int* itrace, int* iflag,
		int* ia, int* ja, double* a, int* nz );

 void wrapdgpadm_(int * ideg,int * m,double * t,double * H,int * ldh,
                        double * wsp,int * lwsp,int * ipiv,int * iexph,int *ns,int *iflag );