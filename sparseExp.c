#include "math.h"
#include "sparseExp.h"
#include <stdio.h>
#include <stdlib.h>
/////////////////////////////////////////////////////  Exponential Stuff

void sparse_exponential(int n, int nnz, int *ia, int *ja, double *a, double t, double *initialArray, double *res){ 
	int i;
	int m = n-1; 
	double w[n];
	int ideg = 6;
	double tol = 1;
	int iflag = 0;
	int lwsp = n*(m+1)+n+pow((m+2.),2)+4*pow((m+2.),2)+ideg+1;
	double  wsp[lwsp];
	int liwsp = m+2;
	int iwsp[liwsp];
	double t1 = t;
	double anorm = 0;
	int itrace = 0;

	for(i=0;i<=n*n;i++)res[i]=0.0;	
	wrapalldmexpv_(&n,&m,&t1,initialArray,w,&tol,&anorm,wsp,&lwsp,iwsp,&liwsp,
			&itrace,&iflag,ia,ja,a,&nnz,res);

}


void sparse_exponential_single_row(int n, int nnz, int *ia, int *ja, double *a, double t, double *initialArray, double *res, int row){ 
	int i;
	int m = n-1; 
	double v[n],w[n];
	int ideg = 16;
	double tol = 0.1;
	int iflag = 0;
	int lwsp = n*(m+1)+n+pow((m+2.),2)+4*pow((m+2.),2)+ideg+1;
	double  *wsp= malloc(sizeof(double) *lwsp);
	int liwsp = m+2;
//	printf("lwsp %d comp %d liwsp %d comp %d\n",lwsp,(n*(m+2))+(5*(m+2)*(m+2))+ideg+1,liwsp,liwsp,m+2);

	int *iwsp = malloc(sizeof(int) *liwsp);
	double t1 = t;
	double anorm = 0;
	int itrace = 0;
	
	for( i=0;i<n;i++) v[i] = 0;
	v[row]= 1;
	for(i=0;i<=n*n;i++)res[i]=0.0;	
	
	wrapsingledmexpv_(&n,&m,&t1,v,w,&tol,&anorm,wsp,&lwsp,iwsp,&liwsp,
			&itrace,&iflag,ja,ia,a,&nnz,res);
	free(wsp);	
	free(iwsp);
}
