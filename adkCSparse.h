struct cs_di_sparse *gslMatrixToCSsparse(gsl_matrix *gm, int numNonZero);
struct cs_di_sparse *arraysToCSsparse(int rank, int numNonZero, int *dim1, int *dim2,double *x);
double cs_max(struct cs_di_sparse *mat);
void cs_scale(struct cs_di_sparse *mat, double c);
void cs_print_adk(struct cs_di_sparse *mat);
