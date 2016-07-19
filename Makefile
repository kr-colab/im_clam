CC = gcc
CCPLUS = g++
FC = gfortran
MY_CFLAGS =  -Wall -lm -lgsl -lcxsparse -llapack -lblas -lgslcblas -I/usr/include/suitesparse
COMMON = AFS.c sites.c multimin.c adkGSL.c adkCSparse.c
TWOLOCCOMMON = twoLocusAFS.c sites.c multimin.c adkGSL.c adkCSparse.c
CTMCCOMMON = AFS_ctmc.c sparseExp.c
HEADERS = AFS.h 
GLIB = `pkg-config --cflags glib-2.0` `pkg-config --libs glib-2.0`

FORTRAN = mataid.f my_expokit.f my_matexp.f blas.f lapack.f clock.f
FOBJECTS =  mataid.o my_expokit.o my_matexp.o clock.o 

#PETSC_DIR = /Users/adk/src/petsc-3.5.3/


include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules
include ${SLEPC_DIR}/conf/slepc_common

all:	fortran cmc_stateSpace cmc_eAFS cmc_im_mle

fortran: mataid.f my_expokit.f my_matexp.f clock.f
	gfortran -c $(FORTRAN) -O3

cmc_stateSpace:	cmc_stateSpace.c AFS.c 
	$(CC) cmc_stateSpace.c $(COMMON) -O3 -lm -o cmc_stateSpace $(GLIB) $(MY_CFLAGS)

cmc_topol:	cmc_topol.c AFS.c AFS_ctmc.c AFS_pthreads.c
	$(CC) cmc_topol.c $(COMMON) -O3 -lm -o cmc_topol $(GLIB) $(MY_CFLAGS) 

im_clam:	im_clam.c AFS.c adkCSparse.c AFS_ctmc.c AFS_ctmc_petsc.c $(FOBJECTS) chkopts
	${CLINKER} im_clam.c AFS_ctmc.c AFS_ctmc_petsc.c -o im_clam ${PETSC_KSP_LIB} ${SLEPC_MFN_LIB} $(MY_CFLAGS) $(COMMON) $(FOBJECTS) -lgfortran -lnlopt -I ${PETSC_DIR}/include/ -I ${SLEPC_DIR}/include/ -I ${SLEPC_DIR}/${PETSC_ARCH}/include/ -O3 -lgsl

#clean:
#	rm -f island cmc_stateSpace cmc_eAFS cmc_im_mle $(FOBJECTS)

test:  test.c adkGSL.c
	$(CC) $(MY_CFLAGS) test.c adkGSL.c -o test

#ex1: ex1.o  chkopts
#	${CLINKER} -o ex1 ex1.o  ${PETSC_KSP_LIB}
#	${RM} ex1.o
