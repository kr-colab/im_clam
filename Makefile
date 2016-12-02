CC = gcc
CCPLUS = g++
FC = gfortran
MY_CFLAGS =  -Wall -lm -lgsl -lcxsparse -llapack -lblas -lgslcblas -I/usr/include/suitesparse -I/usr/include/glib-2.0 -g
COMMON = AFS.c sites.c multimin.c adkGSL.c adkCSparse.c

HEADERS = AFS.h 
GLIB = `pkg-config --cflags glib-2.0 gsl nlopt` `pkg-config --libs glib-2.0 gsl nlopt`


#include ${PETSC_DIR}/conf/variables
#include ${PETSC_DIR}/conf/rules
include ${SLEPC_DIR}/conf/slepc_common

all:	cmc_stateSpace cmc_topol im_clam

cmc_stateSpace:	cmc_stateSpace.c AFS.c 
	$(CC) cmc_stateSpace.c $(COMMON) -O3 -lm -o cmc_stateSpace $(GLIB) $(MY_CFLAGS)

cmc_topol:	cmc_topol.c AFS.c AFS_ctmc.c AFS_pthreads.c
	$(CC) cmc_topol.c $(COMMON) -O3 -lm -o cmc_topol $(GLIB) $(MY_CFLAGS) 

im_clam:	im_clam.c AFS.c adkCSparse.c AFS_ctmc.c AFS_ctmc_petsc.c  chkopts
	${CLINKER} im_clam.c AFS_ctmc.c AFS_ctmc_petsc.c -o im_clam ${PETSC_KSP_LIB} ${SLEPC_MFN_LIB} $(COMMON) $(GLIB) -lnlopt -I ${PETSC_DIR}/include/ -I ${SLEPC_DIR}/include/ -I ${SLEPC_DIR}/${PETSC_ARCH}/include/ $(MY_CFLAGS)   

#clean:	im_clam cmc_stateSpace cmc_topol
#	rm -f im_clam cmc_stateSpace cmc_topol

#test:  test.c adkGSL.c
#	$(CC) $(MY_CFLAGS) test.c adkGSL.c -o test

#ex1: ex1.o  chkopts
#	${CLINKER} -o ex1 ex1.o  ${PETSC_KSP_LIB}
#	${RM} ex1.o
