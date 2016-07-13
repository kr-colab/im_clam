/* siteFile utilities
/
/
/ Andrew Kern
*/

#include "stdlib.h"
#include "stdio.h"
#include "math.h"
#include "assert.h"
#include "sites.h"


//siteFileImport5-- reads a file and stores info into pre-alloc'd data for 5 column
//returns snpNumber
int siteFileImport(char *fileName, struct site *data){
	FILE *infile;
	long pos;
	int  p1, p2, n1, n2,j;

	/* open file, errors? */
	infile = fopen(fileName, "r");
	if (infile == NULL){
		fprintf(stderr,"Error opening infile! ARRRRR!!!!\n");
		exit(1);
	}
	/* go through infile and get snp info; expect pos, i, n */
	j = 0;
	while (fscanf(infile, "%ld %d %d %d %d", &pos, &n1, &n2, &p1, &p2) != EOF){
		data[j].pos = pos;
		data[j].p1 = p1;
		data[j].p2 = p2;
		data[j].n1 = n1;
		data[j].n2 = n2;
		j += 1;
	}
	fclose(infile);
	return(j);
}