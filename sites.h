//sampleConfig-- simple holder to represent a single datum
typedef struct site{
	long pos;
	int n1,n2, p1, p2;
	}site;



int siteFileImport(char *fileName, struct site *data);
