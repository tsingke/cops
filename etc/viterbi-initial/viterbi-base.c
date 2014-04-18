
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<time.h>
#include<ctype.h>
#include<float.h>
#include<math.h>


typedef unsigned int	uint;
typedef unsigned char	uchar;
typedef unsigned char	bool;
typedef unsigned char	byte;

typedef double TYPE;

FILE *outfp;


#define MAX_LINE		2000
#define MAX_MAT_SIZE	200


#define MAX2(a,b)	( (b > a)	? b : a)

#define MAX3(a,b,c)	( (c > b)	? ((c > a)? c : a)	\
								: ((b > a)? b : a) )


char *usage = "Viterbi Decoder Tool\nUsage: program transition-probs-file emission-probs-file";

void process_args(int argc, char** argv, FILE** emissionsfp, FILE** transprobfp)
{
	if (argc < 1)
		exit(fprintf(stderr, "%s\n", usage));
	
	*emissionsfp = fopen(argv[1], "r");
	if (*emissionsfp	== NULL)
		exit(fprintf(stderr, "ERROR file %s not found\n%s\n", argv[1], usage));

	if (transprobfp)
	{
		*transprobfp = fopen(argv[2], "r");
		if (*transprobfp == NULL)
			exit(fprintf(stderr, "ERROR file %s not found\n%s\n", argv[2], usage));
	}
}


// auxiliay
int read_clean_line(FILE *fp, char* buf)
{
	char *ret = fgets(buf, MAX_LINE, fp);
	int j;
	
	if (ret == NULL)
		return -1;

	// trim end
	for (j = strlen(buf)-1; buf[j] <= ' ' && j > 0; j--)
		buf[j] = '\0';
	// trim start
	for (j = 0; buf[j] <= ' ' && buf[j]; j++);

	// remove comments
	if (strlen(buf+j) < 3 )
		return 0; // ignore

	if (buf[j] == '#' || buf[j] == '>' || buf[j] == '/' )
		return 0;

	return 1; // use
}



typedef struct _hmm
{	char tokens[100];
	char transi[10][2];
	int tokenToIdx['Z'+1];

	double* matchprobs[200];
	double* inserprobs[200];
	double* transprobs[200];
	int ntokens;
	int ntrans;
	int nstates;
	
	double** resultsDel;
	double** resultsIns;
	double** resultsMat;
	
	int** tracebacksDel;
	int** tracebacksIns;
	int** tracebacksMat;
} hmm;

#define NEG_INFINITY (DBL_MIN/2)

#define CHECK(val,msg)	\
{	int aaa = 0;	\
	if (isnan(val))	aaa = fprintf(stderr, "NAN in %s val: %f\n", msg, val);	\
	else if (isinf(val)) aaa = fprintf(stderr, "INF in %s val: %f\n", msg, val);	\
}

//	if (aaa != 0) sleep(10000);
//	else if (val > 0.0)	 aaa = fprintf(stderr, "POS in %s val: %f\n", msg, val);	
//	else if (val == 0.0)	fprintf(stderr, "ZERO in %s val: %f | log(val): %f\n", msg, val, log(val));	
//	else if (val < 0.000001)	fprintf(stderr, "SMALL %s val: %f\n", msg, val);	
//	else if (val < 0.0)	 aaa = fprintf(stderr, "NEG in %s val: %f\n", msg, val);	


/* function used in Hmmer
static void printprob(FILE *fp, int fieldwidth, double p)
{
	if      (p == 0.0) fprintf(fp, " %*s",   fieldwidth, "*");
	else if (p == 1.0) fprintf(fp, " %*.5f", fieldwidth, 0.0);
	else               fprintf(fp, " %*.5f", fieldwidth, -logf(p));
}
*/

double convertplog(double plog)
{
	double ret = exp(plog);

//	if (ret > 1.0) printf("plog %f => prob %f\n", plog, ret);
	return plog;
	
	
	CHECK(ret, "viterbi");
	if (ret > 1.1) fprintf(stdout, "IN VITERBI Prob > 1 !!!: %f\n", ret);
	if (isinf(ret))
	{	fprintf(stderr, "IN VITERBI infinity!!!: %f\n", ret);
		return 0.0;
	}
	return ret;
}

double convertprob(char* data)
{
	if (*data == '*')
	{	// printf("Read Neg Infinity\n");
		return NEG_INFINITY;
	}

	double plog = atof(data);
	if (plog < 0.0001)
		return 1.0;

	double ret = exp(-plog);
	if (ret > 1.1)
		fprintf(stderr, "Prob > 1 !!!: %f\n", ret);

	if (ret < 0.0001)
		fprintf(stderr, "Prob 00 !!!: %f\n", ret);
	return ret;
}

hmm* load_model(FILE *matfp)
{
	hmm* hhh = calloc(1, sizeof(hmm));
	char buf[4000], *pp;
	char *delims = " \t\n\r";
	int* tokenToIdx = hhh->tokenToIdx;
	char* tokens = hhh->tokens;
	char (*transi)[2]  = hhh->transi;
	double** matchprobs = hhh->matchprobs;
	double** inserprobs = hhh->inserprobs;
	double** transprobs = hhh->transprobs;
	int Ntokens, Ntrans;
	int i,j;
	
	// process 1st line: with the tokens
	while (read_clean_line(matfp, buf) < 1);
	pp = strtok(buf, delims);
	for (j = 0; pp != NULL; pp = strtok(NULL, delims))
		if (strcmp(pp, "HMM"))	// if diff
		{	tokenToIdx[*pp] = j;
			tokens[j++] = *pp;
		}
	hhh->ntokens = Ntokens = j;
	
	// process 2nd line: with transitions
	while (read_clean_line(matfp, buf) < 1);
	pp = strtok(buf, delims);
	for (j = 0; pp != NULL;  pp = strtok(NULL, delims))
		if (pp[1] == '-' && pp[2] == '>')	// if diff
		{
			transi[j][0] = pp[0];
			transi[j][1] = pp[3];
			j++;
		}
	hhh->ntrans = Ntrans = j;
	
	// process data
	for(i = 0; !feof(matfp); i++)
	{
		double* ppr;
		double total = 0.0;
	
		while ((j = read_clean_line(matfp, buf)) == 0);
		if (j < 0) break;
	
		// line with match probs
		ppr = matchprobs[i] = (double*) malloc((Ntokens+1)*sizeof(double));
		pp = strtok(buf, delims);
		pp = strtok(NULL, delims);	// ignore ID
		for (j = 0; j < Ntokens;  pp = strtok(NULL, delims))
			total += ppr[j++] = convertprob(pp);
		ppr[j] = 0;

		total = 0;

		// line with insert probs
		ppr = inserprobs[i] = (double*) malloc((Ntokens+1)*sizeof(double));
		while (read_clean_line(matfp, buf) == 0);
		pp = strtok(buf, delims);
		for (j = 0; j < Ntokens;  pp = strtok(NULL, delims))
			total += ppr[j++] = convertprob(pp);
		ppr[j] = 0;

		total = 0;

		// line with transition probs
		ppr = transprobs[i] = (double*) malloc((Ntrans+1)*sizeof(double));
		while (read_clean_line(matfp, buf) == 0);
		pp = strtok(buf, delims);
		for (j = 0; j < Ntrans ; pp = strtok(NULL, delims))
			ppr[j++] = convertprob(pp);
		ppr[j] = 0;
	}
	
	hhh->nstates = i-1;
	return hhh;
}

void print_model(hmm* model)
{
	int i, j;
	char* tokens = model->tokens;
	char (*transi)[2] = model->transi;
	double** matchprobs = model->matchprobs;
	double** inserprobs = model->inserprobs;
	double** transprobs = model->transprobs;
	FILE* mmoutfp = stdin;

	fprintf(mmoutfp, "HMM     ");
	for (i = 0; i < model->ntokens; i++)
		fprintf(mmoutfp, "     %c   ", tokens[i]);

	fprintf(mmoutfp, "\n       ");
	for (i = 0; i < model->ntrans; i++)
		fprintf(mmoutfp, "     %c->%c", transi[i][0], transi[i][1]);
	
	for (i = 0; matchprobs[i]; i++)
	{
		double* ppr;
		
		fprintf(mmoutfp, "\n   %4d ", i);
		ppr = matchprobs[i];
		for (j = 0; j < model->ntokens; j++)
		{	if (ppr[j] < 0.001) fprintf(stderr, "MATCH  PROB == 0!!!\n");
			fprintf(mmoutfp, "  %.5f", ppr[j]);
		}			

		fprintf(mmoutfp, "\n        ");
		ppr = inserprobs[i];
		for (j = 0; j < model->ntokens; j++)
		{	if (ppr[j] < 0.01) fprintf(stderr, "INSERT PROB == 0!!!\n");
			fprintf(mmoutfp, "  %.5f", ppr[j]);
		}			
		fprintf(mmoutfp, "\n        ");
		ppr = transprobs[i];
		for (j = 0; j < model->ntrans; j++)
		{	// if (ppr[j] < 0.00001) fprintf(stderr, "TRANS PROB == 0! Log(prob) = %f!!\n", log(ppr[j]));
			fprintf(mmoutfp, "  %.5f", ppr[j]);
		}			
	}
	fprintf(mmoutfp, "\n");
}


void viterbi_decode(hmm* model, char* seq)
{
	int i, j, k;
	char* tokens = model->tokens;
	char (*transi)[2] = model->transi;
	double** matchprobs = model->matchprobs;
	double** inserprobs = model->inserprobs;
	double** transprobs = model->transprobs;
	double **resDel = model->resultsDel;
	double **resIns = model->resultsIns;
	double **resMat = model->resultsMat;
	int **tbDel = model->tracebacksDel;
	int **tbIns = model->tracebacksIns;
	int **tbMat = model->tracebacksMat;
	
//	0		1		2		3		4		5		6
//	m->m	m->i	m->d	i->m	i->i	d->m	d->d
#define EMISSIONPROB(probs,a,b) ((probs == inserprobs)? 1.0 : probs[a][b]/inserprobs[a][b])

	// Initiliaze for state 0
	for (i = 0; i < strlen(seq); i++)	// seq token
		resMat[0][i] = log(1.0);
	for (i = 0; i < strlen(seq); i++)	// seq token
		resDel[0][i] = log(1.0);


	resIns[0][0] = log(1.0);
	for (k = 0, i = 1; i < strlen(seq); i++)	// seq token
	{
		int idx = model->tokenToIdx[seq[i]];
		
		resIns[k][i]= log(EMISSIONPROB(inserprobs,k,idx))
					+ MAX2(	resMat[k][i-1] + log( transprobs[k][1]),	// Mj-1Ij
							resIns[k][i-1] + log( transprobs[k][4]));	// Ij-1Ij
						//	resDel[k][i-1] + log( transprobs[k][xx])	// Dj-1Ij

/*		printf("Em: %f, resMat[k][i-1]: %f, log( transprobs[k][1]): %f\nresIns[k][i-1]: %f, log( transprobs[k][4])): %f || Max2: %f, RESULT: %f\n", 
			log(EMISSIONPROB(inserprobs,k,idx)), resMat[k][i-1], log( transprobs[k][1]), resIns[k][i-1], log( transprobs[k][4]),
			MAX2( resMat[k][i-1] + log( transprobs[k][1]), resIns[k][i-1] + log( transprobs[k][4])), resIns[k][i]);
*/
		CHECK(resIns[k][i], "ins");
	}

	for (k = 1; k < model->nstates; k++)	// state
	{
//		printf("run state %d\n", k);

		i = 0;
		int idx = model->tokenToIdx[seq[0]];
		
		resMat[k][i]= log(EMISSIONPROB(matchprobs,k,idx)) + log(NEG_INFINITY);
		CHECK(resMat[k][i], "mat_init");
		
		resIns[k][i]= log(EMISSIONPROB(inserprobs,k,idx)) + log(NEG_INFINITY);
		
		resDel[k][i]= MAX2(	resMat[k-1][i  ] + log( transprobs[k-1][2]),	// Mj-1Dj
						//	resIns[k-1][i  ] + log( transprobs[k-1][xx]),	// Ij-1Dj
							resDel[k-1][i  ] + log( transprobs[k-1][6]));	// Dj-1Dj
		
		for (i = 1; i < strlen(seq); i++)	// seq token
		{
			int idx = model->tokenToIdx[seq[i]];

//			printf("seq %d. matchprobs: %f, inserprobs: %f, div: %f | logodds: %f => exp(logodds): %f\n",
//					 i, matchprobs[k][idx], inserprobs[k][idx], matchprobs[k][idx]/inserprobs[k][idx], log(matchprobs[k][idx])-log(inserprobs[k][idx]),
//					exp(log(matchprobs[k][idx])-log(inserprobs[k][idx])));
			
			resMat[k][i]= log(EMISSIONPROB(matchprobs,k,idx))
						+ MAX3(	resMat[k-1][i-1] + log( transprobs[k-1][0]),	// Mj-1Mj
								resIns[k-1][i-1] + log( transprobs[k-1][3]),	// Ij-1Mj
								resDel[k-1][i-1] + log( transprobs[k-1][5]));	// Dj-1Mj
			CHECK(resMat[k][i], "mat");
	
			resIns[k][i]= log(EMISSIONPROB(inserprobs,k,idx))
						+ MAX2(	resMat[k  ][i-1] + log( transprobs[k][1]),		// Mj-1Ij
								resIns[k  ][i-1] + log( transprobs[k][4]));		// Ij-1Ij
							//	resDel[k  ][i-1] + log( transprobs[k][xx])		// Dj-1Ij

			resDel[k][i]= MAX2(	resMat[k-1][i  ] + log( transprobs[k-1][2]),	// Mj-1Dj
							//	resIns[k-1][i  ] + log( transprobs[k-1][xx]),	// Ij-1Dj
								resDel[k-1][i  ] + log( transprobs[k-1][6]));	// Dj-1Dj
		}
		
//		printf("run state %d\n", k);
		for (i = 0; i < strlen(seq); i++)	// seq token
		{
			CHECK(resIns[k][i], "ins");
			CHECK(resMat[k][i], "mat");
			CHECK(resDel[k][i], "del");
		}
	}
}



void print_results(double** results, char* seq, FILE* resfp, int nstates)
{
	int i, k;
	
	fprintf(resfp, "     ");
	for (i = 0; i < strlen(seq); i++)	// seq token
		fprintf(resfp, "  %c           ", seq[i]);
	
	for (k = 0; k < nstates; k++)	// state
	{
		fprintf(resfp, "\n%3d  ", k);
		for (i = 0; i < strlen(seq); i++)	// seq token
			fprintf(resfp, " %11.5f", convertplog(results[k][i]));
	}

	fprintf(resfp, "\n");
}


int main(int argc, char** argv)
{
	int i, j;
	int numstates, numsymbols, seqlen=0;
	char *seq;

	FILE* modelfp;
	outfp = stdout;

	if (argc < 3)
		return 0;
	
	process_args(argc, argv, &modelfp, NULL);

	hmm* model = load_model(modelfp);

	print_model(model);

	seq = argv[2];
	seqlen = strlen(seq);
	
	model->resultsDel		= calloc(sizeof(double*), (model->nstates+2));
	model->resultsIns		= calloc(sizeof(double*), (model->nstates+2));
	model->resultsMat		= calloc(sizeof(double*), (model->nstates+2));
	model->tracebacksDel	= calloc(sizeof(int  *), (model->nstates+2));
	model->tracebacksIns	= calloc(sizeof(int  *), (model->nstates+2));
	model->tracebacksMat	= calloc(sizeof(int  *), (model->nstates+2));
	
	for (i = 0; i <= model->nstates; i++)
	{
		model->resultsDel[i]  = calloc(sizeof(double), seqlen+1);
		model->resultsIns[i]  = calloc(sizeof(double), seqlen+1);
		model->resultsMat[i]  = calloc(sizeof(double), seqlen+1);
		
		model->tracebacksDel[i] = calloc(sizeof(int), seqlen+1);
		model->tracebacksIns[i] = calloc(sizeof(int), seqlen+1);
		model->tracebacksMat[i] = calloc(sizeof(int), seqlen+1);
	}
	
	viterbi_decode(model, seq);
	
	print_results(model->resultsMat, seq, fopen("results-matchx.txt", "w"), model->nstates);
	print_results(model->resultsIns, seq, fopen("results-insert.txt", "w"), model->nstates);
	print_results(model->resultsDel, seq, fopen("results-delete.txt", "w"), model->nstates);


	return 0;
}


