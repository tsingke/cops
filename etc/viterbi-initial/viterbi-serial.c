#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<time.h>
#include<ctype.h>
#include<float.h>
#include<math.h>

#include<sys/timeb.h>

typedef unsigned int	uint;
typedef unsigned char	uchar;
typedef unsigned char	bool;
typedef unsigned char	byte;

typedef float TYPE;

FILE *outfp;

#define DEBUG 0


#define MAX_LINE		MAX_LENGTH
#define MAX_MAT_SIZE	MAX_STATES
#define MAX2(a,b)	( (b > a)	? b : a)
#define MAX3(a,b,c)	( (c > b)	? ((c > a)? c : a)	\
								: ((b > a)? b : a) )

char *usage = "Viterbi Decoder Tool\nUsage: program transition-probs-file emission-probs-file";


#define MAX_STATES	200
#define MAX_LENGTH	2000

typedef struct _hmm
{	char tokens[MAX_STATES];
	char transi[10][2];
	int tokenToIdx['Z'+1];
	
	// match emissions probabilities (in log-odds)
	float* matchprobs[MAX_STATES];	// indexed by state
	// transition probabilities (in log-odds)
	float* transprobs[MAX_STATES];	// indexed by state
	// background probabilities (in log-odds)
	float* inserprobs;
	float* inserprobs0;	// for state 0. From HMMER. not used
	
	// probs indexed reordered by sequence element
	float* matchprofs[MAX_STATES];
	float* inserprofs;

	int ntokens;
	int ntrans;
	int nstates;
	
	float** resultsDel;
	float** resultsIns;
	float** resultsMat;
} hmm;

#define NEG_INFINITY (DBL_MIN/2)


int read_clean_line(FILE *fp, char* buf);
void print_model(hmm* model);
void process_args(int argc, char** argv, FILE** emissionsfp, FILE** transprobfp);
void print_results(float** results, char* seq, FILE* resfp, int nstates);
hmm* load_model(FILE *matfp);
void alloc_result_matrixes(hmm* model, int size);


// Convert back to output prob. do nothing
#define convertplog(plog) plog

static float logneginf; // log(NEG_INFINITY)

float convertprob(char* data)
{
	if (*data == '*')	return logneginf;
	else				return -atof(data);
}


void create_profiles(hmm* hhh, char* seq)
{
	float** matchprobs = hhh->matchprobs;
	float** matchprofs = hhh->matchprofs;
	float*  inserprobs = hhh->inserprobs;
	float*  inserprofs = hhh->inserprofs = calloc(MAX_LENGTH, sizeof(float));
	int k, i, seqlen = strlen(seq);
	int *tokens = hhh->tokenToIdx;

	for (i = 0; i < seqlen; i++)	// seq token
		inserprofs[i] = inserprobs[tokens[seq[i]]];

	for (k = 0; k < hhh->nstates; k++)	// seq token
	{
		matchprofs[k] = calloc(MAX_STATES, sizeof(float*));
		for (i = 0; i < seqlen; i++)	// seq token
			matchprofs[k][i] = matchprobs[k][tokens[seq[i]]];
	}
}


//	0		1		2		3		4		5		6
//	m->m	m->i	m->d	i->m	i->i	d->m	d->d

void viterbi_decode(hmm* model, char* seq)
{
	int idx, i, j, k, seqlen = strlen(seq);
	char* tokens = model->tokens;
	char (*transi)[2] = model->transi;
	float** matchprofs = model->matchprofs;
	float** transprobs = model->transprobs;
	float * inserprofs = model->inserprofs;
	float **resDel = model->resultsDel;
	float **resIns = model->resultsIns;
	float **resMat = model->resultsMat;
	
	float vDprev;
	float* vD	= calloc(seqlen+1, sizeof(float));
	float* vI	= calloc(seqlen+1, sizeof(float));
	float* vM1	= calloc(seqlen+1, sizeof(float));
	float* vM2	= calloc(seqlen+1, sizeof(float));
	float* transprobsK1;
	const float loginf = logneginf;

	// Initiliaze for state 0
	memset(vM2, 0, seqlen*sizeof(float));
	memset(vD , 0, seqlen*sizeof(float));

	vI[0] = 0.0; // log(1.0);
	transprobsK1 = transprobs[0];
	for (i = 1; i < strlen(seq); i++)	// seq token
		vI[i] = MAX2(vM2[i-1] + transprobsK1[1], vI [i-1] + transprobsK1[4]);
	vI[0] = loginf;

	if (DEBUG)
		for (i = 0; i < seqlen; i++)
		{	resDel[0][i] = vD [i];
			resMat[0][i] = vM2[i];
			resIns[0][i] = vI[i];
		}

	for (k = 1; k < model->nstates; k++)	// state
	{
		float* matchTemps = matchprofs[k];
		float* transprobsK = transprobs[k];
		transprobsK1 = transprobs[k-1];

		i = 0;
		vDprev = vD[i];
		
		for (i = 1; i < seqlen; i++)	// seq token
		{
			vM1[i] = matchTemps[i] - inserprofs[i]
						+ MAX3( vM2[i-1] + transprobsK1[0],	// Mj-1Mj
								vI [i-1] + transprobsK1[3],	// Ij-1Mj
								vD [i-1] + transprobsK1[5]);// Dj-1Mj
		}


		i = 0;
		vD[i] = MAX2(vM2[i] + transprobsK1[2], vD[i] + transprobsK1[6]);
		for (i = 1; i < seqlen; i++)	// seq token
			vD[i]  = MAX2(vM2[i] + transprobsK1[2],	vD[i] + transprobsK1[6]);

							//idx = model->tokenToIdx[seq[0]];
		vM1[0] = loginf;	// matchprobs[k][idx] - inserprobs[k][idx] + loginf;
		
		// Compute insert values
		vI[0] = loginf;
		for (i = 1; i < seqlen; i++)	// seq token
			vI[i] = MAX2(vM1[i-1] + transprobsK[1], vI[i-1] + transprobsK[4]);
		
	if (DEBUG)		// Debug: copy values to result matrix
		for (i = 0; i < seqlen; i++)
		{	resIns[k][i] = vI[i];
			resDel[k][i] = vD[i];
			resMat[k][i] = vM1[i];
		}

		// swap
		{ float* aa = vM1; vM1 = vM2; vM2 = aa; }
	}

	free(vI); free(vD); free(vM1); free(vM2);
}




int main(int argc, char** argv)
{
	int i, j;
	int numstates, numsymbols, seqlen=0;
	char *seq;
	struct timeb t1, t2;

	FILE* modelfp;
	outfp = stdout;
	logneginf = log(NEG_INFINITY);

	if (argc < 3)
		return 0;
	
	process_args(argc, argv, &modelfp, NULL);

	seq = argv[2];
	seqlen = strlen(seq);

	hmm* model = load_model(modelfp);
	create_profiles(model, seq);
//	print_model(model);

	if (DEBUG) alloc_result_matrixes(model, seqlen+3);

#define NROUNDS 20000
	ftime(&t1);

	// Reset to 0 the SSE arrays
//	memset(model->data, 0, model->totalBytes);
	for (i = 0; i < NROUNDS; i++)
		viterbi_decode(model, seq);

	ftime(&t2);

	printf("Took %f millisecs for %d Rounds\n", (t2.time - t1.time) +  (t2.millitm - t1.millitm)/1000.0, NROUNDS);
	
	if (DEBUG) {
		print_results(model->resultsMat, seq, fopen("results-matchx.txt", "w"), model->nstates);
		print_results(model->resultsIns, seq, fopen("results-insert.txt", "w"), model->nstates);
		print_results(model->resultsDel, seq, fopen("results-delete.txt", "w"), model->nstates);
	}

	return 0;
}



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


void print_model(hmm* model)
{
	int i, j;
	char* tokens = model->tokens;
	char (*transi)[2] = model->transi;
	float** matchprobs = model->matchprobs;
	float* inserprobs = model->inserprobs;
	float** transprobs = model->transprobs;
	FILE* mmoutfp = stdout;

	fprintf(mmoutfp, "HMM     ");
	for (i = 0; i < model->ntokens; i++)
		fprintf(mmoutfp, "     %c   ", tokens[i]);

	fprintf(mmoutfp, "\n       ");
	for (i = 0; i < model->ntrans; i++)
		fprintf(mmoutfp, "     %c->%c", transi[i][0], transi[i][1]);
	
	for (i = 0; matchprobs[i]; i++)
	{
		float* ppr;
		
		fprintf(mmoutfp, "\n   %4d ", i);
		ppr = matchprobs[i];
		for (j = 0; j < model->ntokens; j++)
			fprintf(mmoutfp, "  %.5f", ppr[j]);

		fprintf(mmoutfp, "\n        ");
		ppr = inserprobs;
		for (j = 0; j < model->ntokens; j++)
			fprintf(mmoutfp, "  %.5f", ppr[j]);

		fprintf(mmoutfp, "\n        ");
		ppr = transprobs[i];
		for (j = 0; j < model->ntrans; j++)
			fprintf(mmoutfp, "  %.5f", ppr[j]);
	}
	fprintf(mmoutfp, "\n");
}


hmm* load_model(FILE *matfp)
{
	hmm* hhh = calloc(1, sizeof(hmm));
	char buf[4000], *pp;
	char *delims = " \t\n\r";
	int* tokenToIdx = hhh->tokenToIdx;
	char* tokens = hhh->tokens;
	char (*transi)[2]  = hhh->transi;
	float** matchprobs = hhh->matchprobs;
	float** transprobs = hhh->transprobs;
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
		float* ppr;
	
		while ((j = read_clean_line(matfp, buf)) == 0);
		if (j < 0) break;
	
		// line with match probs
		ppr = matchprobs[i] = (float*) malloc((Ntokens+1)*sizeof(float));
		pp = strtok(buf, delims);
		pp = strtok(NULL, delims);	// ignore ID
		for (j = 0; j < Ntokens;  pp = strtok(NULL, delims))
			ppr[j++] = convertprob(pp);
		ppr[j] = 0;

		while (read_clean_line(matfp, buf) == 0);

		if (i < 2)
		{	// line with insert/background probs
			float* ppr = (float*) malloc((Ntokens+1)*sizeof(float));
			pp = strtok(buf, delims);
			for (j = 0; j < Ntokens;  pp = strtok(NULL, delims))
				ppr[j++] = convertprob(pp);
			ppr[j] = 0;

			if (!hhh->inserprobs0)	hhh->inserprobs0 = ppr;
			else					hhh->inserprobs = ppr;
		}

		// line with transition probs
		ppr = transprobs[i] = (float*) malloc((Ntrans+1)*sizeof(float));
		while (read_clean_line(matfp, buf) == 0);
		pp = strtok(buf, delims);
		for (j = 0; j < Ntrans ; pp = strtok(NULL, delims))
			ppr[j++] = convertprob(pp);
		ppr[j] = 0;
	}
	
	hhh->nstates = i-1;
	return hhh;
}


void print_results(float** results, char* seq, FILE* resfp, int nstates)
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

void alloc_result_matrixes(hmm* model, int size)
{
	int i;
	model->resultsDel = calloc(sizeof(float*), (model->nstates+1));
	model->resultsIns = calloc(sizeof(float*), (model->nstates+1));
	model->resultsMat = calloc(sizeof(float*), (model->nstates+1));
	
	for (i = 0; i <= model->nstates; i++)
	{	model->resultsDel[i]  = calloc(sizeof(float), size+1);
		model->resultsIns[i]  = calloc(sizeof(float), size+1);
		model->resultsMat[i]  = calloc(sizeof(float), size+1);
	}
}

