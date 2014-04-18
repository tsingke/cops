#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<time.h>
#include<ctype.h>
#include<float.h>
#include<math.h>
#include<emmintrin.h>

#include<sys/timeb.h>


typedef unsigned int	uint;
typedef unsigned char	uchar;
typedef unsigned char	bool;
typedef unsigned char	byte;

typedef float TYPE;

FILE *outfp;

#define DEBUG 0

#define MAX_STATES	200

#define MAX_LINE		5000
#define MAX_MAT_SIZE	MAX_STATES
#define MAX2(a,b)	( (b > a)	? b : a)
#define MAX3(a,b,c)	( (c > b)	? ((c > a)? c : a)	\
								: ((b > a)? b : a) )
char *usage = "Viterbi Decoder Tool\nUsage: program transition-probs-file emission-probs-file";

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

	int profileLength;
	int totalBytes;
	byte* data;

	float *vD;
	float *vI;
	float *vM1;
	float *vM2;
	
	float** resultsDel;
	float** resultsIns;
	float** resultsMat;
} hmm;


int read_clean_line(FILE *fp, char* buf);
void print_model(hmm* model, FILE* dfg);
void process_args(int argc, char** argv, FILE** emissionsfp, FILE** transprobfp);
void print_results(float** results, char* seq, FILE* resfp, int nstates);
hmm* load_model(FILE *matfp);
void alloc_result_matrixes(hmm* model, int size);


#define NEG_INFINITY (DBL_MIN/2)
static float logneginf; // log(NEG_INFINITY)

// Convert back to output prob. do nothing
#define convertplog(plog) plog

float convertprob(char* data)
{
	if (*data == '*')	return logneginf;
	else				return -atof(data);
}


size_t align_data(size_t data, int ALIGNMENT)
{
	return ((data + ALIGNMENT-1) & ~(ALIGNMENT-1));
}

byte* alloc_aligned(int size, int alignment)
{
	byte* data = malloc(size+5*alignment+1);
	memset(data, 0, size+5*alignment+1);
	return (byte*) align_data((size_t) data, alignment);
}


void create_profiles(hmm* hhh, char* oseq)
{
	int k, i;
	const int seqlen = strlen(oseq);
	const int segSize = 4;
	const int numSegs = (seqlen+3)/segSize;
	const int totalSize = hhh->profileLength = numSegs*segSize;
	char* seq = calloc(totalSize+1, sizeof(char));

	int *tokens = hhh->tokenToIdx;
	float** matchprobs = hhh->matchprobs;
	float** matchprofs = hhh->matchprofs;
	float*  inserprobs = hhh->inserprobs;	// misalign the arrays by 1 float
	float*  inserprofs = hhh->inserprofs = (float*) (12 + alloc_aligned((totalSize+4)*sizeof(float), 16));
	
	// Alloc the SSE arrays
	const int memBytes = hhh->totalBytes = 4*(totalSize+16)*sizeof(float);
	hhh->data = calloc(memBytes, 1);
	hhh->vD	= (float*) align_data((size_t) hhh->data, 16);
	hhh->vI	= (float*) align_data((size_t) (hhh->vD + totalSize+1), 16);
	hhh->vM1= (float*) align_data((size_t) (hhh->vI + totalSize+1), 16);
	hhh->vM2= (float*) align_data((size_t) (hhh->vM1+ totalSize+1), 16);

	// Pad the sequence
	strcpy(seq, oseq);
	for (i = seqlen; i < totalSize; i++)
		seq[i] = ' ';
	
	// dummy prob for matches and inserts
	tokens[' '] = 0.0;
	
	// Create insert probabilities profile
	for (i = 0; i < totalSize; i++)	// seq token
		inserprofs[i] = inserprobs[tokens[seq[i]]];

	// Create match probabilities profile
	for (k = 0; k < hhh->nstates; k++)	// seq token
	{
		matchprofs[k] = (float*) (12 + alloc_aligned((totalSize+4)*sizeof(float), 16));
		for (i = 0; i < totalSize; i++)	// seq token
			matchprofs[k][i] = matchprobs[k][tokens[seq[i]]];
	}
}


// SSE operations

#define EXPANDV4F(v1,val)		\
	v1 = _mm_load_ss(&val);		\
	v1 = (__m128) _mm_shuffle_epi32((__m128i)v1, 0);

#define SUBV4F(v1,v2,v3)	v1 = _mm_sub_ps(v2, v3)
#define ADDV4F(v1,v2,v3)	v1 = _mm_add_ps(v2, v3)
#define MAXV4F(v1,v2,v3)	v1 = _mm_max_ps(v2, v3)
#define MOVEV4F(v1,v2)		v1 = v2;
#define unit __m128

static int ccc = 0;
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
	float **resDel = 0, **resIns = 0, **resMat = 0;

	const int segSize = 4;
	const int numSegs = (seqlen+3)/segSize;
	const int totalSize = numSegs*segSize;

	float *vD, *vI, *vM1, *vM2, *transprobsK1;
	unit* vDt, *vIt, *vM1t, *vM2t, *vMatch, *vInser;	// iterators
	unit t1, t2, t3, t4;	// auxiliary
	unit tp0, tp1, tp2, tp3, tp4, tp5, tp6;	// transition probabilities
	const float loginf = logneginf;

	vD	= model->vD;	vI	= model->vI;
	vM1	= model->vM1;	vM2	= model->vM2;

	// Initiliaze for state 0
	memset(vM2, 0, totalSize*sizeof(float));
	memset(vD , 0, totalSize*sizeof(float));

	// Not vectorizabble
	vI[0] = 0.0; // log(1.0);
	transprobsK1 = transprobs[0];
	for (i = 0; i < totalSize; i++)	// seq token
		vI[i+1] = MAX2(vM2[i] + transprobsK1[1], vI [i] + transprobsK1[4]);
	vI [0] = loginf;

	if (DEBUG)
	{	alloc_result_matrixes(model, totalSize);
		resDel = model->resultsDel;	resIns = model->resultsIns;	resMat = model->resultsMat;

		for (i = 0; i < totalSize; i++)
		{	resDel[0][i] = vD [i];
			resMat[0][i] = vM2[i];
			resIns[0][i] = vI[i];
		}
	}

	for (k = 1; k < model->nstates; k++)	// state
	{
		float* transprobsK = transprobs[k];
		transprobsK1 = transprobs[k-1];

		EXPANDV4F(tp0,transprobsK1[0]);
		EXPANDV4F(tp1,transprobsK1[1]);
		EXPANDV4F(tp2,transprobsK1[2]);
		EXPANDV4F(tp3,transprobsK1[3]);
		EXPANDV4F(tp4,transprobsK1[4]);
		EXPANDV4F(tp5,transprobsK1[5]);
		EXPANDV4F(tp6,transprobsK1[6]);

		vM2t = (unit*) vM2; vIt = (unit*) vI; vDt = (unit*) vD;
		vM1t = (unit*) (vM1);

		// matchprofs and inserprofs have been misaligned by 1 float
		vMatch = (unit*) (matchprofs[k]+1);
		vInser = (unit*) (inserprofs+1);

		for (i = 0; i < numSegs; i++)	// seq token
		{
			ADDV4F(t1, *vM2t, tp0);
			ADDV4F(t2, *vIt , tp3);
			ADDV4F(t3, *vDt , tp5);

			MAXV4F(t1, t1, t2);
			MAXV4F(t1, t1, t3);

			SUBV4F(t4, *vMatch, *vInser);
			ADDV4F(*vM1t, t4, t1);

			vMatch++; vInser++;
			vM2t++; vDt++; vM1t++; vIt++;
		}

		vM2t = (unit*) vM2; vDt = (unit*) vD;
		for (i = 0; i < numSegs; i++)	// seq token
		{
			ADDV4F(t1, *vM2t, tp2);
			ADDV4F(t2, *vDt , tp6);
			MAXV4F(*vDt, t1, t2);
			vM2t++; vDt++;
		}
		
		// Shift vM1
		memmove(vM1+1, vM1, totalSize*sizeof(float));

		// Compute insert values. Not vectorizable
		vM1[0] = loginf;
		for (j = 0; j < totalSize-1; j++)	// seq token
			vI[j+1] = MAX2(vM1[j] + transprobsK[1], vI[j] + transprobsK[4]);
		
		if (DEBUG)		// Debug: copy values to result matrix
			for (i = 0; i < totalSize; i++)
			{	resIns[k][i] = vI[i];
				resDel[k][i] = vD[i];
				resMat[k][i] = vM1[i];
			}

		// swap
		{ float* aa = vM1; vM1 = vM2; vM2 = aa; }
	}
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
//	print_model(model, fopen("model1.txt", "w"));

#define NROUNDS 20000
	ftime(&t1);

	// Reset to 0 the SSE arrays
//	memset(model->data, 0, model->totalBytes);
	for (i = 0; i < NROUNDS; i++)
		viterbi_decode(model, seq);

	ftime(&t2);

	printf("Took %f millisecs for %d Rounds\n", (t2.time - t1.time) +  (t2.millitm - t1.millitm)/1000.0, NROUNDS);

	if (DEBUG)
	{	print_results(model->resultsMat, seq, fopen("results-matchx.txt", "w"), model->nstates);
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


void print_model(hmm* model, FILE* mmoutfp)
{
	int i, j;
	char* tokens = model->tokens;
	char (*transi)[2] = model->transi;
	float** matchprobs = model->matchprobs;
	float* inserprobs = model->inserprobs;
	float** transprobs = model->transprobs;

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
	{	fprintf(resfp, "\n%3d  ", k);
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

