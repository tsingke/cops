#define _GNU_SOURCE
#include <xmmintrin.h>		/* SSE	*/
#include <emmintrin.h>		/* SSE2 */
#include <string.h>

#include "p7_config.h"
//#include "hmmer.h"
#include "base/p7_profile.h"

#include <pthread.h>
//#include <semaphore.h>
//#include <unistd.h>


typedef unsigned char byte;

void check_range(float vv)
{
	if (isinf(vv))
		return;
	if (vv > 4.0 || vv < -32.0)
		printf("RANGE EXCEEDED: %f\n", vv);
}

int viterbi_serial(ESL_DSQ *dsq, int L, P7_PROFILE *gm,	float *opt_sc, float* cmp)
{
	float *tsc= gm->tsc;
	int M = gm->M;
	int	i,k;
	
	float xmxN, xmxB, xmxE, xmxC, xmxJ;

	float dmx[M+2], imx[M+2], mmx[M+2];
	float loopJ = gm->xsc[p7P_J][p7P_LOOP];
	float loopE = gm->xsc[p7P_E][p7P_LOOP];
	float loopC = gm->xsc[p7P_C][p7P_LOOP];
	float loopN = gm->xsc[p7P_N][p7P_LOOP];
	float moveJ = gm->xsc[p7P_J][p7P_MOVE];
	float moveE = gm->xsc[p7P_E][p7P_MOVE];
	float moveC = gm->xsc[p7P_C][p7P_MOVE];
	float moveN = gm->xsc[p7P_N][p7P_MOVE];
	
	/* Initialization of the zero row.	*/
	xmxN = 0;
	xmxB = gm->xsc[p7P_N][p7P_MOVE];					/* S->N->B, no N-tail	*/
	xmxE = xmxJ = xmxC = -eslINFINITY;					/* need seq to get here */
		
	for (k = 0; k <= gm->M; k++)
		dmx[k] = mmx[k] = imx[k] = -eslINFINITY;		/* need seq to get here */
	
	/* DP recursion */
	for (i = 1; i <= L; i++) 
	{
		float *rsc = gm->rsc[dsq[i]]; // + p7P_NR;
		float sc, ipv, mpv, dpv, dcv;

		dcv = mpv = dpv = ipv = sc = -eslINFINITY;
		imx[0] = dmx[0] = mmx[0] = -eslINFINITY;
		xmxE = -eslINFINITY;
	
		for (k = 1; k <= M; k++)
		{
			dcv	= ESL_MAX(sc + TSC(p7P_MD,k-1), dcv + TSC(p7P_DD,k-1));
			
			/* match state. Sc == new MMX(k) */
			sc	= xmxB + TSC(p7P_LM,k-1);
			sc	= ESL_MAX(sc, mpv + TSC(p7P_MM,k-1));
			sc	= ESL_MAX(sc, ipv + TSC(p7P_IM,k-1));
			sc	= ESL_MAX(sc, dpv + TSC(p7P_DM,k-1));
			sc	= sc + MSC(k); // rsc[p7P_M];	
			xmxE= ESL_MAX(xmxE, sc);	// + esc);

			dpv = dmx[k];
			ipv = imx[k];
			mpv = mmx[k];
			mmx[k] = sc;
			dmx[k] = dcv;

			imx[k]	= ESL_MAX(mpv + TSC(p7P_MI,k), ipv + TSC(p7P_II,k)) + ISC(k);
		}

		/* E state update; transition from M_M scores 0 by def'n */
		xmxE	= ESL_MAX(xmxE, dmx[M]);

		/* Now the special states. E must already be done, and B must follow N,J.
		remember, N, C and J emissions are zero score by definition. */
		xmxJ	= ESL_MAX(xmxJ + loopJ, xmxE + loopE);	/* E->J is E's "loop" */
		xmxC	= ESL_MAX(xmxC + loopC, xmxE + moveE);
		xmxN	= xmxN			 + loopN;
		xmxB	= ESL_MAX(xmxN + moveN, xmxJ + moveJ);	/* J->B is J's move */
		
		if (cmp != NULL) cmp[i] = xmxB;
	}
	
	/* T state (not stored) */
	if (opt_sc != NULL)
		opt_sc[0] = xmxC + moveC;

	return eslOK;
}



// ======================================================================================
//	VITERBI	ROGNES	8-WORD	PARALLEL	CHANNELS
// ======================================================================================

//VER	L1D_PREFETCH L1D_REPL lines in

#define unit	__m128
#define uniti	__m128i

typedef struct _data_interleaved
{	// data TSC
	uniti* tsc_all;
	int Mlinesize;
	int16_t** rsc_isc;
	int16_t** rsc_msc;

	uniti* xmxN;
	uniti* xmxE;
	uniti* psc;
	uniti* pdcv;
	uniti* pmpv;
	uniti* pipv;
	uniti* pdpv;

	uniti moveE;
	uniti loopC;
	uniti moveC;
} dpopt;

#define WORDMAX 32767


#define SETV4(v1,val)	\
{	v1 = _mm_load_ss(&val);		\
	v1 = (__m128) _mm_shuffle_epi32((__m128i)v1, 0);\
}
#define SETV8(v1,val)	v1 = _mm_set1_epi16(val);

#define SET8MOV(d,val)	\
{	int16_t vv = val;	\
	SETV8(*d,vv); d++;	\
}


#define PAGESIZE 4096	//0x1000

// gcc default: 32 bytes alignment

float* aligned_allodffdc(int nfloats)
{
	float *data = malloc((nfloats+4)*sizeof(float));
	size_t idata = (size_t) data;
	if (idata % 16 == 0) return data;
	else	return (float*) ((idata & ~0xf) + 0x10);
}

int16_t* aligned_alloc64(int nwords)
{
	int16_t *data = malloc((nwords+32)*sizeof(int16_t));
	size_t idata = (size_t) data;
	if (idata % 64 == 0) return data;
	else	return (int16_t*) ((idata & ~0x3f) + 0x40);
}



void* align_page_boundary()
{
	void	*data = malloc(1);
	size_t idata	= (size_t) data;
	size_t offset = PAGESIZE - ((idata+32) % PAGESIZE);
	size_t size;

	if (offset < 32)
		size = PAGESIZE;
	else if (offset%32 == 0)
		size = offset-16;
	else
		size = offset - (offset %32) +1;
	
	void *data1 = malloc(size);
//	void *data2 = malloc(1);
//	printf("data: %p | data+48: %p | offset: %d | size to alloc: %d | data1: %p | data2: %p\n",
//			data, idata+48, offset, size, data1, data2);
	return data1;
}

int roundtop(int v, int align)
{
	if (v%align == 0)	return v;
	else			return v + (align - v%align);
}

float undiscretize(int16_t ii)
{
	float ff = (float) (ii);
	return (ff)/1000;
}

int16_t discretize(float ff)
{
	if (isinf(ff))	return -WORDMAX;
//	if (!isfinite(ff))	exit(fprintf(stderr, "NoTA NORMAL FF %f\n", ff));
//	if (ff < 4.0 && ff > -32.0)
		return (int16_t) (roundf (ff * 1000));
//	exit(fprintf(stderr, "Discretize overflow!!! %f\n", ff));
}


dpopt* create_data_structs(P7_PROFILE *gm, int L, ESL_DSQ **ddsq)
{
	float *tsc= gm->tsc;
	uniti *tt;
	int M = gm->M, Mtop, k, i;
//	align_page_boundary();
	dpopt* dr = calloc(1, sizeof(dpopt));

	ESL_DSQ* edata = (ESL_DSQ*) aligned_alloc64(8*(sizeof(ESL_DSQ)*(L+2)));
	ddsq[0] = edata;
	ddsq[1] = edata + 1*(L+2);
	ddsq[2] = edata + 2*(L+2);
	ddsq[3] = edata + 3*(L+2);
	ddsq[4] = edata + 4*(L+2);
	ddsq[5] = edata + 5*(L+2);
	ddsq[6] = edata + 6*(L+2);
	ddsq[7] = edata + 7*(L+2);

	Mtop = dr->Mlinesize = roundtop(M,8); //M + ((M%64 == 0)? 0 : (64 - (M%64)));
//	printf("MMM %d\n", dr->Mlinesize);
	
	tt = (uniti*) aligned_alloc64(8*Mtop*8);

	dr->tsc_all = tt;
	//	#define TSC(s,k) (tsc[(k)p7P_NTRANS + (s)])
	for (k = 0; k < M; k++, tsc += p7P_NTRANS)
	{
		SET8MOV(tt, discretize(tsc[p7P_MD]));
		SET8MOV(tt, discretize(tsc[p7P_DD]));
		SET8MOV(tt, discretize(tsc[p7P_LM])); 
		SET8MOV(tt, discretize(tsc[p7P_MM]));
		SET8MOV(tt, discretize(tsc[p7P_IM]));
		SET8MOV(tt, discretize(tsc[p7P_DM]));
		SET8MOV(tt, discretize(tsc[p7P_NTRANS+p7P_MI]));
		SET8MOV(tt, discretize(tsc[p7P_NTRANS+p7P_II]));
	}

	// fill the left-over places with neutral values
	int16_t neginf = -WORDMAX;
	for (	; k < Mtop; k++, tsc += p7P_NTRANS)
	{
		SET8MOV(tt, neginf);
		SET8MOV(tt, neginf);
		SET8MOV(tt, neginf);
		SET8MOV(tt, neginf);
		SET8MOV(tt, neginf);
		SET8MOV(tt, neginf);
		SET8MOV(tt, neginf);
		SET8MOV(tt, neginf);
	}

	dr->rsc_msc = calloc(21, sizeof(float*));
	int16_t *data = (int16_t*) aligned_alloc64((Mtop)*20);

	for (k = 0; k < 20; k++)
	{
		int16_t *rmsc = dr->rsc_msc[k] = data + (k)*(Mtop);
		float *rsc = gm->rsc[k] + p7P_NR;	// ignore 1st values
		// #define ISC(k)	 (rsc[(k)p7P_NR + p7P_ISC])
		// #define MSC(k)	 (rsc[(k)p7P_NR + p7P_MSC])
		for (i = 0; i< M; i++, rsc += p7P_NR)
			*rmsc++ = discretize(rsc[p7P_M]);

		// fill the left-over places with neutral values
		for (	; i< Mtop; i++)
			*rmsc++ = -WORDMAX;
	}
	return dr;
}


// =====>	Emission scores ISC são agora sempre 0!!!!	<=====

#define NUCLEUS(k,j)	\
{	sc	= _mm_adds_epi16(sc, *tt); tt++;	\
	dcv = _mm_max_epi16(sc, _mm_adds_epi16(dcv, *tt));tt++;	\
	\
	sc = _mm_adds_epi16(xmxB, *tt); tt++;	\
	sc = _mm_max_epi16(sc, _mm_adds_epi16(mpv, *tt)); tt++;	\
	sc = _mm_max_epi16(sc, _mm_adds_epi16(ipv, *tt)); tt++;	\
	sc = _mm_max_epi16(sc, _mm_adds_epi16(dpv, *tt)); tt++;	\
	sc = _mm_adds_epi16(sc,  xmm[j]);		\
	xmxE = _mm_max_epi16(xmxE, sc);			\
	\
	dpv = dmx[k];	\
	ipv = imx[k];	\
	mpv = mmx[k];	\
	mmx[k] = sc;	\
	dmx[k] = dcv;	\
	\
	temp = _mm_adds_epi16(mpv, *tt); tt++;	\
	imx[k] = _mm_max_epi16(temp, _mm_adds_epi16(ipv, *tt)); tt++;\
	k++;			\
}

#define MIX(i,r,range,var)	\
	var[r  ] = _mm_unpacklo_epi##range(var[i], var[i+2]);	\
	var[r+1] = _mm_unpackhi_epi##range(var[i], var[i+2]);

#define MMLOAD(a,b)		\
	xmm[a] = _mm_unpacklo_epi16(*mscore##a, *mscore##b);\
	xmm[b] = _mm_unpackhi_epi16(*mscore##a, *mscore##b);


void pre_compute_sse(P7_PROFILE* gm, int L, dpopt* dd)
{
	uniti* tt = dd->xmxN = (uniti*) aligned_alloc64((L+1)*8);
	uniti loopN, moveN, temp;
	int16_t zero = 0.0; int k;
	dd->xmxE = (uniti*) aligned_alloc64((L+1)*8);

	SETV8(moveN, discretize(gm->xsc[p7P_N][p7P_MOVE]));
	SETV8(loopN, discretize(gm->xsc[p7P_N][p7P_LOOP]));
	SETV8(temp, zero);

	tt[0] = moveN;
	for (k = 1; k < L; k++)
	{	temp = _mm_adds_epi16(temp, loopN);
		tt[k]= _mm_adds_epi16(temp, moveN);
	}

	dd->psc  = (uniti*) aligned_alloc64((L)*8);
	dd->pdcv = (uniti*) aligned_alloc64((L)*8);
	dd->pmpv = (uniti*) aligned_alloc64((L)*8);
	dd->pipv = (uniti*) aligned_alloc64((L)*8);
	dd->pdpv = (uniti*) aligned_alloc64((L)*8);

	SETV8(dd->moveE, discretize(gm->xsc[p7P_E][p7P_MOVE]));
	SETV8(dd->loopC, discretize(gm->xsc[p7P_C][p7P_LOOP]));
	SETV8(dd->moveC, discretize(gm->xsc[p7P_C][p7P_MOVE]));
}


int ccsignals = 0, ccwaits = 0;
pthread_barrier_t barrsynch;

volatile byte *synchflags[1000];

typedef enum _ssynch { START, STOP } ssynch;

volatile ssynch synchcontrol;
#define SYNCH	0
#define DISTURB	0


#define TIME_WAIT() { if (rand() % 10 == 0) usleep(100*(rand()%10)); }
#define YIELD() { if (rand() % 10 == 0) usleep(100*(rand()%10)); }

#define MAX_PARTITION 112	// 112, 96 e 88 parecem ser os melhores
int PARTITION;
int NPARTITIONS;
//#define NTHREADS 1

#define tprintf if(0) printf
#define SSTEST 0

#define dprintf if(0)printf


int viterbi_parallel_channels(	ESL_DSQ** ddsq, int L, P7_PROFILE* gm,
								dpopt* opdata, float* opt_sc, int thrid)
{
	int M = gm->M, i, k, v, t;
	uniti** oprmsc = (uniti**) opdata->rsc_msc;
	uniti* xmxNv = opdata->xmxN;
	uniti* xmxEv = opdata->xmxE;
	uniti xmxB, xmxE, xmxC, loopC, moveE, moveC, Vinf;
	int16_t neginf = -WORDMAX;

#define ARRAY_ALIGNMMENT  
// __attribute__ ((aligned (64)));
	uniti dmx[PARTITION] ARRAY_ALIGNMMENT;
	uniti mmx[PARTITION] ARRAY_ALIGNMMENT;
	uniti imx[PARTITION] ARRAY_ALIGNMMENT;
	uniti xmm[16] ARRAY_ALIGNMMENT;
	uniti* lastdmx = dmx + ((M%PARTITION == 0) ? PARTITION-1 : M%PARTITION-1);

	// Thread 0 (main) must have the last block!
	t = ((NPARTITIONS+thrid)%NTHREADS)*PARTITION;

	tprintf("START viterbi Thr %d in %d\n", thrid, t);

	moveE = opdata->moveE;
	loopC = opdata->loopC;
	moveC = opdata->moveC;
	SETV8(Vinf, neginf);
	xmxC = Vinf;

	if (NTHREADS > 1)
		pthread_barrier_wait(&barrsynch);

	for (	; t < M; t += NTHREADS*PARTITION)
	{
		volatile byte* synchflags1 = synchflags[t/PARTITION-1];
		volatile byte* synchflags2 = synchflags[t/PARTITION];
		int t8 = t/8;

		for (k = 0; k < PARTITION; k++)
			dmx[k] = mmx[k] = imx[k] = Vinf;

		for (i = 1; i <= L; i++)
		{
			uniti sc, dcv, temp, mpv, ipv, dpv;
			uniti *restrict tt = (uniti*) opdata->tsc_all + t*8;
			v = i-1;
			xmxB = xmxNv[v];

			if (t == 0)
				xmxE = mpv = dpv = ipv = sc = dcv = Vinf;
			else {
				if (NTHREADS > 1)
					 while (!synchflags1[v]) sched_yield(); 
				if (SYNCH) __sync_synchronize();
			
			//	printf("Synchflags1: %p\n", synchflags1);
				xmxE = xmxEv[v];
				sc	= opdata->psc [v];
				dcv = opdata->pdcv[v];
				mpv = opdata->pmpv[v];
				ipv = opdata->pipv[v];
				dpv = opdata->pdpv[v];
			}
		
			uniti* mscore0 = oprmsc[ddsq[0][i]] + t8;
			uniti* mscore1 = oprmsc[ddsq[1][i]] + t8;
			uniti* mscore2 = oprmsc[ddsq[2][i]] + t8;
			uniti* mscore3 = oprmsc[ddsq[3][i]] + t8;
			uniti* mscore4 = oprmsc[ddsq[4][i]] + t8;
			uniti* mscore5 = oprmsc[ddsq[5][i]] + t8;
			uniti* mscore6 = oprmsc[ddsq[6][i]] + t8;
			uniti* mscore7 = oprmsc[ddsq[7][i]] + t8;

			for (k = 0; k < PARTITION && t+k < M; )
			{
				MMLOAD(0,1)			MMLOAD(2,3)
				MMLOAD(4,5)			MMLOAD(6,7)
				MIX(0, 8,32,xmm)	MIX(4,10,32,xmm)
				MIX(1,12,32,xmm)	MIX(5,14,32,xmm)
				MIX( 8,0,64,xmm)	MIX( 9,2,64,xmm)
				MIX(12,4,64,xmm)	MIX(13,6,64,xmm)

				NUCLEUS(k,0)	NUCLEUS(k,1)
				NUCLEUS(k,2)	NUCLEUS(k,3)
				NUCLEUS(k,4)	NUCLEUS(k,5)
				NUCLEUS(k,6)	NUCLEUS(k,7)

				mscore0++; mscore1++; mscore2++; mscore3++;
				mscore4++; mscore5++; mscore6++; mscore7++;
			}

			if (t+k < M)
			{	v = i-1;
				xmxEv[v] = xmxE;
				opdata->psc [v] =  sc;
				opdata->pdcv[v] = dcv;
				opdata->pmpv[v] = mpv;
				opdata->pipv[v] = ipv;
				opdata->pdpv[v] = dpv;

				if (SYNCH) __sync_synchronize();
				if (NTHREADS > 1) synchflags2[v] = 1;
			}
			else
			{
				xmxE = _mm_adds_epi16(moveE, _mm_max_epi16(xmxE, *lastdmx)); 
				xmxC = _mm_max_epi16 ( xmxE, _mm_adds_epi16(xmxC, loopC));
			}
		}
	}

	xmxC = _mm_adds_epi16(xmxC, moveC);

	if (opt_sc != NULL && SSTEST)
	{
		tprintf("COPY Result Thr %d - t %d \n", thrid, t);
		int16_t tt[8];
		memmove(tt, &xmxC, sizeof(xmxC));
		for (i = 0; i < 8; i++)
			opt_sc[i] = undiscretize(tt[i]);
	}

	tprintf("END viterbi Thr %d - t %d\n", thrid, t);

//	synchcontrol = STOP;
	return eslOK;
}


typedef struct _thr_args
{	int thrid;
	ESL_DSQ	**ddsq;
	int L;
	int N;
	P7_PROFILE	*gm;
	dpopt		*opdata;
} pthr_args;

void* viterbi_thr_loop(void* argst)
{
	pthr_args* args = (pthr_args*) argst;
	int i, N = args->N;

	for (i = 0; i < N; i++)
	{
		tprintf("wAIT worker thread %d\n", args->thrid);
		// Wait for the main thread to finish and reset flags
//		while (synchcontrol == STOP) sched_yield();

		tprintf("AWOKE worker thread %d\n", args->thrid);
		viterbi_parallel_channels(args->ddsq, args->L, args->gm, args->opdata, NULL, args->thrid);			
	}

	tprintf("FINISHED worker thread %d\n", args->thrid);
	return (void*) 0;
}
 

/*****************************************************************
2. Benchmark driver.
 *****************************************************************/
//#ifdef VITERBI_BENCHMARK
/*
gcc -g -Wall -O3 -std=gnu99 -fstrict-aliasing -Wstrict-aliasing -o benchmark -I. -L. -I../lib/easel -L../lib/easel -DVITERBI_BENCHMARK viterbi-parallel-16bit.c -lhmmer -leasel -lm -lpthread
./benchmark-generic-viterbi <hmmfile>
*/


#include "easel.h"
#include "esl_randomseq.h"
#include "esl_stopwatch.h"
#include "base/p7_hmmfile.h"
#include "search/modelconfig.h"

#include <sys/timeb.h>
#include <sched.h> 

static ESL_OPTIONS options[] = {
	/* name			type		default	env	range toggles reqs incomp	help						docgroup*/
	{ "-h",		eslARG_NONE	, FALSE	, NULL, NULL , NULL, NULL, NULL, "show brief help on version and usage",0 },
	{ "-s",		eslARG_INT	, "42"	, NULL, NULL , NULL, NULL, NULL, "set random number seed to <n>",		0 },
	{ "-L",		eslARG_INT	, "400"	, NULL, "n>0", NULL, NULL, NULL, "length of random target seqs",		0 },
	{ "-N",		eslARG_INT	, "50000",NULL, "n>0", NULL, NULL, NULL, "number of random target seqs",		0 },
	{ "-M",		eslARG_INT	, "1024", NULL, NULL , NULL, NULL, NULL, "Max Partition", 0 },
	{	0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage [] = "[-options] <hmmfile>";
static char banner[] = "benchmark driver for generic Viterbi";

int main(int argc, char **argv)
{
	ESL_GETOPTS	*go		= esl_getopts_CreateDefaultApp(options, 1, argc, argv, banner, usage);
	ESL_STOPWATCH *w	= esl_stopwatch_Create();
	ESL_RANDOMNESS *r	= esl_randomness_CreateFast(esl_opt_GetInteger(go, "-s"));
	ESL_ALPHABET*abc	= NULL;
	P7_HMMFILE	*hfp	= NULL;
	P7_HMM		*hmm	= NULL;
	P7_BG		*bg		= NULL;
	P7_PROFILE	*gm		= NULL;
	char		*hmmfile= esl_opt_GetArg(go, 1);
	int			L		= 2000;	// esl_opt_GetInteger(go, "-L");
	int			N		= 2000;	// esl_opt_GetInteger(go, "-N");
	ESL_DSQ		**ddsq	= calloc(8, sizeof(ESL_DSQ*));
	unit		resdata[10];
	float		*sc1 	= (float*)	(resdata+0);	// uses 2 units
	int			i, j;
	dpopt		*opdata;
	struct timeb tbstart, tbend;

	if (p7_hmmfile_OpenE(hmmfile, NULL, &hfp, NULL)	!= eslOK) p7_Fail("Failed to open HMM file %s", hmmfile);
	if (p7_hmmfile_Read(hfp, &abc, &hmm)	 		!= eslOK) p7_Fail("Failed to read HMM");
//	srand(time(NULL));

	bg = p7_bg_Create(abc);
	p7_bg_SetLength(bg, L);

	gm = p7_profile_Create(hmm->M, abc);
//	p7_ProfileConfig		(hmm, bg, gm, L, p7_LOCAL); // podia ser GLOCAL, GLOBAL, e UNI-qq para só 1 hit
	p7_profile_ConfigUnilocal(gm, hmm, bg, L);
//	MAX_PARTITION = esl_opt_GetInteger(go, "-M");

	opdata = create_data_structs(gm, L, ddsq);
	pre_compute_sse(gm, L, opdata);
	
	int M = gm->M;
	int nblocks = NTHREADS;
	while (M / nblocks > MAX_PARTITION)
		nblocks += NTHREADS;
	PARTITION	= roundtop(M/nblocks, 8);
	NPARTITIONS = M/PARTITION + (M%PARTITION != 0);

	int numflags = roundtop(L+1, 8);
	for (i = 0; i < NPARTITIONS; i++)
		synchflags[i] = (byte*) aligned_alloc64(numflags);
	synchflags[NPARTITIONS] = NULL;

	printf("M: %d | Q: %d | L: %d | N: %d | NThrs: %d | Part: %d | NParts: %d\n", 
			gm->M, (int) ceil(gm->M/4.0), L, N, NTHREADS, PARTITION, NPARTITIONS);

	ESL_DSQ	**seqsdb = calloc(8*(N+1), sizeof(ESL_DSQ*));
	for (i = 0; i < N*8; i++)
	{	seqsdb[i] = malloc((L+2)*sizeof(ESL_DSQ));
		esl_rsq_xfIID(r, bg->f, abc->K, L, seqsdb[i]);
	}

	/* Get a baseline time: how long it takes just to generate the sequences */
/*	ftime(&tbstart);
	for (i = 0; i < N; i++)
		esl_rsq_xfIID(r, bg->f, abc->K, L, ddsq[0]);
	ftime(&tbend);
	basesecs = ((tbend.time - tbstart.time) + (tbend.millitm - tbstart.millitm)*0.001)*8;
	dprintf("Basetime: %f\n", basesecs);
*/

	if(NTHREADS > 1)
	{
		pthread_t threads[NTHREADS];
		pthr_args args;
		args.ddsq = ddsq; args.L = L; args.N = N;
		args.gm = gm; args.opdata = opdata; args.thrid = 0;

		pthread_barrier_init(&barrsynch, NULL, NTHREADS);
	//	synchcontrol = STOP; // to prevent worker thread from starting immediately
		__sync_synchronize();

		int mainthrid = NTHREADS-1;
		cpu_set_t cpuset;
		CPU_ZERO(&cpuset);
		CPU_SET(mainthrid, &cpuset);
		threads[mainthrid] = pthread_self();
		pthread_setaffinity_np(threads[mainthrid], sizeof(cpu_set_t), &cpuset);
	
		for (i = 0; i < NTHREADS-1; i++)
		{
			pthr_args *argscopy = calloc(1, sizeof(pthr_args));
			memcpy(argscopy, &args, sizeof(pthr_args));
			argscopy->thrid	= i;
		
			CPU_ZERO(&cpuset);
			int htcore = (i%2!=0)*4 + (i/2);
			printf("SET HTthread %d in logical core %d\n", i, htcore);
			CPU_SET(htcore, &cpuset);

			pthread_attr_t attr;
			pthread_attr_init(&attr);
			pthread_attr_setaffinity_np(&attr, sizeof(cpu_set_t), &cpuset);

			if (pthread_create(&threads[i], &attr, viterbi_thr_loop, argscopy))
				return fprintf(stderr, "ERROR could not create worker thread\n");
	
		// Alternative method: change affinity after creation
		// pthread_setaffinity_np(threads[i], sizeof(cpu_set_t), &cpuset);
		}

		for (i = 0; i < NTHREADS; i++)
		{
			pthread_getaffinity_np(threads[i], sizeof(cpu_set_t), &cpuset);	
			for (j = 0; j < 8; j++)
				if (CPU_ISSET(j, &cpuset))
					break;
			printf("THR %d running on cpu %d\n", i, j);
		}
	}

	ftime(&tbstart);
	for (i = 0; i < N; i++)
	{
		for (j = 0; j < 8; j++)
		//	esl_rsq_xfIID(r, bg->f, abc->K, L, ddsq[j]);
		//	memcpy(ddsq[j], seqsdb[i+j], L*sizeof(ESL_DSQ));
			ddsq[j] = seqsdb[i+j];

		for (j = 0; j < gm->M/PARTITION; j++)
			memset((void*) synchflags[j], 0, numflags);

		viterbi_parallel_channels(ddsq, L, gm, opdata, sc1, NTHREADS-1);
		
#if SSTEST	// 51segs
	//	fscanf(results, "%f %f %f %f\n", &sc1[0], &sc1[1], &sc1[2], &sc1[3]);
		float *sc2 = (float*) (resdata+2);

		for (j = 0; j < 8; j++)
		{
			viterbi_serial (ddsq[j], L, gm, sc2+j, NULL);
			if (fabs(sc1[j] - sc2[j]) > 0.04)
			{	fprintf(stderr, "WRONG entry %d-%d: %f %f\n", i, j, sc1[j], sc2[j]); getc(stdin); }
		}
#endif
	}
	ftime(&tbend);
	
	double secs = (tbend.time - tbstart.time) + (tbend.millitm - tbstart.millitm)*0.001;
	w->elapsed = w->user = secs;
	esl_stopwatch_Display(stdout, w, "# Opt CPU time: ");
	printf("# M %d: %.1f Mc/s\n", gm->M, N * L * (double) gm->M * 1e-6 / (double) secs);

//	printf("Secs: %f | Millis: %f\n", secs, secs*1000);
	return 0;
}
//#endif 


