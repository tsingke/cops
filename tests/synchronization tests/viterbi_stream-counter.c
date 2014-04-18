#define _GNU_SOURCE
#include <xmmintrin.h>		/* SSE	*/
#include <emmintrin.h>		/* SSE2 */

#include "hmmer.h"
#include "p7_config.h"

#include <pthread.h>
#include <semaphore.h>
#include <string.h>

typedef unsigned char uchar;
typedef unsigned char byte;

extern int p7_ViterbiFilter(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_OMX *ox, float *ret_sc);


#include <viterbi_serial.h>


/***************************************************************************
 *	Viterbi Stream with 8-channel integers and threaded dstream->partitions
 ***************************************************************************/

// comment to use volatile
#define volatile

typedef struct _DATA_STREAM
{	// data TSC
	__m128i* tsc_all;
	int16_t** rsc_msc;

	int Mlinesize;
	int nflags;	
	int partition;
	int Npartitions;

	double scale; 		// configurable
	int16_t wordoffset;	// configurable
	
	// Mutable args. Master thread must set them for workers
	int thrid;
	P7_PROFILE *gm;
	ESL_DSQ	*seqs[8];	// seqs to align
	int allocL;		// size of allocated sequence buffers
	int L;			// max length of seqs to align

	volatile uchar *synchflags[200];	// max threads
	// counter with # of runs to synchronize the threads starting each run
	int synccontrol;
	pthread_barrier_t barrier;

	__m128i* xmxE;
	__m128i* pdcv;
	__m128i* psc;
} DATA_STREAM;

#define MAX_PARTITION	112	// 112, 96 or 88 seem to work best
//#define NTHREADS 4

#define WORDMAX 32767	// constant
#define tprintf if(0) printf
#define dprintf if(0) printf
#define SSTEST 0


__m128i* alloc_m128_aligned64(int nvecs)
{
	void *mem;
	posix_memalign(&mem, 64, (nvecs)*sizeof(__m128i));
	return (__m128i*) mem;
}

void* alloc_aligned64(int nbytes)
{
	void *mem;
	posix_memalign(&mem, 64, nbytes);
	return mem;
}

int roundtop(int v, int align)
{
	if (v%align==0)	return v;
	else			return v + (align - v%align);
}

int roundfloor(int v, int align)
{
	if (v%align==0)	return v;
	else			return v - (v%align);
}

int16_t discretize(float scale, float ff)
{
	if (isinf(ff))	return -WORDMAX;
	else			return (int16_t) (roundf (ff * scale));
}


#if 0
// to see the annotated source in profiling
#define _mm_max_epi16(__A,__B)		(__m128i)__builtin_ia32_pmaxsw128	((__v8hi)__A, (__v8hi)__B)
#define _mm_and_si128(__A,__B)		(__m128i)__builtin_ia32_pand128		((__v2di)__A, (__v2di)__B)
#define _mm_adds_epi16(__A,__B)		(__m128i)__builtin_ia32_paddsw128	((__v8hi)__A, (__v8hi)__B)
#define _mm_cmplt_epi16(__A,__B)	(__m128i)__builtin_ia32_pcmpgtw128	((__v8hi)__B, (__v8hi)__A)
#define _mm_andnot_si128(__A,__B)	(__m128i)__builtin_ia32_pandn128	((__v2di)__A, (__v2di)__B)
#define _mm_unpacklo_epi16(__A,__B)	(__m128i)__builtin_ia32_punpcklwd128((__v8hi)__A, (__v8hi)__B)
#define _mm_unpackhi_epi16(__A,__B)	(__m128i)__builtin_ia32_punpckhwd128((__v8hi)__A, (__v8hi)__B)
#define _mm_unpacklo_epi32(__A,__B)	(__m128i)__builtin_ia32_punpckldq128((__v4si)__A, (__v4si)__B)
#define _mm_unpackhi_epi32(__A,__B)	(__m128i)__builtin_ia32_punpckhdq128((__v4si)__A, (__v4si)__B)
#define _mm_unpacklo_epi64(__A,__B)	(__m128i)__builtin_ia32_punpcklqdq128((__v2di)__A, (__v2di)__B)
#define _mm_unpackhi_epi64(__A,__B)	(__m128i)__builtin_ia32_punpckhqdq128((__v2di)__A, (__v2di)__B)

void printreg(char* msg, __m128i* reg)
{
	int16_t tt[8] __attribute__ ((aligned (16)));
	memmove(tt, reg, 16);
	printf("%s: %6d %6d %6d %6d %6d %6d %6d %6d\n", msg,
			tt[0], tt[1], tt[2], tt[3], tt[4], tt[5], tt[6], tt[7]);
}
#endif

int viterbi_stream_word_partitioned(DATA_STREAM* dstream, float* opt_res, int thrid)
{
//	if (NTHREADS > 1)		pthread_barrier_wait(&dstream->barrier);

	return 0;

	int L = dstream->L;
	P7_PROFILE* gm = dstream->gm;
	ESL_DSQ** ddsq = dstream->seqs;
	int M = gm->M, i, k, v, t, j;
	const int PARTITION = dstream->partition;
	__m128i** oprmsc = (__m128i**) dstream->rsc_msc;
	__m128i* xmxEv = dstream->xmxE;
	__m128i xmxB, xmxE, xmxC, moveC, Vinf = _mm_set1_epi16(-WORDMAX);

	__m128i dmx[PARTITION];
	__m128i mmx[PARTITION];
	__m128i imx[PARTITION];
	__m128i xmm[24];
	__m128i *mscore[8];
	__m128i overflowlimit, overflows;
	overflowlimit = overflows = Vinf;

	if (thrid == NTHREADS-1) 
	{	overflowlimit = _mm_set1_epi16(WORDMAX-1);
		overflows= _mm_xor_si128(overflows,overflows);	// zero out
	}

	t = ((dstream->Npartitions+thrid)%NTHREADS)*PARTITION;
	tprintf("START viterbiThr %d in %d L %d | Seq %d\n", thrid, t, L, 0); // ccount[thrid]++);

	xmxC  = Vinf;
	moveC = _mm_set1_epi16(discretize(dstream->scale, gm->xsc[p7P_C][p7P_MOVE]));
	xmxB  = _mm_set1_epi16(dstream->wordoffset + discretize(dstream->scale, gm->xsc[p7P_N][p7P_MOVE]));

	for (	; t < M; t += NTHREADS*PARTITION)
	{
		volatile uchar* synchflags1 = dstream->synchflags[t/PARTITION];
		volatile uchar* synchflags2 = dstream->synchflags[t/PARTITION+1];
		int t8 = t/8;

		for (k = 0; k < PARTITION; k++)
			dmx[k] = mmx[k] = imx[k] = Vinf;

		for (i = 1; i <= L; i++)
		{
		//	tprintf("Iter Thr %d t %d: I %d\n", thrid, t, i);
			__m128i sc, dcv, temp, mpv, ipv, dpv;
			__m128i *ttsc = dstream->tsc_all + t*8;
			v = i-1;
			ttsc += 3;

			if (t == 0)
				xmxE = mpv = dpv = ipv = sc = dcv = Vinf;
			else {
				if (NTHREADS > 1)
					 while (!synchflags1[v]) sched_yield();
				xmxE = xmxEv[v];
				dcv = dstream->pdcv[v];
				sc  = dstream->psc[v];
			}

			for (j = 0; j < 8; j++)
				mscore[j] = oprmsc[ddsq[j][i]] + t8;

			for (k = 0; k < PARTITION && t+k < M; )
			{
#if 0
#define EMLOAD(i)	xmm[i+24] = _mm_load_si128(mscore[i]); 
				EMLOAD(0) 	EMLOAD(1) 
				EMLOAD(2) 	EMLOAD(3) 
				EMLOAD(4) 	EMLOAD(5) 
				EMLOAD(6) 	EMLOAD(7) 

#define MIX16(i,r,range)	\
	xmm[r  ] = _mm_unpacklo_epi##range(xmm[24+i], xmm[24+i+1]);	\
	xmm[r+1] = _mm_unpackhi_epi##range(xmm[24+i], xmm[24+i+1]);
				MIX16(0,0,16)	MIX16(2,2,16)
				MIX16(4,4,16)	MIX16(6,6,16)
#else

#define MMLOAD(a,b)		\
	xmm[a] = _mm_unpacklo_epi16(*mscore[a], *mscore[b]);	\
	xmm[b] = _mm_unpackhi_epi16(*mscore[a], *mscore[b]);

				MMLOAD(0,1)	MMLOAD(2,3)
				MMLOAD(4,5)	MMLOAD(6,7)
#endif

#define MIX(i,r,range)	\
	xmm[r  ] = _mm_unpacklo_epi##range(xmm[i], xmm[i+2]);	\
	xmm[r+1] = _mm_unpackhi_epi##range(xmm[i], xmm[i+2]);

				MIX(0, 8,32)	MIX(1,12,32)
				MIX(4,10,32)	MIX(5,14,32)

				MIX( 8,16,64)	MIX( 9,18,64)
				MIX(12,20,64)	MIX(13,22,64)


#define TRIPLETCOMPUTE(k,j)	\
{	/* Calculate new M(k), delay store */	\
	sc = _mm_max_epi16(sc, _mm_adds_epi16(xmxB, *ttsc)); ttsc++;	\
	sc = _mm_adds_epi16(sc,  xmm[j]);		\
	/* Update E */							\
	xmxE = _mm_max_epi16(xmxE, sc);			\
	\
	/* Pre-emptive load of M, D, I */		\
	dpv = dmx[k];	\
	ipv = imx[k];	\
	mpv = mmx[k];	\
	\
	/* Calculate current I(k) */			\
	temp = _mm_adds_epi16(mpv, *ttsc); ttsc++;	\
	imx[k] = _mm_max_epi16(temp, _mm_adds_epi16(ipv, *ttsc)); ttsc++;\
	\
	/* Delayed stores of M and D */			\
	mmx[k] = sc;	\
	dmx[k] = dcv;	\
	\
	/* Calculate next D, D(k+1) */			\
	sc	= _mm_adds_epi16(sc, *ttsc); ttsc++;	\
	dcv = _mm_max_epi16(sc, _mm_adds_epi16(dcv, *ttsc));ttsc++;	\
	\
	/* Pre-emptive partial calculation of M(k+1) */	\
	sc = _mm_adds_epi16(mpv, *ttsc); ttsc++;	\
	sc = _mm_max_epi16(sc, _mm_adds_epi16(ipv, *ttsc)); ttsc++;	\
	sc = _mm_max_epi16(sc, _mm_adds_epi16(dpv, *ttsc)); ttsc++;	\
	k++;			\
}
				TRIPLETCOMPUTE(k,16+0)	TRIPLETCOMPUTE(k,16+1)
				TRIPLETCOMPUTE(k,16+2)	TRIPLETCOMPUTE(k,16+3)
				TRIPLETCOMPUTE(k,16+4)	TRIPLETCOMPUTE(k,16+5)
				TRIPLETCOMPUTE(k,16+6)	TRIPLETCOMPUTE(k,16+7)

				mscore[0]++; mscore[1]++; mscore[2]++; mscore[3]++;
				mscore[4]++; mscore[5]++; mscore[6]++; mscore[7]++;
			}

			if (t+k < M)
			{	v = i-1;
				xmxEv[v] = xmxE;
				dstream->pdcv[v] = dcv;
				dstream->psc [v] = sc;

				if (NTHREADS > 1) synchflags2[v] = 1;
			}
			else	// executed only by main thread (NTHRS-1)
			{
				__m128i overfs = _mm_cmpgt_epi16(xmxE, overflowlimit);
				overflows = _mm_or_si128(overflows, overfs);	// select the overflowed channels
				xmxC	= _mm_max_epi16(xmxC, xmxE);
			}
		}
	}

	xmxC = _mm_adds_epi16(xmxC, moveC);

	if (opt_res != NULL)
	{
		float offset = (float) dstream->wordoffset;
		int16_t res[8] __attribute__ ((aligned (16)));
		int16_t ovs[8] __attribute__ ((aligned (16)));

		memmove(res, &xmxC, sizeof(xmxC));
		memmove(ovs, &overflows, sizeof(overflows));

		for (i = 0; i < 8; i++)
			if (ovs[i])	opt_res[i] = eslINFINITY;	// signal overflow
			else		opt_res[i] = ((float) res[i] - offset) / dstream->scale - 2.0;
												// 2.0 nat approximation, UNILOCAL mode
	}

	tprintf("END viterbi Thr %d - t %d\n", thrid, t);
	return eslOK;
}





/***************************************************************************
 *	Thread Management
 ***************************************************************************/

typedef struct _thr_info_t
{	int thrid;
	DATA_STREAM *dstream;
} pthr_info_t;


// variable to synchronize the threads starting each run of the algo.
// has a counter with # of runs
sem_t semsynch[10];
volatile int syncflags[10];

void* viterbi_stream_thread_loop(void* argst)
{
	int i, execcount = 0;
	pthr_info_t* args = (pthr_info_t*) argst;
	DATA_STREAM *dstream = args->dstream;
	
#ifdef _GNU_SOURCE
	printf("THR %d running on cpu %d\n", args->thrid, sched_getcpu());
#endif

	for (i = 0; 1; i++)
	{	
		execcount++; while (dstream->synccontrol != execcount) sched_yield();

//		while (syncflags[args->thrid] == 0) sched_yield(); 	syncflags[args->thrid] = 0;
		
//		sem_wait(&semsynch[args->thrid]);

		tprintf("THR %d entering\n", args->thrid);

		viterbi_stream_word_partitioned(dstream, NULL, args->thrid);
	}

#ifdef _GNU_SOURCE
	printf("THR %d running on cpu %d\n", args->thrid, sched_getcpu());
#endif
	return (void*) 0;
}


#define MAPIDCPU(id) (id % 4)	// ((id)*4 +2)

void viterbi_stream_create_threads(DATA_STREAM* dstream)
{
	int i;
	pthread_barrier_init(&dstream->barrier, NULL, NTHREADS);

	for (i = 0; i < NTHREADS-1; i++)
		sem_init(&semsynch[i], 0, 0);

	if(NTHREADS == 1)
		return;

	pthread_t threads[NTHREADS];
	pthr_info_t args;
	args.dstream = dstream; args.thrid = 0;

#ifdef _GNU_SOURCE
	cpu_set_t cpuset;
	CPU_ZERO(&cpuset);
	CPU_SET(MAPIDCPU(NTHREADS-1), &cpuset);
	threads[NTHREADS-1] = pthread_self();
	pthread_setaffinity_np(threads[NTHREADS-1], sizeof(cpu_set_t), &cpuset);
#endif

	for (i = 0; i < NTHREADS-1; i++)
	{
		pthr_info_t *argscopy = calloc(1, sizeof(pthr_info_t));
		memcpy(argscopy, &args, sizeof(pthr_info_t));
		argscopy->thrid	= i;

		pthread_attr_t attr;
		pthread_attr_init(&attr);
#ifdef _GNU_SOURCE
		CPU_ZERO(&cpuset);
		CPU_SET(MAPIDCPU(i), &cpuset);
		pthread_attr_setaffinity_np(&attr, sizeof(cpu_set_t), &cpuset);
#endif
		if (pthread_create(&threads[i], &attr, viterbi_stream_thread_loop, argscopy))
			exit(fprintf(stderr, "ERROR could not create worker thread\n"));
	}

#ifdef _GNU_SOURCE
	printf("THR %d running on cpu %d\n", NTHREADS-1, sched_getcpu());
#endif
}



/***************************************************************************
 *	Prepare Viterbi Stream
 ***************************************************************************/

void free_stream_buffers(DATA_STREAM* dstream)
{
	int j; 
	for (j = 0; j < 8; j++)
		free(dstream->seqs[j]);
	for (j = 1; j <= dstream->Npartitions; j++)
		free(dstream->synchflags[j]);

	free(dstream->xmxE);
	free(dstream->pdcv);
	free(dstream->psc);
}

void alloc_stream_buffers(DATA_STREAM* dstream, int maxL)
{
	int j; 
	for (j = 0; j < 8; j++)
		dstream->seqs[j] = calloc(dstream->allocL+16, sizeof(ESL_DSQ));

	dstream->nflags = maxL;
	for (j = 1; j <= dstream->Npartitions; j++)
		dstream->synchflags[j] = alloc_aligned64(dstream->nflags);
	dstream->synchflags[0] = dstream->synchflags[dstream->Npartitions+1] = NULL;

	dstream->psc  = alloc_m128_aligned64((maxL));
	dstream->pdcv = alloc_m128_aligned64((maxL));
	dstream->xmxE = alloc_m128_aligned64((maxL+1));
}

DATA_STREAM* p7_ViterbiStream_Create(P7_PROFILE *gm)
{
	float *tsc= gm->tsc;
	__m128i *tt;
	int M = gm->M, Mtop, k, i;
	DATA_STREAM* dr = calloc(1, sizeof(DATA_STREAM));
	Mtop = dr->Mlinesize = roundtop(M,8);
	dr->tsc_all = tt = alloc_m128_aligned64((Mtop+1)*8);
	dr->tsc_all = tt;
	dr->scale	= 1000.0;
	dr->wordoffset = 12000;
	dr->gm = gm;

#define SETSCORE(d,val)	*d++ = _mm_set1_epi16(val);
	for (k = 0; k < M; k++, tsc += p7P_NTRANS)
	{
		SETSCORE(tt, discretize(dr->scale, tsc[p7P_MM]));
		SETSCORE(tt, discretize(dr->scale, tsc[p7P_IM]));
		SETSCORE(tt, discretize(dr->scale, tsc[p7P_DM]));
		SETSCORE(tt, discretize(dr->scale, tsc[p7P_BM]));
		if (k < M-1) 
		{	SETSCORE(tt, discretize(dr->scale, tsc[p7P_NTRANS+p7P_MI]));
			SETSCORE(tt, discretize(dr->scale, tsc[p7P_NTRANS+p7P_II]));
			SETSCORE(tt, discretize(dr->scale, tsc[p7P_NTRANS+p7P_MD]));
			SETSCORE(tt, discretize(dr->scale, tsc[p7P_NTRANS+p7P_DD]));
		}
		else
		{	SETSCORE(tt, -WORDMAX);	// para MI
			SETSCORE(tt, -WORDMAX);	// para II
			SETSCORE(tt, -WORDMAX);	// para MD
			SETSCORE(tt, -WORDMAX);	// para DD
		}
	}

	// fill the left-over places with dummy values
	for (	; k < Mtop+1; k++, tsc += p7P_NTRANS)
		for (i = 0; i < p7P_NTRANS; i++)
			SETSCORE(tt, -WORDMAX);

	const int AlphaSize = gm->abc->Kp;
	dr->rsc_msc = calloc(AlphaSize+2, sizeof(int16_t*));
	int16_t *data = alloc_aligned64((Mtop)*(AlphaSize+1)*sizeof(int16_t));

	for (k = 0; k < AlphaSize; k++)
	{
		int16_t *rmsc = dr->rsc_msc[k] = data + (k)*(Mtop);
		float *rsc = gm->rsc[k] + p7P_NR + p7P_MSC;	// ignore 1st values
		// #define MSC(k)	 (rsc[(k)p7P_NR + p7P_MSC])
		for (i = 0; i< M; i++, rsc += p7P_NR)
			*rmsc++ = discretize(dr->scale, *rsc);
		// fill the left-over places with neutral values
		for (	; i < Mtop; i++)
			*rmsc++ = -WORDMAX;
	}

	// create dummy values for k = AlphaSize, to use the padded sequences
	dr->rsc_msc[k] = data + (k)*(Mtop);
	for (i = 0; i< Mtop; i++)
		dr->rsc_msc[k][i] = -WORDMAX;

	return dr;
}


void p7_ViterbiStream_Setup(DATA_STREAM* dstream, int maxL, int max_partition)
{
	int M = dstream->gm->M, nblocks;

	// choose highest possible no. of partitions, multiple of no. threads
	// use (double) to force an non-integer max over MaxPart
//	for (nblocks = NTHREADS; M / (double) nblocks > max_partition; nblocks += NTHREADS);
	// same as:
	nblocks = roundtop(ceil(((M*NTHREADS)/(double)(NTHREADS*max_partition))),NTHREADS);

	dstream->partition	= roundtop(ceil(M/(double)nblocks), 8);
	dstream->Npartitions= M/dstream->partition + (M % dstream->partition != 0);

	dstream->L = 0;
	dstream->allocL = maxL;
	alloc_stream_buffers(dstream, maxL);

	viterbi_stream_create_threads(dstream);
}


typedef struct _dsq_cmp
{
	ESL_DSQ *seq;
	int length;
} dsq_cmp_t;

int p7_ViterbiStream(DATA_STREAM* dstream, dsq_cmp_t **seqsdb, float* results)
{
	int j, maxL = 0;

#if 0
	for (j = 0; j < 8; j++)
		if (seqsdb[j]->length > maxL)
			maxL = seqsdb[j]->length;

//	maxL = 2000;
//	for (j = 0; j < 8; j++)	printf("%d ", seqsdb[j]->length); 	printf("\n"); 

	if (maxL > dstream->allocL-10)
	{	// realloc
		dstream->allocL = ESL_MAX(dstream->allocL*2, roundtop(maxL,1024));
		free_stream_buffers(dstream);
		alloc_stream_buffers(dstream, dstream->allocL);
	}

	if (dstream->L != maxL)
	{   dstream->L  = maxL;
		p7_ReconfigLength(dstream->gm, maxL);
	}

	for (j = 0; j < 8; j++)
		memcpy(dstream->seqs[j], seqsdb[j]->seq, (seqsdb[j]->length+1)*sizeof(ESL_DSQ));
	
	// pad sequences with k = Alphasize
	for (j = 0; j < 8; j++)
		memset(dstream->seqs[j]+seqsdb[j]->length+1, dstream->gm->abc->Kp, maxL-seqsdb[j]->length);

	for (j = 1; j <= dstream->gm->M/dstream->partition; j++)
		memset((void*) dstream->synchflags[j], 0, dstream->nflags);

	for (j = 1; j <= dstream->gm->M/dstream->partition; j++)
		memset((void*) dstream->synchflags[j], 0, dstream->nflags);
#endif

if(1)
	for (j = 0; j < NTHREADS-1; j++)
	{	sem_post(&semsynch[j]);
		syncflags[j] = 1;
	}
	
	dstream->synccontrol++;
	viterbi_stream_word_partitioned(dstream, results, NTHREADS-1	);
	return maxL;
}

int compare_seqs(const void *vs1, const void *vs2)
{
	dsq_cmp_t *s1 = *((dsq_cmp_t**) vs1);
	dsq_cmp_t *s2 = *((dsq_cmp_t**) vs2);
	if (s1->length == s2->length)
		return 0;
	return (s1->length > s2->length);
}





/*****************************************************************
 *	Benchmark driver.
 *****************************************************************/

/*
gcc -g -Wall -O3 -std=gnu99 -o viterbistream -I. -I.. -L.. -I../impl_sse -L../impl_sse -I../../easel -L../../easel viterbi_stream.c viterbi_serial.c -lhmmer -lhmmerimpl -leasel -lm -lpthread -DNTHREADS=1
*/

#include "esl_stopwatch.h"
#include "esl_sqio.h"
#include "esl_randomseq.h"

#include <sys/timeb.h>
#include <sched.h>

static ESL_OPTIONS options[] = {
	/* name		type	default		env	range toggles reqs incomp	help						docgroup*/
	{ "-h",	eslARG_NONE	, FALSE	, NULL, NULL , NULL, NULL, NULL, "show brief help on version and usage",0 },
	{ "-s",	eslARG_INT	, "42"	, NULL, NULL , NULL, NULL, NULL, "set random number seed to <n>",		0 },
	{ "-L",	eslARG_INT	, "400"	, NULL, "n>0", NULL, NULL, NULL, "length of random target seqs",		0 },
	{ "-N",	eslARG_INT	, "50000",NULL, "n>0", NULL, NULL, NULL, "number of random target seqs",		0 },
	{ "-M", eslARG_INT	, "112" , NULL, "n>0", NULL, NULL, NULL, "maximum size for the partitions",		0 },
	{	 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage [] = "[-options] <hmmfile>";
static char banner[] = "benchmark driver for generic Viterbi";


#define TIMEDIFF(t1,t2) ((t2.time - t1.time) + (t2.millitm - t1.millitm)*0.001)

int main(int argc, char **argv)
{
	ESL_GETOPTS   *go	= esl_getopts_CreateDefaultApp(options, 2, argc, argv, banner, usage);
	char	*hmmfile	= esl_opt_GetArg(go, 1);
	char	*seqfile	= esl_opt_GetArg(go, 2);
	ESL_STOPWATCH *w	= esl_stopwatch_Create();
	ESL_RANDOMNESS*r	= esl_randomness_Create(esl_opt_GetInteger(go, "-s"));
	ESL_ALPHABET*abc	= NULL;
	P7_HMMFILE	*hfp	= NULL;
	P7_HMM		*hmm	= NULL;
	P7_BG		*bg		= NULL;
	P7_PROFILE	*gm1, *gm2;
	int			L		= 2000;		// esl_opt_GetInteger(go, "-L");
	int			N		= esl_opt_GetInteger(go, "-N");
	int			MaxPart	= esl_opt_GetInteger(go, "-M");
	__m128		resdata[10];
	float		*sc1 	= (float*)	(resdata+0);	// uses 1 __m128s
	ESL_DSQ		*dsq	= NULL;
	int			i, j;
	ESL_SQFILE   *sqfp	= NULL;
	DATA_STREAM *dstream;
	struct timeb tbstart, tbend;
	int sumlengths = 0;

	if (p7_hmmfile_Open(hmmfile, NULL, &hfp)!= eslOK) p7_Fail("Failed to open HMM file %s", hmmfile);
	if (p7_hmmfile_Read(hfp, &abc, &hmm)	!= eslOK) p7_Fail("Failed to read HMM");

	bg = p7_bg_Create(abc);
	p7_bg_SetLength(bg, L);

	gm1 = p7_profile_Create(hmm->M, abc);
	gm2 = p7_profile_Create(hmm->M, abc);
	p7_ProfileConfig(hmm, bg, gm1, L, p7_UNILOCAL);
	p7_ProfileConfig(hmm, bg, gm2, L, p7_UNILOCAL);

	dstream = p7_ViterbiStream_Create(gm1);
	p7_ViterbiStream_Setup(dstream, L+100, MaxPart); // use max L
	dstream->L = L;

	// No. of partitions computed without full parallelism ( == no. of threads active while some are idle)
	int Niters_part	= dstream->Npartitions % NTHREADS;
	// No. of Model lines that could be computed but are wasted by idle threads waiting on the end
	int Nwasted_threads	= dstream->partition * ((NTHREADS-Niters_part) % NTHREADS);
	// No. of lines in the last partition that go beyond M. It's wasted comp time by a single thread
	int Nwasted_leftover= (dstream->partition - gm1->M % dstream->partition) % dstream->partition;
	// Total number of wasted lines
	int wastedcomp = Nwasted_threads + Nwasted_leftover;
	// Total number of lines computed and waited for
	int totalcomp = wastedcomp + gm1->M; // same as: roundtop(gm1->M, dstream->partition * NTHREADS);

	printf("Viterbi Stream Word with %d Threads, model %s: Modelsize %d, #Segms: %d, SeqL: %d, Nseqs %d, Part %d, #Parts %d\n",
			NTHREADS, hmmfile, gm1->M, (int) ceil(gm1->M/8.0), L, 8*N, dstream->partition, dstream->Npartitions);
	printf("Total Comp Lines: %d | Wasted Comp Lines: %d\n", totalcomp, wastedcomp);

	// for ViterbiFilter
	P7_OPROFILE *om = p7_oprofile_Create(hmm->M, gm1->abc);
	p7_oprofile_Convert(gm1, om);
	P7_OMX		*ox	= p7_omx_Create(hmm->M, 0, 0);

	dsq_cmp_t **seqsdb= calloc(8*N+64, sizeof(dsq_cmp_t*));
if(0)
{	ESL_SQ* sq = esl_sq_CreateDigital(abc);
	if (esl_sqfile_OpenDigital(abc, seqfile, eslSQFILE_FASTA, NULL, &sqfp) != eslOK)
	{	p7_Fail("Failed to open sequence file\n");	return -1; }

	for (j = 0; j < 8*N; j++)
	{
		int res = esl_sqio_Read(sqfp, sq);
		if (res != eslOK) { printf("ATENCAO: faltam sequencias\n"); break; }
		int len = sq->n;
		dsq = sq->dsq;
		seqsdb[j] = malloc(sizeof(dsq_cmp_t));
		seqsdb[j]->length = len;
		seqsdb[j]->seq = malloc((len+4)*sizeof(ESL_DSQ));
		memcpy(seqsdb[j]->seq, dsq, len+2);
		sumlengths += len;
		esl_sq_Reuse(sq);
	}

	ftime(&tbstart);
	N = j/8;
	printf("N = %d\n", N);
		// Sort sequences by length
	qsort(seqsdb, N*8, sizeof(dsq_cmp_t*), compare_seqs);
}
else if(0)
	for (i = 0; i < N; i++)
	{
		for (j = 0; j < 8; j++)
		{
			int len = L - rand()%1000;
			seqsdb[i*8+j] = malloc(sizeof(dsq_cmp_t));
			seqsdb[i*8+j]->seq = malloc(len+4);
			seqsdb[i*8+j]->length = len;
		    esl_rsq_xfIID(r, bg->f, abc->K, len, seqsdb[i*8+j]->seq);
			sumlengths += len;
		}
	}

//	double sumerrors = 0;
	float* results = (float*) alloc_m128_aligned64(N*2+2);	
	ftime(&tbstart);

	for (j = 0; j < N; j++)
	for (i = 0; i < N; i++)
	{
	//	if (i % 10000 == 0) printf("START %d\n", i);

		p7_ViterbiStream(dstream, seqsdb+i*8, sc1);

	//	memcpy(results+i*8, sc1, 32);
	}
	ftime(&tbend);

	double secs = TIMEDIFF(tbstart,tbend);
//	printf("Qsort time: %6.3f | Viterbi time: %6.3f\n", TIMEDIFF(tbqsort,tbstart), secs);
	w->elapsed = w->user = secs;
	esl_stopwatch_Display(stdout, w, "# Opt CPU time: ");
	printf("# %.0fM cells in %.1f Mc/s\n", (sumlengths * (double) gm1->M) / 1e6, (sumlengths * (double) gm1->M * 1e-6) / secs);

if(0)	// compare results against base version
	for (i = 0; i < 1000 && i < N; i++)
	{
		int maxll = 0; float sc2;
		for (j = 0; j < 8; j++)
			if (maxll < seqsdb[i*8+j]->length)
				maxll = seqsdb[i*8+j]->length;

//		for (j = 0; j < 8; j++)	printf("%d ", seqsdb[i*8+j]->length); printf("\n");
//		if (i % 10 == 0) printf("i %d\n", i);

		p7_oprofile_ReconfigRestLength(om, maxll);
		p7_ReconfigLength(gm2, maxll);	// fazer Reconfig aqui para emular compl o VitStream

		for (j = 0; j < 8; j++)
		{
//			p7_ReconfigLength(gm2, seqsdb[i*8+j]->length);
//			p7_Viterbi_unilocal(seqsdb[i*8+j]->seq, seqsdb[i*8+j]->length, gm2, &sc3);
//			p7_Viterbi_unilocal_word(seqsdb[i*8+j]->seq, seqsdb[i*8+j]->length, gm2, &sc2);

//			p7_oprofile_ReconfigLength(om, seqsdb[i*8+j]->length);
			p7_ViterbiFilter(seqsdb[i*8+j]->seq, seqsdb[i*8+j]->length, om, ox, &sc2);
			
			//sumerrors += fabs(sc1[j]- sc2);
			if (fabs(results[i*8+j] - sc2) > 0.00001)
			{	printf("%3d-%d L %4d: VS %f d %f\t| Base SerI %f\n", i, j, seqsdb[i*8+j]->length, 
						results[i*8+j], fabs(results[i*8+j] - sc2), sc2);
				getc(stdin);
			}
		}
	}

	return 0;
}

