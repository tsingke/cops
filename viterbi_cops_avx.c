/* 
This file is part of the COPS (Cache-Oblivious Parallel SIMD Viterbi) prototype.

 COPS by Miguel M. A. Ferreira, Luis M. S. Russo and Nuno V. Roma / KDBIO & SIPS /
INESC-ID is licensed under a Creative Commons Attribution 3.0 Unported
License.  Permissions beyond the scope of this license may be
available at #email.

 The Software is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

 Any published media that is related with to use of the distributed
software, or derived software, must contain a reference to "Cache-Oblivious
Parallel SIMD Viterbi Decoding for Sequence Search in HMMER, Miguel Ferreira, 
Nuno Roma, Luis Russo (submitted to bmc bioinformatics)", but not in any way that
suggests that they endorse you or your use of the work.

 License available at http://creativecommons.org/licenses/by/3.0/
*/

/*
 * Stub for a AVX port. Incomplete and still untested
 */

#define _GNU_SOURCE
#include <xmmintrin.h>		/* SSE	*/
#include <emmintrin.h>		/* SSE2 */

#include <hmmer.h>
#include <p7_config.h>

#include <pthread.h>
#include <string.h>

typedef unsigned char uchar;
typedef unsigned char byte;

extern int p7_ViterbiFilter(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_OMX *ox, float *ret_sc);


#include "viterbi_serial.h"


/***************************************************************************
 *	Viterbi COPS with 16-channel 16-bit integers and threaded PARTITIONS
 ***************************************************************************/

// comment to use volatile
#define volatile

typedef struct _DATA_COPSAVX
{	// data TSC
	__m256i* tsc_all;
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
	ESL_DSQ	*seqs[16];	// seqs to align
	int allocL;		// size of allocated sequence buffers
	int L;			// max length of seqs to align

	volatile uchar *synchflags[200];	// max threads
	// counter with # of runs to synchronize the threads starting each run
	volatile int synccontrol;

	__m256i* xmxE;
	__m256i* pdcv;
	__m256i* psc;
} DATA_COPSAVX;

#define MAX_PARTITION	112	// 112, 96 or 88 seem to work best
//#define NTHREADS 4

#define WORDMAX 32767	// constant
#define tprintf if(0) printf
#define dprintf if(0) printf


__m256i* alloc_m128_aligned64(int nvecs)
{
	void *mem;
	posix_memalign(&mem, 64, (nvecs)*sizeof(__m256i));
	return (__m256i*) mem;
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
#define _mm256_max_epi16(__A,__B)		(__m256i)__builtin_ia32_pmaxsw128	((__v8hi)__A, (__v8hi)__B)
#define _mm256_and_si128(__A,__B)		(__m256i)__builtin_ia32_pand128		((__v2di)__A, (__v2di)__B)
#define _mm256_adds_epi16(__A,__B)		(__m256i)__builtin_ia32_paddsw128	((__v8hi)__A, (__v8hi)__B)
#define _mm256_cmplt_epi16(__A,__B)	(__m256i)__builtin_ia32_pcmpgtw128	((__v8hi)__B, (__v8hi)__A)
#define _mm256_andnot_si128(__A,__B)	(__m256i)__builtin_ia32_pandn128	((__v2di)__A, (__v2di)__B)
#define _mm256_unpacklo_epi16(__A,__B)	(__m256i)__builtin_ia32_punpcklwd128((__v8hi)__A, (__v8hi)__B)
#define _mm256_unpackhi_epi16(__A,__B)	(__m256i)__builtin_ia32_punpckhwd128((__v8hi)__A, (__v8hi)__B)
#define _mm256_unpacklo_epi32(__A,__B)	(__m256i)__builtin_ia32_punpckldq128((__v4si)__A, (__v4si)__B)
#define _mm256_unpackhi_epi32(__A,__B)	(__m256i)__builtin_ia32_punpckhdq128((__v4si)__A, (__v4si)__B)
#define _mm256_unpacklo_epi64(__A,__B)	(__m256i)__builtin_ia32_punpcklqdq128((__v2di)__A, (__v2di)__B)
#define _mm256_unpackhi_epi64(__A,__B)	(__m256i)__builtin_ia32_punpckhqdq128((__v2di)__A, (__v2di)__B)

void printreg(char* msg, __m256i* reg)
{
	int16_t tt[8] __attribute__ ((aligned (16)));
	memmove(tt, reg, 16);
	printf("%s: %6d %6d %6d %6d %6d %6d %6d %6d\n", msg,
			tt[0], tt[1], tt[2], tt[3], tt[4], tt[5], tt[6], tt[7]);
}

_mm256_add_epi16 
_mm256_sub_epi16 
_mm256_max_epi16 
__m256i _mm256_unpackhi_epi8(__m256i m1, __m256i m2)
__m256i _mm256_unpackhi_epi16(__m256i m1,__m256i m2)
__m256i _mm256_unpackhi_epi32(__m256i m1, __m256i m2)
__m256i _mm256_unpackhi_epi64( __m256i a, __m256i b)

#endif


#define SET1epi16(val) \
 _mm_broadcastsi128_si256(SET1epi16(val))

int viterbi_cops_word_partitioned(DATA_COPSAVX* dcops, float* opt_res, int thrid)
{
	int L = dcops->L;
	P7_PROFILE* gm = dcops->gm;
	ESL_DSQ** ddsq = dcops->seqs;
	int M = gm->M, i, k, v, t, j;
	const int PARTITION = dcops->partition;
	__m256i** oprmsc = (__m256i**) dcops->rsc_msc;
	__m256i* xmxEv = dcops->xmxE;
	__m256i xmxB, xmxE, xmxC, moveC, Vinf = SET1epi16(-WORDMAX);

	__m256i dmx[PARTITION];
	__m256i mmx[PARTITION];
	__m256i imx[PARTITION];
	__m256i xmm[60];
	__m256i *mscore[16];
	__m256i overflowlimit, overflows;
	overflowlimit = overflows = Vinf;

	__m256i masklows, maskhighs;
	masklows = _mm256_xor_si256( masklows,masklows );	// zero out
	maskhighs= _mm256_xor_si256(maskhighs,maskhighs);	// zero out
	// load masks
	memset(&masklows, 0xff, 128);
	memset(((byte*) &maskhighs)+128, 0xff, 128);

	if (thrid == NTHREADS-1) 
	{	overflowlimit = SET1epi16(WORDMAX-1);
		overflows= _mm256_xor_si256(overflows,overflows);	// zero out
	}

	t = ((dcops->Npartitions+thrid)%NTHREADS)*PARTITION;
	tprintf("START viterbiThr %d in %d L %d | Seq %d\n", thrid, t, L, 0); // ccount[thrid]++);

	xmxC  = Vinf;
	moveC = SET1epi16(discretize(dcops->scale, gm->xsc[p7P_C][p7P_MOVE]));
	xmxB  = SET1epi16(dcops->wordoffset + discretize(dcops->scale, gm->xsc[p7P_N][p7P_MOVE]));

	for (	; t < M; t += NTHREADS*PARTITION)
	{
		volatile uchar* synchflags1 = dcops->synchflags[t/PARTITION];
		volatile uchar* synchflags2 = dcops->synchflags[t/PARTITION+1];
		int t16 = t/16;

		for (k = 0; k < PARTITION; k++)
			dmx[k] = mmx[k] = imx[k] = Vinf;

		for (i = 1; i <= L; i++)
		{
		//	tprintf("Iter Thr %d t %d: I %d\n", thrid, t, i);
			__m256i sc, dcv, temp, mpv, ipv, dpv;
			__m256i *ttsc = dcops->tsc_all + t*16;
			v = i-1;
			ttsc += 3;

			if (t == 0)
				xmxE = mpv = dpv = ipv = sc = dcv = Vinf;
			else {
				if (NTHREADS > 1)
					 while (!synchflags1[v]) sched_yield();
				xmxE = xmxEv[v];
				dcv = dcops->pdcv[v];
				sc  = dcops->psc[v];
			}

			for (j = 0; j < 16; j++)
				mscore[j] = oprmsc[ddsq[j][i]] + t16;

			for (k = 0; k < PARTITION && t+k < M; )
			{
#define MMLOAD(a,b)		\
	xmm[a] = _mm256_unpacklo_epi16(*mscore[a], *mscore[b]);	\
	xmm[b] = _mm256_unpackhi_epi16(*mscore[a], *mscore[b]);

				MMLOAD( 0, 1)	MMLOAD( 2, 3)
				MMLOAD( 4, 5)	MMLOAD( 6, 7)
				MMLOAD( 8, 9)	MMLOAD(10,11)
				MMLOAD(12,13)	MMLOAD(14,15)
// lows  0 2 4 ...
// highs 1 3 5

#define MIX(i,r,range)	\
	xmm[r  ] = _mm256_unpacklo_epi##range(xmm[i], xmm[i+2]);	\
	xmm[r+1] = _mm256_unpackhi_epi##range(xmm[i], xmm[i+2]);

				MIX( 0,16,32)	MIX( 1,18,32)
				MIX( 4,20,32)	MIX( 5,22,32)
				MIX( 8,24,32)	MIX( 7,26,32)
				MIX(12,28,32)	MIX(13,30,32)
// lows  16 18 20 22 ...
// highs 17 19 21 23
				MIX(16,32,64)	MIX(17,34,64)
				MIX(20,36,64)	MIX(21,38,64)
				MIX(24,40,64)	MIX(25,42,64)
				MIX(28,44,64)	MIX(29,46,64)
// lows  32 34 36 38 ...
// highs 33 35 37 39

#define MIX(i,r,range)		\
{	__mm256i temp1, temp2;	\
	temp1 = _mm256_and_si256(xmm[i  ], masklows);	\
	temp2 = _mm256_and_si256(xmm[i+2], masklows);	\
	xmm[r]= _mm256_or_si256(temp1, temp2);			\
\
	temp1 = _mm256_and_si256(xmm[i  ], maskhighs);	\
	temp2 = _mm256_and_si256(xmm[i+2], maskhighs);	\
	xmm[r+1] = _mm256_or_si256(temp1, temp2);		\
}
				UMIX(32,0,128)	UMIX(33,2,128)
				UMIX(36,4,128)	UMIX(37,6,128)
				UMIX(40,8,128)	UMIX(41,10,128)
				UMIX(44,12,128)	UMIX(45,14,128)


#define TRIPLETCOMPUTE(k,j)	\
{	/* Calculate new M(k), delay store */	\
	sc = _mm256_max_epi16(sc, _mm256_adds_epi16(xmxB, *ttsc)); ttsc++;	\
	sc = _mm256_adds_epi16(sc,  xmm[j]);		\
	/* Update E */							\
	xmxE = _mm256_max_epi16(xmxE, sc);			\
	\
	/* Pre-emptive load of M, D, I */		\
	dpv = dmx[k];	\
	ipv = imx[k];	\
	mpv = mmx[k];	\
	\
	/* Calculate current I(k) */			\
	temp = _mm256_adds_epi16(mpv, *ttsc); ttsc++;	\
	imx[k] = _mm256_max_epi16(temp, _mm256_adds_epi16(ipv, *ttsc)); ttsc++;\
	\
	/* Delayed stores of M and D */			\
	mmx[k] = sc;	\
	dmx[k] = dcv;	\
	\
	/* Calculate next D, D(k+1) */			\
	sc	= _mm256_adds_epi16(sc, *ttsc); ttsc++;	\
	dcv = _mm256_max_epi16(sc, _mm256_adds_epi16(dcv, *ttsc));ttsc++;	\
	\
	/* Pre-emptive partial calculation of M(k+1) */	\
	sc = _mm256_adds_epi16(mpv, *ttsc); ttsc++;	\
	sc = _mm256_max_epi16(sc, _mm256_adds_epi16(ipv, *ttsc)); ttsc++;	\
	sc = _mm256_max_epi16(sc, _mm256_adds_epi16(dpv, *ttsc)); ttsc++;	\
	k++;			\
}
				TRIPLETCOMPUTE(k, 0)	TRIPLETCOMPUTE(k, 1)
				TRIPLETCOMPUTE(k, 2)	TRIPLETCOMPUTE(k, 3)
				TRIPLETCOMPUTE(k, 4)	TRIPLETCOMPUTE(k, 5)
				TRIPLETCOMPUTE(k, 6)	TRIPLETCOMPUTE(k, 7)
				TRIPLETCOMPUTE(k, 8)	TRIPLETCOMPUTE(k, 9)
				TRIPLETCOMPUTE(k,10)	TRIPLETCOMPUTE(k,11)
				TRIPLETCOMPUTE(k,12)	TRIPLETCOMPUTE(k,13)
				TRIPLETCOMPUTE(k,14)	TRIPLETCOMPUTE(k,15)

				mscore[0]++; mscore[1]++; mscore[2]++; mscore[3]++;
				mscore[4]++; mscore[5]++; mscore[6]++; mscore[7]++;
				mscore[8]++; mscore[9]++; mscore[10]++; mscore[11]++;
				mscore[12]++; mscore[13]++; mscore[14]++; mscore[15]++;
			}

			if (t+k < M)
			{	v = i-1;
				xmxEv[v] = xmxE;
				dcops->pdcv[v] = dcv;
				dcops->psc [v] = sc;

				if (NTHREADS > 1) synchflags2[v] = 1;
			}
			else	// executed only by main thread (NTHRS-1)
			{
				__m256i overfs = _mm256_cmpgt_epi16(xmxE, overflowlimit);
				overflows = _mm256_or_si256(overflows, overfs);	// select the overflowed channels
				xmxC	= _mm256_max_epi16(xmxC, xmxE);
			}
		}
	}

	xmxC = _mm256_adds_epi16(xmxC, moveC);

	if (opt_res != NULL)
	{
		float offset = (float) dcops->wordoffset;
		int16_t res[16] __attribute__ ((aligned (16)));
		int16_t ovs[16] __attribute__ ((aligned (16)));

		memmove(res, &xmxC, sizeof(xmxC));
		memmove(ovs, &overflows, sizeof(overflows));

		for (i = 0; i < 16; i++)
			if (ovs[i])	opt_res[i] = eslINFINITY;	// signal overflow
			else		opt_res[i] = ((float) res[i] - offset) / dcops->scale - 2.0;
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
	DATA_COPSAVX *dcops;
} pthr_info_t;

void* viterbi_cops_thread_loop(void* argst)
{
	int i, execcount = 0;
	pthr_info_t* args = (pthr_info_t*) argst;
	DATA_COPSAVX *dcops = args->dcops;
	
#ifdef _GNU_SOURCE
	printf("THR %d running on cpu %d\n", args->thrid, sched_getcpu());
#endif

	for (i = 0; 1; i++)
	{	execcount++;
		while (dcops->synccontrol != execcount) sched_yield();
		viterbi_cops_word_partitioned(dcops, NULL, args->thrid);
	}

#ifdef _GNU_SOURCE
	printf("THR %d running on cpu %d\n", args->thrid, sched_getcpu());
#endif
	return (void*) 0;
}


// for UMA architectures
//#define MAPIDCPU(id) ((id) % 4)

// for NUMA architectures
#define MAPIDCPU(id) (((id)*4 + (id)/8)%32)

// for Aleph
//#define MAPIDCPU(id) (((id)*2 + 16)


void viterbi_cops_create_threads(DATA_COPSAVX* dcops)
{
	int i;

	dcops->synccontrol = 0;
	pthread_t threads[NTHREADS];
	pthr_info_t args;
	args.dcops = dcops; args.thrid = 0;

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
		if (pthread_create(&threads[i], &attr, viterbi_cops_thread_loop, argscopy))
			exit(fprintf(stderr, "ERROR could not create worker thread\n"));
	}

#ifdef _GNU_SOURCE
	printf("THR %d running on cpu %d\n", NTHREADS-1, sched_getcpu());
#endif
}



/***************************************************************************
 *	Prepare Viterbi COPS
 ***************************************************************************/

void free_cops_buffers(DATA_COPSAVX* dcops)
{
	int j; 
	for (j = 0; j < 16; j++)
		free(dcops->seqs[j]);
	for (j = 1; j <= dcops->Npartitions; j++)
		free(dcops->synchflags[j]);

	free(dcops->xmxE);
	free(dcops->pdcv);
	free(dcops->psc);
}

void alloc_cops_buffers(DATA_COPSAVX* dcops, int maxL)
{
	int j; 
	for (j = 0; j < 16; j++)
		dcops->seqs[j] = calloc(maxL+16, sizeof(ESL_DSQ));

	dcops->nflags = maxL;
	for (j = 1; j <= dcops->Npartitions; j++)
		dcops->synchflags[j] = alloc_aligned64(dcops->nflags);
	dcops->synchflags[0] = dcops->synchflags[dcops->Npartitions+1] = NULL;

	dcops->psc  = alloc_m128_aligned64((maxL));
	dcops->pdcv = alloc_m128_aligned64((maxL));
	dcops->xmxE = alloc_m128_aligned64((maxL+1));
}

DATA_COPSAVX* p7_ViterbiCOPS_Create(P7_PROFILE *gm)
{
	float *tsc= gm->tsc;
	__m256i *tt;
	int M = gm->M, Mtop, k, i;
	DATA_COPSAVX* dr = calloc(1, sizeof(DATA_COPSAVX));
	Mtop = dr->Mlinesize = roundtop(M,16);
	dr->tsc_all = tt = alloc_m128_aligned64((Mtop+1)*16);
	dr->tsc_all = tt;
	dr->scale	= 500.0 / eslCONST_LOG2;;
	dr->wordoffset = 12000;
	dr->gm = gm;

#define SETSCORE(d,val)	*d++ = SET1epi16(val);
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


void p7_ViterbiCOPS_Setup(DATA_COPSAVX* dcops, int maxL, int max_partition)
{
	int M = dcops->gm->M, nblocks;

	// choose highest possible no. of partitions, multiple of no. threads
	// use (double) to force an non-integer max over MaxPart
//	for (nblocks = NTHREADS; M / (double) nblocks > max_partition; nblocks += NTHREADS);
	// same as:
	nblocks = roundtop(ceil(M / (double) max_partition),NTHREADS);

	dcops->partition	= roundtop(ceil(M/(double)nblocks), 16);
	dcops->Npartitions= M/dcops->partition + (M % dcops->partition != 0);

	dcops->L = 0;
	dcops->allocL = maxL+16;
	alloc_cops_buffers(dcops, dcops->allocL);

	viterbi_cops_create_threads(dcops);
}


typedef struct _dsq_cmp
{
	ESL_DSQ *seq;
	int length;
} dsq_cmp_t;

int p7_ViterbiCOPS(DATA_COPSAVX* dcops, dsq_cmp_t **seqsdb, float* results)
{
	int j, maxL = 0;
	for (j = 0; j < 16; j++)
		if (seqsdb[j]->length > maxL)
			maxL = seqsdb[j]->length;

//	maxL = 2000;
//	for (j = 0; j < 8; j++)	printf("%d ", ilengths[j]); 	printf("\n");

	if (maxL > dcops->allocL-10)
	{	// realloc
		dcops->allocL = ESL_MAX(dcops->allocL*2, roundtop(maxL,1024));
		free_cops_buffers(dcops);
		alloc_cops_buffers(dcops, dcops->allocL);
	}

	if (dcops->L != maxL)
	{   dcops->L  = maxL;
		p7_ReconfigLength(dcops->gm, maxL);
	}

	for (j = 0; j < 16; j++)
		// must be copied since we run N seqs in parallel up to the max length of the longest
		memcpy(dcops->seqs[j], seqsdb[j]->seq, (seqsdb[j]->length+1)*sizeof(ESL_DSQ));
	
	// pad sequences with k = Alphasize
	byte alphSize = dcops->gm->abc->Kp;
	for (j = 0; j < 16; j++)
		memset(dcops->seqs[j]+seqsdb[j]->length+1, alphSize, maxL-seqsdb[j]->length);

	for (j = 1; j <= dcops->gm->M/dcops->partition; j++)
		memset((void*) dcops->synchflags[j], 0, dcops->nflags);

	dcops->synccontrol++;
	viterbi_cops_word_partitioned(dcops, results, NTHREADS-1	);
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
gcc -g -Wall -O3 -std=gnu99 -o viteravx -I. -I.. -L.. -I../impl_sse -L../impl_sse -I../../easel -L../../easel viterbi_avx.c viterbi_serial.c -lhmmer -lhmmerimpl -leasel -lm -lpthread -DNTHREADS=1 -DNROUNDS=2
*/

#include <esl_stopwatch.h>
#include <esl_sqio.h>
#include <esl_randomseq.h>

#include <sys/timeb.h>
#include <sched.h>

static ESL_OPTIONS options[] = {
	/* name		type	default		env	range toggles reqs incomp	help						docgroup*/
	{ "-h",	eslARG_NONE	, FALSE	, NULL, NULL , NULL, NULL, NULL, "show brief help on version and usage",0 },
	{ "-s",	eslARG_INT	, "42"	, NULL, NULL , NULL, NULL, NULL, "set random number seed to <n>",		0 },
	{ "-L",	eslARG_INT	, "2000"	, NULL, "n>0", NULL, NULL, NULL, "length of random target seqs",		0 },
	{ "-N",	eslARG_INT	, "500000",NULL, "n>0", NULL, NULL, NULL, "number of random target seqs",		0 },
	{ "-M", eslARG_INT	, "112" , NULL, "n>0", NULL, NULL, NULL, "maximum size for the partitions",		0 },
	{	 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage [] = "[-options] <hmmfile>";
static char banner[] = "benchmark driver for generic Viterbi";

#define SSTEST 0

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
	int			L		= esl_opt_GetInteger(go, "-L");
	int			N		= esl_opt_GetInteger(go, "-N") / 16;
	int			MaxPart	= esl_opt_GetInteger(go, "-M");
	__m128		resdata[10];
	float		*sc1 	= (float*)	(resdata+0);	// uses 1 __m128s
	ESL_DSQ		*dsq	= NULL;
	int			i, j;
	ESL_SQFILE   *sqfp	= NULL;
	DATA_COPSAVX *dcops;
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

	dcops = p7_ViterbiCOPS_Create(gm1);
	p7_ViterbiCOPS_Setup(dcops, L+100, MaxPart); // use max L
	dcops->L = L;

	// No. of partitions computed without full parallelism ( == no. of threads active while some are idle)
	int Niters_part	= dcops->Npartitions % NTHREADS;
	// No. of Model lines that could be computed but are wasted by idle threads waiting on the end
	int Nwasted_threads	= dcops->partition * ((NTHREADS-Niters_part) % NTHREADS);
	// No. of lines in the last partition that go beyond M. It's wasted comp time by a single thread
	int Nwasted_leftover= (dcops->partition - gm1->M % dcops->partition) % dcops->partition;
	// Total number of wasted lines
	int wastedcomp = Nwasted_threads + Nwasted_leftover;
	// Total number of lines computed and waited for
	int totalcomp = wastedcomp + gm1->M; // same as: roundtop(gm1->M, dcops->partition * NTHREADS);

	// for ViterbiFilter
#if SSTEST
	P7_OPROFILE *om = p7_oprofile_Create(hmm->M, gm1->abc);
	p7_oprofile_Convert(gm1, om);
	P7_OMX		*ox	= p7_omx_Create(hmm->M, 0, 0);
	double sumerrors = 0;
	float* results = (float*) alloc_m128_aligned64(N*4+2);
#endif

	dsq_cmp_t **seqsdb= calloc(16*N+64, sizeof(dsq_cmp_t*));

	if (esl_sqfile_OpenDigital(abc, seqfile, eslSQFILE_FASTA, NULL, &sqfp) == eslOK)
//	{	p7_Fail("Failed to open sequence file\n");	return -1; }
	{	// Use Sequence file
		ESL_SQ* sq = esl_sq_CreateDigital(abc);
		for (j = 0; j < 16*N; j++)
		{
			int res = esl_sqio_Read(sqfp, sq);
			if (res != eslOK) { /* printf("ATENCAO: faltam sequencias\n"); */ break; }
			int len = sq->n;
			dsq = sq->dsq;
			seqsdb[j] = malloc(sizeof(dsq_cmp_t));
			seqsdb[j]->length = len;
			seqsdb[j]->seq = malloc((len+4)*sizeof(ESL_DSQ));
			memcpy(seqsdb[j]->seq, dsq, len+2);
			sumlengths += len;
			esl_sq_Reuse(sq);
		}

		N = j/16;
		tprintf("N = %d\n", N);

		ftime(&tbstart);
		// Sort sequences by length
		qsort(seqsdb, N*16, sizeof(dsq_cmp_t*), compare_seqs);
	}
	else	// Generate random sequences
		for (i = 0; i < N; i++)
		{
			for (j = 0; j < 16; j++)
			{
				int len = L; // - rand()%1000;
				seqsdb[i*16+j] = malloc(sizeof(dsq_cmp_t));
				seqsdb[i*16+j]->seq = malloc(len+4);
				seqsdb[i*16+j]->length = len;
				esl_rsq_xfIID(r, bg->f, abc->K, len, seqsdb[i*16+j]->seq);
				sumlengths += len;
			}
			ftime(&tbstart);
		}

	printf("Viterbi COPS Word with %d Threads, model %s: Modelsize %d, #Segms: %d, SeqL: %d, Nseqs %d, Part %d, #Parts %d\n",
			NTHREADS, hmmfile, gm1->M, (int) ceil(gm1->M/16.0), L, 16*N, dcops->partition, dcops->Npartitions);
	printf("Total Comp Lines: %d | Wasted Comp Lines: %d\n", totalcomp, wastedcomp);


	for (j = 0; j < NROUNDS; j++) 
	for (i = 0; i < N; i++)
	{
	//	if (i % 10000 == 0) printf("START %d\n", i);

		p7_ViterbiCOPS(dcops, seqsdb+i*16, sc1);

		if (SSTEST) memcpy(results+i*16, sc1, 64);	// 32 bytes indeed! 16 floats
	}
	ftime(&tbend);

	double secs = TIMEDIFF(tbstart,tbend);
//	printf("Qsort time: %6.3f | Viterbi time: %6.3f\n", TIMEDIFF(tbqsort,tbstart), secs);
	w->elapsed = w->user = secs;
	esl_stopwatch_Display(stdout, w, "# Opt CPU time: ");
	double compmillioncells = NROUNDS * (double) sumlengths * (double) hmm->M * 1e-6;
	printf("# %.0fM cells in %.1f Mc/s\n", compmillioncells, compmillioncells / secs);

	printf("XXXXX\tViterbiCOPS\t%dthrs\t%s\t%.1f\n\n", NTHREADS, hmmfile, compmillioncells / secs);


#if SSTEST	// compare results against base version
	for (i = 0; i < 1000 && i < N; i++)
	{
		int maxll = 0; float sc2;
		for (j = 0; j < 16; j++)
			if (maxll < seqsdb[i*16+j]->length)
				maxll = seqsdb[i*16+j]->length;

//		for (j = 0; j < 16; j++)	printf("%d ", seqsdb[i*16+j]->length); printf("\n");
//		if (i % 10 == 0) printf("i %d\n", i);

		p7_oprofile_ReconfigRestLength(om, maxll);
		p7_ReconfigLength(gm2, maxll);	// fazer Reconfig aqui para emular compl o VitStream

		for (j = 0; j < 16; j++)
		{
//			p7_ReconfigLength(gm2, seqsdb[i*16+j]->length);
//			p7_Viterbi_unilocal(seqsdb[i*16+j]->seq, seqsdb[i*16+j]->length, gm2, &sc3);
//			p7_Viterbi_unilocal_word(seqsdb[i*16+j]->seq, seqsdb[i*16+j]->length, gm2, &sc2);

//			p7_oprofile_ReconfigLength(om, seqsdb[i*16+j]->length);
			p7_ViterbiFilter(seqsdb[i*16+j]->seq, seqsdb[i*16+j]->length, om, ox, &sc2);
			
			//sumerrors += fabs(sc1[j]- sc2);
			if (fabs(results[i*16+j] - sc2) > 0.00001)
			{	printf("%3d-%d L %4d: VS %f d %f\t| Base SerI %f\n", i, j, seqsdb[i*16+j]->length, 
						results[i*16+j], fabs(results[i*16+j] - sc2), sc2);
				getc(stdin);
			}
		}
	}
#endif
	return 0;
}

