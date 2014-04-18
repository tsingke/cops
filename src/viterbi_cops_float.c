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

#define _GNU_SOURCE
#include <xmmintrin.h>		/* SSE	*/
#include <emmintrin.h>		/* SSE2 */

#include <pthread.h>
#include <string.h>

typedef unsigned char byte;

#include "viterbi_cops.h"
#include "viterbi_serial.h"

extern int p7_ViterbiFilter(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_OMX *ox, float *ret_sc);

/***************************************************************************/

/***************************************************************************/

#ifndef MAX_PARTITION
	#define MAX_PARTITION	120
#endif

#define tprintf if(0) printf
#define dprintf if(0) printf


/***************************************************************************
 *	Auxiliary functions
 ***************************************************************************/

__m128* alloc_m128_aligned64(int nvecs)
{
	void *mem;
	posix_memalign(&mem, 64, (nvecs)*sizeof(__m128));
	return (__m128*) mem;
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

#if 0	// to see the annotated source in profiling
#define _mm_add_ps(__A,__B)		(__m128) __builtin_ia32_addps ((__v4sf)__A, (__v4sf)__B)
#define _mm_max_ps(__A,__B)		(__m128) __builtin_ia32_maxps ((__v4sf)__A, (__v4sf)__B)

#define _mm_and_si128(__A,__B)	(__m128)__builtin_ia32_pand128	((__v2di)__A, (__v2di)__B)
#define _mm_andnot_si128(__A,__B)(__m128)__builtin_ia32_pandn128((__v2di)__A, (__v2di)__B)
#define _mm_or_si128(__A,__B)	(__m128)__builtin_ia32_por128	((__v2di)__A, (__v2di)__B)
#define _mm_xor_si128(__A,__B)	(__m128)__builtin_ia32_pxor128	((__v2di)__A, (__v2di)__B)

#define _mm_and_ps(__A,__B)		__builtin_ia32_andps	(__A, __B)
#define _mm_andnot_ps(__A,__B)	__builtin_ia32_andnps	(__A, __B)
#define _mm_or_ps(__A,__B)		__builtin_ia32_orps	(__A, __B)
#define _mm_xor_ps(__A,__B)		__builtin_ia32_xorps	(__A, __B)

#define _mm_and_pd(__A,__B)		(__m128d)__builtin_ia32_andpd	((__v2df)__A, (__v2df)__B)
#define _mm_andnot_pd(__A,__B)	(__m128d)__builtin_ia32_andnpd	((__v2df)__A, (__v2df)__B)
#define _mm_or_pd(__A,__B)		(__m128d)__builtin_ia32_orpd	((__v2df)__A, (__v2df)__B)
#define _mm_xor_pd(__A,__B)		(__m128d)__builtin_ia32_xorpd	((__v2df)__A, (__v2df)__B)

#define _mm_unpacklo_ps(__A,__B)(__m128)__builtin_ia32_unpcklps((__v4sf)__A, (__v4sf)__B)
#define _mm_unpackhi_ps(__A,__B)(__m128)__builtin_ia32_unpckhps((__v4sf)__A, (__v4sf)__B)
#define _mm_unpacklo_pd(__A,__B)(__m128)__builtin_ia32_unpcklpd((__v2df)__A, (__v2df)__B)
#define _mm_unpackhi_pd(__A,__B)(__m128)__builtin_ia32_unpckhpd((__v2df)__A, (__v2df)__B)

void printreg(char* msg, __m128* reg)
{
	float tt[SSE32_NVALS] __attribute__ ((aligned (16)));
	memmove(tt, reg, 16);
	printf("%s: %6f %6f %6f %6f\n", msg, tt[0], tt[1], tt[2], tt[3]);
}

#endif


/***************************************************************************
 *	Viterbi algorithm
***************************************************************************/

int p7_ViterbiCOPSf_partitioned_threaded(DATA_COPS32* dcops, float* opt_res, int thrid, __m128* maxlengthsptr)
{
	int L = dcops->L;
	P7_PROFILE* gm = dcops->gm;
	ESL_DSQ** ddsq = dcops->seqs;
	int M = gm->M, i, k, v, t, j;
	const int PARTITION = dcops->partition;
	float** oprmsc = dcops->rsc_msc;
	__m128* xmxEv = dcops->xmxE;
	__m128 xmxB, xmxN, xmxE, xmxC, loopN, moveN, moveC, Vinf  = _mm_set1_ps(-eslINFINITY);

	__m128  dmx[PARTITION];
	__m128  mmx[PARTITION];
	__m128  imx[PARTITION];
	__m128	xmm[12];
	float *mscore[SSE32_NVALS];
	__m128 loopC, maxlengths;
	maxlengths = loopC = Vinf;

	if (thrid == NTHREADS-1) 
	{	
		loopC	= _mm_set1_ps(gm->xsc[p7P_C][p7P_LOOP]);	
		// lengths -> lengths -1, para usar GT em vez de GT+Eq
		maxlengths = _mm_add_ps(*maxlengthsptr, _mm_set1_ps(1));
	}

	t = ((dcops->Npartitions+thrid)%NTHREADS)*PARTITION;
	tprintf("START Viterbi Thr %d in %d L %d\n", thrid, t, L);

	moveC = _mm_set1_ps(gm->xsc[p7P_C][p7P_MOVE]);
	moveN = _mm_set1_ps(gm->xsc[p7P_N][p7P_MOVE]);
	loopN = _mm_set1_ps(gm->xsc[p7P_N][p7P_LOOP]);
	xmxC  = xmxN = Vinf;

	for (	; t < M; t += NTHREADS*PARTITION)
	{
		volatile byte* synchflags1 = dcops->synchflags[t/PARTITION];
		volatile byte* synchflags2 = dcops->synchflags[t/PARTITION+1];
		int t4 = t;	// do *not* divide by SSE32_NVALS, t/SSE32_NVALS, since we are using a float* pointer
	
		xmxN = _mm_xor_ps(xmxN,xmxN);
		xmxB = moveN;
		for (k = 0; k < PARTITION; k++)
			dmx[k] = mmx[k] = imx[k] = Vinf;

		for (i = 1; i <= L; i++)
		{
		//	tprintf("Iter Thr %d t %d: I %d\n", thrid, t, i);
			__m128 sc, dcv, temp, mpv, ipv, dpv;
			__m128 *ttsc = dcops->tsc_all + t*8;
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
			
			for (j = 0; j < SSE32_NVALS; j++)
				mscore[j] =  oprmsc[ddsq[j][i]] + t4;
			
			for (k = 0; k < PARTITION && t+k < M; )
			{
#define EMLOAD(i)	xmm[i+8] = _mm_load_ps(mscore[i]);
				EMLOAD(0) 	EMLOAD(1) 
				EMLOAD(2) 	EMLOAD(3) 

#define MIX32(i,r)	\
	xmm[r  ] = _mm_unpacklo_ps(xmm[i], xmm[i+1]);	\
	xmm[r+1] = _mm_unpackhi_ps(xmm[i], xmm[i+1]);	
#define MIX64(i,r)	\
	xmm[r  ] = (__m128) _mm_unpacklo_pd((__m128d) xmm[i], (__m128d) xmm[i+2]);	\
	xmm[r+1] = (__m128) _mm_unpackhi_pd((__m128d) xmm[i], (__m128d) xmm[i+2]);

				MIX32(8,0)	MIX32(10,2)
				MIX64(0,4)	MIX64(1,6)

#define TRIPLETCOMPUTE(k,j)			\
{	/* Calculate new M(k), delay store */			\
	sc = _mm_max_ps(sc, _mm_add_ps(xmxB,*ttsc++));	\
	sc = _mm_add_ps(sc, xmm[j]);	\
	/* Update E */					\
	xmxE = _mm_max_ps(xmxE, sc);	\
	\
	/* Pre-emptive load of M,D,I */	\
	dpv = dmx[k];	\
	ipv = imx[k];	\
	mpv = mmx[k];	\
	\
	/* Calculate current I(k) */	\
	temp = _mm_add_ps(mpv, *ttsc++);	\
	imx[k] = _mm_max_ps(temp, _mm_add_ps(ipv, *ttsc++));	\
	\
	/* Delayed stores of M and D */	\
	mmx[k] = sc;	\
	dmx[k] = dcv;	\
	\
	/* Calculate next D, D(k+1) */	\
	sc	= _mm_add_ps(sc, *ttsc++);	\
	dcv = _mm_max_ps(sc, _mm_add_ps(dcv, *ttsc++));	\
	\
	/* Pre-emptive partial calculation of M(k+1) */	\
	sc = _mm_add_ps(mpv, *ttsc++);\
	sc = _mm_max_ps(sc, _mm_add_ps(ipv, *ttsc++));\
	sc = _mm_max_ps(sc, _mm_add_ps(dpv, *ttsc++));\
	k++;			\
}
				TRIPLETCOMPUTE(k,4)
				TRIPLETCOMPUTE(k,5)
				TRIPLETCOMPUTE(k,6)
				TRIPLETCOMPUTE(k,7)

				mscore[0] += SSE32_NVALS; mscore[1] += SSE32_NVALS;
				mscore[2] += SSE32_NVALS; mscore[3] += SSE32_NVALS;
			}

			xmxN = _mm_add_ps(xmxN, loopN);
			xmxB = _mm_add_ps(xmxN, moveN);

			if (t+k < M)
			{	v = i-1;
				xmxEv[v] = xmxE;
				dcops->pdcv[v] = dcv;
				dcops->psc [v] = sc;

				if (NTHREADS > 1) synchflags2[v] = 1;
			}
			else	// executed only by main thread (NTHRS-1)
			{
				__m128 lengths = _mm_set1_ps(i);
				// Greater than or Equal, GtEq. Impl as (maxL +1) > L
				__m128 cmpsgteq= _mm_cmpgt_ps(maxlengths, lengths);
				__m128 tloopC	= _mm_and_ps(loopC, cmpsgteq); 	// nullify add of loopC if after the end

				xmxC	= _mm_add_ps(xmxC, tloopC); 
				xmxC	= _mm_max_ps(xmxC, xmxE);
			}
		}
	}

	xmxC = _mm_add_ps(xmxC, moveC);
	if (opt_res != NULL)
		memmove(opt_res, &xmxC, sizeof(xmxC));

	tprintf("END viterbi Thr %d - t %d\n", thrid, t);
	return eslOK;
}



/***************************************************************************
 *	Thread Management
 ***************************************************************************/

typedef struct _thr_info_t
{	int thrid;
	DATA_COPS32 *dcops;
} pthr_info_t;

void* p7_ViterbiCOPSf_thread_loop(void* argst)
{
	int i, execcount = 0;
	pthr_info_t* args = (pthr_info_t*) argst;
	DATA_COPS32 *dcops = args->dcops;

#ifdef _GNU_SOURCE
	printf("THR %d running on cpu %d\n", args->thrid, sched_getcpu());
#endif

	for (i = 0; 1; i++)
	{	execcount++;
		while (dcops->synccontrol != execcount) sched_yield();
		p7_ViterbiCOPSf_partitioned_threaded(dcops, NULL, args->thrid, NULL);
	}

#ifdef _GNU_SOURCE
	printf("THR %d running on cpu %d\n", args->thrid, sched_getcpu());
#endif
	return (void*) 0;
}


// Control the best mapping of threads for pinning

#define UMA		100
#define NUMA	101
#define ALEPH	102

#if	  ARCH == UMA
	#define MAPIDCPU(id) ((id) % 8)
#elif ARCH == NUMA
	#define MAPIDCPU(id) (((id)*4 + (id)/8) % 32)
#elif ARCH == ALEPH
	#define MAPIDCPU(id) ((id)*2 + 16)
#else
    #define MAPIDCPU(id) (id)
	#error "WRONG ARCH OR NO ARCH SPECIFIED"
#endif


void p7_ViterbiCOPSf_create_threads(DATA_COPS32* dcops)
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
		if (pthread_create(&threads[i], &attr, p7_ViterbiCOPSf_thread_loop, argscopy))
			exit(fprintf(stderr, "ERROR could not create worker thread\n"));
	}

#ifdef _GNU_SOURCE
	if (NTHREADS > 1) printf("THR %d running on cpu %d\n", NTHREADS-1, sched_getcpu());
#endif
}



/***************************************************************************
 *	Setup COPS
 ***************************************************************************/

void free_cops_buffers(DATA_COPS32* dcops)
{
	int j;
	for (j = 0; j < SSE32_NVALS; j++)
		free(dcops->seqs[j]);
    for (j = 1; j <= dcops->Npartitions; j++)
		free(dcops->synchflags[j]);

	free(dcops->xmxE);
	free(dcops->pdcv);
	free(dcops->psc);
}

void alloc_cops_buffers(DATA_COPS32* dcops, int maxL)
{
	int j;
	for (j = 0; j < SSE32_NVALS; j++)
		dcops->seqs[j] = calloc(maxL+16, sizeof(ESL_DSQ));

	dcops->nflags = maxL;
	for (j = 1; j <= dcops->Npartitions; j++)
		dcops->synchflags[j] = alloc_aligned64(dcops->nflags);
	dcops->synchflags[0] = dcops->synchflags[dcops->Npartitions+1] = NULL;

	dcops->psc  = alloc_m128_aligned64((maxL));
	dcops->pdcv = alloc_m128_aligned64((maxL));
	dcops->xmxE = alloc_m128_aligned64((maxL+1));
}

DATA_COPS32* p7_ViterbiCOPSf_Create(P7_PROFILE *gm)
{
	float *tsc= gm->tsc;
	__m128 *tt;
	int M = gm->M, Mtop, k, i;
	DATA_COPS32* dr = calloc(1, sizeof(DATA_COPS32));
	Mtop = roundtop(M,SSE32_NVALS);
	dr->tsc_all = tt = alloc_m128_aligned64((Mtop+1)*8);
	dr->gm = gm;

#define SETSCORE(d,val)	*d++ = _mm_set1_ps(val);
	for (k = 0; k < M; k++, tsc += p7P_NTRANS)
	{
		SETSCORE(tt, tsc[p7P_MM]);
		SETSCORE(tt, tsc[p7P_IM]);
		SETSCORE(tt, tsc[p7P_DM]); 
		SETSCORE(tt, tsc[p7P_BM]); 
		if (k < M-1) 
		{	SETSCORE(tt, tsc[p7P_NTRANS+p7P_MI]);
			SETSCORE(tt, tsc[p7P_NTRANS+p7P_II]);
			SETSCORE(tt, tsc[p7P_NTRANS+p7P_MD]);
			SETSCORE(tt, tsc[p7P_NTRANS+p7P_DD]);
		}
		else
		{	SETSCORE(tt, -eslINFINITY);	// para MI
			SETSCORE(tt, -eslINFINITY);	// para II
			SETSCORE(tt, -eslINFINITY);	// para MD
			SETSCORE(tt, -eslINFINITY);	// para DD
		}
	}

	// fill the left-over places with neutral values
	// +SSE32_NVALS to account for the spurious pre-emptive last loads of MM, IM, DM
	for (	; k < Mtop + 1; k++, tsc += p7P_NTRANS)
		for (i = 0; i < p7P_NTRANS; i++)
			SETSCORE(tt, -eslINFINITY);

	const int AlphaSize = gm->abc->Kp;
	dr->rsc_msc = calloc(AlphaSize+2, sizeof(int16_t*));
	float *data = alloc_aligned64(Mtop*(AlphaSize+1)*sizeof(float));

	for (k = 0; k < AlphaSize; k++)
	{
		float *rmsc = dr->rsc_msc[k] = data + k*Mtop;
		float *rsc = gm->rsc[k] + p7P_NR + p7P_MSC;	// ignore 1st values
		// #define MSC(k)   (rsc[(k)p7P_NR + p7P_MSC])
		for (i = 0; i< M; i++, rsc += p7P_NR)
			*rmsc++ = *rsc;
		// fill the left-over places with neutral values
		for (	; i < Mtop; i++)
			*rmsc++ = -eslINFINITY;
	}

	// create dummy values for k = AlphaSize, to use the padded sequences
	dr->rsc_msc[k] = data + (k)*(Mtop);
	for (i = 0; i< Mtop; i++)
		dr->rsc_msc[k][i] = -eslINFINITY;

	return dr;
}

int choose_partition_length(int M, int maxpart, int Nthreads)
{
	int nblocks;
    //	Method of computing partition lengths which gives always a multiple of Nthreads, until the minimum 8
/*	int i;
	for (i = maxpart; i > 8; i -= 8)
		if (((int) ceil(M / (double) i)) % Nthreads == 0)
			break;
*/
	// approximate method for the highest possible no. of partitions, in most cases multiple of no. threads
	// use (double) to force an non-integer max over MaxPart
//	for (nblocks = NTHREADS; M / (double) nblocks > maxpart; nblocks += NTHREADS);
	// same as:
	nblocks = roundtop(ceil(M / (double) maxpart), Nthreads);
	return roundtop(ceil(M/(double)nblocks), SSE32_NVALS);
}

void p7_ViterbiCOPSf_Setup(DATA_COPS32* dcops, int maxL, int max_partition)
{
	int M = dcops->gm->M;
    
	if (maxL <= 10) maxL = 1000;
	if (max_partition <= 8 || max_partition % 8 != 0)
		max_partition = MAX_PARTITION;

	dcops->partition	= choose_partition_length(M, max_partition, NTHREADS);
	dcops->Npartitions= M/dcops->partition + (M % dcops->partition != 0);

	dcops->L = 0;
	dcops->allocL = maxL+16;
	alloc_cops_buffers(dcops, dcops->allocL);

	p7_ViterbiCOPSf_create_threads(dcops);
}

int p7_ViterbiCOPSf_run(DATA_COPS32* dcops, SEQ **seqsdb, float* results)
{
	int j, maxL = 0;
	__m128 lengths;
	float *ilengths = (float*) &lengths;

	for (j = 0; j < SSE32_NVALS; j++)
	{	ilengths[j] = seqsdb[j]->length;
		if (ilengths[j] > maxL)
			maxL = ilengths[j];
	}

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

	for (j = 0; j < SSE32_NVALS; j++)
		// must be copied since we run N seqs in parallel up to the max length of the longest
		memcpy(dcops->seqs[j], seqsdb[j]->seq, (seqsdb[j]->length+1)*sizeof(ESL_DSQ));

	// pad sequences with k = Alphasize
	byte alphSize = dcops->gm->abc->Kp;
	for (j = 0; j < SSE32_NVALS; j++)
		memset(dcops->seqs[j]+seqsdb[j]->length+1, alphSize, maxL-seqsdb[j]->length);

	for (j = 1; j <= dcops->gm->M/dcops->partition; j++)
		memset((void*) dcops->synchflags[j], 0, dcops->nflags);

	dcops->synccontrol++;
	p7_ViterbiCOPSf_partitioned_threaded(dcops, results, NTHREADS-1, &lengths);
	return maxL;
}

int compare_seqs(const void *vs1, const void *vs2)
{
	SEQ *s1 = *((SEQ**) vs1);
	SEQ *s2 = *((SEQ**) vs2);
	if (s1->length == s2->length)
		return 0;
	return (s1->length > s2->length);
}




/****************************************************************
 *	Benchmark driver											*
 ****************************************************************/

/*
gcc -g -Wall -Wextra -O3 -std=gnu99 -o viterbicops-float -I. -I.. -L.. -I../../easel -L../../easel p7_ViterbiCOPSf_float.c viterbi_serial.c -lhmmer -leasel -lm -lpthread -DNTHREADS=1 -DARCH=UMA -DMAX_PARTITION=120 -DCOPS_FLOAT_BENCHMARK
*/

#ifdef COPS_FLOAT_BENCHMARK

#include <sys/timeb.h>
#include <limits.h>
#include <easel.h>
#include <esl_randomseq.h>
#include <esl_sqio.h>

#define STRHELPER(x) #x
#define STR(x) STRHELPER(x)

static ESL_OPTIONS options[] = {
	/* name			type		default	env	range toggles reqs incomp	help						docgroup*/
	{ "-h",	eslARG_NONE	, FALSE	, NULL, NULL , NULL, NULL, NULL, "show brief help on version and usage",0 },
  	{ "-c", eslARG_NONE	, FALSE	, NULL, NULL,  NULL, NULL, NULL, "compare scores to trusted serial implementation (debug)", 0 }, 
	{ "-s",	eslARG_INT	, "42"	, NULL, NULL , NULL, NULL, NULL, "set random number seed to <n>",		0 },
	{ "-L",	eslARG_INT	, "400"	, NULL, "n>0", NULL, NULL, NULL, "length of random target seqs",		0 },
	{ "-N",	eslARG_INT	, "8800", NULL, "n>0", NULL, NULL, NULL, "number of random target seqs",		0 },
	{ "-M", eslARG_INT	, STR(MAX_PARTITION), NULL, "n>0", NULL, NULL, NULL, "maximum size for the partitions",		0 },
	{ "-R", eslARG_INT	, "1"	, NULL, "n>0", NULL, NULL, NULL, "number of test rounds",				0 },
	{	0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage [] = "[-options] <hmmfile>";
static char banner[] = "benchmark driver for generic Viterbi";

#define TIMEDIFF(t1,t2) ((t2.time - t1.time) + (t2.millitm - t1.millitm)*0.001)

int main(int argc, char **argv)
{
	ESL_GETOPTS	*go		= esl_getopts_CreateDefaultApp(options, 2, argc, argv, banner, usage);
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
	int			N		= esl_opt_GetInteger(go, "-N") / SSE32_NVALS;
	int			MaxPart	= esl_opt_GetInteger(go, "-M");
	int			NROUNDS	= esl_opt_GetInteger(go, "-R");
	int			check	= esl_opt_GetBoolean(go, "-c");
	int			i, j;
	__m128		resdata[10];
	float 		*sc1	= (float*) resdata;
	ESL_SQFILE   *sqfp	= NULL;
	DATA_COPS32	*dcops;
	struct timeb tbstart, tbend;
	int sumlengths = 0;
	float* results = NULL;

	srand(time(NULL));
	if (p7_hmmfile_Open(hmmfile, NULL, &hfp) != eslOK) p7_Fail("Failed to open HMM file %s", hmmfile);
	if (p7_hmmfile_Read(hfp, &abc, &hmm)	 != eslOK) p7_Fail("Failed to read HMM");

	bg = p7_bg_Create(abc);
	p7_bg_SetLength(bg, L);
	gm1 = p7_profile_Create(hmm->M, abc);
	gm2 = p7_profile_Create(hmm->M, abc);
	p7_ProfileConfig(hmm, bg, gm1, L, p7_UNILOCAL);
	p7_ProfileConfig(hmm, bg, gm2, L, p7_UNILOCAL);

	dcops = p7_ViterbiCOPSf_Create(gm1);
	p7_ViterbiCOPSf_Setup(dcops, L+100, MaxPart); // use max L
	dcops->L = L;

    int dbsize = SSE32_NVALS*N;
	SEQ **seqsdb= calloc(dbsize+1, sizeof(SEQ*));
	int equallength = 1;

	if (esl_sqfile_OpenDigital(abc, seqfile, eslSQFILE_FASTA, NULL, &sqfp) == eslOK)
	{	// Use Sequence file
		ESL_SQ* sq = esl_sq_CreateDigital(abc);
		int maxseqs, len=0;

        if (esl_opt_IsDefault(go, "-N"))    // N not specified in cmdline
            maxseqs = INT_MAX;   // no limit
        else
            maxseqs = SSE32_NVALS*N;      // use cmdline limit

		for (j = 0; j < maxseqs && esl_sqio_Read(sqfp, sq) == eslOK; j++)
		{
		 	if (equallength && sq->n != len && j > 0)
                equallength = 0;
               
			len = sq->n;
			if (j > dbsize)
			{	seqsdb = realloc(seqsdb, 2*(dbsize+1)*sizeof(SEQ*));
				dbsize *= 2;
			}
            
			ESL_DSQ* dsq = sq->dsq;
			seqsdb[j] = malloc(sizeof(SEQ));
			seqsdb[j]->length = len;
			seqsdb[j]->seq = malloc((len+4)*sizeof(ESL_DSQ));
			memcpy(seqsdb[j]->seq, dsq, len+2);
			sumlengths += len;
			esl_sq_Reuse(sq);
		}
		N = j/SSE32_NVALS;
	}
    else	// Not found database. Generate random sequences
		for (i = 0; i < N; i++)
			for (j = 0; j < SSE32_NVALS; j++)
			{
				int len = L; // - rand()%1000;
				seqsdb[i*SSE32_NVALS+j] = malloc(sizeof(SEQ));
				seqsdb[i*SSE32_NVALS+j]->seq = malloc(len+4);
				seqsdb[i*SSE32_NVALS+j]->length = len;
				esl_rsq_xfIID(r, bg->f, abc->K, len, seqsdb[i*SSE32_NVALS+j]->seq);
				sumlengths += len;
			}

   	printf("Viterbi COPS float with %d threads, model %s. ModelLen: %d, #Segms: %d, SeqL.: %d, #seqs: %d, Partition: %d, #parts: %d\n",
			NTHREADS, hmmfile, gm1->M, (int) ceil(gm1->M/SSE32_NVALS), L, SSE32_NVALS*N*NROUNDS, dcops->partition, dcops->Npartitions);
            
/*	// No. of partitions computed without full parallelism ( == no. of threads active while some are idle)
	int Niters_part	= dcops->Npartitions % NTHREADS;
	// No. of Model lines that could be computed but are wasted by idle threads waiting on the end
	int Nwasted_threads	= dcops->partition * ((NTHREADS-Niters_part) % NTHREADS);
	// No. of lines in the last partition that go beyond M. It's wasted comp time by a single thread
	int Nwasted_leftover= (dcops->partition - gm1->M % dcops->partition) % dcops->partition;
	// Total number of wasted lines
	int wastedcomp = Nwasted_threads + Nwasted_leftover;
	// Total number of lines computed and waited for
	int totalcomp = wastedcomp + gm1->M; // same as: roundtop(gm1->M, dcops->partition * NTHREADS);
	printf("Total Comp Lines: %d | Wasted Comp Lines: %d\n", totalcomp, wastedcomp);
*/

	if (check) results = (float*) alloc_m128_aligned64((N+1)*2);

	ftime(&tbstart);
    
	if (!equallength)
	{	// Sort sequences by length
		qsort(seqsdb, N*SSE32_NVALS, sizeof(SEQ*), compare_seqs);
	}

	for (j = 0; j < NROUNDS; j++)
        for (i = 0; i < N; i++)
        {
        //	if (i % 1000 == 0) printf("Seq %d\n", i);

            p7_ViterbiCOPSf_run(dcops, seqsdb+i*SSE32_NVALS, sc1);
            
			if (check) memcpy(results+i*SSE32_NVALS, sc1, 16);	// 16 bytes, SSE32_NVALS floats
        }
        
	ftime(&tbend);

	double secs = TIMEDIFF(tbstart,tbend);
	w->elapsed = w->user = secs;
	esl_stopwatch_Display(stdout, w, "# Opt CPU time: ");
	double compmillioncells = NROUNDS * (double) sumlengths * (double) hmm->M * 1e-6;
	printf("# %.0fM cells in %.1f Mc/s\n", compmillioncells, compmillioncells / secs);

	if (check)
	{	printf("Compare results against base version\n");

        for (i = 0; i < N; i++)
        {
            float sc2;
            int m1 = ESL_MAX(seqsdb[i*SSE32_NVALS+0]->length, seqsdb[i*SSE32_NVALS+1]->length);
            int m2 = ESL_MAX(seqsdb[i*SSE32_NVALS+2]->length, seqsdb[i*SSE32_NVALS+3]->length);
            int maxll = ESL_MAX(m1,m2);
            p7_ReconfigLength(gm2, maxll);

            for (j = 0; j < SSE32_NVALS; j++)
            {
    //			p7_ReconfigLength(gm2, seqsdb[i*SSE32_NVALS+j]->length);
                p7_Viterbi_unilocal(seqsdb[i*SSE32_NVALS+j]->seq, seqsdb[i*SSE32_NVALS+j]->length, gm2, &sc2);
        
                if (results[i*SSE32_NVALS+j] != sc2)
                {	printf("Seq %d Len %4d: %f - %f\tdiff: %f\n", i*SSE32_NVALS+j, seqsdb[i*SSE32_NVALS+j]->length, 
							results[i*SSE32_NVALS+j], sc2, fabs(results[i*SSE32_NVALS+j] - sc2));
                }
            }
        }
    }

	return 0;
}

#endif

