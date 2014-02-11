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

#include <string.h>

typedef unsigned char byte;

#include "viterbi_cops.h"
#include "viterbi_serial.h"

extern int p7_ViterbiFilter(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_OMX *ox, float *ret_sc);


/***************************************************************************/

#ifndef MAX_PARTITION
	#define MAX_PARTITION	120
#endif

#define WORDMAX 32767	// SHRT_MAX
#define tprintf if(0) printf
#define dprintf if(0) printf


/***************************************************************************
 *	Auxiliary functions
 ***************************************************************************/

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

// To show up in the annotated source in profiling
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
	int16_t tt[SSE16_NVALS] __attribute__ ((aligned (16)));
	memmove(tt, reg, 16);
	printf("%s: %6d %6d %6d %6d %6d %6d %6d %6d\n", msg,
			tt[0], tt[1], tt[2], tt[3], tt[4], tt[5], tt[6], tt[7]);
}
#endif



/***************************************************************************
 *	Viterbi algorithm
***************************************************************************/

int p7_ViterbiCOPSw_inlined(DATA_COPS16* dcops, float* opt_res)
{
	int L = dcops->L;
	P7_PROFILE* gm = dcops->gm;
	ESL_DSQ** ddsq = dcops->seqs;
	int i, k, j;
	const int Mtop = roundtop(gm->M,SSE16_NVALS);
	__m128i** oprmsc = (__m128i**) dcops->rsc_msc;
	__m128i xmxB, xmxE, xmxC, moveC, Vinf = _mm_set1_epi16(-WORDMAX);

	__m128i dmx[Mtop];
	__m128i mmx[Mtop];
	__m128i imx[Mtop];
	__m128i xmm[24];
	__m128i *mscore[SSE16_NVALS];
	__m128i overflowlimit, overflows;

	overflowlimit = _mm_set1_epi16(WORDMAX-1);
	overflows= _mm_xor_si128(overflows,overflows);	// zero out

	xmxC  = Vinf;
	moveC = _mm_set1_epi16(discretize(dcops->scale, gm->xsc[p7P_C][p7P_MOVE]));
	xmxB  = _mm_set1_epi16(dcops->wordoffset + discretize(dcops->scale, gm->xsc[p7P_N][p7P_MOVE]));

	for (k = 0; k < Mtop; k++)
		dmx[k] = mmx[k] = imx[k] = Vinf;

	for (i = 1; i <= L; i++)
	{
		__m128i sc, dcv, temp, mpv, ipv, dpv;
		__m128i *ttsc = dcops->tsc_all;
		ttsc += 3;

		xmxE = mpv = dpv = ipv = sc = dcv = Vinf;

		for (j = 0; j < SSE16_NVALS; j++)
			mscore[j] = oprmsc[ddsq[j][i]];

		for (k = 0; k < Mtop; )
		{
#define MMLOAD(a,b)		\
xmm[a] = _mm_unpacklo_epi16(*mscore[a], *mscore[b]);	\
xmm[b] = _mm_unpackhi_epi16(*mscore[a], *mscore[b]);

			MMLOAD(0,1)	MMLOAD(2,3)
			MMLOAD(4,5)	MMLOAD(6,7)

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

		__m128i overfs = _mm_cmpgt_epi16(xmxE, overflowlimit);
		overflows = _mm_or_si128(overflows, overfs);	// select the overflowed channels
		xmxC	= _mm_max_epi16(xmxC, xmxE);
	}

	xmxC = _mm_adds_epi16(xmxC, moveC);

	if (opt_res != NULL)
	{
		float offset = (float) dcops->wordoffset;
		int16_t res[SSE16_NVALS] __attribute__ ((aligned (16)));
		int16_t ovs[SSE16_NVALS] __attribute__ ((aligned (16)));

		memmove(res, &xmxC, sizeof(xmxC));
		memmove(ovs, &overflows, sizeof(overflows));

		for (i = 0; i < SSE16_NVALS; i++)
			if (ovs[i])	opt_res[i] = eslINFINITY;	// signal overflow
			else		opt_res[i] = ((float) res[i] - offset) / dcops->scale - 2.0;
												// 2.0 nat approximation, UNILOCAL mode
	}

	return eslOK;
}




/***************************************************************************
 *	Setup COPS
 ***************************************************************************/

void free_cops_buffers(DATA_COPS16* dcops)
{
	int j; 
	for (j = 0; j < SSE16_NVALS; j++)
		free(dcops->seqs[j]);
}

void alloc_cops_buffers(DATA_COPS16* dcops, int maxL)
{
	int j; 
	for (j = 0; j < SSE16_NVALS; j++)
		dcops->seqs[j] = calloc(maxL+16, sizeof(ESL_DSQ));
}

DATA_COPS16* p7_ViterbiCOPSw_Create(P7_PROFILE *gm)
{
	float *tsc= gm->tsc;
	__m128i *tt;
	int M = gm->M, Mtop, k, i;
	DATA_COPS16* dr = calloc(1, sizeof(DATA_COPS16));
	Mtop = roundtop(M,SSE16_NVALS);
	dr->tsc_all = tt = alloc_m128_aligned64((Mtop+1)*SSE16_NVALS);
	dr->tsc_all = tt;
	dr->scale	= 500.0 / eslCONST_LOG2;;
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
   	// +4 to account for the spurious pre-emptive last loads of MM, IM, DM
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

void p7_ViterbiCOPSW_Setup(DATA_COPS16* dcops, int maxL)
{
	dcops->L = 0;
	dcops->allocL = maxL+16;
	alloc_cops_buffers(dcops, dcops->allocL);
}

int p7_ViterbiCOPSw_run(DATA_COPS16* dcops, SEQ **seqsdb, float* results)
{
	int j, maxL = 0;
	for (j = 0; j < SSE16_NVALS; j++)
		if (seqsdb[j]->length > maxL)
			maxL = seqsdb[j]->length;

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

	for (j = 0; j < SSE16_NVALS; j++)
		// must be copied since we run N seqs in parallel up to the max length of the longest
		memcpy(dcops->seqs[j], seqsdb[j]->seq, (seqsdb[j]->length+1)*sizeof(ESL_DSQ));
	
	// pad sequences with k = Alphasize
	byte alphSize = dcops->gm->abc->Kp;
	for (j = 0; j < SSE16_NVALS; j++)
		memset(dcops->seqs[j]+seqsdb[j]->length+1, alphSize, maxL-seqsdb[j]->length);

	p7_ViterbiCOPSw_inlined(dcops, results);
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
gcc -g -Wall -Wextra -O3 -std=gnu99 -o viterbicops -I. -I.. -L.. -I../../easel -L../../easel viterbi_cops.c viterbi_serial.c -lhmmer -leasel -lm -DARCH=UMA -DNTHREADS=1 -DMAX_PARTITION=120 
*/

#include <string.h>
#include <sys/timeb.h>
#include <easel.h>
#include <esl_randomseq.h>
#include <esl_sqio.h>

#define STRHELPER(x) #x
#define STR(x) STRHELPER(x)

static ESL_OPTIONS options[] = {
	/* name		type	default		env	range toggles reqs incomp	help						docgroup*/
	{ "-h",	eslARG_NONE	, FALSE	, NULL, NULL , NULL, NULL, NULL, "show brief help on version and usage",0 },
    { "-c", eslARG_NONE	, FALSE	, NULL, NULL , NULL, NULL, NULL, "compare scores to HMMER's VitFilter implementation (debug)", 0 }, 
    { "-s",	eslARG_INT	, "42"	, NULL, NULL , NULL, NULL, NULL, "set random number seed to <n>",		0 },
	{ "-L",	eslARG_INT	, "400" , NULL, "n>0", NULL, NULL, NULL, "length of random target seqs",	    0 },
	{ "-N",	eslARG_INT	, "8800", NULL, "n>0", NULL, NULL, NULL, "number of random target seqs",		0 },
	{ "-R", eslARG_INT	, "1"	, NULL, "n>0", NULL, NULL, NULL, "number of test rounds",				0 },
	{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
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
	int			L		= esl_opt_GetInteger(go, "-L");
	int			N		= esl_opt_GetInteger(go, "-N") / SSE16_NVALS;
	int			NROUNDS	= esl_opt_GetInteger(go, "-R");
   	int			check	= esl_opt_GetBoolean(go, "-c");
	__m128		resdata[10];
	int			i, j;
	float		*sc1 	= (float*) resdata;
	ESL_SQFILE   *sqfp	= NULL;
	DATA_COPS16 *dcops;
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

	dcops = p7_ViterbiCOPSw_Create(gm1);
	p7_ViterbiCOPSW_Setup(dcops, L+100); // use max L
	dcops->L = L;

    int dbsize = SSE16_NVALS*N;
	SEQ **seqsdb= calloc(dbsize+1, sizeof(SEQ*));
	int equallength = 1;

	if (esl_sqfile_OpenDigital(abc, seqfile, eslSQFILE_FASTA, NULL, &sqfp) == eslOK)
	{	// Use Sequence file
		ESL_SQ* sq = esl_sq_CreateDigital(abc);
        int maxseqs, len=0;
        
        if (esl_opt_IsDefault(go, "-N"))    // N not specified in cmdline
            maxseqs = (1<<31);   // no limit
        else
            maxseqs = SSE16_NVALS*N;      // use cmdline limit

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
		N = j/SSE16_NVALS;
	}
    else	// Not found database. Generate random sequences
        for (i = 0; i < N; i++)
			for (j = 0; j < SSE16_NVALS; j++)
			{
				int len = L; // - rand()%1000;
				seqsdb[i*SSE16_NVALS+j] = malloc(sizeof(SEQ));
				seqsdb[i*SSE16_NVALS+j]->seq = malloc(len+4);
				seqsdb[i*SSE16_NVALS+j]->length = len;
				esl_rsq_xfIID(r, bg->f, abc->K, len, seqsdb[i*SSE16_NVALS+j]->seq);
				sumlengths += len;
			}

   	printf("Viterbi COPS Word inlined, model %s. ModelLen: %d, #Segms: %d, SeqL.: %d, #seqs: %d\n",
			hmmfile, gm1->M, (int) ceil(gm1->M/SSE16_NVALS), L, SSE16_NVALS*N*NROUNDS);
            
   	if (check) results = (float*) alloc_m128_aligned64((N+1)*2);

	ftime(&tbstart);
	
    if (!equallength)
	{	// Sort sequences by length
		qsort(seqsdb, N*SSE16_NVALS, sizeof(SEQ*), compare_seqs);
	}

	for (j = 0; j < NROUNDS; j++) 
        for (i = 0; i < N; i++)
        {
        //	if (i % 1000 == 0) printf("Seq %d\n", i);

            p7_ViterbiCOPSw_run(dcops, seqsdb+i*SSE16_NVALS, sc1);

         	if (check) memcpy(results+i*SSE16_NVALS, sc1, 32);	// 32 bytes indeed! SSE16_NVALS floats
        }
    
	ftime(&tbend);

	double secs = TIMEDIFF(tbstart,tbend);
	w->elapsed = w->user = secs;
	esl_stopwatch_Display(stdout, w, "# Opt CPU time: ");
	double compmillioncells = NROUNDS * (double) sumlengths * (double) hmm->M * 1e-6;
	printf("# %.0fM cells in %.1f Mc/s\n", compmillioncells, compmillioncells / secs);

	if (check)
    {   P7_OPROFILE *om = p7_oprofile_Create(hmm->M, gm1->abc);
		p7_oprofile_Convert(gm1, om);
		P7_OMX		*ox	= p7_omx_Create(hmm->M, 0, 0);
		printf("Compare results against base version\n");

        for (i = 0; i < N; i++)
        {
            int maxll = 0; float sc2;
            for (j = 0; j < SSE16_NVALS; j++)
                if (maxll < seqsdb[i*SSE16_NVALS+j]->length)
                    maxll = seqsdb[i*SSE16_NVALS+j]->length;

            p7_oprofile_ReconfigRestLength(om, maxll);
    //      p7_ReconfigLength(gm2, maxll);	// emulate the lock-step inter-sequence reconfigs

            for (j = 0; j < SSE16_NVALS; j++)
            {
    //			p7_ReconfigLength(gm2, seqsdb[i*SSE16_NVALS+j]->length);
    //			p7_Viterbi_unilocal(seqsdb[i*SSE16_NVALS+j]->seq, seqsdb[i*SSE16_NVALS+j]->length, gm2, &sc3);
    //			p7_Viterbi_unilocal_word(seqsdb[i*SSE16_NVALS+j]->seq, seqsdb[i*SSE16_NVALS+j]->length, gm2, &sc2);

    //			p7_oprofile_ReconfigLength(om, seqsdb[i*SSE16_NVALS+j]->length);
                p7_ViterbiFilter(seqsdb[i*SSE16_NVALS+j]->seq, seqsdb[i*SSE16_NVALS+j]->length, om, ox, &sc2);
   				sc2 += 1.0;	// -2.0nat optimization, Local to Unilocal mode		
  
                if (fabs(results[i*SSE16_NVALS+j] - sc2) > 0.0001)
                {	printf("Seq %d Len %4d: %f - %f\tdiff: %f\n", i*SSE16_NVALS+j, seqsdb[i*SSE16_NVALS+j]->length, 
							results[i*SSE16_NVALS+j], sc2, fabs(results[i*SSE16_NVALS+j] - sc2));
                }
            }
        }
    }
    
	return 0;
}


