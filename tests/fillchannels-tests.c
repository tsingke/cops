#include <xmmintrin.h>		/* SSE	*/
#include <emmintrin.h>		/* SSE2 */

#include "hmmer.h"

#include <../viterbi_stream/viterbi_serial.h>

//static int ntotal = 0, nlazy = 0, nloops= 0;


int viterbi_serial_opt(ESL_DSQ *dsq, int L, P7_PROFILE *gm,  float *opt_sc)
{
	float *tsc= gm->tsc;
	int M = gm->M;
	int	i,k;
	
	// float esc = p7_profile_IsLocal(gm) ? 0 : -eslINFINITY;
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
		float *rsc = gm->rsc[dsq[i]];
		float sc, ipv, mpv, dpv, dcv;

		dcv = mpv = dpv = ipv = sc = -eslINFINITY;
		imx[0] = dmx[0] = mmx[0] = -eslINFINITY;
		xmxE = -eslINFINITY;
	
		for (k = 1; k <= M; k++)
		{
			dcv	= ESL_MAX(sc + TSC(p7P_MD,k-1), dcv + TSC(p7P_DD,k-1));

			/* match state. Sc == new MMX(k) */
			sc	= xmxB + TSC(p7P_BM,k-1);
			sc	= ESL_MAX(sc, mpv + TSC(p7P_MM,k-1));
			sc	= ESL_MAX(sc, ipv + TSC(p7P_IM,k-1));
			sc	= ESL_MAX(sc, dpv + TSC(p7P_DM,k-1));
			sc	= sc + MSC(k);

			xmxE	= ESL_MAX(xmxE, sc); // + esc);

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
		 * remember, N, C and J emissions are zero score by definition. */
		xmxJ	= ESL_MAX(xmxJ + loopJ, xmxE + loopE);	/* E->J is E's "loop" */
		xmxC	= ESL_MAX(xmxC + loopC, xmxE + moveE);
		xmxN	= xmxN		   + loopN;
		xmxB	= ESL_MAX(xmxN + moveN, xmxJ + moveJ);	/* J->B is J's move */
	}
	
	/* T state (not stored) */
	if (opt_sc != NULL)
		opt_sc[0] = xmxC + moveC;

	return eslOK;
}


#define unit	__m128
#define uniti	__m128i


typedef struct _data_interleaved
{	// data TSC
	float* tsc_mm;
	float* tsc_im;
	float* tsc_dm;
	float* tsc_bm;
	float* tsc_md;
	float* tsc_dd;
	float* tsc_mi;
	float* tsc_ii;

	unit* tsc_all;
	unit* buffer;	
	float** rsc_isc;
	float** rsc_msc;
	unit* risc;
	unit* rmsc;
} dpopt;


#define SET4MOV(dest,val)	{ dest[0] = dest[1] = dest[2] = dest[3] = val; dest += 4; }

float* aligned_alloc(int nfloats)
{
	float *data = malloc((nfloats+4)*sizeof(float));
	size_t idata = (size_t) data;
	if (idata % 16 == 0) return data;
	else	return (float*) ((idata & ~0xf) + 0x10);
}

float* aligned_alloc64(int nfloats)
{
	float *data = malloc((nfloats+16)*sizeof(float));
	size_t idata = (size_t) data;
	if (idata % 64 == 0) return data;
	else	return (float*) ((idata & ~0x3f) + 0x40);
}

dpopt* create_data_structs(P7_PROFILE *gm)
{
	float *tsc= gm->tsc, *tt;
	int M = gm->M, k, i;
//	float* mm, *im, *dm, *bm, *md, *dd, *mi, *ii;
	dpopt* dr = calloc(1, sizeof(dpopt));
	
	dr->buffer = (unit*) aligned_alloc64((M)*4*3);

	tt = aligned_alloc64(4*M*8);
	dr->tsc_all = (unit*) tt;
	//	#define TSC(s,k) (tsc[(k) * p7P_NTRANS + (s)])
	for (k = 0; k < M; k++, tsc += p7P_NTRANS)
	{
		SET4MOV(tt, tsc[p7P_MD]);
		SET4MOV(tt, tsc[p7P_DD]);
		SET4MOV(tt, tsc[p7P_BM]); 
		SET4MOV(tt, tsc[p7P_MM]);
		SET4MOV(tt, tsc[p7P_IM]);
		SET4MOV(tt, tsc[p7P_DM]); 
		SET4MOV(tt, tsc[p7P_NTRANS+p7P_MI]);
		SET4MOV(tt, tsc[p7P_NTRANS+p7P_II]);
	}

	dr->risc = (unit*) aligned_alloc64((2*M+8)*4);
	dr->rmsc = (unit*) aligned_alloc64((2*M+8)*4);
	dr->rsc_isc = (float**) aligned_alloc64(4*6*2);
	dr->rsc_msc = (float**) aligned_alloc64(4*6*2);

	for (k = 0; k < 20; k++)
	{
		float *risc = dr->rsc_isc[k] = aligned_alloc64((M+1));
		float *rmsc = dr->rsc_msc[k] = aligned_alloc64((M+1));
		float *rsc = gm->rsc[k] + p7P_NR;	// ignore 1st values

		for (i = 0; i<= M; i++, rsc += p7P_NR)
		{	*risc = rsc[p7P_ISC];
			*rmsc = rsc[p7P_MSC];
			risc++; rmsc++;
		}
	}
	return dr;
}


void fill_channels_separate(ESL_DSQ *ddsq[4], int M, int i, dpopt* opdata)
{
	int k;
	uniti xmm0, xmm1, xmm2, xmm3, xmm4, xmm5, xmm6, xmm7, xmm8, xmm9, xmm10, xmm11;
	uniti xii0, xii1, xii2, xii3, xii4, xii5, xii6, xii7, xii8, xii9, xii10, xii11;
	uniti *vrisci = (uniti*) opdata->risc;
	uniti *vrmsci = (uniti*) opdata->rmsc;
		
	for (k = 0; k*4 < M; k++, vrmsci += 4,	vrisci += 4)
	{
		xmm0 = * (((uniti*) opdata->rsc_msc[ddsq[0][i]])+k);
		xmm1 = * (((uniti*) opdata->rsc_msc[ddsq[1][i]])+k);
		xmm2 = * (((uniti*) opdata->rsc_msc[ddsq[2][i]])+k);
		xmm3 = * (((uniti*) opdata->rsc_msc[ddsq[3][i]])+k);
		xii0 = * (((uniti*) opdata->rsc_isc[ddsq[0][i]])+k);
		xii1 = * (((uniti*) opdata->rsc_isc[ddsq[1][i]])+k);
		xii2 = * (((uniti*) opdata->rsc_isc[ddsq[2][i]])+k);
		xii3 = * (((uniti*) opdata->rsc_isc[ddsq[3][i]])+k);

		xmm4 = _mm_unpacklo_epi32(xmm0, xmm1);
		xmm5 = _mm_unpackhi_epi32(xmm0, xmm1);
		xmm6 = _mm_unpacklo_epi32(xmm2, xmm3);
		xmm7 = _mm_unpackhi_epi32(xmm2, xmm3);
		xii4 = _mm_unpacklo_epi32(xii0, xii1);
		xii5 = _mm_unpackhi_epi32(xii0, xii1);
		xii6 = _mm_unpacklo_epi32(xii2, xii3);
		xii7 = _mm_unpackhi_epi32(xii2, xii3);

		xmm8 = _mm_unpacklo_epi64(xmm4, xmm6);
		xmm9 = _mm_unpackhi_epi64(xmm4, xmm6);
		xmm10= _mm_unpacklo_epi64(xmm5, xmm7);
		xmm11= _mm_unpackhi_epi64(xmm5, xmm7);
		xii8 = _mm_unpacklo_epi64(xii4, xii6);
		xii9 = _mm_unpackhi_epi64(xii4, xii6);
		xii10= _mm_unpacklo_epi64(xii5, xii7);
		xii11= _mm_unpackhi_epi64(xii5, xii7);

		vrmsci[0] = xmm8;	vrmsci[1] = xmm9;
		vrmsci[2] = xmm10;	vrmsci[3] = xmm11;
		vrisci[0] = xii8;	vrisci[1] = xii9;
		vrisci[2] = xii10;	vrisci[3] = xii11;
	}
}


void fill_channels_4shifts(ESL_DSQ *ddsq[4], int M, int i, dpopt* opdata)
{
	int k;
	uniti xmm0, xmm1, xmm2, xmm3; // xmm4, xmm5, xmm6, xmm7, xmm8, xmm9, xmm10, xmm11;
	uniti xii0, xii1, xii2, xii3; // xii4, xii5, xii6, xii7, xii8, xii9, xii10, xii11;
	uniti *vrisci = (uniti*) (opdata->risc); // = aligned_alloc(8*(M+1)));
	uniti *vrmsci = (uniti*) opdata->rmsc;
		
	for (k = 0; k < M; k++, vrisci++, vrmsci++)
	{
		xmm0 = * (((uniti*) opdata->rsc_msc[ddsq[0][i]])+k);
		xmm1 = * (((uniti*) opdata->rsc_msc[ddsq[1][i]])+k);
		xmm2 = * (((uniti*) opdata->rsc_msc[ddsq[2][i]])+k);
		xmm3 = * (((uniti*) opdata->rsc_msc[ddsq[3][i]])+k);
		xii0 = * (((uniti*) opdata->rsc_isc[ddsq[0][i]])+k);
		xii1 = * (((uniti*) opdata->rsc_isc[ddsq[1][i]])+k);
		xii2 = * (((uniti*) opdata->rsc_isc[ddsq[2][i]])+k);
		xii3 = * (((uniti*) opdata->rsc_isc[ddsq[3][i]])+k);

		xmm1 = _mm_slli_si128(xmm1, 4 );
		xmm2 = _mm_slli_si128(xmm2, 8 );
		xmm3 = _mm_slli_si128(xmm3, 12);
		xii1 = _mm_slli_si128(xii1, 4 );
		xii2 = _mm_slli_si128(xii2, 8 );
		xii3 = _mm_slli_si128(xii3, 12);

		xmm0 = _mm_or_si128(xmm0, xmm1);
		xmm2 = _mm_or_si128(xmm2, xmm3);
		xii0 = _mm_or_si128(xii0, xii1);
		xii2 = _mm_or_si128(xii2, xii3);
		vrmsci[0] = _mm_or_si128(xmm0, xmm2);
		vrisci[0] = _mm_or_si128(xii0, xii2);
	}
}



void fill_channels_ps(ESL_DSQ *seqs[4], int M, int iResidue, dpopt* opdata)
{
	int k;
	unit xmm0, xmm1, xmm2, xmm3, xmm4, xmm5, xmm6, xmm7, xmm8, xmm9, xmm10, xmm11;
	uniti *Vrmsci = (__m128i*) opdata->rmsc;
	float **rsc = opdata->rsc_msc;

	for (k = 0; k < M; k += 4, Vrmsci += 4)
	{
		xmm0 = _mm_load_ps(rsc[seqs[0][iResidue]]+k);
		xmm1 = _mm_load_ps(rsc[seqs[1][iResidue]]+k);
		xmm2 = _mm_load_ps(rsc[seqs[2][iResidue]]+k);
		xmm3 = _mm_load_ps(rsc[seqs[3][iResidue]]+k);

		xmm4 = _mm_unpacklo_ps(xmm0, xmm1);
		xmm5 = _mm_unpackhi_ps(xmm0, xmm1);
		xmm6 = _mm_unpacklo_ps(xmm2, xmm3);
		xmm7 = _mm_unpackhi_ps(xmm2, xmm3);

		xmm8  = (__m128) _mm_unpacklo_pd((__m128d) xmm4, (__m128d) xmm6);
		xmm9  = (__m128) _mm_unpackhi_pd((__m128d) xmm4, (__m128d) xmm6);
		xmm10 = (__m128) _mm_unpacklo_pd((__m128d) xmm5, (__m128d) xmm7);
		xmm11 = (__m128) _mm_unpackhi_pd((__m128d) xmm5, (__m128d) xmm7);

		_mm_store_si128(Vrmsci+0, (__m128i) xmm8);
		_mm_store_si128(Vrmsci+1, (__m128i) xmm9);
		_mm_store_si128(Vrmsci+2, (__m128i) xmm10);
		_mm_store_si128(Vrmsci+3, (__m128i) xmm11);
	}
}

void fill_channels_epi32(ESL_DSQ *seqs[4], int M, int iResidue, dpopt* opdata)
{
	int k;
	uniti xmm0, xmm1, xmm2, xmm3, xmm4, xmm5, xmm6, xmm7, xmm8, xmm9, xmm10, xmm11;
	uniti *Vrmsci = (__m128i*) opdata->rmsc;
	float **rsc = opdata->rsc_msc;

	for (k = 0; k < M; k += 4, Vrmsci += 4)
	{
		xmm0 = (__m128i) _mm_load_ps(rsc[seqs[0][iResidue]]+k);
		xmm1 = (__m128i) _mm_load_ps(rsc[seqs[1][iResidue]]+k);
		xmm2 = (__m128i) _mm_load_ps(rsc[seqs[2][iResidue]]+k);
		xmm3 = (__m128i) _mm_load_ps(rsc[seqs[3][iResidue]]+k);

		xmm4 = _mm_unpacklo_epi32(xmm0, xmm1);
		xmm5 = _mm_unpackhi_epi32(xmm0, xmm1);
		xmm6 = _mm_unpacklo_epi32(xmm2, xmm3);
		xmm7 = _mm_unpackhi_epi32(xmm2, xmm3);

		xmm8 = _mm_unpacklo_epi64(xmm4, xmm6);
		xmm9 = _mm_unpackhi_epi64(xmm4, xmm6);
		xmm10= _mm_unpacklo_epi64(xmm5, xmm7);
		xmm11= _mm_unpackhi_epi64(xmm5, xmm7);

		_mm_store_si128(Vrmsci+0, xmm8);
		_mm_store_si128(Vrmsci+1, xmm9);
		_mm_store_si128(Vrmsci+2, xmm10);
		_mm_store_si128(Vrmsci+3, xmm11);
	}
}	


void fill_channels_aux(ESL_DSQ *seqs[4], int M, int iResidue, dpopt* opdata)
{
	int k;
	uniti xmm0, xmm1, xmm2, xmm3, xmm4, xmm5; // xmm6, xmm7, xmm8, xmm9, xmm10, xmm11;
	uniti *Vrmsci = (__m128i*) opdata->rmsc;
	float **rsc = opdata->rsc_msc;

	for (k = 0; k < M; k += 4, Vrmsci += 4)
	{
		xmm0 = (__m128i) _mm_load_ps(rsc[seqs[0][iResidue]]+k);
		xmm1 = (__m128i) _mm_load_ps(rsc[seqs[1][iResidue]]+k);
		xmm2 = (__m128i) _mm_load_ps(rsc[seqs[2][iResidue]]+k);
		xmm3 = (__m128i) _mm_load_ps(rsc[seqs[3][iResidue]]+k);

		xmm4 = xmm0;
		xmm0 = _mm_unpacklo_epi32(xmm0, xmm1);
		xmm4 = _mm_unpackhi_epi32(xmm4, xmm1);
		xmm5 = xmm2;
		xmm2 = _mm_unpacklo_epi32(xmm2, xmm3);
		xmm5 = _mm_unpackhi_epi32(xmm5, xmm3);

		xmm1 = xmm0;
		xmm0 = _mm_unpacklo_epi64(xmm0, xmm2);
		xmm1 = _mm_unpackhi_epi64(xmm1, xmm2);		
		xmm3 = xmm4;
		xmm4 = _mm_unpacklo_epi64(xmm4, xmm5);
		xmm3 = _mm_unpackhi_epi64(xmm3, xmm5);

		_mm_store_si128(Vrmsci+0, xmm0);
		_mm_store_si128(Vrmsci+1, xmm1);
		_mm_store_si128(Vrmsci+2, xmm4);
		_mm_store_si128(Vrmsci+3, xmm3);
	}
}	

__inline void fill_channels(ESL_DSQ *ddsq[4], int M, int i, dpopt* opdata)
{
//	fill_channels_4shifts(ddsq, M, i, opdata);
	fill_channels_epi32(ddsq, M, i, opdata);
}

//#define EXPANDV4F(v1,val)	
#define SETV4(v1,val)	\
	v1 = _mm_load_ss(&val);		\
	v1 = (__m128) _mm_shuffle_epi32((__m128i)v1, 0);

#define SSTEST 0

int viterbi_parallel_channels(ESL_DSQ *ddsq[4], int L, P7_PROFILE *gm, dpopt* opdata, float *opt_sc)
{
	float (*xsc)[p7P_NXTRANS] = gm->xsc;
	int M = gm->M, i, k;
	float neginf = -eslINFINITY;
	float zero = 0;
//	float esc1 = p7_profile_IsLocal(gm) ? 0 : -eslINFINITY;
	unit Vinf, xmxN, xmxB, xmxE, xmxC, xmxJ;
	unit *dmx = opdata->buffer, *imx = dmx+M, *mmx = imx+M;
	unit loopJ, loopE, loopC, loopN, moveJ, moveE, moveC, moveN;

	SETV4(loopJ, xsc[p7P_J][p7P_LOOP]);
	SETV4(loopE, xsc[p7P_E][p7P_LOOP]);
	SETV4(loopC, xsc[p7P_C][p7P_LOOP]);
	SETV4(loopN, xsc[p7P_N][p7P_LOOP]);
	SETV4(moveJ, xsc[p7P_J][p7P_MOVE]);
	SETV4(moveE, xsc[p7P_E][p7P_MOVE]);
	SETV4(moveC, xsc[p7P_C][p7P_MOVE]);
	SETV4(moveN, xsc[p7P_N][p7P_MOVE]);

	/* Initialization of the zero row.	*/
	SETV4(Vinf, neginf);
	SETV4(xmxN, zero);
	xmxE = xmxJ = xmxC = Vinf;
	xmxB = moveN;	// S->N->B, no N-tail

	for (k = 0; k < gm->M; k++)
		dmx[k] = mmx[k] = imx[k] = Vinf;

	/* DP recursion */
	for (i = 1; i <= L; i++) 
	{
		unit sc, dcv, temp, *tt = (unit*) opdata->tsc_all;
		unit *vrisc = (unit*) opdata->risc;
		unit *vrmsc = (unit*) opdata->rmsc;
		unit mpv, ipv, dpv;
		
		// prepare SSE vectors
		fill_channels(ddsq, gm->M, i, opdata);
		xmxE = mpv = dpv = ipv = sc = dcv = Vinf;
		continue;
		for (k = 0; k < M; k++)
		{
			sc	= _mm_add_ps(sc, *tt); tt++;
			dcv = _mm_max_ps(sc, _mm_add_ps(dcv, *tt));tt++;

			sc = _mm_add_ps(xmxB, *tt); tt++;
			sc = _mm_max_ps(sc, _mm_add_ps(mpv, *tt)); tt++;
			sc = _mm_max_ps(sc, _mm_add_ps(ipv, *tt)); tt++;
			sc = _mm_max_ps(sc, _mm_add_ps(dpv, *tt)); tt++;
			sc = _mm_add_ps(sc, *vrmsc); vrmsc++;
			xmxE = _mm_max_ps(xmxE, sc);

			dpv = dmx[k];
			ipv = imx[k];
			mpv = mmx[k];
			mmx[k] = sc;
			dmx[k] = dcv;
			
			temp = _mm_add_ps(mpv, *tt); tt++;
			temp = _mm_max_ps(temp, _mm_add_ps(ipv, *tt)); tt++;
			imx[k] = _mm_add_ps(temp, *vrisc); vrisc++;
		}
	
		xmxE = _mm_max_ps(xmxE, dmx[M-1]);
		xmxJ = _mm_max_ps(_mm_add_ps(xmxJ, loopJ), _mm_add_ps(xmxE, loopE));
		xmxC = _mm_max_ps(_mm_add_ps(xmxC, loopC), _mm_add_ps(xmxE, moveE));
		xmxN = _mm_add_ps(xmxN, loopN);
		xmxB = _mm_max_ps(_mm_add_ps(xmxN, moveN), _mm_add_ps(xmxJ, moveJ));
	}
	
	if (opt_sc != NULL)
		*((unit*) opt_sc) = _mm_add_ps(xmxC, moveC);

	return eslOK;
}




/*
gcc -g -Wall -O3 -std=gnu99 -o vitstream-float -I. -I.. -L.. -I../impl_sse -L../impl_sse -I../../easel -L../../easel zz*.c viterbi_serial.c -lhmmer -lhmmerimpl -leasel -lm -lpthread -DNTHREADS=1
*/

#include <string.h>
#include "easel.h"
#include "esl_randomseq.h"
#include <sys/timeb.h>

static ESL_OPTIONS options[] = {
	/* name			type		default	env	range toggles reqs incomp	help						docgroup*/
	{ "-h",		eslARG_NONE	, FALSE	, NULL, NULL , NULL, NULL, NULL, "show brief help on version and usage",0 },
	{ "-s",		eslARG_INT	, "42"	, NULL, NULL , NULL, NULL, NULL, "set random number seed to <n>",		0 },
	{ "-L",		eslARG_INT	, "400"	, NULL, "n>0", NULL, NULL, NULL, "length of random target seqs",		0 },
	{ "-N",		eslARG_INT	, "50000",NULL, "n>0", NULL, NULL, NULL, "number of random target seqs",		0 },
	{	0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage [] = "[-options] <hmmfile>";
static char banner[] = "benchmark driver for generic Viterbi";

int main(int argc, char **argv)
{
	ESL_GETOPTS	*go		= esl_getopts_CreateDefaultApp(options, 1, argc, argv, banner, usage);
	char		*hmmfile= esl_opt_GetArg(go, 1);
	ESL_STOPWATCH*w		= esl_stopwatch_Create();
	ESL_RANDOMNESS*r	= esl_randomness_Create(esl_opt_GetInteger(go, "-s"));
	ESL_ALPHABET*abc	= NULL;
	P7_HMMFILE	*hfp	= NULL;
	P7_HMM		*hmm	= NULL;
	P7_BG		*bg		= NULL;
	P7_PROFILE	*gm		= NULL;
	int			L		= 2000;	// esl_opt_GetInteger(go, "-L");
	int			N		= 20000; //esl_opt_GetInteger(go, "-N");
	ESL_DSQ		*dsq1, *dsq2, *dsq3, *dsq4; 
	ESL_DSQ		**ddsq = calloc(4, sizeof(ESL_DSQ*));
	int			i, j;
	__m128		resdata[10];
	float		*sc2 = (float*) (resdata+1);
	double		basesecs, secs;
	dpopt		*opdata;
	struct timeb tbstart, tbend;

	if (p7_hmmfile_Open(hmmfile, NULL, &hfp) != eslOK) p7_Fail("Failed to open HMM file %s", hmmfile);
	if (p7_hmmfile_Read(hfp, &abc, &hmm)	 != eslOK) p7_Fail("Failed to read HMM");
//	srand(time(NULL));

	bg = p7_bg_Create(abc);
	p7_bg_SetLength(bg, L);
	gm = p7_profile_Create(hmm->M, abc);
	p7_ProfileConfig(hmm, bg, gm, L, p7_UNILOCAL); // podia ser GLOCAL, GLOBAL, e UNI-qq para só 1 hit
	opdata = create_data_structs(gm);

	int M = gm->M;
	
	dsq1 = ddsq[0];
	dsq2 = ddsq[1];
	dsq3 = ddsq[2];
	dsq4 = ddsq[3];

	printf("M: %d | Q: %d | L: %d | N: %d\n", gm->M, (int) ceil(gm->M/4.0), L, N);
	/* Get a baseline time: how long it takes just to generate the sequences */

	ESL_DSQ	**seqs = calloc(N*4+100, sizeof(ESL_DSQ	*));
	for (i = 0; i < N*4; i++)
	{	seqs[i] = malloc((L+2)*sizeof(ESL_DSQ));
		esl_rsq_xfIID(r, bg->f, abc->K, L, seqs[i]);
	}
	basesecs =	0; //((tbend.time - tbstart.time) + (tbend.millitm - tbstart.millitm)*0.001)*4;

//	FILE* results = fopen("N-5-results.txt", "r");

	ftime(&tbstart);
	for (i = 0; i < N; i++)
	{
	//	if (i % 100 == 0) printf("RUN %d\n", i);
	
		for (j = 0; j < 4; j++)
			ddsq[j] = seqs[i*4 + j];
	
		viterbi_parallel_channels(ddsq, L, gm, opdata, sc2);
		
#if 0
		{
		float *sc1 = (float*) resdata;

	//	fscanf(results, "%f %f %f %f\n", &sc1[0], &sc1[1], &sc1[2], &sc1[3]);

		for (j = 0; j < 4; j++)
			p7_Viterbi_unilocal (ddsq[j], L, gm, sc1+j);

if(1)	for (j = 0; j < 4; j++)
			if (sc1[j] != sc2[j])
			{	printf("WRONG entry %d-%d: %f %f\n", i, j, sc1[j], sc2[j]);  getc(stdin); }
		}
#endif
	}
	ftime(&tbend);
	
	secs = (tbend.time - tbstart.time) + (tbend.millitm - tbstart.millitm)*0.001 - basesecs;
	w->user = secs;
	esl_stopwatch_Display(stdout, w, "# Opt CPU time: ");
	return 0;
}








#if 0
// ANTIGO original do svn
//#ifdef VITERBI_BENCHMARKdgshfh
/*	gcc -g -O2 -o benchmark -I. -L. -I../easel -L../easel -DVITERBI_BENCHMARK viterbi-parallel.c -lhmmer -leasel -lm
	./benchmark-generic-viterbi <hmmfile>
*/

#include "easel.h"
#include "esl_randomseq.h"

static ESL_OPTIONS options[] = {
	/* name			type		default	env	range toggles reqs incomp	help						docgroup*/
	{ "-h",		eslARG_NONE	, FALSE	, NULL, NULL , NULL, NULL, NULL, "show brief help on version and usage",0 },
	{ "-s",		eslARG_INT	, "42"	, NULL, NULL , NULL, NULL, NULL, "set random number seed to <n>",		0 },
	{ "-L",		eslARG_INT	, "400"	, NULL, "n>0", NULL, NULL, NULL, "length of random target seqs",		0 },
	{ "-N",		eslARG_INT	, "50000",NULL, "n>0", NULL, NULL, NULL, "number of random target seqs",		0 },
	{	0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]	= "[-options] <hmmfile>";
static char banner[] = "benchmark driver for generic Viterbi";

int main(int argc, char **argv)
{
	ESL_GETOPTS	*go		= esl_getopts_CreateDefaultApp(options, 1, argc, argv, banner, usage);
	char		*hmmfile= esl_opt_GetArg(go, 1);
	ESL_STOPWATCH	*w	= esl_stopwatch_Create();
	ESL_RANDOMNESS	*r	= esl_randomness_Create(esl_opt_GetInteger(go, "-s"));
	ESL_ALPHABET*abc	= NULL;
	P7_HMMFILE	*hfp	= NULL;
	P7_HMM		*hmm	= NULL;
	P7_BG		*bg		= NULL;
	P7_PROFILE	*gm		= NULL;
	P7_GMX		*gx		= NULL;
	int			L		= 0; //esl_opt_GetInteger(go, "-L");
	int			N		= esl_opt_GetInteger(go, "-N");
	ESL_DSQ		*dsq1, *dsq2, *dsq3, *dsq4; 
	ESL_DSQ		*ddsq[4];
	int			i, j;
	unit		resdata[10];
	float		*sc1 = (float*) resdata, *sc2 = (float*) (resdata+1);
	double		base_time, bench_time, Mcs;
	dpopt		*opdata;	

	if (p7_hmmfile_Open(hmmfile, NULL, &hfp) != eslOK) p7_Fail("Failed to open HMM file %s", hmmfile);
	if (p7_hmmfile_Read(hfp, &abc, &hmm)	 != eslOK) p7_Fail("Failed to read HMM");

	L = 200; // (hmm->M & ~0x7) +8;
	ddsq[0] = dsq1	= malloc(sizeof(ESL_DSQ) * (L+2));
	ddsq[1] = dsq2	= malloc(sizeof(ESL_DSQ) * (L+2));
	ddsq[2] = dsq3	= malloc(sizeof(ESL_DSQ) * (L+2));
	ddsq[3] = dsq4	= malloc(sizeof(ESL_DSQ) * (L+2));

	bg = p7_bg_Create(abc);
	p7_bg_SetLength(bg, L);
	gm = p7_profile_Create(hmm->M, abc);
	p7_ProfileConfig(hmm, bg, gm, L, p7_UNILOCAL);
	gx = p7_gmx_Create(gm->M, L);
	opdata = create_data_structs(gm);

	// Baseline time. 
	esl_stopwatch_Start(w);
	for (i = 0; i < (N/5); i++)
		esl_rsq_xfIID(r, bg->f, abc->K, L, dsq1);
	esl_stopwatch_Stop(w);
	base_time = w->user*4;

	printf("M = %d | Q = %d | L = %d | N = %d\n", gm->M, (int) ceil(gm->M/4.0), L, N);

	/* Benchmark time. */
	esl_stopwatch_Start(w);
	for (i = 0; i < 2*N/5; i++)
	{
		esl_rsq_xfIID(r, bg->f, abc->K, L, dsq1);
		esl_rsq_xfIID(r, bg->f, abc->K, L, dsq2);
		esl_rsq_xfIID(r, bg->f, abc->K, L, dsq3);
		esl_rsq_xfIID(r, bg->f, abc->K, L, dsq4);

		// original SSE com Farrar: 4.3 segs

		// tempo gasto em preparaçoes, encher canais etc: 1.8segs
		viterbi_parallel_channels(ddsq, L, gm, opdata, sc2);		// 6.0segs

#if SSTEST	// 51segs
		viterbi_serial_opt (dsq1, L, gm, sc1+0);
		viterbi_serial_opt (dsq2, L, gm, sc1+1);
		viterbi_serial_opt (dsq3, L, gm, sc1+2);
		viterbi_serial_opt (dsq4, L, gm, sc1+3);
		for (j = 0; j < 4; j++)
		{
			if (sc1[j] != sc2[j])
			{	printf("WRONG entry %d-%d: %f %f\n", i, j, sc1[j], sc2[j]);
				getc(stdin);
			}
		}
#endif
	}
	esl_stopwatch_Stop(w);
	
	bench_time = w->user - base_time+0.02;
	w->user = bench_time;
	Mcs		= (double) N * (double) L * (double) gm->M * 1e-6 / (double) bench_time;
	esl_stopwatch_Display(stdout, w, "# Opt CPU time: ");
	printf("# M	= %d # %.1f Mc/s\n", gm->M, Mcs);

//	printf("Total: %d | Lazy: %d\n", ntotal, nlazy);

	free(dsq1);	free(dsq2);	free(dsq3);	free(dsq4);
	p7_gmx_Destroy(gx);
	p7_profile_Destroy(gm);
	p7_bg_Destroy(bg);
	p7_hmm_Destroy(hmm);
	p7_hmmfile_Close(hfp);
	esl_alphabet_Destroy(abc);
	esl_stopwatch_Destroy(w);
	esl_randomness_Destroy(r);
	esl_getopts_Destroy(go);
	return 0;
}
#endif 



