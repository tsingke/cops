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


#include <string.h>
#include <hmmer.h>
#include <p7_config.h>



/*****************************************************************
 *  Viterbi Serial optimized general
 *****************************************************************/

int p7_Viterbi_general(ESL_DSQ *dsq, int L, P7_PROFILE *gm,  float *opt_res)
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
	xmxN = 0.0;
	xmxB = xmxN + gm->xsc[p7P_N][p7P_MOVE];				/* S->N->B, no N-tail	*/
	xmxE = xmxJ = xmxC = -eslINFINITY;					/* need seq to get here */
		
	for (k = 0; k <= gm->M; k++)
		dmx[k] = mmx[k] = imx[k] = -eslINFINITY;		/* need seq to get here */
	
	/* DP recursion */
	for (i = 1; i <= L; i++) 
	{
		float *rsc = gm->rsc[dsq[i]];
		float sc, ipv, mpv, dpv, dcv;

		dcv = mpv = dpv = ipv = sc = -eslINFINITY;
		// nao faz nada aqui, e estraga no caso de usarmos os indices 0!
//		imx[0] = dmx[0] = mmx[0] = -eslINFINITY;
		xmxE = -eslINFINITY;

		for (k = 1; k <= M; k++)
		{
			dcv	= ESL_MAX(sc + TSC(p7P_MD,k-1), dcv + TSC(p7P_DD,k-1));	// para k=0, sempre -infinito

			/* match state. Sc == new MMX(k) */
			sc	= xmxB + TSC(p7P_BM,k-1);	// LM == antigo BM
			sc	= ESL_MAX(sc, mpv + TSC(p7P_MM,k-1));
			sc	= ESL_MAX(sc, ipv + TSC(p7P_IM,k-1));
			sc	= ESL_MAX(sc, dpv + TSC(p7P_DM,k-1));
			sc	= sc + MSC(k); // rsc[p7P_M];
			xmxE= ESL_MAX(xmxE, sc);	// + esc);

			ipv = imx[k];
			mpv = mmx[k];
			dpv = dmx[k];

			mmx[k] = sc;
			dmx[k] = dcv;
		
			if (k < M)
			{	imx[k] = ESL_MAX(mpv + TSC(p7P_MI,k), ipv + TSC(p7P_II,k)) + ISC(k);	// II e MI vao de 1 a M-1
			//	dcv = ESL_MAX(sc + TSC(p7P_MD,k), dcv + TSC(p7P_DD,k));	// MD e MD vao de 0 a M-1, com o [0] sempre a -infinito
			}
		}

		/* E state update; transition from M_M scores 0 by def'n */
//		xmxE	= ESL_MAX(xmxE, dmx[M]);	// 100% always useless

		/* Now the special states. E must already be done, and B must follow N,J.
		 * remember, N, C and J emissions are zero score by definition. */
		xmxJ	= ESL_MAX(xmxJ + loopJ, xmxE + loopE);	/* E->J is E's "loop" */
		xmxC	= ESL_MAX(xmxC + loopC, xmxE + moveE);
		xmxN	= xmxN		   + loopN;
		xmxB	= ESL_MAX(xmxN + moveN, xmxJ + moveJ);	/* J->B is J's move */
	}
	
	/* T state (not stored) */
	if (opt_res != NULL)
		opt_res[0] = xmxC + moveC;

	return eslOK;
}




/*****************************************************************
 *  Viterbi Serial optimized for Unilocal alignments
 *****************************************************************/

static float* eslInf = NULL;
float max = -eslINFINITY;
float min = eslINFINITY;

int p7_Viterbi_unilocal(ESL_DSQ *dsq, int L, P7_PROFILE *gm,  float *opt_res)
{
	int M = gm->M;
	int	i,k;

	// float esc = p7_profile_IsLocal(gm) ? 0 : -eslINFINITY;
	float xmxN, xmxB, xmxE, xmxC, xmxJ;
	float dmx[M], imx[M], mmx[M];
	float moveC = gm->xsc[p7P_C][p7P_MOVE];
	float moveN = gm->xsc[p7P_N][p7P_MOVE];
	float loopC = gm->xsc[p7P_C][p7P_LOOP];
	float loopN = gm->xsc[p7P_N][p7P_LOOP];

	xmxN = 0.0;
	xmxB = xmxN + moveN;								/* S->N->B, no N-tail	*/
	xmxE = xmxJ = xmxC = -eslINFINITY;					/* need seq to get here */		

	if (eslInf == NULL)
	{	eslInf = malloc((gm->M+1)*sizeof(float));
		for (k = 0; k <= gm->M; k++)
			eslInf[k] = -eslINFINITY;
	}

	memcpy(dmx, eslInf, (M)*sizeof(float));
	memcpy(mmx, eslInf, (M)*sizeof(float));
	memcpy(imx, eslInf, (M)*sizeof(float));

	for (i = 1; i <= L; i++) 
	{
		float *tsc= gm->tsc;
		float *rsc = gm->rsc[dsq[i]];
		float sc, ipv, mpv, dpv, dcv;

		dcv = mpv = dpv = ipv = sc = xmxE = -eslINFINITY;

//	#define TSC(s,k) (tsc[(k) * p7P_NTRANS + (s)])
//	#define MSC(k)   (rsc[(k) * p7P_NR     + p7P_MSC])

		for (k = 0; k < M-1; k++)
		{
		//	dcv	= ESL_MAX(sc + TSC(p7P_MD,k), dcv + TSC(p7P_DD,k));

			sc	= xmxB + tsc[p7P_BM];
			sc	= ESL_MAX(sc, mpv + tsc[p7P_MM]);
			sc	= ESL_MAX(sc, ipv + tsc[p7P_IM]);
			sc	= ESL_MAX(sc, dpv + tsc[p7P_DM]);
			sc	= sc + MSC(k+1);
			xmxE	= ESL_MAX(xmxE, sc);

			dpv = dmx[k];
			ipv = imx[k];
			mpv = mmx[k];
			mmx[k] = sc;
			dmx[k] = dcv;

			tsc += p7P_NTRANS;
			imx[k]	= ESL_MAX(mpv + tsc[p7P_MI], ipv + tsc[p7P_II]);
			dcv	= ESL_MAX(sc + tsc[p7P_MD], dcv + tsc[p7P_DD]);
		}

		// unroll last iteration, M-1, to remove computation of imx[M-1] and new dcv
		sc	= xmxB + tsc[p7P_BM];
		sc	= ESL_MAX(sc, mpv + tsc[p7P_MM]);
		sc	= ESL_MAX(sc, ipv + tsc[p7P_IM]);
		sc	= ESL_MAX(sc, dpv + tsc[p7P_DM]);
		sc	= sc + MSC(k+1);
		xmxE	= ESL_MAX(xmxE, sc);

//		xmxE	= ESL_MAX(xmxE, dcv /* == dmx[M-1] */);	// Wing-retraction. Redundant in unilocal mode
		xmxC	= ESL_MAX(xmxC + loopC, xmxE /* moveE == 0*/);
		xmxN	= xmxN + loopN;
		xmxB	= xmxN + moveN;
	}
	
	if (opt_res != NULL)		opt_res[0] = xmxC + moveC;

	return eslOK;
}





/*****************************************************************
 *		Viterbi Serial optimized for Unilocal alignments
 *			and discretized for 16bit integers
 *****************************************************************/

//	A precisao dos resultados e' controlado so' pela escala

// configurar
#define OFFSET	12000
#define SCALE	(500.0/eslCONST_LOG2)

#define WORDMAX 32767


	// simple and slow
static int16_t adds16(int16_t a, int16_t b)
{
	int sum = ((int) a) + (int) b;
	if (sum < -WORDMAX)	return -WORDMAX;
	else if (sum > WORDMAX)	return WORDMAX;
	else return sum;
}


#if 0
// from http://locklessinc.com/articles/sat_arithmetic/
// In Core2 arch, it's much slower than the simple branching version!!!
int16_t donotuse_aadds16(int16_t x, int16_t y)
{
	uint16_t ux = x;
	uint16_t uy = y;
	uint16_t res = ux + uy;
	// Calculate overflowed result. (Don't change the sign bit of ux)
	ux = (ux >> 15) + WORDMAX;
	// Force compiler to use cmovns instruction
	if ((int16_t) ((ux ^ uy) | ~(uy ^ res)) >= 0)
		res = ux;		
	return res;
}
#endif

//#undef adds16
//#define adds16(a,b)  (a+b)

static int16_t discretize_word(float ff)
{
	if (isinf(ff))	return -WORDMAX;
	else			return (int16_t) (roundf (ff * SCALE));
}


static 	int16_t* tsci = NULL;
static 	int16_t** rsci = NULL;
static 	int16_t* infM = NULL;

static void viterbi_unilocal_word_prepare(P7_PROFILE *gm)
{
	if (tsci != NULL) return;

	float *tsc= gm->tsc;
	int16_t* tt;
	int k, i;

	tsci = calloc((gm->M+1)*8, sizeof(int16_t));
	infM = calloc((gm->M+8), sizeof(int16_t));
	
	for (k = 0; k <= gm->M; k++)
		infM[k] = -WORDMAX;

#undef TSC
#define TSC(s,k) discretize_word(tsc[(k) * p7P_NTRANS + (s)])

	for (k = 0, tt = tsci; k < gm->M-1; k++, tt += p7P_NTRANS, tsc += p7P_NTRANS)
	{	tt[0] = discretize_word(tsc[p7P_BM]);
		tt[1] = discretize_word(tsc[p7P_MM]);
		tt[2] = discretize_word(tsc[p7P_IM]);
		tt[3] = discretize_word(tsc[p7P_DM]);
		tt[4] = discretize_word(tsc[p7P_NTRANS+p7P_MI]);
		tt[5] = discretize_word(tsc[p7P_NTRANS+p7P_II]);
		tt[6] = discretize_word(tsc[p7P_NTRANS+p7P_MD]);
		tt[7] = discretize_word(tsc[p7P_NTRANS+p7P_DD]);
	}
	tt[0] = discretize_word(tsc[p7P_BM]);
	tt[1] = discretize_word(tsc[p7P_MM]);
	tt[2] = discretize_word(tsc[p7P_IM]);
	tt[3] = discretize_word(tsc[p7P_DM]);
	tt[4] = -WORDMAX;	// para MI
	tt[5] = -WORDMAX;	// para II
	tt[6] = -WORDMAX;	// para MD
	tt[7] = -WORDMAX;	// para DD

	const int AlphaSize = gm->abc->Kp;
	rsci = calloc(AlphaSize+1, sizeof(int16_t*));
	int16_t *data = calloc((gm->M+1)*AlphaSize, sizeof(int16_t));

	for (k = 0; k < AlphaSize; k++)
	{
		int16_t *rmsc = rsci[k] = data + (k)*(gm->M);
		float *rsc = gm->rsc[k] + p7P_NR + p7P_MSC;	// ignore 1st values
//#undef MSC
//#define MSC(k)   discretize_word(rsc[(k) * p7P_NR     + p7P_MSC])
		for (i = 0; i< gm->M; i++, rsc += p7P_NR)
			*rmsc++ = discretize_word(*rsc);
	}
}

int p7_Viterbi_unilocal_word(ESL_DSQ *dsq, int L, P7_PROFILE *gm,  float *opt_ret)
{
	int M = gm->M, i, k;
	int16_t* tt;
	int16_t xmxB, xmxE, xmxC;
	int16_t moveN, moveC;
	int16_t dmx[M+2], imx[M+2], mmx[M+2];

	if (tsci == NULL)
		viterbi_unilocal_word_prepare(gm);
		
	moveN = discretize_word(gm->xsc[p7P_N][p7P_MOVE]);
	moveC = discretize_word(gm->xsc[p7P_C][p7P_MOVE]);

//	xmxN = 0;
	xmxB = moveN;
	xmxE = xmxC = -WORDMAX;

	memcpy(dmx, infM, (M+1)*sizeof(int16_t));
	memcpy(mmx, infM, (M+1)*sizeof(int16_t));
	memcpy(imx, infM, (M+1)*sizeof(int16_t));

	/* DP recursion */
	for (i = 1; i <= L; i++) 
	{
		int16_t *rsc = rsci[dsq[i]];
		int16_t sc, ipv, mpv, dpv, dcv;
		dcv = mpv = dpv = ipv = sc = xmxE = -WORDMAX;

		// 2nat optimizacao:
		xmxB   = adds16(OFFSET, moveN);
	
		for (k = 0, tt = tsci; k < M-1; k++, tt += p7P_NTRANS)
		{
			sc	= adds16(xmxB, tt[0]);
			sc	= ESL_MAX(sc, adds16(mpv, tt[1]));
			sc	= ESL_MAX(sc, adds16(ipv, tt[2]));
			sc	= ESL_MAX(sc, adds16(dpv, tt[3]));
			sc	= adds16(sc, *rsc++);
			xmxE= ESL_MAX(xmxE, sc);

			dpv = dmx[k];
			ipv = imx[k];
			mpv = mmx[k];
			mmx[k] = sc;
			dmx[k] = dcv;
			
			imx[k] = ESL_MAX(adds16(mpv, tt[4]), adds16(ipv, tt[5]));
			dcv    = ESL_MAX(adds16(sc , tt[6]), adds16(dcv, tt[7]));
		}

		// unroll last iteration, M-1, to remove computation of imx[M-1] and new dcv
		sc	= adds16(xmxB, tt[0]);
		sc	= ESL_MAX(sc, adds16(mpv, tt[1]));
		sc	= ESL_MAX(sc, adds16(ipv, tt[2]));
		sc	= ESL_MAX(sc, adds16(dpv, tt[3]));
		sc	= adds16(sc, *rsc);
		xmxE= ESL_MAX(xmxE, sc);

//		xmxE = ESL_MAX(xmxE, dmx[M-1]);	// Wing-retraction. Redundant in unilocal mode

		if (xmxE > WORDMAX-1) { *opt_ret = eslINFINITY; return eslERANGE; }
		xmxC = ESL_MAX(	xmxC, // adds16(xmxC, loopC), // removido pela opt 2nat
						xmxE);
	//	removidos pela opt 3nat
	//	xmxN = adds16(xmxN, loopN);
	//	xmxB = adds16(xmxN, moveN);
	}

	xmxC = adds16(xmxC, moveC);
	
	if (opt_ret != NULL)
		opt_ret[0] = (xmxC-OFFSET)/SCALE  -2.0;	// 2nat optimizacao

	return eslOK;
}





/*****************************************************************
 *	Benchmark driver.
 *****************************************************************/

/*
gcc -g -Wall -O3 -std=gnu99 -o viterbiserial -I. -I.. -L.. -I../impl_sse -L../impl_sse -I../../easel -L../../easel viterbi_serial.c -lhmmer -leasel -lm  -DSERIAL
*/
#ifdef SERIAL

#include <sys/timeb.h>
#include <limits.h>
#include <easel.h>
#include <esl_randomseq.h>
#include <esl_sqio.h>

static ESL_OPTIONS options[] = {
	/* name		type	default		env	range toggles reqs incomp	help						docgroup*/
	{ "-h",	eslARG_NONE	, FALSE	, NULL, NULL , NULL, NULL, NULL, "show brief help on version and usage",0 },
	{ "-s",	eslARG_INT	, "42"	, NULL, NULL , NULL, NULL, NULL, "set random number seed to <n>",		0 },
	{ "-L",	eslARG_INT	, "400"	, NULL, "n>0", NULL, NULL, NULL, "length of random target seqs",	0 },
	{ "-N",	eslARG_INT	, "8800",NULL, "n>0", NULL, NULL, NULL, "number of random target seqs",		0 },
	{ "-R", eslARG_INT	, "1"	, NULL, "n>0", NULL, NULL, NULL, "number of test rounds",				0 },
	{	 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage [] = "[-options] <hmmfile>";
static char banner[] = "benchmark driver for generic Viterbi";


typedef struct _dsq_cmp
{
	ESL_DSQ *seq;
	int length;
} SEQ;


#define TIMEDIFF(t1,t2) ((t2.time - t1.time) + (t2.millitm - t1.millitm)*0.001)

int main(int argc, char **argv)
{
	ESL_GETOPTS *go     = esl_getopts_CreateDefaultApp(options, 2, argc, argv, banner, usage);
	char	*hmmfile	= esl_opt_GetArg(go, 1);
	char	*seqfile	= esl_opt_GetArg(go, 2);
	ESL_STOPWATCH *w	= esl_stopwatch_Create();
	ESL_RANDOMNESS*r	= esl_randomness_Create(esl_opt_GetInteger(go, "-s"));
	ESL_ALPHABET*abc	= NULL;
	P7_HMMFILE	*hfp	= NULL;
	P7_HMM		*hmm	= NULL;
	P7_BG		*bg		= NULL;
	P7_PROFILE	*gm1;
	int			L		= esl_opt_GetInteger(go, "-L");
	int			N		= esl_opt_GetInteger(go, "-N");
   	int			NROUNDS	= esl_opt_GetInteger(go, "-R");
	__m128		resdata[10];
	float		*sc1 	= (float*)	(resdata+0);	// uses 1 __m128s
	int			i, j;
	ESL_SQFILE   *sqfp	= NULL;
	struct timeb tbstart, tbend;
	int sumlengths = 0;

   	srand(time(NULL));
	if (p7_hmmfile_Open(hmmfile, NULL, &hfp)!= eslOK) p7_Fail("Failed to open HMM file %s", hmmfile);
	if (p7_hmmfile_Read(hfp, &abc, &hmm)	!= eslOK) p7_Fail("Failed to read HMM");

	bg = p7_bg_Create(abc);
	p7_bg_SetLength(bg, L);

	gm1 = p7_profile_Create(hmm->M, abc);
	p7_ProfileConfig(hmm, bg, gm1, L, p7_UNILOCAL);

    int dbsize = N;
	SEQ **seqsdb= calloc(dbsize+1, sizeof(SEQ*));

	if (esl_sqfile_OpenDigital(abc, seqfile, eslSQFILE_FASTA, NULL, &sqfp) == eslOK)
	{	// Use Sequence file
		ESL_SQ* sq = esl_sq_CreateDigital(abc);
		int maxseqs, len=0;

        if (esl_opt_IsDefault(go, "-N"))    // N not specified in cmdline
            maxseqs = INT_MAX;   // no limit
        else
            maxseqs = N;      // use cmdline limit
            
		for (j = 0; j < maxseqs && esl_sqio_Read(sqfp, sq) == eslOK; j++)
		{
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
		N = j;
	}
	else	// Not found database. Generate random sequences
        for (i = 0; i < N; i++)
		{
			int len = L; // - rand()%1000;
			seqsdb[i] = malloc(sizeof(SEQ));
			seqsdb[i]->seq = malloc(len+4);
			seqsdb[i]->length = len;
			esl_rsq_xfIID(r, bg->f, abc->K, len, seqsdb[i]->seq);
			sumlengths += len;
		}


    printf("Viterbi Serial, model %s. ModelLen: %d, SeqL.: %d, #seqs: %d\n", hmmfile, gm1->M, L, N*NROUNDS);
          
    ftime(&tbstart);

	for (j = 0; j < NROUNDS; j++) 
        for (i = 0; i < N; i++)
        {
        //	if (i % 10000 == 0) printf("Seq %d\n", i);

            p7_Viterbi_unilocal(seqsdb[i]->seq, seqsdb[i]->length, gm1, sc1);
        }
        
	ftime(&tbend);

	double secs = TIMEDIFF(tbstart,tbend);
	w->elapsed = w->user = secs;
	esl_stopwatch_Display(stdout, w, "# Opt CPU time: ");
	double compmillioncells = NROUNDS * (double) sumlengths * (double) hmm->M * 1e-6;
	printf("# %.0fM cells in %.1f Mc/s\n", compmillioncells, compmillioncells / secs);

	return 0;
}

#endif


