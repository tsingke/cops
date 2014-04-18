#include "hmmer.h"
#include "p7_config.h"
#include "string.h"



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
#define SCALE	1000.0

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

		if (xmxE > WORDMAX-1) { *opt_ret = eslINFINITY; printf("Positive Overflow: %i\n", xmxE); return eslERANGE; }
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

