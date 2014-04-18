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

#include <p7_config.h>

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <sys/timeb.h>

#include <easel.h>
#include <hmmer.h>


// Viterbi Filter SSE implementation, with Farrar scheme, 16bit integer scores. vitfilter.c
extern int p7_ViterbiFilter(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_OMX *ox, float *ret_sc);

// Viterbi Filter SSE implementation, with Farrar scheme, 32bit float scores. vitscore.c
extern int p7_ViterbiScore (const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_OMX *ox, float *ret_sc);

#include "viterbi_serial.h"


/*****************************************************************
 *  ViterbiFilter threaded Benchmark
 *****************************************************************/

#include <pthread.h>

typedef struct _dsq_cmp
{
	ESL_DSQ *seq;
	int length;
} SEQ;

typedef struct _thr_info
{	int thrid;
	int N;
	int L;
    int tocheck;
    int nrounds;
	P7_PROFILE	*gm;
	P7_OPROFILE	*om;
	SEQ **seqs;
	int* lastseq;
	pthread_mutex_t* dbmutex;
	pthread_barrier_t* barrier;
} pthr_info_t;

#define TIMEDIFF(t1,t2) ((t2.time - t1.time) + (t2.millitm - t1.millitm)*0.001)

void* vitfilter_thr_loop(void* args)
{
	pthr_info_t* info = (pthr_info_t*) args;
	int i, N = info->N, L = 0, M = info->gm->M;
	float sc1;
	ESL_DSQ* dsq;
	SEQ *dsqcmp;
	P7_OPROFILE *om = info->om;
	P7_OMX		*ox	= p7_omx_Create(M, 0, 0);
	int nseqs = 0;
    float sumerrors = 0.0;

#ifdef _GNU_SOURCE
	printf("THR %d running on cpu %d\n", info->thrid, sched_getcpu());
#endif

    int j;
    for (j = 0; j < info->nrounds; j++) 
    {
        pthread_barrier_wait(info->barrier);
        __sync_fetch_and_and(info->lastseq, 0);
        // Wait for all threads to be ready
        pthread_barrier_wait(info->barrier);
        
        for (i = 0; 1; i++)
        {
            int lseq, newlen;
            // fetch new sequence
            pthread_mutex_lock(info->dbmutex);
            // full memory barrier. make sure other threads see updated value of lastseq
            lseq = __sync_fetch_and_add(info->lastseq, 1);
            dsqcmp = info->seqs[lseq];	// padded with dummies
            // info->seqs[lseq] = NULL;
            pthread_mutex_unlock(info->dbmutex);

            if (lseq >= N) break;

            dsq = dsqcmp->seq;
            newlen = dsqcmp->length;

            if (L != newlen)	// resize if L changes
            {	L = newlen;
                p7_ReconfigLength(info->gm, L);
                p7_oprofile_ReconfigRestLength(om, L);
            }

            nseqs++;
        //	p7_GViterbi(dsq, L, info->gm, gx, &sc1);
            p7_ViterbiFilter(dsq, L, om, ox, &sc1);		// 8x16bit integers
        //	p7_ViterbiScore (dsq, L, om, ox, &sc1);		// 4x32bit floats

            if (info->tocheck)
            {
                float sc2;

                p7_Viterbi_unilocal(dsq, L, info->gm, &sc2);
            //	p7_Viterbi_unilocal_word(dsq, L, info->gm, &sc2);
                
                sc1 = sc1 + 1.0;	// -2.0nat optimization, Local to Unilocal mode
                sumerrors += fabs(sc1 - sc2);
                if (fabs(sc1 - sc2) >  0.7) 
                    printf("[T%d-%d] L %d |  %.4f %.4f => diff %f\n", info->thrid, i, L, sc1, sc2, fabs(sc1 - sc2));
            }
        }
    }

//	printf("Finished Thread %d. Computed %d seqs. Sum errors: %f\n", info->thrid, nseqs, sumerrors);
	return (void*) 0;
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

/* -c, -x are used for debugging, testing; see msvfilter.c for explanation
gcc vitfilter_test_threaded.c viterbi_serial.c -O3 -o vitfilterthr -std=gnu99 -g -Wall -msse2  -I.. -L..  -I../impl_sse -L../impl_sse -I../../easel -L../../easel -Dp7VITFILTER_BENCHMARK -lhmmer -lm -leasel -lpthread -DARCH=UMA -DNTHREADS=1
*/


#define UMA     100
#define NUMA    101
#define ALEPH   102

#if   ARCH == UMA
    #define MAPIDCPU(id) ((id) % 8)
#elif ARCH == NUMA
    #define MAPIDCPU(id) (((id)*4 + (id)/8) % 32)
#elif ARCH == ALEPH
	#define MAPIDCPU(id) ((id)*2 + 16)
#else
    #define MAPIDCPU(id) (id)
    #error "WRONG ARCH OR NO ARCH SPECIFIED"
#endif


#include <limits.h>
#include <esl_sqio.h>
#include <esl_randomseq.h>


static ESL_OPTIONS options[] = {
	/* name	type		default	env	range toggles reqs incomp	help							docgroup*/
	{ "-h",	eslARG_NONE,FALSE, NULL, NULL,	NULL,	NULL, NULL, "show brief help",					0 },
	{ "-c",	eslARG_NONE,FALSE, NULL, NULL,	NULL,	NULL, "-x", "compare scores to generic impl",	0 }, 
	{ "-s",	eslARG_INT,	"42",  NULL, NULL,	NULL,	NULL, NULL, "set random number seed to <n>",	0 },
	{ "-x",	eslARG_NONE,FALSE, NULL, NULL,	NULL,	NULL, "-c", "equate scores to trusted impl",	0 },
	{ "-L",	eslARG_INT,	"400", NULL, "n>0", NULL,	NULL, NULL, "length of random target seqs",		0 },
	{ "-N",	eslARG_INT,	"8800",NULL, "n>0", NULL,	NULL, NULL, "number of random target seqs",		0 },
	{ "-R", eslARG_INT, "1"	 , NULL, "n>0",	NULL,	NULL, NULL, "number of test rounds",			0 },
	{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]	= "[-options] <hmmfile>";
static char banner[] = "benchmark driver for Viterbi filter";

#define TIMEDIFF(t1,t2) ((t2.time - t1.time) + (t2.millitm - t1.millitm)*0.001)

int main(int argc, char **argv)
{
	ESL_GETOPTS	*go		= esl_getopts_CreateDefaultApp(options, 2, argc, argv, banner, usage);
	ESL_STOPWATCH*w		= esl_stopwatch_Create();
	ESL_RANDOMNESS*r	= esl_randomness_Create(esl_opt_GetInteger(go, "-s"));
	ESL_ALPHABET*abc	= NULL;
	P7_HMMFILE	*hfp	= NULL;
	P7_HMM		*hmm	= NULL;
	P7_BG		*bg		= NULL;
	P7_PROFILE	*gm		= NULL;
	P7_OPROFILE	*om		= NULL;
	int			L		= esl_opt_GetInteger(go, "-L");
	int			N   	= esl_opt_GetInteger(go, "-N");
	int 		NROUNDS	= esl_opt_GetInteger(go, "-R");
   	int			check	= esl_opt_GetBoolean(go, "-c");	
	ESL_DSQ		*dsq	= NULL;	
	ESL_SQFILE   *sqfp	= NULL;
	char	*hmmfile	= esl_opt_GetArg(go, 1);
	char	*seqfile	= esl_opt_GetArg(go, 2);
	int		i, j;
	int		sumlengths	= 0;	
	struct timeb tbstart, tbend;

	if (p7_hmmfile_Open(hmmfile, NULL, &hfp)!= eslOK) p7_Fail("Failed to open HMM file %s", hmmfile);
	if (p7_hmmfile_Read (hfp, &abc, &hmm)	!= eslOK) p7_Fail("Failed to read HMM");

	bg = p7_bg_Create(abc);
	p7_bg_SetLength(bg, L);
	gm = p7_profile_Create(hmm->M, abc);
	p7_ProfileConfig(hmm, bg, gm, L, p7_UNILOCAL);
	om = p7_oprofile_Create(hmm->M, abc);
	p7_oprofile_Convert(gm, om);

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
	else	// Generate random sequences
		for (j = 0; j < N; j++)
		{
			int len = L; // - rand()%1000;
			seqsdb[j] = malloc(sizeof(SEQ));
			seqsdb[j]->seq = malloc(len+4);
			seqsdb[j]->length = len;
			esl_rsq_xfIID(r, bg->f, abc->K, len, seqsdb[j]->seq);
			sumlengths += len;
		}

//	qsort(seqsdb, N, sizeof(SEQ*), compare_seqs);

	printf("ViterbiFilter HMMER with %d threads, %s. ModelL.: %d, #Segms: %d, SeqL.: %d, #seqs: %d\n",
			NTHREADS, hmmfile, gm->M, (int) ceil(gm->M/8.0), L, N*NROUNDS);

	pthread_mutex_t dbmutex;
	pthread_mutex_init(&dbmutex, NULL);
	pthread_barrier_t barrier;
	pthread_barrier_init(&barrier, NULL, NTHREADS+1);

	pthread_t *threads[NTHREADS];
	pthr_info_t args;
	int lastseqidx = 00000;
	args.L = L; args.N = N;	args.gm = gm; args.thrid = 0; args.om = NULL;
	args.seqs = seqsdb; args.lastseq = &lastseqidx;
	args.dbmutex = &dbmutex; args.barrier = &barrier;
    args.tocheck = check; args.nrounds = NROUNDS;

	for (i = 0; i < NTHREADS; i++)
	{
		pthr_info_t *argscopy = calloc(1, sizeof(pthr_info_t));
		memcpy(argscopy, &args, sizeof(pthr_info_t));
		argscopy->thrid	= i;
		argscopy->om	= p7_oprofile_Clone(om);

		pthread_attr_t attr;
		pthread_attr_init(&attr);
#ifdef _GNU_SOURCE
		cpu_set_t cpuset;
		CPU_ZERO(&cpuset);
		CPU_SET(MAPIDCPU(i), &cpuset);
		pthread_attr_setaffinity_np(&attr, sizeof(cpu_set_t), &cpuset);
#endif
		threads[i] = calloc(1, sizeof(pthread_t));
		if (pthread_create(threads[i], &attr, vitfilter_thr_loop, argscopy))
			return fprintf(stderr, "ERROR could not create worker thread %d\n", i);
	}

	sleep(1);

	ftime(&tbstart);

	for (j = 0; j < NROUNDS; j++) 
	{	pthread_barrier_wait(&barrier);
		pthread_barrier_wait(&barrier);
	}

	void* retval;
	for (i = 0; i < NTHREADS; i++)
		pthread_join(*threads[i], &retval);

	ftime(&tbend);

	double secs = TIMEDIFF(tbstart,tbend);
	w->elapsed = w->user = secs;
	esl_stopwatch_Display(stdout, w, "# Opt CPU time: ");
	double compmillioncells = NROUNDS * (double) sumlengths * (double) hmm->M * 1e-6;
	printf("# %.0fM cells in %.1f Mc/s\n", compmillioncells, compmillioncells / secs);

	free(dsq);
	p7_oprofile_Destroy(om);
	p7_profile_Destroy(gm);
	p7_bg_Destroy(bg);
	p7_hmm_Destroy(hmm);
	p7_hmmfile_Close(hfp);
	esl_sqfile_Close(sqfp);
	esl_alphabet_Destroy(abc);
	esl_stopwatch_Destroy(w);
	esl_getopts_Destroy(go);
	return 0;
}



