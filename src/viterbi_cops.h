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


/** @file 
 * ViterbiCOPS public interface
 */

#include <hmmer.h>

typedef struct _dsq_cmp
{	ESL_DSQ *seq;
	int length;
} SEQ;


/***************************************************************************
 *	Viterbi COPS with 8-Word channels and multi-threaded partitions
 ***************************************************************************/

#define SSE16_NVALS 8
 
typedef struct _DATA_COPS16
{	// data TSC
	__m128i* tsc_all;
	int16_t** rsc_msc;

	int nflags;	
	int partition;
	int Npartitions;

	double scale; 		// configurable
	int16_t wordoffset;	// configurable
	
	// Mutable args. Master thread must set them for workers
	int thrid;
	P7_PROFILE *gm;
	ESL_DSQ	*seqs[SSE16_NVALS];	// seqs to align
	int allocL;		// size of allocated sequence buffers
	int L;			// max length of seqs to align

	byte *synchflags[200];	// max threads
	// counter with # of runs to synchronize the threads starting each run
	int synccontrol;

	__m128i* xmxE;
	__m128i* pdcv;
	__m128i* psc;
} DATA_COPS16;

/** 
 * COPS processing function. Executes Viterbi algorithm on the 8-sequence batch
 * 	\arg COPS data
 *	\arg Array with 8 sequences to score
 *	\arg Result array with 8 float values [out]
 */
int p7_ViterbiCOPSw_run(DATA_COPS16* dcops, SEQ **seqsdb, float* results);

/**
 *	Setup COPS. Needs to be called before calling the processing function
 *	\arg COPS data
 *	\arg Estimate of the maximum Sequences' length. Used to optimize the memory allocation. Defaults to 1000
 *	\arg Maximum partition length allowed. Must be a multiple of 8. Defaults to 120
 */
void p7_ViterbiCOPSw_Setup(DATA_COPS16* dcops, int maxL, int max_partition);

/**
 *	Create COPS data structure for the 16-bit integer implementation
 *	\arg A valid P7_PROFILE
 */
DATA_COPS16* p7_ViterbiCOPSw_Create(P7_PROFILE *gm);



/***************************************************************************
 *	Viterbi COPS with 4-Float channels and multi-threaded partitions
 ***************************************************************************/

 #define SSE32_NVALS 4

typedef struct _DATA_COPS32
{
	__m128* tsc_all;
	float** rsc_msc;

	int nflags;	
	int partition;
	int Npartitions;
	
	// Mutable args. Master thread must set them for workers
	int thrid;
	P7_PROFILE *gm;
	ESL_DSQ	*seqs[SSE32_NVALS];	// seqs to align
	int allocL;		// size of allocated sequence buffers
	int L;			// max length of seqs to align

	byte *synchflags[200];	// max threads
	// counter with # of runs to synchronize the threads starting each run
	int synccontrol;

	__m128* xmxE;
	__m128* psc;
	__m128* pdcv;
} DATA_COPS32;

/** 
 * COPS processing function. Executes Viterbi algorithm on the 4-sequence batch
 * 	\arg COPS data
 *	\arg Array with 4 sequences to score
 *	\arg Result array with 4 float values [out]
 */
int p7_ViterbiCOPSf_run(DATA_COPS32* dcops, SEQ **seqsdb, float* results);

/**
 *	Setup COPS. Needs to be called before calling the processing function
 *	\arg COPS data
 *	\arg Estimate of the maximum Sequences' length. Used to optimize the memory allocation. Defaults to 1000
 *	\arg Maximum partition length allowed. Must be a multiple of 4. Defaults to 120
 */
void p7_ViterbiCOPSf_Setup(DATA_COPS32* dcops, int maxL, int max_partition);

/**
 *	Create COPS data structure for the 32-bit float implementation
 *	\arg A valid P7_PROFILE
 */
DATA_COPS32* p7_ViterbiCOPSf_Create(P7_PROFILE *gm);


