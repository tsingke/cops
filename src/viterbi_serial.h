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
 * Serial implementation of Viterbi decoding
 */

/**
 * Viterbi Serial, for any alignment mode
 */
int p7_Viterbi_general(ESL_DSQ *dsq, int L, P7_PROFILE *gm,  float *opt_res);

/**
 * Viterbi Serial optimized for Unilocal alignments
 */
int p7_Viterbi_unilocal(ESL_DSQ *dsq, int L, P7_PROFILE *gm,  float *opt_res);

/**
 * Viterbi Serial optimized for Unilocal alignments and discretized with 16bit integers
 * Emulates saturated arithmetic. Very slow.
 */
int p7_Viterbi_unilocal_word(ESL_DSQ *dsq, int L, P7_PROFILE *gm,  float *opt_ret);


