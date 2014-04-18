
//	Viterbi Serial optimized general
int p7_Viterbi_general(ESL_DSQ *dsq, int L, P7_PROFILE *gm,  float *opt_res);

//  Viterbi Serial optimized for Unilocal alignments
int p7_Viterbi_unilocal(ESL_DSQ *dsq, int L, P7_PROFILE *gm,  float *opt_res);

//	Viterbi Serial optimized for Unilocal alignments and Discretized for 16bit integers
int p7_Viterbi_unilocal_word(ESL_DSQ *dsq, int L, P7_PROFILE *gm,  float *opt_ret);


