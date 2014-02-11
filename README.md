COPS: Cache Oblivious Parallel SIMD Viterbi
=====================================

A cache-oblivious SSE-based implementation of the Viterbi algorithm for Hidden Markov Models (HMMER), using inter-task parallelism on the SSE units of x86 processors, and multi-threading. Developed for the latest version (3.1b1) of the HMMER suite: http://hmmer.janelia.org/

COPS is mainly a research project, developed as the final application of my master's thesis work. It was the subject of a paper submitted to the BMC journal of Bioinformatics, one the most prestigious worldwide. Currently (February 2014) the paper is in its final peer-review.

A simplified, single-threaded, version of the project, that was presented in the research paper, is available here:
https://kdbio.inesc-id.pt/~lsr/COPS/#home



**Build & Install:**

1) Download HMMER 3.1b1 from here:
http://selab.janelia.org/software/hmmer3/3.1b1/hmmer-3.1b1-linux-intel-ia32.tar.gz

2) Download the present COPS code.

3) Build HMMER and install

4) Move the folder and header 'impl_sse/impl_sse.h' to the system include dir (i.e. /usr/local/include). 
This is a small bug in HMMER's install script, in which it does nto create the 'impl_sse' in the includes directory as it should.

5) Type 'make' to build COPS

6) Unzip the datasets.zip archive, and use the script 'fetch-databases.sh' to obtain some testing instances


The programs may be compiled with some specific options to optimize their performance. These options can be enabled by exporting as environment variables:

    export OPTION=VALUE
    
The options are:

    COPS_NTHREADS   number of threads. default is number of cores available
    
    COPS_PARTITION  max partition length. default is 120
    
    COPS_ARCH       type of multi-core architecture. Default is UMA (unified memory access), NUMA (non-uniform memory access) is also available
    


**Run:**

    ./ < prog >  [options]  < HMM model file >  < DB sequence file >

The tools may also run against a set of randomly generated sequences. To enable this, specify a non-existent sequence file.

Options allowed:

    -c             Check results against trusted implementation

    -N < int >     No. of seqs to generate

    -L < int >     Sequence length of generated sequences

    -M < int >     Maximum partition length. Only for COPS implementations. Defaults to 112

    -R < int >     No. of rounds to search the whole sequence set. Useful for benchmarking speeds. Default is 1



