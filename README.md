COPS: Cache Oblivious Parallel SIMD Viterbi
=====================================

A cache-oblivious SSE-based implementation of the Viterbi algorithm for Hidden Markov Models, using inter-task parallelism on the SSE units of x86 processors, and additionally multi-threading. Developed for the latest version (3.1b1) of the HMMER suite: http://hmmer.janelia.org/

COPS is mainly a research project, developed as the final application of my mas bter's thesis work. It was the subject of a paper submitted to the BMC journal of Bioinformatics, one the most prestigious bioinformatics papers worldwide. Currently (February 2014) the paper is in its final peer-review.

A simplified, single-threaded, version of the project, the version presented in the research paper, is available here:
https://kdbio.inesc-id.pt/~lsr/COPS/#home


The version published on the website is the one submitted to BMC Bioinformatics, and does not include the multi-threaded version (since this version is not so clearly universally competitive vs non-multi-threaded implementations)



**Build & Install:**

1) Download [HMMER latest version (3.1b1)](http://selab.janelia.org/software/hmmer3/3.1b1/hmmer-3.1b1-linux-intel-ia32.tar.gz)

2) Download the present COPS code.

3) Build HMMER and install

4) Copy the folder and header `impl_sse/impl_sse.h` from the HMMMER distribution, over to a global include dir (i.e. /usr/local/include).   
(This is a bug in HMMER's install script, which does not create the 'impl_sse' in the include directory as it should.)

5) Type `make` to build COPS

6) Unzip the datasets.zip archive, and use the script `fetch-databases.sh` to obtain some testing instances


The programs may be compiled with some specific options to optimize their performance. These options can be enabled by exporting as environment variables:

   >export OPTION=VALUE
    
The options are:

	COPS_NTHREADS   Number of threads. Defaults to the number of available cores   
   
	COPS_PARTITION  Maximum partition length. Default is 120

	COPS_ARCH       Type of multi-core architecture. Default is UMA (unified memory access), NUMA (non-uniform memory access) is also available



**Run:**

> < prog >  [options]  < HMM model file >  < DB sequence file >

The tools may also run against a set of randomly generated sequences. To enable this, specify a non-existent sequence file.

Options allowed:

    -c             Check results against trusted implementation

    -N < int >     No. of seqs to generate

    -L < int >     Sequence length of generated sequences

    -M < int >     Maximum partition length. Only for COPS implementations. Defaults to 112

    -R < int >     No. of rounds to search the whole sequence set. Useful for benchmarking speeds. Default is 1



### References


**SIMD for sequence alignment:**

[Introduction to SIMD computing in the Cell processor](https://www.kernel.org/pub/linux/kernel/people/geoff/cell/ps3-linux-docs/CellProgrammingTutorial/BasicsOfSIMDProgramming.html)

[Sequence Alignment](https://en.wikipedia.org/wiki/Sequence_alignment)

[Striped SIMD alignment](http://bioinformatics.oxfordjournals.org/content/23/2/156.abstract)

[Inter-sequence SIMD alignment](http://dna.uio.no/swipe)


**Hidden Markov Models:**

[Introduction to Hidden Markov Models](https://en.wikipedia.org/wiki/Hidden_Markov_model)

[Probabilistic Models in Biological Sequence Analysis](http://selab.janelia.org/cupbook.html)

[HMMER suite website](http://hmmer.janelia.org)

[HMMER3: Accelerated Profile HMM Searches (2011)](http://www.ploscompbiol.org/article/info%3Adoi%2F10.1371%2Fjournal.pcbi.1002195)



**Intel's SIMD (SSE/AVX):**

[Intel's Developer Manual](http://www.intel.com/content/www/us/en/processors/architectures-software-developer-manuals.html)

[Intel SSE Tutorial : An Introduction to the SSE Instruction Set](http://neilkemp.us/src/sse_tutorial/sse_tutorial.html)

[SSE Primer](http://tommesani.com/index.php/component/content/article/2-simd/35-sse-primer.html)

[Intel Compiler Intrinsics](https://software.intel.com/sites/products/documentation/doclib/iss/2013/compiler/cpp-lin/index.htm#GUID-7478B278-2240-44D8-B396-1DC508E3656E.htm)

[Intel Intrinsics Guide](https://software.intel.com/sites/landingpage/IntrinsicsGuide)


[MSDN Compiler Intrinsics](http://msdn.microsoft.com/en-us/library/26td21ds%28v=vs.100%29.aspx)

[Intel ISA Extensions Programming Reference](https://software.intel.com/sites/default/files/managed/68/8b/319433-019.pdf)












