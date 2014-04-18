The program striped is a test framework to time different Smith-Waterman implementations.


BUILDING:

The makefile support three targets:

    x86 - build an SSE2 enabled executable for an x86 processor
    sim - build an executable for the Cell B.E. simulator
    cell - build an executable for the Cell B.E. hardware

For the Cell B.E. builds, the gcc and xlc compilers are supported with the COMPILER flag.

    make sim COMPILER=xlc
    make cell COMPILER=gcc
    make x86


RUNNING:

striped [-h] [-i num] [-e num] [-t num] [-c num] [-T num] [-M num] matrix query db
    -h       : this help message
    -i num   : gap init penalty (default -10)
    -e num   : gap extension penalty (default -2)
    -t num   : minimum score threshold (default 20)
    -c num   : number of scores to be displayed (default 250)
    -T num   : number of calculation threads
    -M num   : maximum length of database sequence

    matrix   : scoring matrix file
    query    : query sequence file (fasta format)
    db       : sequence database file (fasta format) 

EXAMPLES:

To build for an x86 and run using two threads and displaying the top 25 matches:

   make x86
   striped -T 2 -c 25 blosum50.mat qseq.fasta dbseq.fasta 

To build for a PS3 and run using all six SPUs and displaying the top 10 matches:

   make cell
   striped -T 6 -c 10 blosum50.mat qseq.fasta dbseq.fasta 
