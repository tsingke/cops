.PHONY: example

# from shell environnment variable
LIBS=-lhmmer -leasel -lm -lpthread
DIRS=
#-I. -I.. -L.. -I../../easel -L../../easel
FLAGS=-Wall -Wextra -O3 -std=gnu99 -msse2
NPROCS=`nproc`

all: compile

build: compile

compile: setdefaultvals
	gcc ${FLAGS} -o viterbicops  viterbi_cops.c viterbi_serial.c ${DIRS} ${LIBS} -DCOPS_WORD_BENCHMARK -DMAX_PARTITION=${COPS_PARTITION} -DARCH=${COPS_ARCH} -DNTHREADS=${COPS_NTHREADS}
	gcc ${FLAGS} -o viterbicops-float  viterbi_cops_float.c viterbi_serial.c ${DIRS} ${LIBS} -DCOPS_FLOAT_BENCHMARK -DMAX_PARTITION=${COPS_PARTITION} -DARCH=${COPS_ARCH} -DNTHREADS=${COPS_NTHREADS}
	gcc ${FLAGS} -o viterbihmmer  vitfilter_test_threaded.c viterbi_serial.c ${DIRS} ${LIBS} -DARCH=${COPS_ARCH} -DNTHREADS=${COPS_NTHREADS}
	gcc ${FLAGS} -o viterbiserial  viterbi_serial.c ${DIRS} ${LIBS} -DSERIAL

compile-incr-versions: setdefaultvals
	gcc ${FLAGS} -o viterbicops-partitioned  viterbi_cops-partitioned.c viterbi_serial.c ${DIRS} ${LIBS} -DMAX_PARTITION=${COPS_PARTITION}
	gcc ${FLAGS} -o viterbicops-inlined  viterbi_cops-inlined.c viterbi_serial.c ${DIRS} ${LIBS}
	gcc ${FLAGS} -o viterbicops-initial  viterbi_cops-initial.c viterbi_serial.c ${DIRS} ${LIBS}


setdefaultvals:
ifeq (${COPS_PARTITION}, )
  COPS_PARTITION=120
endif
ifeq (${COPS_ARCH}, )
  COPS_ARCH=UMA
endif
ifeq (${COPS_NTHREADS}, )
  COPS_NTHREADS=${NPROCS}
endif


test-cops:	
	./viterbicops  -N 100000 -L 1000  datasets/hmms-dna/M1000-L1MEg2_5end.hmm dummy

test-hmmer:	
	./viterbihmmer -N 100000 -L 1000  datasets/hmms-dna/M1000-L1MEg2_5end.hmm dummy

check-cops:	
	./viterbicops  -N 10000 -L 1000 -c datasets/hmms-dna/M1000-L1MEg2_5end.hmm dummy

check-hmmer:	
	./viterbihmmer -N 10000 -L 1000 -c datasets/hmms-dna/M1000-L1MEg2_5end.hmm dummy

example: compile
	./viterbicops datasets/hmms-proteins/M0300-Aldose_epim.hmm  datasets/nrdb90
	./viterbicops-float -c -N 10000 datasets/hmms-dna/M1000-L1MEg2_5end.hmm  datasets/rna.fa
	./viterbihmmer -N 10000 -L 1000  datasets/hmms-proteins/M2203-P27732.hmm  dummyfile

tar-source:
	tar -czf sourceCOPS.tgz *.c *.h *.sh Makefile README Doxygen datasets/hmms* 

tar-bins:
	tar -czf binsCOPS.tgz  viterbicops  viterbicops-float  viterbihmmer

clean:
	rm -f viterbicops  viterbicops-float  viterbihmmer viterbiserial viterbicops-partitioned viterbicops-inlined viterbicops-initial  


