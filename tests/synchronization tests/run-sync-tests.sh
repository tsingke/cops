NN=2000	#2000*2000

gcc -g -Wall -O3 -std=gnu99 -o viterbistream -I. -I.. -L.. -I../impl_sse -L../impl_sse -I../../easel -L../../easel viterbi_stream-barrier.c viterbi_serial.c -lhmmer -lhmmerimpl -leasel -lm -lpthread -DNTHREADS=4

echo Barrier
./viterbistream -N 5000 BEX.hmm sghfgh

gcc -g -Wall -O3 -std=gnu99 -o viterbistream -I. -I.. -L.. -I../impl_sse -L../impl_sse -I../../easel -L../../easel viterbi_stream-sems.c viterbi_serial.c -lhmmer -lhmmerimpl -leasel -lm -lpthread -DNTHREADS=4

echo Semaphores
./viterbistream -N 5000 BEX.hmm sghfgh

gcc -g -Wall -O3 -std=gnu99 -o viterbistream -I. -I.. -L.. -I../impl_sse -L../impl_sse -I../../easel -L../../easel viterbi_stream-syncflags.c viterbi_serial.c -lhmmer -lhmmerimpl -leasel -lm -lpthread -DNTHREADS=4

echo Syncflags
./viterbistream -N 5000 BEX.hmm sghfgh

gcc -g -Wall -O3 -std=gnu99 -o viterbistream -I. -I.. -L.. -I../impl_sse -L../impl_sse -I../../easel -L../../easel viterbi_stream-counter.c viterbi_serial.c -lhmmer -lhmmerimpl -leasel -lm -lpthread -DNTHREADS=4

echo Counter
./viterbistream -N 5000 BEX.hmm sghfgh


