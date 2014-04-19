#!/bin/bash

#gcc vitfilter_test_threaded.c -O3 -o vitfilterthr -std=gnu99 -g -Wall -msse2  -I.. -L..  -I../impl_sse -L../impl_sse -I../../easel -L../../easel -Dp7VITFILTER_BENCHMARK -lhmmer -lhmmerimpl -lm -leasel -lpthread -DNTHREADS=1

gcc -g -Wall -O3 -std=gnu99 -o viterbicops -I. -I.. -L.. -I../impl_sse -L../impl_sse -I../../easel -L../../easel viterbi_cops.c viterbi_serial.c -lhmmer -lhmmerimpl -leasel -lm -lpthread -DNTHREADS=1

gcc -g -Wall -Wextra -O3 -std=gnu99 -o viterbi_cops-float -I. -I.. -L.. -I../impl_sse -L../impl_sse -I../../easel -L../../easel viterbi_cops_float.c viterbi_serial.c -lhmmer -lhmmerimpl -leasel -lm -lpthread -DNTHREADS=1 

exec=viterbicops
startM=4000
endM=5000
totalN=1250000

aaa=$startM
while [ $aaa -le 5000 ]
do
	sudo opcontrol --reset
	sudo opcontrol --start

	if   [ $aaa -lt 1000 ]; then dirname='0000-0999';
	elif [ $aaa -lt 3000  ]; then dirname='1000-2999';
	else dirname='3000-9999';
	fi
	
	n=$(($totalN/$aaa))

	./$exec -N $n /home/miguel/dbs/conserv/conserv-$dirname/conserv-$aaa.hmm dgfh
	if [[ $? -ne 0 ]]; then exit; fi;

	sudo opcontrol --shutdown
#	sudo opcontrol --stop
#	sudo opcontrol --dump
	sudo opreport -n
	
	if   [ $aaa -lt 200  ]; then aaa=$(($aaa+10)) ;
	elif [ $aaa -lt 500  ]; then aaa=$(($aaa+20)) ;
	elif [ $aaa -lt 1000 ]; then aaa=$(($aaa+50)) ;
	else aaa=$(($aaa+100)) ;
	fi
done

