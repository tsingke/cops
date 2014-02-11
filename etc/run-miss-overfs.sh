seqfile=~/dbs/human-dna.fasta
#/home/miguelzf/dbs/chimp-all.fasta

hmmspath=~/dbs/dfam-hmms

nn=25000000

for seqfile in ~/dbs/human-dna.fasta ~/dbs/macaca-dna.fasta
do
echo "RUN with seqfile $seqfile"

for tt in 1 2 3 4; do
	
	rounds=$((5 * $tt))
	
	gcc vitfilter_test_threaded.c viterbi_serial.c -O3 -o vitfilterthr -std=gnu99 -g -Wall -msse2  -I.. -L.. -I../../easel -L../../easel -Dp7VITFILTER_BENCHMARK -lhmmer  -lm -leasel -lpthread -DNTHREADS=$tt -DNROUNDS=$rounds

	gcc -g -Wall -O3 -std=gnu99 -o viterbicops -I. -I.. -L..  -I../../easel -L../../easel viterbi_cops.c viterbi_serial.c -lhmmer -leasel -lm -lpthread -DNTHREADS=$tt -DNROUNDS=$rounds

#	gcc -g -Wall -Wextra -O3 -std=gnu99 -o vitfloats -I. -I.. -L.. -I../../easel -L../../easel viterbi_cops_float.c viterbi_serial.c -lhmmer -leasel -lm -lpthread -DNTHREADS=$tt 

	gcc -g -Wall -O3 -std=gnu99 -o viterbi_cops-initial -I. -I.. -L.. -I../../easel -L../../easel viterbi_cops-initial.c -lhmmer -leasel -lm -DNROUNDS=$rounds
	gcc -g -Wall -O3 -std=gnu99 -o viterbi_cops-inlined -I. -I.. -L.. -I../../easel -L../../easel viterbi_cops-inlined.c -lhmmer -leasel -lm -DNROUNDS=$rounds
	gcc -g -Wall -O3 -std=gnu99 -o viterbi_cops-partitions -I. -I.. -L..  -L../impl_sse -I../../easel -L../../easel viterbi_cops-partitions.c -lhmmer -leasel -lm -DNROUNDS=$rounds

	echo "Running with $tt threads in $rounds rounds"

	pathlen=${#hmmspath}
	pathlen=$(($pathlen + 2))

	for ff in M0200* M0301*
#$hmmspath/*
	do
		hmmlen=${ff:1:4}
		hmmlen="$(echo $hmmlen | sed 's/0*//')"
		nuse=$(($nn / $hmmlen))
		
		if [ $tt = 1 ]
		then
		./viterbicops-initial	-N $nuse	$ff  $seqfile
		./viterbicops-inlined	-N $nuse	$ff  $seqfile
		./viterbicops-partitions	-N $nuse	$ff  $seqfile
		fi

		./vitfilterthr	-N $nuse $ff  $seqfile
		./viterbicops -N $nuse $ff  $seqfile
	#	./vitfloats	-N $nn	$ff  $seqfile
	done
done
done


