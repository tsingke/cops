seqfile=~/dbs/human-dna.fasta
hmmspath=~/dbs/dfam-hmms
arch=UMA
partition=112

nn=12300
nreal=2000
compsize=5	# factor to control total computation load/size/time. higher => less comp.


export COPS_ARCH=$arch
export COPS_PARTITION=$partition

echo "RUN with seqfile $seqfile"

for tt in 1 2 4 6 8; do

	echo "Running with $tt threads"

	export COPS_NTHREADS=$tt
	make compile
	make compile-incr-versions

	pathlen=${#hmmspath}
	pathlen=$(($pathlen + 2))

	for ff in $hmmspath/*
	do
		hmmlen=${ff:$pathlen:4}
		hmmlen="$(echo $hmmlen | sed 's/0*//')"
	#	nuse=$((($nn * $tt) / $hmmlen ))
	
		ngood=$(( ($nn * 3000)/($hmmlen * $compsize ) ))
		check1=$(( $ngood * $hmmlen ))
		
		rounds=$(( $ngood / $nreal ))
		check2=$(( $rounds * $nreal * $hmmlen ))
		echo Rounds $rounds
		#echo Check2 $check2

		nuse=$(( $nreal * $tt ))
		
		if [ $rounds != 1 ]; then
			rounds=$(( $rounds / 2 ))
		fi

	if [ $tt = 1 ]; then	
		./viterbicops-initial		-N $nuse	-R $rounds	$ff  $seqfile
		./viterbicops-inlined		-N $nuse	-R $rounds	$ff  $seqfile
		./viterbicops-partitioned	-N $nuse	-R $rounds	$ff  $seqfile
	fi

	#	./viterbiserial 	    -N $nuse -R $rounds $ff  $seqfile
        ./viterbihmmer -c       -N $nuse -R $rounds $ff  $seqfile
		./viterbicops  -c       -N $nuse -R $rounds $ff  $seqfile
		./viterbicopsfloats -c  -N $nuse -R $rounds $ff  $seqfile
	done
done


