#!/bin/sh

###########################################################
#	Configurable parameters
###########################################################

dbfile=$1
#~/datasets/nrdb90.
hmmspath=$2
#~/datasets/hmms-proteins

# Max Partition length for Viterbi COPS
partition=120
# Desired cumulative length of all sequences to evaluate
runsize=5000000
# Max expected model length. Controls the number of rounds
maxhmmlen=3000
# Source directory
srcdir=src/


###########################################################
###########################################################

if [ $# -ne 2 ]; then
	echo "Must specify Database file and HMMs directory."
elif [ ! -f "$dbfile" ]; then
	echo "Database file $dbfile not found."
elif [ ! -d "$hmmspath" ]; then
	echo "HMMs directory $hmmspath not found."
fi

if  [ $# -ne 2 ] || [ ! -f "$dbfile" ] || [ ! -d "$hmmspath" ]; then
	echo "Usage ./run-benchmarks  <DB fasta file>  <directory with HMMs>"
	exit 1
fi
	
echo "=== Run Benchmarks with Database $dbfile and HMMs $hmmspath"

dbsize=`wc -c < $dbfile`
dbnseqs=`grep "^>" "$dbfile" | wc -l`
if [ $dbnseqs -eq 0 ]; then
	echo "No sequences found in database."
	exit 1
fi
avgseqlen=$(( $dbsize / $dbnseqs ))
runseqs=$(( $runsize / $avgseqlen ))

echo "=== Database with size $dbsize, $dbnseqs sequences of average length $avgseqlen. Use $runseqs sequences"

export COPS_PARTITION=$partition
make -C $srcdir

echo "=== All programs compiled"

for ff in "$hmmspath"/*.hmm
do

	hmmlen=`sed -n 's/LENG *\([0-9]*\)/\1/p' "$ff"`
	
	if [ -z $hmmlen ] || [ "$hmmlen" -ne "$hmmlen" ]; then	
		# not an integer, invalid length
		echo "File $ff not a valid HMM-3.1"
		continue
	fi
		
	# Control number of rounds to even out the longer runtimes of larger models
	nrounds=$(( $maxhmmlen / $hmmlen ))
	if [ $nrounds -lt 1 ]; then
	   nrounds=1
	fi
	
	echo "Call with $nrounds rounds"
		
	$srcdir/viterbihmmer		 -N $runseqs -R $nrounds $ff $dbfile
	$srcdir/viterbicops		 -N $runseqs -R $nrounds $ff $dbfile
#	./viterbicops-floats -N $runseqs -R $nrounds $ff $dbfile
done


