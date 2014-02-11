#!/bin/bash 


ii=1930
inc=0

while [ $ii -lt 1960 ]; do
for ii in 4576 5491; do

	~/hmmer-execs/hmmbuild conserv-$ii.hmm temp-$ii.sto

	for inc in {-5..5}
	do
		val=$(($ii + $inc))
		val=$1
		echo Generate $val

		#gcc -Wall gen.c -o gen-data -DTOTALSIZE=$ii
		./gen-data $val > temp.sto 2> log.txt
		cat log.txt | sed -n 's/.*empties\. Real size: \([0-9]\+\).*/ mv temp\.sto conserv-\1\.sto/p' | sh
		cat log.txt
		hmmbuild temp.hmm temp.sto | sed -n 's/1 \+temp \+[0-9]\+ \+[0-9]\+ \+\([0-9]\+\)/mv temp.hmm conserv-\1.hmm/p' | bash
		ii=$(($ii+1))
	done

	ii=$(($ii+1))
done



