#!/usr/local/bin/perl5
# 
##############################################################################
#
# nrdb90.pl generates a representative set at 90 % identity threshold
#
# Copyright (C) 1997  Liisa Holm
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
# 
##############################################################################
#
# reference: Removing near-neighbour redundancy from large protein
#	sequence collections.  L. Holm & C. Sander, Bioinformatics 14:423-429 (1998)
#
#
# input: 
#	1st argument [trivial.rdb] - empty file or result from previous run
#	2nd argument [courant_nrdb] - fasta (pearson) format sequence database
#		(it is recommended to remove duplicates beforehand using NCBI's nrdb program)
# output: 
#	clean_nrdb - fasta format sequence database 
#	clusters - representatives + groupies
#	trivial.rdb - relational database
#
# work files:
#	new.rdb - sequences with undefined rep
#	old.rdb - sequences with defined rep
#	deleted.rdb - sequences not in courant_nrdb
#

$[=1;
$|=1;

# parameters should satisfy
# 	idecutoff>=1-1/ktup1
# 	fcutoff<=1-(1-idecutoff)*ktup2
#
# For example (ktup1=10, ktup2=5):
# - nrdb90 (default, 90 % identity threshold): idecutoff=0.90 fcutoff=0.50
# - nrdb95 (95 % identity): idecutoff=0.95 fcutoff=0.75
# - nrdb98 (98 % identity): idecutoff=0.98 fcutoff=0.90
# - nrdb99 (99 % identity): idecutoff=0.99 fcutoff=0.95
#

$ktup1=10;
$ktup2=5;
$fcutoff=0.5;
$idecutoff=0.9;
$gapfrac=0.05; # percentage of gaps allowed for diagonal hashing
$mintupfrac=0.02; # required tuple matches to check diagonal
$maxnkey=1500000; # memory limit
$maxislesize=2000; # island size limit
$rdbheader="\# Perl-RDB\n\# $0 $rdbfile $courant_nrdb \nseqno\trepno\tseqlen\tname\tseq\n8N\t8N\t8N\t40S\t40S\n";
$quiet=1;

chop($jobstart=`date`); 

if($#ARGV!=2) { die "Usage: $0 trivial.rdb '/data/research/bin/nrdb -l10 /data/research/db/nrdb \"//\" |' \n"; }
($rdbfile,$courant_nrdb)=@ARGV;

# in: rdbfile+fastafile; out: old.rdb new.rdb deleted.rdb [renumbered sequences]
&diffnrdb; 

$nsw=0; $nal=0;

open(IN,'<new.rdb');
while(<IN>) { last if(!/^\#/); } $_=<IN>; # skip header
while(<IN>) { 
	chop; 
	($iprot,$dummy1,$dummy2,$dummy3,$seq)=split(/\t/); 
	$seq{$iprot}=$seq; 
	$new{$iprot}=1; 
	push(@outpool,$iprot); 
}
close(IN);
print "$#outpool new sequences from new.rdb\n";
open(IN,'<old.rdb');
while(<IN>) { last if(!/^\#/); } $_=<IN>; # skip header
while(<IN>) { 
	($iprot,$rep)=split(/\t/); 
	next if($iprot!=$rep);
	push(@old,$iprot); 
	$rev{$iprot}=$iprot; 
}
close(IN);
print "$#old reps from old.rdb\n"; $n=0;
open(IN,'<old.rdb');
while(<IN>) { last if(!/^\#/); } $_=<IN>; # skip header
while(<IN>) { 
	chop; 
	($iprot,$rep,$dummy2,$dummy3,$seq)=split(/\t/); 
	$seq{$iprot}=$seq; 
	if($iprot!=$rep) { $n++; $rev{$rep}.=":$iprot"; } 
}
close(IN);
print "$n groupies from old.rdb\n";

@x=@outpool; @outpool=sort { $a <=> $b } @x; undef (@x);
while($outpool[$[]) {
	$iprot=shift(@outpool);
	$seq=$seq{$iprot}; $lenseq=length($seq); 
	&cluster($iprot,$lenseq,$seq);
	if($nkey>$maxnkey||$maxsize>$maxislesize||!$outpool[$[]) {
		print "$#outpool sequences left in outpool; $#old old representatives\n";
		# map oldreps to clusters
		$n=0; undef(%oldmates);
		foreach $jprot (@old) {
			$n++; if($n=~/000$/) { print "oldmates ... $n\n" unless $quiet;}
			next if(!$rev{$jprot}); ## splice replaced oldreps from @old!
			undef(%mate);
			(@x)=split(//,$seq{$jprot}); (@w)=splice(@x,1,$ktup1-1);
			while($x[$[]) {
				push(@w,shift(@x)); $key=join('',@w); shift(@w);
				$x=$dekahash{$key};
				if($x) { $mate{$current{$x}}++; } 
			}
			foreach(keys(%mate)) { $oldmates{$_}.="$jprot "; }
			undef(%mate);
		}
		undef(%dekahash); undef(%current); undef(%tnerruc); undef(%cluster);
		undef(%size);
		foreach $r (keys(%protclus)) { 
			print "protclus $r -> $protclus{$r} oldmates -> $oldmates{$r}\n" unless $quiet;
			# pentahash oldmates
			foreach $q (split(/\s+/,$oldmates{$q})) {
				&load_m($seq{$q},$ktup2,$q,*hash);
			}
			(@pool)=split(/\s+/,$protclus{$r}); 
			undef(@rep2); undef(@newreps);
			while($newreps[$[]||$pool[$[]) {
				if($pool[$[]) { 
					(@newreps)=@rep2; undef(@rep2);
					&pool_loop; 
				}
				@x=@newreps; @newreps=sort { $a <=> $b } @x; undef(@x);
				while($newreps[$[]) { 
					$rep2=shift(@newreps);
					&mates_m($rep2);
					if(!$skip) { push(@rep2,$rep2); }
				}
			}
			undef(%hash);
			foreach $q (@rep2) { delete $new{$q}; }
			push(@old,@rep2); 
			# print "cleared reps: @rep2\n";
		}
		undef(%oldmates); undef(%protclus);
		$iclus=0; $nkey=0; $maxsize=0;
		@x=@outpool; @outpool= sort { $a <=> $b } @x; undef(@x);
	}
}

open(OUT,">clusters"); open(OUT1,">clean_nrdb");
&rdb2clusters; 
close(OUT); close(OUT1);

print "$nsw calls to diagonal; $nal calls to align\n";

print "writing to $rdbfile and clean_nrdb\n";
open(OUT,">$rdbfile");
print OUT $rdbheader;
foreach $rep (keys(%rev)) { foreach(split(/\:/,$rev{$rep})) { $rep{$_}=$rep; } }
foreach $iprot (sort { $a <=> $b } keys(%seq)) {
	$l=length($seq{$iprot}); 
	print OUT "$iprot\t$rep{$iprot}\t$l\t$name{$iprot}\t$seq{$iprot}\n";
}
close(OUT);

chop($jobend=`date`);
print "timestamp: $jobstart ... $jobend\n";

exit(1);

###############################################################################

sub getused {
	local ($seq,$ktup,*used)=@_;
	# hash seq
	(@x)=split(//,$seq);
	(@w)=splice(@x,1,$ktup-1);
	while($x[$[]) { 
		push(@w,shift(@x)); $key=join('',@w); shift(@w);
		if(!defined($used{$key})) { 
			$used{$key}=1; 
		} else {
			$used{$key}++; 
		}
	}
}

sub load_m {
	local ($seq,$ktup,$label,*hash)=@_;
	# print "this is load_m $label $ktup ",length($seq),"\n";
	undef(%used);
	&getused($seq,$ktup,*used);
	foreach $key (keys(%used)) {
		$hash{$key}.="$label:$used{$key} ";
	}
}

##############################################################################

sub mates_m {
	local($iprot)=@_; # returns ilen minmatch window minwinmatch %mate $seq
	$seq=$seq{$iprot};
	$ilen=length($seq); 
	$minmatch=$ilen*$fcutoff; # require 50 % match of tetrapeptides -> 90%ide
	undef(%mate); 
	undef(%used); 
	&getused($seq,$ktup2,*used);
	foreach $key (keys(%used)) {
		$n0=$used{$key};
		# print "$key ($n0) -> $hash{$key}\n";
		foreach(split(/\s+/,$hash{$key})) {
			($jprot,$n)=split(/\:/);
			if($n0<$n) { $n=$n0; }
			$mate{$jprot}+=$n; 
		}
	}
	# only keep mates that are above threshold! compare sorted by matecount
	undef(@array);
	foreach $jprot (keys(%mate)) {
		$n=$mate{$jprot};
		if($iprot>$jprot&&$n>=$minmatch) {
			push(@array,$jprot); 
		} elsif ($iprot<$jprot) { # use 0.5*jlen-4 exact formula
			$x=$fcutoff*length($seq{$jprot})-$ktup2+1;
			if($n>=$x) { push(@array,$jprot); }
		}
	}
	$skip=0;
	foreach $jprot (sort { $mate{$b} <=> $mate{$a} } @array) {
		# print "testing $iprot jprot $jprot [$rev{$jprot}] ($new{$jprot})\n";
		next if(!$rev{$jprot}); # has been merged before
		$skipthis=0;
		$domflag=0;
		$n=$mate{$jprot};
		$jlen=length($seq{$jprot}); 
		if($ilen<=$jlen) { ## default: query is shorter than match
			$domflag=1;
			$skipthis=&diagonal($seq{$iprot},$seq{$jprot},$ktup2,$gapfrac,$fcutoff,$mintupfrac,$idecutoff); 
		} else { ## compare longer query to shorter match
			$domflag=2;
			$skipthis=&diagonal($seq{$jprot},$seq{$iprot},$ktup2,$gapfrac,$fcutoff,$mintupfrac,$idecutoff); 
		}
		# print "matched $iprot=$rev{$iprot} $jprot=$rev{$jprot} $n $ilen $jlen $skipthis\n";
		if($skipthis) { # means merger
			if($domflag==1) {
				&dissolve($iprot,$jprot,$jlen);		
				$skip=1;
			} else {
				if(!defined($rev{$iprot})) { $rev{$iprot}="$iprot"; } # replacement inside iterated pool_loop
				&dissolve($jprot,$iprot,$ilen);
			}
		}
		last if($skip);
	}
}

sub dissolve {
	local($iprot,$jprot,$jlen)=@_;
	# merge rep2==iprot with rep1==jprot
	$rev{$jprot}.=":$iprot";
	(@groupies)=split(/\:/,$rev{$iprot}); shift(@groupies);
	print "merge $iprot ($new{$iprot}) -> $jprot ($new{$jprot}) $jlen $#groupies groupies (domflag $domflag)\n" unless $quiet;
	delete $rev{$iprot};
	# dissolve groupies of iprot -> follow jprot else ->@pool
	foreach $groupie (@groupies) {
		# if($groupie<$jprot) { ## never possible since iprot is shorter than jprot
		#	$skipthis=0; 
		#} else {
			$skipthis=&diagonal($seq{$groupie},$seq{$jprot},$ktup2,$gapfrac,$fcutoff,$mintupfrac,$idecutoff);
		#}
		if($skipthis) { 
			$rev{$jprot}.=":$groupie"; 
			print "groupie $groupie followed $iprot -> $jprot\n" unless $quiet;
		} else { 
			if($new{$groupie}) { 
				push(@pool,$groupie); 
				print "groupie $groupie of $iprot pooled ($#pool)\n" unless $quiet;
			} else { 
				push(@outpool,$groupie);
				print "groupie $groupie of $iprot outpooled ($#outpool)\n" unless $quiet;
				$new{$groupie}=1;
			} 
		}
	}
}

##############################################################################
#              diagonal identity check                                       #
##############################################################################

sub diagonal {
	local ($seq1,$seq2,$ktup,$gapfrac,$fcutoff,$mintupfrac,$idecutoff)=@_;
	local ($ilen,$jlen,$minsumtup,$minsumres,$ires,@x,@w,%hash,%diag,$d,
		@d,@diags,$maxdiag,%mat,%best,$score);
	$ilen=length($seq1);
	$jlen=length($seq2);
	$mintup=$mintupfrac*$ilen; # diagonal must match at least 2 tups per 100 aa
	$minsumtup=$fcutoff*$ilen; # seq1 is shorter one!
	$minsumres=$idecutoff*$ilen; # seq1 is shorter one!
	# print "this is subroutine diagonal($ktup,$gapfrac,$fcutoff,$mintupfrac)  $ilen $jlen $minsumtup $minsumres\n";
	$nsw++;
	# step 1: find best gapfrac diagonals by ktup hashing
	(@x)=split(//,$seq1); (@w)=splice(@x,1,$ktup-1);
	$ires=0;
	while($x[$[]) {
		$ires++;
		push(@w,shift(@x)); $key=join('',@w); shift(@w);
#		ambiguous letters always mismatch!
		next if($key=~/[BXZ]/);
		if(!defined($hash{$key})) { $hash{$key}=$ires; } 
		else { $hash{$key}.=":$ires"; }
	}
	
	(@x)=split(//,$seq2); (@w)=splice(@x,1,$ktup-1);
	$jres=0;
	while($x[$[]) {
		$jres++;
		push(@w,shift(@x)); $key=join('',@w); shift(@w); 
#		ambiguous letters always mismatch!
		next if($key=~/[BXZ]/);
		foreach(split(/\:/,$hash{$key})) {
			$d=$_-$jres;
			$diag{$d}+=1;
		}
	}
	
	(@d)=sort { $diag{$b} <=> $diag{$a} } keys(%diag);
	$maxdiag=1+$gapfrac*$jlen;
	splice(@d,$maxdiag+1);
	$s=0; $above=0; undef(@diags); 
	foreach(@d) { 
		last if($diag{$_}<$mintup);
		if($diag{$_}>$minsumres) { # single diagonal
			# print "single diagonal $diag{$_} > $minsumres\n";
			return(1); 
		} 
		$s+=$diag{$_}; $f=$s/$minsumtup;
		# print "best diagonals: $_\t$diag{$_}\t$s\t$f\t$minsumtup\n"; 
		# build up alignment matrix
		push(@diags,$_);
	}
	if($#diags>0) { 
		($above)=&align_sw; 
	}
	# print "return $above from diagonal $#diags\n";
	return($above);
}	

sub align_sw {
	local($i,$ires,$jres,$xres);
	local(@mat,@prev,@best);
	# align few best diagonals -> return score>minsum
	$nal++;
	# print "aligned diagonals @diags\n";
	(@seq1)=split(//,$seq1);
	(@seq2)=split(//,$seq2);
	foreach $ires (1..$ilen) {
		$xres=$seq1[$ires];
		foreach $shift (1..$#diags) { 
			$jres=$ires-$diags[$shift];
			next if($jres<1 || $jres>$jlen);
			if($xres ne $seq2[$jres]) {
				$mat[$ires][$shift]=-0.15; # mismatch
			} else {
				$mat[$ires][$shift]=1; # identity
			}
			# print "$mat[$ires][$shift]\t$seq2[$jres]\t";
		}
		# print "$xres\n";
	}
	$globscore=0; $globcell="0:0";
	foreach $shift (1..$#diags) { 
		$prev[1][$shift]="0:0"; 
		$best[1][$shift]=$mat[1][$shift]; 
	}
	foreach $ires (2..$ilen) {
		foreach $shift (1..$#diags) {
			$s=$mat[$ires][$shift];
			$curdiag=$diags[$shift];
			$jres=$ires-$curdiag;
			# determine best place to come from
			# 	(ires-1,[1..jres-1])  \
			# 	(ires-1,jres-1)	       ---> (ires,jres)
			# 	([1..ires-1],jres-1)  /
			$bestscore=0;
			$bestcell="0:0";
			foreach $diag (1..$#diags) {
				if($diag==$shift) { # continue match
					$xres=$ires-1;
					$score=$s+$best[$xres][$diag];
				} elsif($diags[$diag]>$curdiag) { # deletion in query
					$xres=$ires-1;
					$score=$s+$best[$xres][$diag]-1;
				} else { # insert in query
					$xres=$jres+$diags[$diag];
					next if($xres<1 || $xres>$ilen);
					$score=$s+$best[$xres][$diag]-($ires-$xres);
				}
				if($score>$bestscore) {
					$bestcell="$xres:$diag"; 
					$bestscore=$score;
				}
			}
			$prev[$ires][$shift]=$bestcell;
			$best[$ires][$shift]=$bestscore;
			if($bestscore>$globscore) {
				$globscore=$bestscore;
				$globcell="$ires:$shift";
			}
			# print "$diags[$shift]:$best[$ires][$shift]\t";
		}
		# print "$ires\t$bestcell\n";
	}
	# traceback to calculate nide
	# print "globcell: $globcell\n";
	$nide=0;
	($ires,$shift)=split(/\:/,$globcell);
	while($ires>0) {
		if($mat[$ires][$shift]>0) { $nide++; }
		# print "$ires:$shift\t$nide\n";
		($ires,$shift)=split(/\:/,$prev[$ires][$shift]);
	}

	if($nide>$minsumres) { $x=1; } else { $x=0; }
	# print "score = $nide / $minsumres\treturn $x from align\n";
	return ($x);
}

##############################################################################
##############################################################################

sub cluster {
	local($iprot,$lenseq,$seq)=@_;
	#
	# dekapeptide dekahash
	#
        # ktups of query sequence
        (@x)=split(//,$seq); (@w)=splice(@x,1,$ktup1-1);
	undef(%merge);
	$iclus++;
        while($x[$[]) {
                push(@w,shift(@x)); ($key)=join('',@w); shift(@w);
                next if($key=~/[BXZ]/); # ignore ambiguity codes
		# lookup islands linked by query protein
		$xclus=$dekahash{$key};
		if(!$xclus) {
			$dekahash{$key}=$iclus;
			$nkey++;
		} else {
			next if($xclus==$iclus);
			$merge{$current{$xclus}}=1;
		}
        }
	(@merge)=sort {$a <=> $b} keys(%merge);
	# print "$nkey words; to merge clusters @merge with $iclus\n";
	$kclus=$merge[$[];
	if($kclus) { 
		# merge linked clusters
		foreach $jclus (@merge[$[+1..$#merge]) {
			# print "* * * jclus->current: $jclus, $current{$jclus}, $tnerruc{$jclus} kclus: $kclus <- $tnerruc{$kclus}\n";
			foreach(split(/\s+/,$protclus{$jclus})) { $cluster{$_}=$kclus; }
			$protclus{$kclus}.=$protclus{$jclus};
			$size{$kclus}+=$size{$jclus};
			delete $protclus{$jclus};
			# update short-link to current cluster
			foreach(split(/\s+/,$tnerruc{$jclus})) { 
				$current{$_}=$kclus; 
				# print "* * * current $_ -> $kclus\n";
			}
			$tnerruc{$kclus}.=$tnerruc{$jclus};
			delete $tnerruc{$jclus};
			$size{$jclus}=0;
		}
		# append iprot
		# print "* append $iprot -> $kclus\n";
		$cluster{$iprot}=$kclus;
		$protclus{$kclus}.="$iprot ";
		$current{$iclus}=$kclus;
		$tnerruc{$kclus}.="$iclus ";
		$size{$kclus}+=1;
		if($size{$kclus}>$maxsize) { $maxsize=$size{$kclus}; }
	} else  { # iprot is new unique
		# print "* $kclus * new $iprot -> $iclus\n";
	        $cluster{$iprot}=$iclus; 
	     	$protclus{$iclus}="$iprot ";
		$current{$iclus}=$iclus;
		$tnerruc{$iclus}="$iclus ";
		$size{$iclus}=1;
	}
}

##############################################################################
##############################################################################

sub diffnrdb {
	local(%ix,%xi,%rep,%hash,%newix,$oldix,%newrep,%oldrep,%name,%lenseq);
	# read in new nrdb -> %name %seq %lenseq %hash
	print "reading new: $courant_nrdb\n" unless $quiet;
	open(IN,"$courant_nrdb");
	$nprot=0; $t=0; $m=0;
	while(<IN>) {
	        chop;
	        if(/^>/) { 
			# save old entry
			if($nprot>0) { &save_entry; }
	                ($name)=/^>(.*)/; 
			$name=~s/\t/ /g; # tab is field separator, not allowed inside field
	                $seq=''; 
	                $nprot++; 
			$ix{$nprot}=$nprot;
	                if($nprot=~/0000$/) { print "new...$nprot\n" unless $quiet; }
	        }
	        else { s/\s+//mg; $seq.=$_; }
	}
	close(IN);
	# save last entry
	&save_entry;
	print "$nprot sequences saved; $t keys; $m multiple hash key occupancies\n";

	# read in old nrdb -> %newix{oldseqno}=new-seqno %oldix{new-seqno}=oldseqno
	if(-e $rdbfile) {
	  $ndel=0;
	  print "reading old: $rdbfile\n" unless $quiet;
	  open(IN,"<$rdbfile");
	  open(OUT,">deleted.rdb"); print OUT $rdbheader;
	  while(<IN>) { last if(!/^\#/); } $_=<IN>; # skip rdb header
	  $n=0;
	  while(<IN>) {
		chop($line=$_);
		$n++;
		($i,$rep,$l,$name,$seq)=split(/\t+/,$line);
		$rep{$i}=$rep;
		if($n=~/0000$/) { print "old...$n\n" unless $quiet; }
		($checksum)=&do_checksum($l,$seq);
		(@x)=split(/\:/,$hash{$checksum});
		# print "old-new $i $l $checksum vs. @x \n";
		$keep=0;
		foreach(@x) {
			if($seq eq $xseq{$_}) {
				$newix{$i}=$_;
				$oldix{$_}=$i;
				$keep=1;
				# print "old $i ==  $_ ($l) <=> ($lenseq{$_})\n";
				last;
			}
		}
		if(!$keep) { 
			$ndel++;
			print OUT "$i\t$rep\t$l\t$name\t$seq\n";
		}
	  }
	  close(IN); close(OUT);
	  print "$ndel deleted entries in deleted.rdb\n";
	}

	# transfer reps
	foreach $ix (keys(%lenseq)) {# ix keys=courant-seqno
		$oldix=$oldix{$ix};  
		if($oldix) { $oldrep=$rep{$oldix}; } else { $oldrep=0; }
		if($oldrep) { $newrep{$ix}=$newix{$oldrep}; } else { $newrep{$ix}=0; }
		# print "transfer $ix ($lenseq{$ix})-> $oldix -> $oldrep -> $newix{$oldrep} \n";
	}
	
	# renumber sorted by length -> %ix{courant-seqno}=sorted-seqno
	foreach $i (keys(%lenseq)) {
		if($newrep{$i}==$i) { $lenseq{$i}+=0.2; }
		elsif($newrep{$i}>0) { $lenseq{$i}+=0.1; }
	}
	$i=0;
	foreach (sort { $lenseq{$b} <=> $lenseq{$a} } keys(%lenseq) ) {
		$i++;
		$ix{$_}=$i;
		$xi{$i}=$_;
		# print "sorted $_ -> $i $lenseq{$_} \n";
	}

	# write renumbered sequences/reps to new.rdb, old.rdb
	print "writing renumbered sequences/reps to new.rdb, old.rdb\n";
	open(OUT,">old.rdb");
	open(OUT1,">new.rdb");
	$nnew=0; $nold=0;
	print OUT $rdbheader; print OUT1 $rdbheader;
	foreach $i (sort { $a <=> $b } keys(%ix) ) {
		$ix=$xi{$i};
		$_=$lenseq{$ix}; ($l)=/^(\d+)/;
		if($newrep{$ix}) {
			$nold++;
			print OUT "$i\t$ix{$newrep{$ix}}\t$l\t$name{$ix}\t$xseq{$ix}\n";
		} else {
			$nnew++;
			print OUT1 "$i\t$ix{$newrep{$ix}}\t$l\t$name{$ix}\t$xseq{$ix}\n";
		}
	}
	close(OUT); close(OUT1);
	print "$nnew new sequences in new.rdb\n$nold old sequences in old.rdb\n";
	undef(%ix); undef(%xi);
	undef(%rep); undef(%hash);
	undef(%newix); undef(%oldix);
	undef(%newrep); undef(%oldrep);
	undef(%xseq); undef(%name); undef(%lenseq);
}


sub save_entry {
	$xseq{$nprot}=$seq; 
	$name{$nprot}=$name;
	$l=length($seq);
	$lenseq{$nprot}=$l;
	($checksum)=&do_checksum($l,$seq);
	if(!defined($hash{$checksum})) { $hash{$checksum}="$nprot"; $t++; } 
	else { $hash{$checksum}.=":$nprot"; $m++; }
	# print "saved $nprot\t$checksum\t$hash{$checksum}\t$l\t$name\n";
}

sub do_checksum {
	local ($l,$seq)=@_;
	local($checksum,$k);
	$checksum=0;
	$k=$l/20; if($k<2) { $k=2; }
	while($seq) {
		$checksum+=unpack("%32C*",$seq);
		$seq=substr($seq,$k);
	}
	$checksum%=429496796; # 2^32;
	$checksum.=":$l";
	return($checksum);
}

###############################################################################
###############################################################################

sub pool_loop { # redefine representatives for dissolved groupies in set2
	# print "this is pool_loop dealing with $#pool sequences and $#rep2 rep2s: @rep2\n";
	# pentahash rep2s within island2
	local(%hash);
	undef(%hash);
	foreach $rep2 (@rep2) { # "clean" reps passed through rep_loop
		next if(!$rev{$rep2}); # has been merged before
		&load_m($seq{$rep2},$ktup2,$rep2,*hash);
	}
	# sort pool by length
	@x=@pool; @pool=sort { $a <=> $b } @x;
	while($pool[$[]) {
		$iprot=shift(@pool);
		&mates_m($iprot);
		# unique iprot within island2
		if(!$skip) {  # "clean" rep2 -> add to hash table!
			push(@newreps,$iprot);
			$rev{$iprot}="$iprot"; 
			# print "new rep $iprot from pool_loop $#newreps newreps\n" unless $quiet;
			&load_m($seq{$iprot},$ktup2,$iprot,*hash);
		}
	}
	undef(%hash);
	return;
}

##############################################################################
##############################################################################

sub rdb2clusters {
	open(IN,'<new.rdb');
	while(<IN>) { last if(!/^\#/); } $_=<IN>; # skip rdb header
	while(<IN>) {
		chop;
		($iprot,$dummy1,$l,$name)=split(/\t/);
		$lenseq{$iprot}=$l;
		$name{$iprot}=$name;
	}
	close(IN);
	open(IN,'<old.rdb');
	while(<IN>) { last if(!/^\#/); } $_=<IN>; # skip rdb header
	while(<IN>) {
		chop;
		($iprot,$dummy1,$l,$name)=split(/\t/);
		$lenseq{$iprot}=$l;
		$name{$iprot}=$name;
	}
	close(IN);
	
	$nclus=0;
	foreach $rep (sort { $rep{$a} <=> $rep{$b} } keys(%rev)) {
		(@prot)=sort { $a <=> $b } split(/\:/,$rev{$rep});
		$nclus++;
		# summarize information from names
		$info=&expert(@prot);
		print OUT "Cluster $nclus with $#prot members $rep $lenseq{$rep} $name{$rep} $info\n";
		print OUT1 "\>$name{$rep} $info\n$seq{$rep}\n";
		foreach(@prot) {
			print OUT "$nclus\t$_\t$lenseq{$_}\t$name{$_}\n";
		}
	}
	
}


sub expert {
	local(@prot)=@_;
	$size=$#prot;
	shift(@prot); # remove rep
	# aim: select informative annotation from cluster groupies
	$info='';
	# report first swissprot or pdb
	foreach(@prot) { if($name{$_}=~/\:swiss\|/) { $info.=$name{$_}; last; } }
	foreach(@prot) { if($name{$_}=~/\:pdb\|/) { $info.=$name{$_}; last; } }
	# exclude junk
	if(!$info) {
		foreach(@prot) {
			$x=$name{$_}; ($x)=~tr/[A-Z]/[a-z]/;
			next if($x=~/hypothetical/);
			next if($x=~/homolog/);
			next if($x=~/ orf /);
			next if($x=~/putative/);
			next if($x=~/open reading frame/);
			next if($x=~/ polyprotein /);
			next if($x=~/ cosmid /);
			next if($x=~/product\: \"[a-z]+[0-9]+\"/);
			$info.=$name{$_}; last;
		}
	}
	if($info) { return("\/\/info(N\=$size): $info"); } else { return(''); }
}

##############################################################################
##############################################################################

__END__

GNU GENERAL PUBLIC LICENSE

Version 2, June 1991 

Copyright (C) 1989, 1991 Free Software Foundation, Inc.  
59 Temple Place - Suite 330, Boston, MA  02111-1307, USA

Everyone is permitted to copy and distribute verbatim copies
of this license document, but changing it is not allowed.

Preamble

The licenses for most software are designed to take away your freedom to share and change it. By contrast, the
GNU General Public License is intended to guarantee your freedom to share and change free software--to
make sure the software is free for all its users. This General Public License applies to most of the Free
Software Foundation's software and to any other program whose authors commit to using it. (Some other Free
Software Foundation software is covered by the GNU Library General Public License instead.) You can apply
it to your programs, too. 

When we speak of free software, we are referring to freedom, not price. Our General Public Licenses are
designed to make sure that you have the freedom to distribute copies of free software (and charge for this
service if you wish), that you receive source code or can get it if you want it, that you can change the software
or use pieces of it in new free programs; and that you know you can do these things. 

To protect your rights, we need to make restrictions that forbid anyone to deny you these rights or to ask you to
surrender the rights. These restrictions translate to certain responsibilities for you if you distribute copies of the
software, or if you modify it. 

For example, if you distribute copies of such a program, whether gratis or for a fee, you must give the
recipients all the rights that you have. You must make sure that they, too, receive or can get the source code.
And you must show them these terms so they know their rights. 

We protect your rights with two steps: (1) copyright the software, and (2) offer you this license which gives you
legal permission to copy, distribute and/or modify the software. 

Also, for each author's protection and ours, we want to make certain that everyone understands that there is no
warranty for this free software. If the software is modified by someone else and passed on, we want its
recipients to know that what they have is not the original, so that any problems introduced by others will not
reflect on the original authors' reputations. 

Finally, any free program is threatened constantly by software patents. We wish to avoid the danger that
redistributors of a free program will individually obtain patent licenses, in effect making the program
proprietary. To prevent this, we have made it clear that any patent must be licensed for everyone's free use or
not licensed at all. 

The precise terms and conditions for copying, distribution and modification follow. 

TERMS AND CONDITIONS FOR COPYING, DISTRIBUTION
AND MODIFICATION

0. This License applies to any program or other work which contains a notice placed by the copyright holder
saying it may be distributed under the terms of this General Public License. The "Program", below, refers to
any such program or work, and a "work based on the Program" means either the Program or any derivative
work under copyright law: that is to say, a work containing the Program or a portion of it, either verbatim or
with modifications and/or translated into another language. (Hereinafter, translation is included without
limitation in the term "modification".) Each licensee is addressed as "you". 

Activities other than copying, distribution and modification are not covered by this License; they are outside its
scope. The act of running the Program is not restricted, and the output from the Program is covered only if its
contents constitute a work based on the Program (independent of having been made by running the Program).
Whether that is true depends on what the Program does. 

1. You may copy and distribute verbatim copies of the Program's source code as you receive it, in any medium,
provided that you conspicuously and appropriately publish on each copy an appropriate copyright notice and
disclaimer of warranty; keep intact all the notices that refer to this License and to the absence of any warranty;
and give any other recipients of the Program a copy of this License along with the Program. 

You may charge a fee for the physical act of transferring a copy, and you may at your option offer warranty
protection in exchange for a fee. 

2. You may modify your copy or copies of the Program or any portion of it, thus forming a work based on the
Program, and copy and distribute such modifications or work under the terms of Section 1 above, provided that
you also meet all of these conditions: 

     a) You must cause the modified files to carry prominent notices stating that you changed the files and the
     date of any change. 

     b) You must cause any work that you distribute or publish, that in whole or in part contains or is derived
     from the Program or any part thereof, to be licensed as a whole at no charge to all third parties under the
     terms of this License. 

     c) If the modified program normally reads commands interactively when run, you must cause it, when
     started running for such interactive use in the most ordinary way, to print or display an announcement
     including an appropriate copyright notice and a notice that there is no warranty (or else, saying that you
     provide a warranty) and that users may redistribute the program under these conditions, and telling the
     user how to view a copy of this License. (Exception: if the Program itself is interactive but does not
     normally print such an announcement, your work based on the Program is not required to print an
     announcement.) 

These requirements apply to the modified work as a whole. If identifiable sections of that work are not derived
from the Program, and can be reasonably considered independent and separate works in themselves, then this
License, and its terms, do not apply to those sections when you distribute them as separate works. But when you
distribute the same sections as part of a whole which is a work based on the Program, the distribution of the
whole must be on the terms of this License, whose permissions for other licensees extend to the entire whole,
and thus to each and every part regardless of who wrote it. 

Thus, it is not the intent of this section to claim rights or contest your rights to work written entirely by you;
rather, the intent is to exercise the right to control the distribution of derivative or collective works based on
the Program. 

In addition, mere aggregation of another work not based on the Program with the Program (or with a work
based on the Program) on a volume of a storage or distribution medium does not bring the other work under the
scope of this License. 

3. You may copy and distribute the Program (or a work based on it, under Section 2) in object code or
executable form under the terms of Sections 1 and 2 above provided that you also do one of the following: 

     a) Accompany it with the complete corresponding machine-readable source code, which must be
     distributed under the terms of Sections 1 and 2 above on a medium customarily used for software
     interchange; or, 

     b) Accompany it with a written offer, valid for at least three years, to give any third party, for a charge
     no more than your cost of physically performing source distribution, a complete machine-readable copy
     of the corresponding source code, to be distributed under the terms of Sections 1 and 2 above on a
     medium customarily used for software interchange; or, 

     c) Accompany it with the information you received as to the offer to distribute corresponding source
     code. (This alternative is allowed only for noncommercial distribution and only if you received the
     program in object code or executable form with such an offer, in accord with Subsection b above.) 

The source code for a work means the preferred form of the work for making modifications to it. For an
executable work, complete source code means all the source code for all modules it contains, plus any
associated interface definition files, plus the scripts used to control compilation and installation of the
executable. However, as a special exception, the source code distributed need not include anything that is
normally distributed (in either source or binary form) with the major components (compiler, kernel, and so on)
of the operating system on which the executable runs, unless that component itself accompanies the executable. 

If distribution of executable or object code is made by offering access to copy from a designated place, then
offering equivalent access to copy the source code from the same place counts as distribution of the source code,
even though third parties are not compelled to copy the source along with the object code. 

4. You may not copy, modify, sublicense, or distribute the Program except as expressly provided under this
License. Any attempt otherwise to copy, modify, sublicense or distribute the Program is void, and will
automatically terminate your rights under this License. However, parties who have received copies, or rights,
from you under this License will not have their licenses terminated so long as such parties remain in full
compliance. 

5. You are not required to accept this License, since you have not signed it. However, nothing else grants you
permission to modify or distribute the Program or its derivative works. These actions are prohibited by law if
you do not accept this License. Therefore, by modifying or distributing the Program (or any work based on the
Program), you indicate your acceptance of this License to do so, and all its terms and conditions for copying,
distributing or modifying the Program or works based on it. 

6. Each time you redistribute the Program (or any work based on the Program), the recipient automatically
receives a license from the original licensor to copy, distribute or modify the Program subject to these terms
and conditions. You may not impose any further restrictions on the recipients' exercise of the rights granted
herein. You are not responsible for enforcing compliance by third parties to this License. 

7. If, as a consequence of a court judgment or allegation of patent infringement or for any other reason (not
limited to patent issues), conditions are imposed on you (whether by court order, agreement or otherwise) that
contradict the conditions of this License, they do not excuse you from the conditions of this License. If you
cannot distribute so as to satisfy simultaneously your obligations under this License and any other pertinent
obligations, then as a consequence you may not distribute the Program at all. For example, if a patent license
would not permit royalty-free redistribution of the Program by all those who receive copies directly or
indirectly through you, then the only way you could satisfy both it and this License would be to refrain entirely
from distribution of the Program. 

If any portion of this section is held invalid or unenforceable under any particular circumstance, the balance of
the section is intended to apply and the section as a whole is intended to apply in other circumstances. 

It is not the purpose of this section to induce you to infringe any patents or other property right claims or to
contest validity of any such claims; this section has the sole purpose of protecting the integrity of the free
software distribution system, which is implemented by public license practices. Many people have made
generous contributions to the wide range of software distributed through that system in reliance on consistent
application of that system; it is up to the author/donor to decide if he or she is willing to distribute software
through any other system and a licensee cannot impose that choice. 

This section is intended to make thoroughly clear what is believed to be a consequence of the rest of this
License. 

8. If the distribution and/or use of the Program is restricted in certain countries either by patents or by
copyrighted interfaces, the original copyright holder who places the Program under this License may add an
explicit geographical distribution limitation excluding those countries, so that distribution is permitted only in
or among countries not thus excluded. In such case, this License incorporates the limitation as if written in the
body of this License. 

9. The Free Software Foundation may publish revised and/or new versions of the General Public License from
time to time. Such new versions will be similar in spirit to the present version, but may differ in detail to
address new problems or concerns. 

Each version is given a distinguishing version number. If the Program specifies a version number of this
License which applies to it and "any later version", you have the option of following the terms and conditions
either of that version or of any later version published by the Free Software Foundation. If the Program does
not specify a version number of this License, you may choose any version ever published by the Free Software
Foundation. 

10. If you wish to incorporate parts of the Program into other free programs whose distribution conditions are
different, write to the author to ask for permission. For software which is copyrighted by the Free Software
Foundation, write to the Free Software Foundation; we sometimes make exceptions for this. Our decision will
be guided by the two goals of preserving the free status of all derivatives of our free software and of promoting
the sharing and reuse of software generally. 

NO WARRANTY

11. BECAUSE THE PROGRAM IS LICENSED FREE OF CHARGE, THERE IS NO WARRANTY FOR
THE PROGRAM, TO THE EXTENT PERMITTED BY APPLICABLE LAW. EXCEPT WHEN
OTHERWISE STATED IN WRITING THE COPYRIGHT HOLDERS AND/OR OTHER PARTIES
PROVIDE THE PROGRAM "AS IS" WITHOUT WARRANTY OF ANY KIND, EITHER EXPRESSED OR
IMPLIED, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE ENTIRE RISK AS TO THE
QUALITY AND PERFORMANCE OF THE PROGRAM IS WITH YOU. SHOULD THE PROGRAM
PROVE DEFECTIVE, YOU ASSUME THE COST OF ALL NECESSARY SERVICING, REPAIR OR
CORRECTION. 

12. IN NO EVENT UNLESS REQUIRED BY APPLICABLE LAW OR AGREED TO IN WRITING WILL
ANY COPYRIGHT HOLDER, OR ANY OTHER PARTY WHO MAY MODIFY AND/OR
REDISTRIBUTE THE PROGRAM AS PERMITTED ABOVE, BE LIABLE TO YOU FOR DAMAGES,
INCLUDING ANY GENERAL, SPECIAL, INCIDENTAL OR CONSEQUENTIAL DAMAGES ARISING
OUT OF THE USE OR INABILITY TO USE THE PROGRAM (INCLUDING BUT NOT LIMITED TO
LOSS OF DATA OR DATA BEING RENDERED INACCURATE OR LOSSES SUSTAINED BY YOU OR
THIRD PARTIES OR A FAILURE OF THE PROGRAM TO OPERATE WITH ANY OTHER
PROGRAMS), EVEN IF SUCH HOLDER OR OTHER PARTY HAS BEEN ADVISED OF THE
POSSIBILITY OF SUCH DAMAGES. 

END OF TERMS AND CONDITIONS

