#score_ppt_tracts.pl

#This script accepts fasta file of introns as input and identifies/scores polypyrimidine tracts

#Script adapted from rules in: PMID11854178
#Searches for PPT in 3' terminal 40 nucleotides
#Nucleotides assigned values of -2 for A/G, +2 for C, +3 for U
#Positional weighting of 1 for +/- 2 positions, 2 for +/- positions and 3 for position 0
#Runs with >=2 are considered potential PPT
#Edges pruned/extended by one nt to ensure that it starts/ends with pyrimidine
#Runs of >9 nt are candidate PPT
#Runs of 5-9 nt with >4 uridines are candidate PPT

#Script outputs all candidate PPT in window and PPT scores

use warnings;

$/=">";

my %ppt;

#demo file - 100 introns in chromosome 1
open (FH, "hg19_introns_demo.fa") or die $!;

while (<FH>){
	if ($_ =~ /chr/){
	$line = $_;
	$line =~ s/\s//g; chomp $line;

	$line =~ /(chr\w+\:\d+\-\d+\([+-]\))([ACGTNacgtn]+)/;
	$coord = $1;
	
	
	$seq = $2;
	
	#convert string to uppercase
	$seq =~ tr/[a-z]/[A-Z]/;

	#only use U2 introns
	if ($seq =~ /^GT.+AG$/){

		$length = length($seq);

		#filter introns <40 nt in length
		if ($length >=40){

			$query = substr $seq,-40;
			my @query = split('',$query);
			my @perntscores;
			my @scores;

			#tile through 40 nucleotide window, identify central position scores in each position
			for (my $j = 0; $j<40;$j++){
				if (($query[$j]eq"A")||($query[$j]eq"G")){
					$perntscores[$j]= 0-2;
				}
				elsif ($query[$j] eq "C"){
					$perntscores[$j] = 2;
				}
				elsif ($query[$j] eq "T"){
					$perntscores[$j] = 3;
				}
				elsif ($query[$j] eq "N"){
					$perntscores[$j]=0;
				}
				else {print "error: $query[$j]\n";}
			}

			#add flanking scores to each position
			for (my $j = 0; $j<40;$j++){
				$min2 = $j-2;
				$min1 = $j-1;
				$plus1 = $j+1;
				$plus2 = $j+2;
				$runningscore = 0;
				if ($min2>=0){
					$runningscore = $runningscore+$perntscores[$min2];
				}
				if ($min1>=0){
					$runningscore = $runningscore+(2*($perntscores[$min1]));
				}
				if ($plus1<40){
					$runningscore = $runningscore+(2*($perntscores[$plus1]));
				}
				if ($plus2<40){
					$runningscore = $runningscore+$perntscores[$plus2];
				}
				$scores[$j] = $runningscore; 
			}
			$flag = 0;
			my $ppt=""; my $score;
			print $coord,"\t";
			my @ppt;
			
			#loop through positions to identify runs with scores >=2 and build PPT sequences
			for (my $j = 0;$j<40;$j++){
				#if central position has score >=2
				if ($scores[$j]>2){
					if ($flag==1){
						$ppt = $ppt.$query[$j];
						$score=$score+$scores[$j];
					}
					else {
						if (($query[$j] eq "C")||($query[$j] eq "T")){
						$ppt = $query[$j];
						$score=$scores[$j];
						$flag = 1;}
					}
				}

				else {
					#determine end of PPT, subtract edges to end on C/T
					if ($flag==1){
						$flag=0;
						$endflag = 0;
						while ($endflag==0){
							$lastnt = substr $ppt,-1;
							if (($lastnt eq "C")||($lastnt eq "T")){
								$endflag=1;
							}
							else {
								$length2 = length($ppt);
								$newppt = substr $ppt,0,($length2-1);
								$ppt = $newppt;
							}
						}
						#Decide whether candidate PPT passes threshold (either 9nt long, or 5-9 nt long with >=5 Ts)
						if (length($ppt)>=9){
							print $ppt,":",$score,",";
						}
						elsif (length($ppt)>=5){
							@ppt = split ('',$ppt);
							$uridines=0;	
							foreach (@ppt){
								if ($_ eq "T"){$uridines++;}
							}
							if ($uridines>4){
								print $ppt,":",$score,",";
							}
						}
						#Print associated score
						$ppt="";$score="";					
					}
				}
			}
			print "\n";
		}
	}}
	
}		
