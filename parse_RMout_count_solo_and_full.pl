use strict;
use Getopt::Std;
use warnings;
use List::Util qw[min max];
use List::MoreUtils qw(uniq);

my $usage = "perl parse_RMout_count_fulllength_LTR.pl RMoutfile <-l LTR> <-i INT>";
my $file = shift or die $usage;
my %opts;
getopts('i:l:', \%opts); #not really use it.
my @RM_array = ();
open(IN, "<$file");
# all lines in RMout staraged in @RM_array.

while(<IN>){
	chomp;
	my @line = split();
	my %RM_hash; #hash of each line
	if (defined $line[10] && ($line[10] eq "Unspecified" || $line[10] =~ /^LTR\// || $line[10] =~ /^ERV/) && $line[9] =~ m/ltr/i){ #only parse ltr sequence, ignore internal regions.
		my $INT;
		if($line[9] =~ m/int/){
			$INT = "INT";
		}
		else{
			$INT = "LTR";
		}
		if ($line[8] eq "+"){
			%RM_hash = ( #the hash of each line
					"div" => $line[1],
					"chr" => $line[4],
					"genostart" => $line[5],
					"genoend" => $line[6],
					"strand" => $line[8],
					"repname" => $line[9],
					"repfamily" => $line[10],
					"repstart" => $line[11],
					"repend" => $line[12],
					#"LI" => $hash{$line[10]},	if LI information comes from an additional file
					"LI" => $INT,	#if LI information obtained from repname directly
			);
		}
		if ($line[8] eq "-" or $line[8] eq "C"){
			%RM_hash = ( #the hash of each line
					"div" => $line[1],
					"chr" => $line[4],
					"genostart" => $line[5],
					"genoend" => $line[6],
					"strand" => "-",   #$line[8],
					"repname" => $line[9],
					"repfamily" => $line[10],
					"repstart" => $line[13],
					"repend" => $line[12],
					#"LI" => $hash{$line[10]},	if LI information comes from an additional file
					"LI" => $INT,	#if LI information obtained from repname directly
			);
		}
		push @RM_array, \%RM_hash;	#storage data in array of hashes
	}
}
close(IN);
@RM_array = sort {	$a->{"repname"} cmp $b->{"repname"} or
			$a->{"chr"} cmp $b->{"chr"} or 
			$a->{"genostart"} <=> $b->{"genostart"} or
			$a->{"genoend"} <=> $b->{"genoend"}
		} @RM_array;
#variable: 
#			$i: the upper internal region line number;
#			$j: the lower internal region line number;
#			$i -$m + 1: upper 20000 block line number;
#			$i + $n -1: lower 20000 block line number;
#			$xx: upper ERV LTR line number;
#			$yy: lower ERV LTR line number;
#parameters:
my $ERV_max_len = 10000;
my $ERV_min_len = 3000;
my $LTR_len = 1000; #so those ERVs with other TE insertion within LTR won't be identified here!!
my $LTR_gap_len = 100;
my @LTR_array = (); #the array for merged LTRs.
my %curr_hash = ();

##########################################################################
my $species = "pangolin"; #change this parameter accordingly!!!!



my %ERV_len = ("MLERV1_1.1" => 475,
		"MLERV1_2" => 615,
		"MLERV1_3.1" => 444,
		"MLERV1_3.2" => 432,
		"MLERV1_1.2" => 486,
		"ERV1-1_FCa-LTR" => 353,
		"manis_ltr" => 302,
		);
#next: merge fragmented soloLTRs (if any)
for (my $i = 0;$i<=$#RM_array;$i++){
	unless (%curr_hash && $RM_array[$i]{"chr"} eq $curr_hash{"chr"} && $RM_array[$i]{"genostart"} < $curr_hash{"genoend"} + $LTR_len && $RM_array[$i]{"strand"} eq $curr_hash{"strand"} && $RM_array[$i]{"repname"} eq $curr_hash{"repname"}){
		my %LTRline = %curr_hash;
		push @LTR_array, \%LTRline if %curr_hash;
		%curr_hash = %{$RM_array[$i]};
	}
	#merge LTR lines:
	elsif ($curr_hash{"strand"} eq "+" && $RM_array[$i]{"genostart"} - $curr_hash{"genoend"} < $LTR_gap_len && $RM_array[$i]{"repstart"} - $curr_hash{"repend"} < $LTR_gap_len){
		$curr_hash{"genoend"} = $RM_array[$i]{"genoend"};
		$curr_hash{"repend"} = $RM_array[$i]{"repend"};
	}
	elsif ($curr_hash{"strand"} eq "-" && $RM_array[$i]{"genostart"} - $curr_hash{"genoend"} < $LTR_gap_len &&  $curr_hash{"repstart"} - $RM_array[$i]{"repend"} < $LTR_gap_len){
		$curr_hash{"genoend"} = $RM_array[$i]{"genoend"};
		$curr_hash{"repstart"} = $RM_array[$i]{"repstart"};
	}
	else{
		print "hmm, check line $i \n";
		my %LTRline = %curr_hash;
		push @LTR_array, \%LTRline; 
		%curr_hash = %{$RM_array[$i]};
	}
}
push @LTR_array, \%curr_hash;

#parse @LTR_array to remove partial hits:
my @intact_array = ();
for (my $i = 0;$i<=$#LTR_array;$i++){
	#get intact LTR length:
	my $intact_len = $ERV_len{$LTR_array[$i]{"repname"}};
	push @intact_array, $LTR_array[$i] if ($LTR_array[$i]{"repend"} > $intact_len - 10 and $LTR_array[$i]{"repstart"} < 150); #5' end is generally more divergenced, so 100bp missing in 5' end is accepted; while >10bp missing in the 3' end is considered as incomplete LTR.
	#100bp in 5'end for MB,
	#150bp in 5'end for EF
}

open (BED, ">$file.sololtr.bed");
#print bed file of all soloLTR for phylogeny:
for (my $i =0;$i<=$#intact_array;$i++){
	my $start = $intact_array[$i]{'genostart'}-1; #transform to 0 based coordinate for bedtools.
	print BED "$intact_array[$i]{'chr'}\t$start\t$intact_array[$i]{'genoend'}\t$species"."_$intact_array[$i]{'repname'}_$i\t$intact_array[$i]{'div'}\t$intact_array[$i]{'strand'}\t$intact_array[$i]{'repstart'}\t$intact_array[$i]{'repend'}\n";
}
close BED;

open(OUT, ">$file".".all.bed");
#parse @intact_array to annotate full length ERVs:
for (my $i =0;$i<$#intact_array;$i++){
	my $j = $i+1;
	while($j<=$#intact_array and $intact_array[$j]{"chr"} eq $intact_array[$i]{"chr"} and $intact_array[$j]{"genoend"} < $intact_array[$i]{"genostart"} + $ERV_max_len){
		$j++;
	}
	if($j > $i+1){
		for ($i+1..$j-1){
			if($intact_array[$_]{"strand"} eq $intact_array[$i]{"strand"} and $intact_array[$_]{"repname"} eq $intact_array[$i]{"repname"} and $intact_array[$_]{"genoend"} > $intact_array[$i]{"genostart"} + $ERV_min_len){
				#print every pair of full length ERVs, 0 based bed format
				my $start = $intact_array[$i]{'genostart'}-1; #transform to 0 based coordinate for bedtools.
				print OUT "$intact_array[$i]{'chr'}\t$start\t$intact_array[$_]{'genoend'}\t$species"."_$intact_array[$i]{'repname'}_$i\t$intact_array[$i]{'div'}\t$intact_array[$i]{'strand'}\t$intact_array[$i]{'repstart'}\t$intact_array[$_]{'repend'}\tfull\n" if $intact_array[$i]{'strand'} eq "+";
				print OUT "$intact_array[$i]{'chr'}\t$start\t$intact_array[$_]{'genoend'}\t$species"."_$intact_array[$i]{'repname'}_$i\t$intact_array[$i]{'div'}\t$intact_array[$i]{'strand'}\t$intact_array[$_]{'repstart'}\t$intact_array[$i]{'repend'}\tfull\n" if $intact_array[$i]{'strand'} eq "-";
				$intact_array[$i]{"LI"} = $intact_array[$i]{"LI"}."full5";
				$intact_array[$_]{"LI"} = $intact_array[$_]{"LI"}."full3";
			}
		}
	}
}
#print chr\tstart\end\n for all soloLTRs and full length ERVs:
for (my $i =0;$i<=$#intact_array;$i++){
	if ($intact_array[$i]{"LI"} eq "LTR"){
		#0 based bed format
		my $start = $intact_array[$i]{'genostart'}-1; #transform to 0 based coordinate for bedtools.
		print OUT "$intact_array[$i]{'chr'}\t$start\t$intact_array[$i]{'genoend'}\t$species"."_$intact_array[$i]{'repname'}_$i\t$intact_array[$i]{'div'}\t$intact_array[$i]{'strand'}\t$intact_array[$i]{'repstart'}\t$intact_array[$i]{'repend'}\tsolo\n";
	}
}
close OUT;
exit;
