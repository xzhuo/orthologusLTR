#to annotate the present/absent of all orthologus ERV/soloLTR in the given genome

use strict;
use warnings;
my $file = shift;
my $flank_len = 300;
my $int_len = 100; #flanking and internal length are consistent with previous setting.
open IN, "<$file";
my $gap = 50; #the acceptable missing length at both ends.
my $ERV_min_len = 4000;
my $LTR_min_len = 100;
my $LTR_max_len = 800;
my $min_hit = 150;
my $space = 10;
my %ortholog = ();
<IN>; #ignore headline
while(<IN>){
	chomp;
	my @pos = split();
	my ($erv, $start, $end, $L_start, $L_end, $R_start, $R_end) = ($pos[5],$pos[2],$pos[3],$pos[6],$pos[7],$pos[8],$pos[9]); #$start is approximately real start + $int_len, $end is real end -$int_len for full length ERV.
	my $length = abs($end-$start);
	if($L_start < $gap && $R_end > $flank_len + $int_len - $gap){
		if($L_end > $flank_len+$int_len -$gap && $R_start < $gap){ #means ortholog identified
			if ($length > $ERV_min_len){ #full length ERV
				#annotate as fulllength
				if (exists $ortholog{$erv}{"present"}){
					$ortholog{$erv}{"present"} +=1;
				}
				else{
					$ortholog{$erv}{"present"} = 1;
				}
			}
			elsif($length > $LTR_min_len && $length < $LTR_max_len){ #soloLTR
				#annotate as soloLTR
				if (exists $ortholog{$erv}{"present"}){
					$ortholog{$erv}{"present"} +=1;
				}
				else{
					$ortholog{$erv}{"present"} = 1;
				}
			}
		}
		elsif($L_end > $min_hit && $L_end < $flank_len + $space && $R_start < $flank_len+$int_len - $min_hit && $R_start > $int_len - $space){ #should be an empty site
			if($length < 20){
				#annotate as absent site
				if (exists $ortholog{$erv}{"absent"}){
					$ortholog{$erv}{"absent"} +=1;
				}
				else{
					$ortholog{$erv}{"absent"} = 1;
				}
			}
		}
	}
}
close IN;
open OUT, ">$file.ortholog.txt";
for my $LTR (keys %ortholog){
	if (exists $ortholog{$LTR}{"present"} and $ortholog{$LTR}{"present"} == 1){
		print OUT "$LTR\tpresent\n";
	}
	elsif(exists $ortholog{$LTR}{"absent"} and $ortholog{$LTR}{"absent"} == 1){
		print OUT "$LTR\tabsent\n";
	}
}
close OUT;
exit;
