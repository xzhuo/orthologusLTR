use strict;
use warnings;
my $usage = "perl script.pl csv";
my $csv = shift or die $usage;
my %total;
open (IN, "<$csv");
open (OUT, ">$csv.out");
print OUT "chr\tchr_start\tchr_end\tchr_start\tchr_end\terv\terv_start\t_erv_end\terv_start\terv_end\n";
while (<IN>){
	next if $_ =~ /^#/ or $_ eq "\n";;
	chomp;
	my @position = split /,/, $_;
	my ($name,$chr,$rep_start,$rep_end,$chr_start,$chr_end,$score) = ($position[0],$position[1],$position[6],$position[7],$position[8],$position[9],$position[11]);
	my ($erv,$flank) = ($name =~ /_(\d+)_(\d+)$/); #$flank can be start or end.
	my $entry_ref;
	$entry_ref->{'chr'} = $chr;
	$entry_ref->{'chr_start'} = $chr_start;
	$entry_ref->{'chr_end'} = $chr_end;
	$entry_ref->{'rep_start'} = $rep_start;
	$entry_ref->{'rep_end'} = $rep_end;
	$entry_ref->{'score'} = $score;
	if (exists $total{$erv}{$flank}){
		push @{$total{$erv}{$flank}}, $entry_ref;
	}
	else{
		$total{$erv}{$flank} = [$entry_ref];
	}
	#the hash table structure:
	#1:{start:[entry1,entry2,...],end:[entry1,entry2...]}
	#2:{...} 3:{...}...
}

close IN;
foreach my $erv (sort keys %total){
	my @start = ();
	my @end = ();
	@start = @{$total{$erv}{'1'}} if defined $total{$erv}{'1'};
	@end = @{$total{$erv}{'2'}} if defined $total{$erv}{'2'};
	foreach my $hit (@start){
		foreach my $target (@end){
			if ($target->{'chr'} eq $hit->{'chr'} and abs($target->{'chr_start'} - $hit->{'chr_end'}) < 10000){
				print OUT "$hit->{'chr'}\t$hit->{'chr_start'}\t$hit->{'chr_end'}\t$target->{'chr_start'}\t$target->{'chr_end'}\t";
				print OUT "$erv\t$hit->{'rep_start'}\t$hit->{'rep_end'}\t$target->{'rep_start'}\t$target->{'rep_end'}\n";
			}
		}
	}
}

close OUT;

exit;
