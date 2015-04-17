use strict;
use warnings;
use Bio::DB::Fasta;
sub flanking_start{
	my $a = shift;
	my $len = shift;
	$a>$len?$a-$len:1;
}
sub flanking_end{
	my $a = shift;
	my $b = shift;
	my $len = shift;
	$a < $b - $len?$a + $len:$b;
}
my $len = 100; #the length of flanking region before and after the boundary.
my $flank_len = 300; # the length out side of ERV. final length = $len + $flank_len
my $usage = "perl script.pl coordinates genome";
my $file = shift;
my $fasta = shift;
my $db = Bio::DB::Fasta->new("$fasta", -reindex =>1);
open IN, "<$file" or die "$!";
open OUT, ">$file.out" or die "$!";
while (<IN>){
	chomp;
	my @pos = split /\t/,$_;
	my $chr = $pos[0];
	my $start = $pos[1];
	my $end = $pos[2];
	my $id = "$pos[3]"."_$."; #get line number here
	my $obj = $db->get_Seq_by_id($chr);
	my $chr_len = $obj->length;
	my $flanking_start = flanking_start($start,$flank_len);
	my $flanking_end = flanking_end($end,$chr_len,$flank_len);
	my $seq1 = $obj->subseq($flanking_start => $start+$len);
	my $seq2 = $obj->subseq($end-$len => $flanking_end);
	print OUT ">$id"."_1\n$seq1\n>$id"."_2\n$seq2\n";
}
close IN;
close OUT;
exit;
