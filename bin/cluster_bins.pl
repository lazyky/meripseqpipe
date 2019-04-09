#! /usr/bin/perl -w
use strict;

my $input=shift;

open FH,$input or die;

my $pos_pre=0;
my $chr_pre="chr";
my %peak;
my $idx=0;
while(<FH>){
	chomp;
	my @tmp=split "\t";
	if($tmp[1]-$pos_pre>1){
		$idx++;
		my $peakname="meyer_peak_".$idx;
		$peak{$peakname}{chr}=$tmp[0];
		$peak{$peakname}{start}=$tmp[1];
		$peak{$peakname}{end}=$tmp[2];
		$pos_pre=$tmp[2];
	}else{
		my $peakname="meyer_peak_".$idx;		
		$peak{$peakname}{end}=$tmp[2];
		$pos_pre=$tmp[2];
	}
}

foreach my $k (sort keys %peak){
	print $peak{$k}{chr}."\t".$peak{$k}{start}."\t".$peak{$k}{end}."\t".$k."\n" if $peak{$k}{end}-$peak{$k}{start}>=100;
}
