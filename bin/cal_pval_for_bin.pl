#! /usr/bin/perl 
use strict;
#no warnings;
use Text::NSP::Measures::2D::Fisher::left;
use Text::NSP::Measures::2D::Fisher::right;

my $inputfile=shift;
my $IPrc=shift;
my $INPUTrc=shift;

open FH,$IPrc or die;
$IPrc=<FH>;
chomp($IPrc);
open FH,$INPUTrc or die;
$INPUTrc=<FH>;
chomp($INPUTrc);

#print $IPrc."\t".$INPUTrc."\n";
open FH,$inputfile or die;

while(<FH>){
	chomp;
	my @tmp=split "\t";

	my $input_pval=Text::NSP::Measures::2D::Fisher::left::calculateStatistic(n11 => $tmp[3],n1p => $tmp[3]+$tmp[7],np1 => $IPrc,npp=>$IPrc+$INPUTrc);
	my $ip_pval=Text::NSP::Measures::2D::Fisher::right::calculateStatistic(n11 => $tmp[3],n1p => $tmp[3]+$tmp[7],np1 => $IPrc,npp=>$IPrc+$INPUTrc);
	print $tmp[0]."\t".$tmp[1]."\t".$tmp[2]."\t".$ip_pval."\t".$input_pval."\n";
}






