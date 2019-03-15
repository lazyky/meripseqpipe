#! /usr/bin/perl -w
#perl m6A_annotate_forGTF.pl /data1/database/hg38/GENCODE/gencode.v25.annotation.gtf macs2/merged_Peak.bed macs2/merged_Peak

use strict;
use warnings;
use FindBin qw($Bin);

if (@ARGV < 3) {
	print "
		usage: perl m6A_annotate_forGTF.pl <in ref gene.gtf> <in merged_Peak.bed> <out prefix>
		program will create multi file: <out>.center <out>.anno.txt <out>.unanno.txt ...
		\n";
	exit;
}


my ($ref_gene_gtf, $peak_bed, $outPrefix) = @ARGV;

#make a gene type list
my %GeneType; #refer to ftp://ftp.sanger.ac.uk/pub/gencode/_README_stats.txt of 1) GENES. of https://www.gencodegenes.org/stats/current.html
foreach (qw (protein_coding)) {
	$GeneType{"$_"} = "mRNA"; }
foreach (qw (3prime_overlapping_ncRNA antisense bidirectional_promoter_lncRNA known_ncrna lincRNA macro_lncRNA non_coding nonsense_mediated_decay non_stop_decay processed_transcript retained_intron sense_intronic sense_overlapping)) {
	$GeneType{"$_"} = "Long non-coding RNA"; }
foreach (qw (miRNA misc_RNA Mt_rRNA Mt_tRNA ribozyme rRNA scaRNA scRNA snoRNA snRNA sRNA vaultRNA)) {
	$GeneType{"$_"} = "Others"; }
foreach (qw (pseudogene polymorphic_pseudogene processed_pseudogene transcribed_processed_pseudogene transcribed_unitary_pseudogene transcribed_unprocessed_pseudogene unitary_pseudogene unprocessed_pseudogene)) {
	$GeneType{"$_"} = "Pseudogene"; }
foreach (qw (IG_C_gene IG_D_gene IG_J_gene IG_V_gene IG_pseudogene IG_C_pseudogene IG_J_pseudogene IG_V_pseudogene TR_C_gene TR_D_gene TR_J_gene TR_V_gene TR_J_pseudogene TR_V_pseudogene processed_transcript TEC)) {
	$GeneType{"$_"} = "Others"; }

`awk '{print \$1"\t"int((\$2+\$3)/2)"\t"int((\$2+\$3)/2)+1"\t"\$1":"\$2"-"\$3}' $peak_bed > $outPrefix.peak_bed.center`;

#reads all genes of reference and save transcript's information(exon & CDS).
my (%RefAllGene, %RefPickTran, %RefPickTran_final);
open IN, $ref_gene_gtf || die;
open RBED,">$outPrefix.tmp.refSeq.bed" or die;
while (<IN>) {
	chomp;
	next if (/^\s*$|^\#/);
	my @w = split (/\t/);
	next if (@w < 9);
	my ($chr, $ftype, $start, $end, $strand, $features) = @w[0,2,3,4,6,8];
	next if ($ftype !~ /^exon$|^CDS$|^transcript$/);
	my ($gene_id, $transcript_id) = ("") x 2;
	$gene_id = $1 if ($features =~ /\bgene_id\s+\"([^\"]+)\";/);
	$transcript_id = $1 if ($features =~ /\btranscript_id\s+\"([^\"]+)\";/);
	if ($ftype eq "transcript") {
		my ($gene_name, $transcript_type) = ("") x 2;
		$gene_name = $1 if ($features =~ /\bgene_name\s+\"([^\"]+)\";/);
		$transcript_type = $1 if ($features =~ /\btranscript_type\s+\"([^\"]+)\";/);
		$RefAllGene{$gene_id}{$transcript_id}{chr} = $chr;
		$RefAllGene{$gene_id}{$transcript_id}{start} = $start;
		$RefAllGene{$gene_id}{$transcript_id}{end} = $end;
		$RefAllGene{$gene_id}{$transcript_id}{strand} = $strand;
		$RefAllGene{$gene_id}{$transcript_id}{gene_id} = $gene_id;
		$RefAllGene{$gene_id}{$transcript_id}{gene_name} = $gene_name;
		$RefAllGene{$gene_id}{$transcript_id}{transcript_type} = $transcript_type;
		$RefAllGene{$gene_id}{$transcript_id}{Gene_Type} = (exists $GeneType{$transcript_type}) ? $GeneType{$transcript_type} : "Unknown";
		print RBED "$RefAllGene{$gene_id}{$transcript_id}{chr}\t$RefAllGene{$gene_id}{$transcript_id}{start}\t$RefAllGene{$gene_id}{$transcript_id}{end}\t$transcript_id\t0\t$RefAllGene{$gene_id}{$transcript_id}{strand}\t$gene_id\n";
	}
	else { #$ftype =~ /^exon$|^CDS$/
		push @{$RefAllGene{$gene_id}{$transcript_id}{$ftype}}, [$start, $end, $end-$start+1];
	}
}
close RBED;
close IN;

`intersectBed -a $outPrefix.peak_bed.center -b $outPrefix.tmp.refSeq.bed -wa -wb > $outPrefix.refSeq.all.bed`;

#delete exon if exists CDS. delete shorter transcript if there are multi transcript in a gene.
open FH,"$outPrefix.refSeq.all.bed" or die; 
while (<FH>) {
	chomp;
	my @fields = split "\t";
	my ($peakid,$overlap_chr,$overlap_start,$overlap_end,$overlap_transid,$overlap_strand,$overlap_gene_id,$ppos) = @fields[3,4,5,6,7,9,10,1];
	#$RefPickTran{$peakid}{$overlap_gene_id}{$overlap_transid}{chr} = $RefAllGene{$overlap_gene_id}{$overlap_transid}{chr};
	#$RefPickTran{$peakid}{$overlap_gene_id}{$overlap_transid}{start} = $RefAllGene{$overlap_gene_id}{$overlap_transid}{start};
	#$RefPickTran{$peakid}{$overlap_gene_id}{$overlap_transid}{end} = $RefAllGene{$overlap_gene_id}{$overlap_transid}{end};
	#$RefPickTran{$peakid}{$overlap_gene_id}{$overlap_transid}{strand} = $RefAllGene{$overlap_gene_id}{$overlap_transid}{strand};
	#$RefPickTran{$peakid}{$overlap_gene_id}{$overlap_transid}{gene_id} = $RefAllGene{$overlap_gene_id}{$overlap_transid}{gene_id};
	#$RefPickTran{$peakid}{$overlap_gene_id}{$overlap_transid}{gene_name} = $RefAllGene{$overlap_gene_id}{$overlap_transid}{gene_name};
	#$RefPickTran{$peakid}{$overlap_gene_id}{$overlap_transid}{transcript_type} = $RefAllGene{$overlap_gene_id}{$overlap_transid}{tanscript_type};
	#$RefPickTran{$peakid}{$overlap_gene_id}{$overlap_transid}{Gene_Type} = $RefAllGene{$overlap_gene_id}{$overlap_transid}{Gene_Type};
	%{$RefPickTran{$peakid}{$overlap_gene_id}{$overlap_transid}} = %{$RefAllGene{$overlap_gene_id}{$overlap_transid}};
	$RefPickTran{$peakid}{$overlap_gene_id}{$overlap_transid}{chr} = $overlap_chr;
        $RefPickTran{$peakid}{$overlap_gene_id}{$overlap_transid}{center_start} = $ppos;
        $RefPickTran{$peakid}{$overlap_gene_id}{$overlap_transid}{center_end} = $ppos+1;
        $RefPickTran{$peakid}{$overlap_gene_id}{$overlap_transid}{chr_trans} = $overlap_chr;
        $RefPickTran{$peakid}{$overlap_gene_id}{$overlap_transid}{trans_start} = $overlap_start;
        $RefPickTran{$peakid}{$overlap_gene_id}{$overlap_transid}{trans_end} = $overlap_end;
	$RefPickTran{$peakid}{$overlap_gene_id}{$overlap_transid}{fuck_peakid} = $peakid;
        $RefPickTran{$peakid}{$overlap_gene_id}{$overlap_transid}{zero} = 0;
        $RefPickTran{$peakid}{$overlap_gene_id}{$overlap_transid}{strand} = $overlap_strand;
	$RefPickTran{$peakid}{$overlap_gene_id}{$overlap_transid}{ppos} = $ppos;
}
print "\n";
foreach my $peakid (keys %RefPickTran) {
	my ($cur_array, $cur_total_len, $cur_ftype);
        my ($longest_transcript_id, $longest_sum_len, $ftype, $longest_gene_id) = ("", 0, "","");
	foreach my $gene_id_temp (keys %{$RefPickTran{$peakid}}){
		foreach my $transcript_id_temp (keys %{$RefPickTran{$peakid}{$gene_id_temp}}){
			if (exists $RefPickTran{$peakid}{$gene_id_temp}{$transcript_id_temp}{CDS}) {
   	        		$cur_array = $RefPickTran{$peakid}{$gene_id_temp}{$transcript_id_temp}{CDS};
   	            		$cur_ftype = "CDS";
    			} else {
				$cur_array = $RefPickTran{$peakid}{$gene_id_temp}{$transcript_id_temp}{exon};
				$cur_ftype = "exon";
        		}
       	 		foreach (@{$cur_array}) {
                		$cur_total_len += $_->[2];
       	 		}
       	 		if ($cur_ftype eq $ftype and $cur_total_len > $longest_sum_len) {
                		$longest_transcript_id = $transcript_id_temp;
				$longest_gene_id = $gene_id_temp;
                		$longest_sum_len = $cur_total_len;
                		$ftype = $cur_ftype;
        		} else {
                		if($cur_ftype ne $ftype and $cur_ftype eq "CDS"){
                        		$longest_transcript_id = $transcript_id_temp;
					$longest_gene_id = $gene_id_temp;
                        		$longest_sum_len = $cur_total_len;
                        		$ftype = $cur_ftype;
                		} else {
                        		if($cur_ftype eq "exon" and $ftype eq ""){
                               			$longest_transcript_id = $transcript_id_temp;
						$longest_gene_id = $gene_id_temp;
                               			$longest_sum_len = $cur_total_len;
                               			$ftype = $cur_ftype;
                        		}
                		}
     	   		}
		}
	}
	$RefPickTran_final{$peakid}{$longest_transcript_id} = $RefPickTran{$peakid}{$longest_gene_id}{$longest_transcript_id};
   	$RefPickTran_final{$peakid}{$longest_transcript_id}{gene_type} = $ftype;
	@{$RefPickTran_final{$peakid}{$longest_transcript_id}{exon}} = sort {$a->[0] <=> $b->[0]} @{$RefPickTran_final{$peakid}{$longest_transcript_id}{exon}};
        if ($ftype eq "CDS") {
        	@{$RefPickTran_final{$peakid}{$longest_transcript_id}{CDS}} = sort {$a->[0] <=> $b->[0]} @{$RefPickTran_final{$peakid}{$longest_transcript_id}{CDS}};
                my ($total_tran_len, $total_cds_len, $utr5_len, $utr3_len) = (0, 0, 0, 0);
                my ($cds_start, $cds_end) = ($RefPickTran_final{$peakid}{$longest_transcript_id}{CDS}->[0][0], $RefPickTran_final{$peakid}{$longest_transcript_id}{CDS}->[-1][1]);
                foreach (@{$RefPickTran_final{$peakid}{$longest_transcript_id}{CDS}}) {
                       	$total_cds_len += $_->[2];
                }
                foreach (@{$RefPickTran_final{$peakid}{$longest_transcript_id}{exon}}) {
                       	$total_tran_len += $_->[2];
                       	my ($cur_start, $cur_end, $cur_len) = @{$_};
                       	if ($cds_start > $cur_start) {
                               	if ($cds_start <= $cur_end) { $utr5_len += $cds_start - $cur_start; }
                               	else { $utr5_len += $cur_len; }
                       	}
                       	if ($cds_end < $cur_end) {
                               	if ($cds_end >= $cur_start) { $utr3_len += $cur_end - $cds_end; }
                               	else { $utr3_len += $cur_len; }
                       	}
                }
                $RefPickTran_final{$peakid}{$longest_transcript_id}{total_tran_len} = $total_tran_len;
                $RefPickTran_final{$peakid}{$longest_transcript_id}{total_cds_len} = $total_cds_len;
                $RefPickTran_final{$peakid}{$longest_transcript_id}{utr5_len} = $utr5_len;
		$RefPickTran_final{$peakid}{$longest_transcript_id}{utr3_len} = $utr3_len;
        }
	else { #$ftype eq "exon", that is not coding RNA
   	     my $total_tran_len;
   	     foreach (@{$RefPickTran_final{$peakid}{$longest_transcript_id}{exon}}) { $total_tran_len += $_->[2]; }
   	     $RefPickTran_final{$peakid}{$longest_transcript_id}{total_tran_len} = $total_tran_len;
   	     $RefPickTran_final{$peakid}{$longest_transcript_id}{total_cds_len} = 0;
   	     $RefPickTran_final{$peakid}{$longest_transcript_id}{utr5_len} = 0;
   	     $RefPickTran_final{$peakid}{$longest_transcript_id}{utr3_len} = 0;
    	}
}

sub sort_2dArray {
	my $in_arr = $_[0];
	my @out_arr = sort {$a->[0] <=> $b->[0]} @{$in_arr};
	return \@out_arr;
}


#annotate
my %p2t;
foreach my $peak_id (keys %RefPickTran_final){
	my $tran_id = (keys  %{$RefPickTran_final{$peak_id}})[0];
	my $ppos = $RefPickTran_final{$peak_id}{$tran_id}{ppos};
	if(exists $p2t{$peak_id}){
		if($p2t{$peak_id}{cdslen}==0){
			if($RefPickTran_final{$peak_id}{$tran_id}{total_cds_len}<=$p2t{$peak_id}{cdslen}){
				next;
			}else{
				if($RefPickTran_final{$peak_id}{$tran_id}{total_tran_len}<=$p2t{$peak_id}{tlen}){
					next;
				}
			}
		}else{
			if($RefPickTran_final{$peak_id}{$tran_id}{total_cds_len}<=$p2t{$peak_id}{cdslen}){
				next;
			}
		}
	}
	$p2t{$peak_id}{cdslen}= $RefPickTran_final{$peak_id}{$tran_id}{total_cds_len}; #$cdslen{$tran_id};
	$p2t{$peak_id}{tlen}=$RefPickTran_final{$peak_id}{$tran_id}{total_tran_len};
	$p2t{$peak_id}{gene}=$RefPickTran_final{$peak_id}{$tran_id}{gene_name};
	$p2t{$peak_id}{ts}=$tran_id;
	$p2t{$peak_id}{ppos}=$ppos;
	$p2t{$peak_id}{intersect}= $RefPickTran_final{$peak_id}{$tran_id}{chr}."\t".$RefPickTran_final{$peak_id}{$tran_id}{center_start}."\t".$RefPickTran_final{$peak_id}{$tran_id}{center_end}."\t".$peak_id."\t".$RefPickTran_final{$peak_id}{$tran_id}{chr_trans}."\t".$RefPickTran_final{$peak_id}{$tran_id}{trans_start}."\t".$RefPickTran_final{$peak_id}{$tran_id}{trans_end}."\t".$tran_id."\t".$RefPickTran_final{$peak_id}{$tran_id}{zero}."\t".$RefPickTran_final{$peak_id}{$tran_id}{strand};
}

open OUT, ">$outPrefix.anno.txt" || die;
foreach my $peak_id (keys %p2t){
	my $ppos=$p2t{$peak_id}{ppos};
	my $tran_id=$p2t{$peak_id}{ts};
	my %Tran = %{$RefPickTran_final{$peak_id}{$tran_id}};
	my ($bin, $exon_sum_len, $segtype) = (0, 0, "");
	my $cstatus = ($Tran{gene_type} eq "CDS") ? "coding" : "noncoding";
	my @cur_array= @{$Tran{exon}};
	my ($cds_start, $cds_end) = ($Tran{CDS}->[0][0], $Tran{CDS}->[-1][1]);
	for (my $i=0; $i<=$#cur_array; $i++) {
		my ($exon_start, $exon_end, $exon_len) = @{$cur_array[$i]};
		$exon_sum_len += $exon_len;
		if ($ppos >= $exon_start && $ppos <= $exon_end) {
			if ($cstatus eq "noncoding") {
#				$bin = int (($ppos - $exon_start) / $exon_len * 100); # $bin is the percentage of ppos in each exon.
				$bin = int (($exon_sum_len - ($exon_end - $ppos)) / $Tran{total_tran_len} * 100); # $bin is the percentage of ppos in each transcript.
				$segtype = "exon";
			}
			else { #$cstatus eq "coding"
				if ($ppos < $cds_start) {
					if ($Tran{utr5_len} == 0) {print join ("\t", $peak_id, $cds_start, $cds_end, $ppos, "\n");}
					$bin = int (($exon_sum_len - ($exon_end - $ppos)) / $Tran{utr5_len} * 100);
					$segtype = ($Tran{strand} eq "+") ? "5UTR" : "3UTR";
				}elsif ($ppos > $cds_end) {
					$bin = int (($Tran{total_tran_len} - $exon_sum_len + ($exon_end - $ppos)) / $Tran{utr3_len} * 100);
					$bin = 100 - $bin;
					$segtype = ($Tran{strand} eq "+") ? "3UTR" : "5UTR";
				}else {
					$bin = int (($exon_sum_len - ($exon_end - $ppos) - $Tran{utr5_len}) / $Tran{total_cds_len} * 100);
					$segtype="CDS";
				}
			}
			last;
		}
		else {
			if ($i < $#cur_array) {# isn't the last one
				my $next_exon_start = $cur_array[$i+1]->[0];
				if ($ppos > $exon_end && $ppos < $next_exon_start) {
					$bin = int (($ppos - $exon_end) / ($next_exon_start - $exon_end) * 100);
					$segtype = "intron";
					last;
				}
			}
		}
	} #end: for (my $i=0; $i<=$#cur_array; $i++)
	if ($segtype) { #can find peak in the transcript
		$bin = 100 - $bin if ($Tran{strand} eq "-");
		print OUT $p2t{$peak_id}{intersect}, "\t", join ("\t", $Tran{gene_name}, $cstatus, $segtype, $bin, $Tran{gene_id}, $Tran{transcript_type}, $Tran{Gene_Type}), "\n";
	}
}
close OUT;

print `perl $Bin/intersec.pl -a $outPrefix.peak_bed.center -na 4 -b $outPrefix.anno.txt -nb 4 -t ua > $outPrefix.unanno.txt`;
print `/usr/bin/Rscript $Bin/m6A_annotate.v2.R $outPrefix.anno.txt $outPrefix.unanno.txt $outPrefix`;



__END__

	for(my $i=0;$i<@starts;$i++){
		$exon_len+=$ends[$i]-$starts[$i];
		if($ppos>=$starts[$i]&&$ppos<=$ends[$i]){
			if($cstatus eq "noncoding"){
				$bin=int(($ppos-$starts[$i])/($ends[$i]-$starts[$i])*100);
				$segtype="exon";
				last;
			}			
			if($ppos<$cds_start){
				$bin=int(($exon_len-($ends[$i]-$ppos))/$utr5len{$tname}*100);
			}elsif($ppos>$cds_end){
				$bin=int(($tlen{$tname}-$exon_len+($ends[$i]-$ppos))/$utr3len{$tname}*100);
				$bin=100-$bin;
			}else{
				$bin=int(($exon_len-($ends[$i]-$ppos)-$utr5len{$tname})/$cdslen{$tname}*100);		
				$segtype="CDS";
			}
		}
		if($ppos<$starts[$i+1]&&$ppos>$ends[$i]){
			$bin=int(($ppos-$ends[$i])/($starts[$i+1]-$ends[$i])*100);
			$segtype="intron";
			last;
		}
	}
	$bin=100-$bin if $strand eq "-";
	print $p2t{$peak_id}{intersect}."\t$gname\t$cstatus\t$segtype\t$bin\n";
		



exit;





my $geneanno=shift;
my $bed=shift;


open G,$geneanno or die;

my %t2g;
my %tlen;
my %cdslen;
my %ref;
my %utr5len;
my %utr3len;

#canonical transcripts: longest cds, otherwise longest transcripts
open RBED,">tmp.refSeq.bed" or die;

while(<G>){
	chomp;
	next if $_=~/^#/;
	my @field=split "\t";
	my $strand=$field[3];
	my $chr=$field[2];
	my $tname=$field[1];
	my $gname=$field[12];
	my $start=$field[4];
	my $end=$field[5];
	my $cds_start=$field[6];
	my $cds_end=$field[7];
	my $starts=$field[9];
    my @starts=split ",",$starts;
    my $ends=$field[10];
    my @ends=split ",",$ends;
	my $tlen=0;
	my $cds_len=0;
	my $utr5_len=0;
	my $utr3_len=0;

	print RBED $chr."\t".$start."\t".$end."\t".$tname."\t0\t$strand\n";
    for(my $i=0;$i<@starts;$i++){
       	$tlen+=$ends[$i]-$starts[$i];
		if($cds_start<=$ends[$i]&&$cds_start>=$starts[$i]){
			$utr5_len=$tlen-($ends[$i]-$cds_start);
		}
		if($cds_end<=$ends[$i]&&$cds_end>=$starts[$i]){
			$utr3_len=$tlen-($ends[$i]-$cds_end);
		}
	}
	$utr3_len=$tlen-$utr3_len;
	$cds_len=$tlen-$utr5_len-$utr3_len;
	
	$ref{$tname}=$_;
	$t2g{$tname}=$gname;
	$tlen{$tname}=$tlen;
	$cdslen{$tname}=$cds_len;
	$utr5len{$tname}=$utr5_len;
	$utr3len{$tname}=$utr3_len;
}
close RBED;

`awk '{print \$1"\t"int((\$2+\$3)/2)"\t"int((\$2+\$3)/2)+1"\t"\$1":"\$2"-"\$3}' $bed > $bed.center`;
`intersectBed -a $bed.center -b tmp.refSeq.bed  -wa -wb > $bed.refSeq.all.bed`;

open FH,"$bed.refSeq.all.bed" or die;


my %p2t;
while(<FH>){
	chomp;
	my @fields=split "\t";
	if(exists $p2t{$fields[3]}){
		if($p2t{$fields[3]}{cdslen}==0){
			if($cdslen{$fields[7]}<=$p2t{$fields[3]}{cdslen}){
				next;
			}else{
				if($tlen{$fields[7]}<=$p2t{$fields[3]}{tlen}){
					next;
				}
			}
		}else{
			if($cdslen{$fields[7]}<=$p2t{$fields[3]}{cdslen}){
				next;
			}
		}
	}
	$p2t{$fields[3]}{cdslen}=$cdslen{$fields[7]};
	$p2t{$fields[3]}{tlen}=$tlen{$fields[7]};
	$p2t{$fields[3]}{gene}=$t2g{$fields[7]};
	$p2t{$fields[3]}{ts}=$fields[7];
	$p2t{$fields[3]}{ppos}=$fields[1];
	$p2t{$fields[3]}{intersect}=$_;
}
close FH;

foreach my $k (keys %p2t){
	my $ppos=$p2t{$k}{ppos};
	my $tname=$p2t{$k}{ts};
	my @field=split "\t",$ref{$tname};
	my $strand=$field[3];
	my $chr=$field[2];
	my $gname=$field[12];
	my $start=$field[4];
	my $end=$field[5];
	my $cds_start=$field[6];
	my $cds_end=$field[7];
	my $starts=$field[9];
    my @starts=split ",",$starts;
    my $ends=$field[10];
    my @ends=split ",",$ends;
	my $exon_len=0;
	my $cstatus="coding";
	$cstatus="noncoding" if $cds_end-$cds_start==0;
	my $bin;
	my $segtype="";
	for(my $i=0;$i<@starts;$i++){
		$exon_len+=$ends[$i]-$starts[$i];
		if($ppos>=$starts[$i]&&$ppos<=$ends[$i]){
			if($cstatus eq "noncoding"){
				$bin=int(($ppos-$starts[$i])/($ends[$i]-$starts[$i])*100);
				$segtype="exon";
				last;
			}			
			if($ppos<$cds_start){
				$bin=int(($exon_len-($ends[$i]-$ppos))/$utr5len{$tname}*100);
				
				$segtype="5UTR" if $strand eq "+";
				$segtype="3UTR" if $strand eq "-";
				
			}elsif($ppos>$cds_end){
				$bin=int(($tlen{$tname}-$exon_len+($ends[$i]-$ppos))/$utr3len{$tname}*100);
				$bin=100-$bin;
				$segtype="3UTR" if $strand eq "+";
				$segtype="5UTR" if $strand eq "-";
				
			}else{
				$bin=int(($exon_len-($ends[$i]-$ppos)-$utr5len{$tname})/$cdslen{$tname}*100);		
				$segtype="CDS";
			}
			last;
		}
		if ($ppos>$ends[-1]) { print "tail, $ppos>$ends[-1]\n"; }
		if ($ppos<$starts[0]) { print "head, $ppos<$starts[0]\n"; }
		if($ppos<$starts[$i+1]&&$ppos>$ends[$i]){
			$bin=int(($ppos-$ends[$i])/($starts[$i+1]-$ends[$i])*100);
			
			$segtype="intron";
			
			last;
		}
		
	}
	$bin=100-$bin if $strand eq "-";
	print $p2t{$k}{intersect}."\t$gname\t$cstatus\t$segtype\t$bin\n";
		
}



