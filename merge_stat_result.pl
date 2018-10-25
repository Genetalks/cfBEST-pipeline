#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
require "$Bin/mylib.pl";
my $BEGIN_TIME=time();
my $version="1.0.0";
#######################################################################################

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($indir, $Uniq_QC, $Bar_QC, $Primer_QC, $Map_QC, $Variant_Stat, $Hotspot_Stat, $fkey,$outdir);
my ($filter, $Variant_Known, $Consens_Stat);
GetOptions(
				"help|?" =>\&USAGE,
				"id:s"=>\$indir,
				"filter:s"=>\$filter,
				"Uniq_QC:s"=>\$Uniq_QC,
				"Bar_QC:s"=>\$Bar_QC,
				"Primer_QC:s"=>\$Primer_QC,
				"Map_QC:s"=>\$Map_QC,
				"Consens_Stat:s"=>\$Consens_Stat,
				"HotSpot_Stat:s"=>\$Hotspot_Stat,
				"Variant_Known:s"=>\$Variant_Known,

				"od:s"=>\$outdir,
				) or &USAGE;
&USAGE unless ($indir);

$outdir||="./";
`mkdir $outdir`	unless (-d $outdir);
$outdir=AbsolutePath("dir",$outdir);
$filter = defined $filter? $filter: 1;

my $out;
#open ($out,">$outdir/$fkey.stat") or die $!;

if (defined $Uniq_QC) {
	my $out;
	open ($out,">$outdir/01.Uniq.stat") or die $!;
	&merge_horizon($out, $indir, ".fquniq.stat.txt", 0, 0); ## is_Total, dir_rank
	close ($out);
}


if (defined $Bar_QC) {
	my $max_num = 20;
	my ($bout,$eout,$dout,$sout);
	open ($eout,">$outdir/02.Error_Bar.stat") or die $!;
	open ($sout,">$outdir/02.Specific_Bar.stat") or die $!;

	## error barcode stat
	my @ttstat=glob("$indir/*.barcode.stat");
	if(@ttstat){  
		&merge_horizon($eout, $indir, ".barcode.stat", 0, 0); ## is_Total, dir_rank
	}else{  
		&merge_horizon($eout, $indir, ".barcode_filter.stat.txt", 0, 0); ## is_Total, dir_rank
	}

	## specific bar
	my @ttinfo=glob("$indir/*.barcode.info");
	if(@ttinfo){  
		&merge_vertical($sout, $indir, ".barcode.info", 1, 0, 0, 1);# is_get_aver, $is_head, $is_subhead, $row_head_column
	}else{  
		&merge_vertical($sout, $indir, "barcode.stat.txt", 1, 0, 0, 1);# is_get_aver, $is_head, $is_subhead, $row_head_column
	}

	close ($eout) ;
	close ($sout) ;
}

if (defined $Primer_QC) {

	my ($dout, $pout, $iout);
#	open ($dout,">$outdir/02.Data.stat") or die $!;
	open ($pout,">$outdir/03.Primer.stat") or die $!;
	open ($iout, ">$outdir/03.Primer.info") or die $!;

	## stat primer 
	&merge_horizon($pout, $indir, ".primer.stat", 0, 0); ## is_Total, dir_rank

	## stat primer info
	&merge_vertical($iout, $indir, ".primer.info", 0, 0, 0, 1);# is_get_aver, $is_head, $is_subhead, $row_head_column

	close ($iout);
#	close ($dout);
	close ($pout);

}


if (defined $Map_QC) {

	## open 
	my ($mout, $iout);
	open ($mout,">$outdir/04.Map.stat") or die $!;
	open ($iout,">$outdir/04.Map.info") or die $!;

	## stat map QC 
	&merge_horizon($mout, $indir, ".mapping.stat.txt", 0, 0); ## is_Total, dir_rank

	## stat primer map QC
	&merge_vertical($iout, $indir, ".mapping.primer.stat.txt", 0, 0, 1, 1);# is_get_aver, $is_head, $is_subhead, $row_head_column

	close($mout);
	close($iout);
}

if (defined $Consens_Stat) {
	my ($dout, $uout, $pout, $iout);
	open ($dout,">$outdir/04.Dedup.stat") or die $!;
	open ($uout,">$outdir/04.Uniq.stat") or die $!;
	open ($iout,">$outdir/04.Uniq.info") or die $!;

	## dedup stat 
	&merge_horizon($dout, $indir, ".dedup.metrics.txt", 0, 0); ## is_Total, dir_rank

	## uniq stat 
	&merge_horizon($uout, $indir, ".consensus.stat", 0, 0); ## is_Total, dir_rank

	## stat uniq info
	&merge_vertical($iout, $indir, ".consensus.info", 0, 0, 1, 1);# is_get_aver, $is_head, $is_subhead, $row_head_column

	close ($dout);
	close ($uout);
	close ($iout);

}


if (defined $Hotspot_Stat) {
	my %stat;

	## stat hotspot.depth
	my @file = glob("$indir/QC/*hotspot.depth"); 
	foreach my $f (@file) {
		my ($sample)=$f=~/QC\/(\S+).hotspot.depth/;
		&read_in_table_stat_vertical($f, \%stat, $sample, 1, 1, 0, 1);	# is_get_aver, $is_head, $is_subhead, $row_head_column
	}

	## output 
	my ($hout, $fout);
	open ($hout,">$outdir/05.Hotspot.stat") or die $!;
	open ($fout,">$outdir/05.Fusion.stat") or die $!;
	&output_vertical($hout, \%stat, 1, 1, 0, 1);
	
	## stat fusion region
	&stat_consens_fusion($fout, $indir);

	close($hout);
	close($fout);

}


if (defined $Variant_Known) {
	my @result = glob("$indir/variants/*variants_known.result");
	my %stat;
	foreach my $f (@result) {
		my $fname = basename($f);
		my ($sample)=$fname=~/(\S+).variants_known.result/; 
		&read_in_table_stat_vertical($f, \%stat, $sample, 0, 1, 0, 11); # is_get_aver, $is_head, $is_subhead, $row_head_column
	}
	
	## merge variant
	{
		my ( $ffout);
		open ($ffout,">$outdir/Total.fusion.final.txt") or die $!;
		&merge_horizon($ffout, $indir, ".fusion.final.txt", 0, 1);

		close ($ffout) ;
	}
}

#close ($out);


#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################

# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------
sub merge_vertical{#
	my ($out, $indir, $fsuffix, $is_get_aver, $is_head, $is_subhead, $row_head_column)=@_;
	my %stat;
	my @file = glob("$indir/*$fsuffix"); 
	foreach my $f (@file) {
		my $fname = basename ($f);
		my ($sample)=$fname=~/(\S+)$fsuffix/;
		&read_in_table_stat_vertical($f, \%stat, $sample, $is_get_aver, $is_head, $is_subhead, $row_head_column);
	}

	&output_vertical($out, \%stat, $is_get_aver, $is_head, $is_subhead, $row_head_column);
}

## is_Total: is there a "Total" row, and only merge "Total" row
## dir_rank: which directory rank is the fsuffix files in, count from 0
sub merge_horizon {#
	my ($out, $indir, $fsuffix, $is_Total, $dir_rank) = @_;
	my $head;

	my %stat;
	my @file;
	if ($dir_rank == 0) {
		@file = glob("$indir/*$fsuffix");
	}elsif($dir_rank == 1){
		@file = glob("$indir/*/*$fsuffix");
	}elsif($dir_rank == 2){
		@file = glob("$indir/*/*/*$fsuffix");
	}
	foreach my $f (@file) {
		my $fname = basename ($f);
		my ($sample)=$fname=~/(\S+)$fsuffix/;

		if ($is_Total) {
			$head = &read_in_row_stat_horizon_total($f, \%stat, $sample);
		}else{
			$head = &read_in_row_stat_horizon($f, \%stat, $sample);
		}
	}

	## write 
	&output_row_stat_horizon($out, \%stat, $head);
}


sub output_row_stat_horizon {#
	my ($out, $astat, $head) = @_;
	print $out "SampleID\t",$head,"\n";
	foreach my $s (sort {$a cmp $b} keys %{$astat}) {
		for (my $i=0;$i<@{$astat->{$s}} ; $i++) {
			print $out $s,"\t",$astat->{$s}->[$i],"\n";
		}
	}
}


sub read_in_row_stat_horizon {#
	my ($f, $astat, $sample) = @_;
	open (I,"$f") or die $!;
	my $head = <I>;
	chomp $head;

	while (<I>) {
		chomp;
		next if(/^$/);
		push @{$astat->{$sample}}, $_;
	}
	close (I) ;
	return($head);
}


sub read_in_row_stat_horizon_total {#
	my ($f, $astat, $sample) = @_;
	open (I,"$f") or die $!;
	my $head = <I>;
	chomp $head;
	my @hunit = split /\t/, $head;
	shift @hunit;
	$head = join("\t",@hunit);

	while (<I>) {
		chomp;
		next if(/^$/ || $_!~/Total/);
		my @unit = split /\t/, $_;
		shift @unit;
		push @{$astat->{$sample}}, join("\t",@unit);
	}
	close (I) ;

	return($head);
}


sub stat_consens_fusion {#
	my ($out, $indir) = @_;

	my %fusion_stat;
	my @fusion = glob("$indir/QC/*fusion.stat");
	my $head_flag = 0;

	foreach my $f (@fusion) {
		my ($sample)=$f=~/QC\/(\S+).fusion.stat/;
		open (I,$f) or die $!;
		if ($head_flag == 0) {
			my $head = <I>;
			print $out $head;
			$head_flag = 1;
		}else{
			<I>;
		}
		my $info = <I>;
		print $out $info;
		close (I) ;
	}
}


sub output_vertical {#
	my ($out, $astat, $is_get_aver, $is_head, $is_subhead, $row_head_column)=@_;

	my @sample = sort {$a cmp $b} keys %{$astat->{"Average"}};
	## head and average

	my $head_info;
	if ($is_head) {
		$head_info = $astat->{"Head"};
	}else{
		$head_info = "#Info"."\t-" x ($row_head_column-1);
	}
	print $out $head_info;
	foreach my $s (@sample) {
		print $out "\t",$s."\t"x($astat->{"Value_c"}-1);
	}
	print $out "\n";

	if ($is_subhead) {
		my ($type)=keys %{$astat->{"SubHead"}};
		print $out $type;
		foreach my $sample (@sample) {
			print $out "\t",$astat->{"SubHead"}{$type}{$sample};
		}
		print $out "\n";
	}

	if ($is_get_aver) {
		print $out "Average"."\t-"x ($row_head_column-1);
		foreach my $sample (@sample) {
			print $out "\t",sprintf("%.2f",$astat->{"Average"}{$sample});
		}
		print $out "\n";
	}
	if (exists $astat->{"Total"}) {
		my $type = "Total";
		print $out $type;
		foreach my $sample (@sample) {
			my $n = exists $astat->{$type}{$sample}? $astat->{$type}{$sample}: "NA"."\tNA"x($astat->{"Value_c"}-1);
			print $out "\t",$n;
		}
		print $out "\n";
	}


	foreach my $type (sort {$a cmp $b} keys %{$astat}) {
		next if($type eq "Average" || $type eq "Value_c" || $type eq "Head" || $type eq "SubHead" || $type eq "Total");
		print $out $type;
		foreach my $sample (@sample) {
			my $n = exists $astat->{$type}{$sample}? $astat->{$type}{$sample}: "NA"."\tNA"x($astat->{"Value_c"}-1);
			print $out "\t",$n;
		}
		print $out "\n";
	}
}

#is_get_aver: get average only when the value column is 1
#is_head: is there a head line contain sample name
#is_subhead: is there a head line which not contain sample name 
#row_head_column: row head columns

### eg
#Info	BRCA31						BRCA405					**head
#Amplicon	TotalUniq	PosUniq	BarUniq	TotalReads	Duplicate	DupRatio	TotalUniq	PosUniq	BarUniq	TotalReads	Duplicate	DupRatio    **subhead
#Total	25666	7668	17998	5743803	5718137	99.55%	11444	2752	8692	4823441	4811997  **total
#AKT1-17-D-TN1	782	205	577	386992	386210	99.80%	704	148	556	533093	532389	99.87%
#AKT1-17-U-TN1	914	174	740	483727	482813	99.81%	NA	NA	NA	NA	NA	NA
#AKT1-17-U-TN2	613	167	446	183150	182537	99.67%	415	88	327	285642	285227	99.85%
sub read_in_table_stat_vertical{
	my ($f, $astat, $sample, $is_get_aver, $is_head, $is_subhead, $row_head_column)=@_;
	open (I, $f) or die $!;
	my $n =0;
	my $s = 0;
	my $value_c;
	if ($is_head) {
		my $head = <I>;
		chomp $head;
		my @unit = split /\t/, $head;
		my $head_type = join("\t",@unit[0..($row_head_column-1)]);
		$astat->{"Head"}=$head_type;
	}
	if ($is_subhead) {
		my $head = <I>;
		chomp $head;
		my @unit = split /\t/, $head;
		my $type = join("\t",@unit[0..($row_head_column-1)]);
		my $value = join("\t",@unit[$row_head_column..$#unit]);
		$astat->{"SubHead"}{$type}{$sample}=$value;
	}
	while (<I>) {
		chomp;
		next if(/^$/ || /#/);
		my @unit = split /\t/, $_;
		my $type = join("\t",@unit[0..($row_head_column-1)]);
		my $value = join("\t",@unit[$row_head_column..$#unit]);
		$value_c = $#unit+1-$row_head_column;

		$n++;
		if ($is_get_aver == 1) {			
			$s+=$value;
		}
		$astat->{$type}{$sample}=$value;
	}
	close (I);
	
	if ($is_get_aver == 1) {
		my $aver = $n>0? sprintf("%.2f",$s/$n): sprintf("%.2f",0);
		$astat->{"Average"}{$sample}=$aver;
	}else{
		$astat->{"Average"}{$sample}="NA";
	}
	$astat->{"Value_c"} = $value_c;
}





sub GetTime {
	my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
	return sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}


sub USAGE {#
	my $usage=<<"USAGE";
Program:
Version: $version
Contact:zeng huaping<huaping.zeng\@genetalks.com> 

Usage:
  Options:
  -id  <dir>   dir of Input file, forced

  --Uniq_QC        merge Uniq QC statistic
  --Bar_QC         merge Barcode QC statisic
  --Primer_QC      merge Primer QC statistic
  --Map_QC         merge Map QC statistic
  --Consens_Stat   merge group info statistic
  --HotSpot_Stat    merge hotspot depth statistic
  --Variant_Known	merge variant known info

  -filter   <int>    1: do flexbar, remove adapter and filter low quality; 0: not do flexbar, [1]
  -od <dir>    outdir of output file, default ./
  -h         Help

USAGE
	print $usage;
	exit;
}
