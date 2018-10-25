#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $BEGIN_TIME=time();
my $version="1.0.0";
my($fq1,$fq2,$fkey);
my $outdir = "./";
my $distance = 2000;
my $SEAsmall = 0;
GetOptions(
	"help|?" =>\&USAGE,
	"fq1:s"=>\$fq1,
	"fq2:s"=>\$fq2,
	"k:s"=>\$fkey,
	"od:s"=>\$outdir,
	"distance:i"=>\$distance,
	"SEAsmall:i"=>\$SEAsmall,
) or &USAGE;
&USAGE unless ($fq1 and $fq2 and $fkey);
`mkdir $outdir` unless(-d $outdir);

my (%uniqstat,%readsstat);
if($fq1 =~ /\.gz$/){
	open F1,"gzip -dc $fq1|" || die "$!\n";
}else{
	open F1,"$fq1" || die "$!\n";
}
if($fq2 =~ /\.gz$/){
	open F2,"gzip -dc $fq2|" || die "$!\n";
}else{
	open F2,"$fq2" || die "$!\n";
}
open O1, ">$outdir/$fkey.consensus.final.R1.fastq" || die "$!\n";
open O2, ">$outdir/$fkey.consensus.final.R2.fastq" || die "$!\n";
open OUT, ">$outdir/$fkey.diffchr.largefragment.txt" || die "$!\n";
open OB1,"|gzip >$outdir/$fkey.consensus.diffchr.largefragment.R1.fastq.gz" || die "$!\n";
open OB2,"|gzip >$outdir/$fkey.consensus.diffchr.largefragment.R2.fastq.gz" || die "$!\n";
while(<F1>){
	chomp;
	my $name1 = $_;
	my $seq1 = <F1>;
	my $or1 = <F1>;
	my $qual1 = <F1>;
	my $name2 = <F2>;chomp $name2;
	my $seq2 = <F2>;
	my $or2 = <F2>;
	my $qual2 = <F2>;
	chomp $name1;
	my @infor = split(/\s+/,$name1);
	my ($chr1, $pos1, $strand1, $chr2, $pos2, $strand2)=$name1=~/YP:Z:(chr\w+)_(\d+)_([+-]),(chr\w+)_(\d+)_([+-])/;
	my $sea = $name1=~/SEA/? "YES":"NO";
	my ($dep) = $name1=~/XF:i:(\d+)/;
	if(!defined $chr1){
		print join("\t",  $name1,($chr1, $pos1, $strand1, $chr2, $pos2, $strand2)),"\n";
		die;
	}
	$uniqstat{'total'}++;
	$readsstat{'total'}+=$dep;
	my $filter = "NO";
	my $size;
	if($chr1 ne $chr2){
		$readsstat{'diffchr'} += $dep;
		$uniqstat{'diffchr'}++;
		$filter = "DiffChr";
	}else{
		if($sea eq "NO"){
			$size = abs($pos1 - $pos2) + 1;
			if($size > $distance){
				$readsstat{'distance'} += $dep;
				$uniqstat{'distance'}++;
				$filter = "Size$size";
			}
		}else{
			if($pos2 - 215395 >10000){
				$size = abs($pos2 - 234700) +1;
				if(($size > $distance) || ( $size < $SEAsmall)){
					$readsstat{'distance'} += $dep;
					$uniqstat{'distance'}++;
					$filter = "Size$size";
				}
			}else{
				$size = abs($pos2 - 215395)+1;
				if(($size > $distance) || ($size < $SEAsmall)){
					$readsstat{'distance'} += $dep;
					$uniqstat{'distance'}++;
					$filter = "Size$size";
				}
			}
		}
	}
	if($filter eq "NO"){
		print O1 "$name1\tSZ:i:$size\n$seq1$or1$qual1";
		print O2 "$name2\tSZ:i:$size\n$seq2$or2$qual2";
	}else{
		print OB1 "$name1\tXR:Z:$filter\n$seq1$or1$qual1";
		print OB2 "$name2\tXR:Z:$filter\n$seq2$or2$qual2";

	}
}
close F1;
close F2;
close O1;
close O2;
close OB1;
close OB2;
my $ratio_diff_u = sprintf("%.2f",$uniqstat{'diffchr'}/$uniqstat{'total'}*100);
my $ratio_diff_r = sprintf("%.2f",$readsstat{'diffchr'}/$readsstat{'total'}*100);
my $ratio_dis_u = sprintf("%.2f",$uniqstat{'distance'}/$uniqstat{'total'}*100);
my $ratio_dis_r = sprintf("%.2f",$readsstat{'distance'}/$readsstat{'total'}*100);
print OUT "Total_unique\tTotal_reads\tDiffChr_unique\tDiffChr_unique_ratio(%)\tDiffChr_reads\tDiffChr_reads_ratio(%)\tLarge_distance_unique\tLarge_distance_unique_ratio(%)\tLarge_distance_reads\tLarge_distance_reads_ratio(%)\n";
print OUT "$uniqstat{'total'}\t$readsstat{'total'}\t$uniqstat{'diffchr'}\t$ratio_diff_u\t$readsstat{'diffchr'}\t$ratio_diff_r\t$uniqstat{'distance'}\t$ratio_dis_u\t$readsstat{'distance'}\t$ratio_dis_r\n";

sub USAGE {
	my $usage=<<"USAGE";
Program:
Version: $version
Contact:Du Yuanping<yuanping.du\@genetalks.com>

Usage:
Options:
-fq1   <file>  consensus.remain.R1.fastq file, forced
-fq2   <file>  consensus.remain.R2.fastq file, forced
-k     <str>   key of output file, forced
-od    <dir>   dir of output file, [$outdir]
-distance  <int>  maximum cutoff value of read1 and read2 alignment position in same chr, [$distance]
-SEAsmall  <int>  minimum cutoff value of read1 fnd read2 alignment position of SEA, [$SEAsmall]
-h Help
USAGE
	print $usage;
	exit;
}
