#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Cwd qw(abs_path);
require "$Bin/../mylib.pl";
my $BEGIN_TIME=time();
my $version="1.0.0";
my @Original_ARGV=@ARGV;
#######################################################################################
# ==============================================================
# Get Options
# ==============================================================
my $fsample;
my $indir;
my $dOut = "./"; 
my $min_depth = 100;
my $infer_genptype_by_cffDNA;
my $fadjust_aver_freq = "$Bin/adjust_average_freq.txt";
my $fconfig = "$Bin/../config/config.txt";
my ($fkey, $configID);
GetOptions(
	"help|?" =>\&USAGE,
	"i:s"=>\$fsample,
	"ic:s"=>\$fconfig,
	"cid:s"=>\$configID,
	"id:s"=>\$indir,
	"od:s"=>\$dOut,
	"k:s"=>\$fkey,
	"d:i"=>\$min_depth,
	"e"=>\$infer_genptype_by_cffDNA,
) or &USAGE;
&USAGE unless ($fsample and $indir and $fkey and $configID);

mkdir $dOut if (defined $dOut && not -d $dOut);
$dOut = Cwd::abs_path($dOut);

#===============================================================
# Main
#===============================================================

my %config;
&read_data_sheet($fsample,$fconfig, \%config, $configID); ## @{$adata->{$sample}->{"sample_info"}}=@sample_info;



# write shell for fetal genotyping and fetal fraction estimation using EM algorithm
my $sh = "$dOut/cffdna_genotype.sh";
open (SH,">$sh") or die $!;
##SampleID	Index	BarcodeGroup	SampleType	Description	GenotypeFetal	GenotypeMaternal	GenotypePaternal
#demo	GATGCCA	B3,B4,B5	gravida-plasma-thalassemia	sample descriptions	CD41-42/N	-28/N	CD41-42/N
#demo2	CTTCGTT	B6,B7,B8	gravida-plasma-thalassemia	sample descriptions	-28/N	-28/N	CD41-42/N
foreach my $s (sort {$a cmp $b} keys %config) {
	my @saminfo = @{$config{$s}{"sample_info"}};
	my $fin = "$indir/$s.hotspots.result";
	my $fout = "$dOut/$s.genotype.txt";
	my $cmd;
	next if($config{$s}{"sample_info"}->[0] ne "gravida-plasma-thalassemia" && $config{$s}{"sample_info"}->[0] ne "standards-thalassemia");

	## adjust by average freq
	&adjust_by_aver($fin, "$dOut/$s.hotspots.result", $fadjust_aver_freq);
	$fin = "$dOut/$s.hotspots.result";

	## genotyping
	if($saminfo[0] eq "gravida-plasma-thalassemia"){
		$cmd = "Rscript $Bin/refine_fetal_genotyping_EM.r -i $fin -o $fout -d $min_depth -s 3 -a 2 ";
		$cmd .= "-e " if (defined $infer_genptype_by_cffDNA);
		$cmd .= ">$fout.log 2>&1";
	}elsif($saminfo[0] eq "standards-thalassemia"){
		my $mix_ratio = $saminfo[3];
		$cmd = "Rscript $Bin/cffDNA_to_genotype.r -i $fin -c $mix_ratio -o $fout >$fout.log 2>&1";
	}

	if($#saminfo >1 && $saminfo[2] ne ""){
		my $genotype;
		if($saminfo[0] eq "gravida-plasma-thalassemia"){
			$genotype = join(":", @saminfo[2..4]);
		}elsif($saminfo[0] eq "standards-thalassemia"){
			my ($spot, @geno)=$saminfo[4]=~/(\S+)\(([-N]+)\/([-N]+)\)/;
			for(my $i=0; $i<@geno; $i++){
				if($geno[$i] eq "-"){
					$geno[$i] = $spot."/N";
				}elsif($geno[$i] eq "--"){
					$geno[$i] = $spot."/".$spot;
				}elsif($geno[$i] eq "N"){
					$geno[$i] = "N/N";
				}else{
					die "Wrong info $saminfo[2] of $s\n";
				}
			}
			$genotype = join(":",$geno[1],$geno[0],"");
		}
		if($genotype=~/^-/){
			$genotype ="'".$genotype;
		}

		$cmd .= " && $Bin/TF_evaluation -i $fout -g \"$genotype\" -o $fout >>$fout.log 2>&1";
	}
	print SH $cmd, "\n"; 
}
close (SH) ;

my $r = system("parallel -j 10 < $sh ");
die "infer fetal genotypes failed\n" if ($r);


## merge result
open(O,">$dOut/$fkey.genotype.result") or die $!;
print O "Sample\tMutID\tMutTotal\tDepTotal\tFreqTotal\tfetal.EM\tMaternalGenotype\tFetalGenotype\tRefined.genotype.EM\tp\tAverDepth\tAverDepthHotspot\tCheck\tFetalGiving\tMotherGiving\tFatherGiving\tTrueOrFalse\tRatio_E\tGenotype_E\n";
foreach my $s (sort {$a cmp $b} keys %config) {
	next if($config{$s}{"sample_info"}->[0] ne "gravida-plasma-thalassemia" && $config{$s}{"sample_info"}->[0] ne "standards-thalassemia");
	my $file = "$dOut/$s.genotype.txt";
	open(I, $file) or die $!;
	<I>;
	while(<I>){
		chomp;
		next if(/rs/);
		my @unit = split /\t/, $_;
		my ($cffdna, $genoM, $genoF, $p, $depa, $deph)=@unit[5,6,7,9,11,12];
		if(scalar @unit==19 || ($genoM ne "N/N" && $genoM ne "N" && $genoM ne "NA") || ($genoF ne "N/N" && $genoF ne "N" && $genoF ne "NA") ){
			## filter check
			my $is_filter = 0;
			my @finfo;
			if($cffdna < 0.05){
				$is_filter = 1;
				push @finfo, "LowCffdna";
			}
			if($depa < 1000){
				$is_filter = 1;
				push @finfo, "LowDep";
			}
			if($p < 0.75){
				$is_filter = 1;
				push @finfo, "Lowp";
			}
			my $check = scalar @finfo>0? join(",",@finfo): "Normal";
			if($#unit>13){
				print O join("\t",@unit[0..9], @unit[11,12], $check, @unit[13..$#unit]),"\n";
			}else{
				print O join("\t",@unit[0..9],@unit[11,12], $check),"\n";
			}
		}
	}
	close(I);
}
close(O);

# gather sample cffdna result
open (O,">$dOut/cffdna.txt") or die $!;
print O join("\t","Sample","SNPNum","AverDepth","AverDepthHotspot","%cffdna_EM",),"\n";
foreach my $s (sort {$a cmp $b} keys %config) {
	next if($config{$s}{"sample_info"}->[0] ne "gravida-plasma-thalassemia" && $config{$s}{"sample_info"}->[0] ne "standards-thalassemia");
	my $file = "$dOut/$s.genotype.txt";
	my ($cffdna, $nsnp, $dep, $dephot)=("NA","NA","NA","NA");
	if (-f $file) {
		my $line = `cat $file|head -2|tail -1`;
		chomp $line;
		($cffdna, $nsnp, $dep, $dephot)=(split /\t/, $line)[5,10..12];
	}

	print O join("\t",$s,$nsnp,$dep,$dephot,$cffdna),"\n";
}
close (O) ;


#

######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################
# ==============================================================
# sub function
# ==============================================================
sub adjust_by_aver{
	my ($fin, $fin_ad, $fadjust)=@_;
	my %adjust;
	open(L, $fadjust) or die $!;
	while(<L>){
		chomp;
		next if(/^$/ || /^#/);
		my ($id, $freq)=split;
		$adjust{$id}=$freq;
	}
	close(L);

	open(I, $fin) or die $!;
	open(O, ">$fin_ad") or die $!;
	while(<I>){
		chomp;
		next if(/^$/);
		my ($s, $id, $mut, $total, $freq)=split;
		if(exists $adjust{$id}){
			my $freq_ad = $freq*0.5/$adjust{$id};
			$mut = sprintf("%.0f",$total*$freq_ad);
			$freq = $mut/$total;
		}
		print O join("\t", $s,  $id, $mut, $total, $freq),"\n";
	}
	close(I);
	close(O);
}




sub GetTime {
	my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
	return sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

sub USAGE {#
	my $usage=<<"USAGE";
Program: $0
Version: $version
Contact: Ma chouxian <chouxian.ma\@genetalks.com> 
Discription:
	
Usage:
  Options:
  -i		<file>	input samplesheet file, forced
  -ic  		<file>   Input config file, [$fconfig]
  -cid 		<str>	config ID in config file, forced
  -id		<dir>	dir of input file(xxx.hotspots.result),xxx/06.hotspots_result/, forced
  -k		<str>	key of output file, forced
  -od		<dir>	output dir, optional, default [$dOut]
  -d        <int>   min depth, default [$min_depth]
  -e                infer genotype by cffdna
  -h			Help

USAGE
	print $usage;
	exit;
}

