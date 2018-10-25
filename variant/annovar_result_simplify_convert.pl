#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $BEGIN_TIME=time();
my $version="1.0.0";
#######################################################################################

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($fIn,$ftarget, $sample,$outdir,$thalassemia);
my $min_depth = 50;
my $min_alt_depth = 1;
GetOptions(
				"help|?" =>\&USAGE,
				"i1:s"=>\$fIn,
				"i2:s"=>\$ftarget,
				"s:s"=>\$sample,
				"thalassemia:s"=>\$thalassemia,
				"d:s"=>\$min_depth,
				"a:s"=>\$min_alt_depth,
				"od:s"=>\$outdir,
				) or &USAGE;
&USAGE unless ($fIn and $sample);

$outdir||="./";
`mkdir $outdir`	unless (-d $outdir);
$outdir=AbsolutePath("dir",$outdir);

my %detect;
my %index_mult;
my %panle_gene;
my $panle_name = (split(/\_/,$ftarget))[1];
open F,"$ftarget" || die "$!\n";
while(<F>){
	chomp;
	my @aa = split(/\t/,$_);
	if($aa[-1]=~/GENE=/){
		my $gene = (split(/\;/,(split(/GENE\=/,$aa[-1]))[1]))[0];
		$panle_gene{$gene} = 1;
	}else{
		my $gene = (split(/\_/,(split(/\=/,$aa[-1]))[1]))[0];
		$panle_gene{$gene} = 1;
	}
}
close F;
open (O,">$outdir/$sample.variant_total.txt") or die $!;
open (I,"$fIn") or die $!;
my $head = <I>;
chomp $head;
my @hunit = split /\t/, $head;
my $anno_num = scalar @hunit - 1;
$hunit[$anno_num-1] = "Target";
print O join("\t",@hunit[0..($anno_num-1)]),"\tVarDepth\tTotalDepth\tFreq\n"; 
while (<I>) {
	chomp;
	next if(/^$/ || /^#/);
	my @unit = split /\t/, $_;
	my $exists=0;
	if(defined $thalassemia){
#	if($panle_name =~ /thalassemia/){
		$exists = 1;
	}else{
		if($unit[6]=~/,|;/){
			my @gene_out = split(/\,|\;/,$unit[6]);
			for(my $ii=0;$ii<=$#gene_out;$ii++){
				if(exists $panle_gene{$gene_out[$ii]}){
					$exists += 1;
				}
			}
		}else{
			if(exists $panle_gene{$unit[6]}){
				$exists += 1;
			}
		}
	}
	next if($exists == 0);
	my ($chr, $pos, $id, $ref, $alt, $qual, $filter, $info, $format, @sample)=@unit[($anno_num+3)..$#unit];

	my ($dpr)=$info=~/AD=([\d,]+)/;
	my ($depr, @dep) = split /,/, $dpr;
	my $depth_total = $depr;
	for (my $i=0; $i<@dep ;$i++) {
		$depth_total+=$dep[$i];
	}
	my @alt = split /,/, $alt;
	my $index = 0;
	if (scalar @alt > 1) {
		$index_mult{$chr.".".$pos.".".$alt}++;
		$index = $index_mult{$chr.".".$pos.".".$alt} - 1;
	}

#	next if($depth_total < $min_depth);
#	next if($dep[$index] < $min_alt_depth);
	if(!defined $dep[$index]){
		print $_,"\n";
		die;
	}
	my $det = $dep[$index].",".$depth_total."_".sprintf("%.2f",$dep[$index]/($depth_total)*100)."%";
	if ($unit[$anno_num-1] ne ".") {
		push @{$detect{$unit[$anno_num-1]}},$det;
	}

	print O join("\t",@unit[0..($anno_num-1)]),"\t", join("\t",$dep[$index], $depth_total, sprintf("%.2f",$dep[$index]/($depth_total)*100)."%"),"\n";
	
}
close (I) ;
close (O);

open (O,">$outdir/$sample.variants_known.result") or die $!;
#1       27087413        27087413        A       T       ID=COSM4143737;GENE=ARID1A;STRAND=+;CDS=c.1987A>T;AA=p.T663S;CNT=1
#1       27087417        27087417        C       G       ID=COSM79135;GENE=ARID1A;STRAND=+;CDS=c.1991C>G;AA=p.S664*;CNT=2
#1       27087418        27087420        AGG     A       ID=COSM3358384;GENE=ARID1A;STRAND=+;CDS=c.1993_1994delGG;AA=p.G665fs*10;CNT=1
#11      5247904 5247904 -       A       ID=HBB:CD71-72;GENE=UNKOWN;STRAND=-;CDS=c.216_217insA/T;AA=UNKOWN
#11      5247904 5247904 -       T       ID=HBB:CD71-72;GENE=UNKOWN;STRAND=-;CDS=c.216_217insA/T;AA=UNKOWN
#11      5247993 5247996 AAAG    -       ID=HBB:CD41-42;GENE=UNKOWN;STRAND=-;CDS=c.126_129delCTTT;AA=UNKOWN
#11      5247994 5247997 AAGA    -       ID=HBB:CD41-42;GENE=UNKOWN;STRAND=-;CDS=c.126_129delCTTT;AA=UNKOWN
#
print O "Chr\tStart\tEnd\tRef\tAlt\tGene\tStrand\tMut-CDS\tMut-AA\tMutID\tMut-description\t$sample\n";
my %info;
my %result;
open (V,"$ftarget") or die $!;
while (<V>) {
	chomp;
	my ($chr, $start, $end, $ref, $alt, $info)=split /\t/, $_;
	if(exists $info{$info}){
		my ($ref0, $alt0)=@{$info{$info}}[3,4];
		my @refs = split /\//, $ref0;
		my @alts = split /\//, $alt0;
		my %refs;
		my %alts;
		map($refs{$_}=1, @refs, $ref);
		map($alts{$_}=1, @alts, $alt);
		$info{$info}->[3]=join("/", sort {$a cmp $b}keys %refs);
		$info{$info}->[4]=join("/", sort {$a cmp $b}keys %alts);
		next;
	}
	my @unit=split /;/, $info;
	my %value;
	for (my $i=0;$i<@unit ;$i++) {
		my ($key,$value)=split /=/, $unit[$i];
		$value{$key}=$value;
		if ($key eq "CNT" || $key eq "SNP" || $key eq "TARGET") {
			push @{$value{"DES"}}, $unit[$i];
		}
	}
	
	my @des=exists $value{"DES"}? @{$value{"DES"}}: (".");
	my $gene = exists $value{"GENE"}? $value{"GENE"}: ".";
	my $strand = exists $value{"STRAND"}? $value{"STRAND"}: ".";
	my $cds = exists  $value{"CDS"}?  $value{"CDS"}: ".";
	my $aa = exists $value{"AA"}? $value{"AA"}: ".";
	my $id = exists $value{"ID"}? $value{"ID"}: ".";
	@{$info{$info}}=("chr".$chr, $start,$end,$ref,$alt, $gene, $strand, $cds, $aa, $id, join(";",@des));
	
	my $detect;
	if(!exists $detect{$info}){
		$detect = "-";
	}else{
		my ($mutant, $total)=(0,0);
		foreach my $det(@{$detect{$info}}){
			my ($m, $t, $f)=$det=~/(\d+),(\d+)_([\d\.]+)%/;
			$mutant+=$m;
			$total+=$t;
		}
		my $freq = $mutant/$total;
		$detect = $mutant.",".$total."_".sprintf("%.2f",$freq*100)."%";
	}
	$result{$info}=$detect;
}

foreach my $info(sort {$a cmp $b} keys %result){
	print O join("\t",  @{$info{$info}}, $result{$info}),"\n";
}
close (V) ;
close (O);

#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################

# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------
sub AbsolutePath
{		#获取指定目录或文件的决定路径
        my ($type,$input) = @_;

        my $return;
		$/="\n";

        if ($type eq 'dir')
        {
                my $pwd = `pwd`;
                chomp $pwd;
                chdir($input);
                $return = `pwd`;
                chomp $return;
                chdir($pwd);
        }
        elsif($type eq 'file')
        {
                my $pwd = `pwd`;
                chomp $pwd;

                my $dir=dirname($input);
                my $file=basename($input);
                chdir($dir);
                $return = `pwd`;
                chomp $return;
                $return .="\/".$file;
                chdir($pwd);
        }
        return $return;
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
  -i1  <file>   Input annovar result file, forced
  -i2  <file>   Input target variants file, forced
  --thalassemia  optional
  -s   <str>    sample name, forced
  
  -od <dir>    Dir of output file, default ./
  -h         Help

USAGE
	print $usage;
	exit;
}
