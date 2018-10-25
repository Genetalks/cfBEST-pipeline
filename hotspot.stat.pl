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
my ($fbam, $fvariants, $fusion_pos, $fkey,$extend,$outdir);
my $draw;
GetOptions(
				"help|?" =>\&USAGE,
				"i1:s"=>\$fbam,
				"i2:s"=>\$fvariants,
				"i3:s"=>\$fusion_pos,
				"k:s"=>\$fkey,
				"draw:s"=>\$draw,
				"e:s"=>\$extend,
				"od:s"=>\$outdir,
				) or &USAGE;
&USAGE unless ($fbam and $fvariants and $fkey);

$outdir||="./";
`mkdir $outdir`	unless (-d $outdir);
$outdir=AbsolutePath("dir",$outdir);
$extend = defined $extend? $extend: 200;

my $cmd;
open (SH,">$outdir/$fkey.stat.sh") or die $!;

#==================================================================
# get depth file
#==================================================================

my $fdepth = "$outdir/$fkey.depth";
$cmd = "perl $Bin/depth.stat.pl -i $fbam -o $fdepth";
print SH $cmd,"\n";
`$cmd`;

#==================================================================
# get hotspot depth
#==================================================================
my %depth;
my %ID;
open (V, $fvariants) or die $!;
while (<V>) {
	chomp;
	next if(/^$/ || /^Chr/);
	my ($c, $p, $p2, $ref, $alt, $info)=split /\s+/, $_;
	$c=~s/chr//;
	my ($id)=split /;/, $info;
	$id=~s/ID=//;
	$ID{$c}{$id}=$p;
}
close (V);
if(defined $fusion_pos){
	open (I,"$fusion_pos") or die $!;
	while (<I>) {
		chomp;
		next if(/^$/);
		my ($type, $id, $c, $p, @pos)=split;
		$c=~s/chr//;
		if($type eq "Break"){
			$pos[0]=~s/chr//;
			$ID{$c}{$id."_1"}=$p;
			$ID{$pos[0]}{$id."_2"}=$pos[1];
			next;
		}
	}
	close(I);
}


open (D, $fdepth) or die $!;
while (<D>) {
	chomp;
	next if(/^$/);
	my ($c, $p, $d)=split /\s+/, $_;
	$c=~s/chr//;
	$depth{$c}{$p}=$d;
}
close (D);

my %stat;
my $num = 0;
open (O, ">$outdir/$fkey.hotspot.depth") or die $!;
print O "#HotSpot\tDepth\n";
foreach my $chr (sort {$a cmp $b} keys %ID) {
	foreach my $id (sort {$a cmp $b} keys %{$ID{$chr}}) {
		my $pos = $ID{$chr}{$id};
		print O $id,"\t", exists $depth{$chr}{$pos}? $depth{$chr}{$pos}: 0,"\n";

		$num ++;
		for (my $i=$pos-$extend; $i<=$pos+$extend ; $i++) {
			my $dis = $i-$pos;
			$stat{$dis}+=exists $depth{$chr}{$i}? $depth{$chr}{$i}: 0;
		}
	}
}
close (O);
open (O, ">$outdir/$fkey.hotspot.depth") or die $!;
print O "#HotSpot\tDepth\n";
foreach my $chr (sort {$a cmp $b} keys %ID) {
	foreach my $id (sort {$a cmp $b} keys %{$ID{$chr}}) {
		my $pos = $ID{$chr}{$id};
		print O $id,"\t", exists $depth{$chr}{$pos}? $depth{$chr}{$pos}: 0,"\n";

		$num ++;
		for (my $i=$pos-$extend; $i<=$pos+$extend ; $i++) {
			my $dis = $i-$pos;
			$stat{$dis}+=exists $depth{$chr}{$i}? $depth{$chr}{$i}: 0;
		}
	}
}
close (O);


#==================================================================
# get fusion region depth
#==================================================================
if (defined $fusion_pos) {
	my $fusion_bed = "$outdir/$fkey.fusion.bed";
	my $fusion_num = 0;
	open (FB,">$fusion_bed") or die $!;
	open (I,"$fusion_pos") or die $!;
	while (<I>) {
		chomp;
		next if(/^$/);
		my ($type, $id, $c, $p, @pos)=split;
		next if($type eq "Break");
		print FB join("\t",$c, $p, $pos[0], $id),"\n";
		$fusion_num++;
		
		if(defined $draw){
			my $nid = $c."_".$id;
			open (B,">$outdir/$nid.bed") or die $!;
			print B join("\t",$c, $p, $pos[0], $id),"\n";
			close (B) ;
			
			## 碱基深度统计
			$cmd="$Bin/third_party/samtools depth -b $outdir/$nid.bed $fbam -m 1000000 |gzip >$outdir/$fkey.$nid.depth.gz";
			print SH $cmd,"\n";
			`$cmd`;
			$cmd = "zcat $outdir/$fkey.$nid.depth.gz|cut -f 2,3 >$outdir/$fkey.$nid.depth.dis";
			print SH $cmd,"\n";
			`$cmd`;
		}
	}
	close(FB);
	close (I);
	if($fusion_num==0){
		`rm $fusion_bed`;
	}else{
		# 绘制均一性图
		#------------------------------------------------------------------
		my $totallen = `less -S $fusion_bed|awk \'{s+=\$3-\$2} END{print s}\'`;
		chomp $totallen;
		my %cover_depth;
		$cmd="$Bin/third_party/samtools depth -b $fusion_bed $fbam -m 1000000 |gzip >$outdir/$fkey.fusion.depth.gz";
		print SH $cmd,"\n";
		`$cmd`;
		my ($coverage, $coverage_100, $coverage_500, $averdepth)=&cover_depth("$outdir/$fkey.fusion.depth.gz",\%cover_depth, $totallen);
	
		open (OUT,">$outdir/$fkey.fusion.stat") or die $!;
		print OUT "SampleID\tAverDepth\tCoverage(>0X)\tCoverage(>100X)\tCoverage(>500X)\n";
		print OUT $fkey,"\t",$averdepth,"\t",$coverage,"\t", $coverage_100,"\t", $coverage_500, "\n";
		close (OUT) ;
		
		if(defined $draw){
			&draw_coverage_rate_cutoff_depth_png(\%cover_depth, $totallen);
		}
	}
}

#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################

# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------
sub draw_coverage_rate_cutoff_depth_png {#
	my ($cover_depth, $totallen) = @_;

	open OUT,">$outdir/$fkey.coverage_rate_cutoff_depth.txt" or die $!;
	my @key=sort {$a<=>$b} keys %{$cover_depth};
	foreach my $k (sort {$a<=>$b} keys %{$cover_depth}) {
		if ($k%10==1 && $k<1000) {
			my $temp_depth=0;
			for (my $i=$k;$i<=$key[-1] ;$i++) {
				if (exists $cover_depth->{$i}) {
					$temp_depth+=$cover_depth->{$i};
				}
			}
			my $per=$temp_depth/$totallen*100;
			print OUT "$k\t$per\n";
		}
	}
	close OUT;
}


sub cover_depth {
	my ($depth, $cover_depth, $totallen)=@_;
	my ($totaldep, $coverlen, $cover_100_len, $cover_500_len)=(0,0,0,0);

	$/="\n";
	open IN,"gzip -dc $depth| " or die $!;
	while (<IN>) {
		chomp;
		next if (/^$/);
		my ($chr,$pos,$dep)=(split /\s+/,$_)[0,1,2];
		$cover_depth->{$dep}++;
		$totaldep+=$dep;
		if ($dep>0) {  
			$coverlen++;
		}
		if ($dep>100) {
			$cover_100_len++;
		}
		if ($dep>500) {
			$cover_500_len++;
		}
	}
	close IN;

	my $coverage=(join "",(sprintf "%.2f",($coverlen/$totallen*100)),"%");
	my $coverage_100=(join "",(sprintf "%.2f",($cover_100_len/$totallen*100)),"%");
	my $coverage_500=(join "",(sprintf "%.2f",($cover_500_len/$totallen*100)),"%");
	my $averdepth=$coverlen>0? sprintf "%.2f",($totaldep/$coverlen): 0;

	return ($coverage, $coverage_100, $coverage_500, $averdepth);
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
	v2: hotspots be extended

Usage:
  Options:
  -i1  <file>   Input bam file, forced
  -i2  <file>   Input variants known file, forced
  -i3  <file>   Input fusion bed file, optional
  -k  <str>    key of Output file, forced

  --draw 		draw png
  -e  <int>    extend length near hotspot, [200]
  -od <dir>    outdir of output file, default ./
  -h         Help

USAGE
	print $usage;
	exit;
}
