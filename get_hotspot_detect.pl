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
my (@indir,$fsamsheet, $prefix,$fkey,$outdir, $date);
my ($fconfig,$configID);
GetOptions(
				"help|?" =>\&USAGE,
				"id:s{,}"=>\@indir,
				"is:s"=>\$fsamsheet,
				"ic:s"=>\$fconfig,
				"cid:s"=>\$configID,
				"d:s"=>\$date,
				"p:s"=>\$prefix,
				"k:s"=>\$fkey,
				"od:s"=>\$outdir,
				) or &USAGE;
&USAGE unless (@indir and $fsamsheet and $configID and $fkey);
$fconfig = defined $fconfig? $fconfig: "$Bin/config/config.txt";


$outdir||="./";
`mkdir $outdir`	unless (-d $outdir);
$outdir=AbsolutePath("dir",$outdir);

my $Det;
open($Det,">$outdir/$fkey.hotspots.variants") or die $!;
print $Det join("\t","Sample","MutID","MutTotal","DepTotal","FreqTotal","MutU","DepU","FreqU","MutD","DepD","FreqD","PrimerU","PrimerD"),"\n";

my $result_dir = "$outdir/05.hotspots_result";
`mkdir $result_dir` unless(-d $result_dir);


# get variants info
foreach my $indir(@indir){
	my %config;
	&read_data_sheet($fsamsheet, $fconfig, \%config, $configID);
	
	foreach my $s(sort {$a cmp $b}keys %config){
		open(O,">$result_dir/$s.hotspots.result") or die $!;
		print O join("\t","Sample","MutID","MutTotal","DepTotal","FreqTotal"),"\n";
		if(defined $date){
			$s=~s/_$date//;
		}
		my $fvcf = "$indir/$s.vcf";
		my $ftarget = $config{$s}{"ftarget"};
		my $ffusion = "$indir/$s.fusion.final.txt";
		my $ftarget_fusion = $config{$s}{"fusion_pos"};
		my %detect;
		if(!-e $fvcf or !-e $ffusion){
			print "Warn: vcf file or fusion file not found, $fvcf or $ffusion\n";
			next;
		}
		&get_variant_vcf(\%detect, $fvcf, $ftarget, $config{$s}{"fregion_primer"});
		&get_fusion_result(\%detect, $ffusion, $ftarget_fusion, $fvcf);
	
		foreach my $id(sort{$detect{$b}->[2]<=>$detect{$a}->[2]}keys %detect){
			my @det=@{$detect{$id}}; ##$depa,$dept,$ratio, $depaU,$deptU,$ratioU, $depaD,$deptD,$ratioD,$primerU, $primerD
			next if($id=~/rs/ || $det[0]==0);
			print $Det join("\t",$s, $id, @{$detect{$id}}),"\n";
		}
		foreach my $id(sort {$a cmp $b} keys %detect){
			my @det=@{$detect{$id}}; ##$depa,$dept,$ratio, $depaU,$deptU,$ratioU, $depaD,$deptD,$ratioD,$primerU, $primerD
			print O join("\t", $s, $id, @det[0..2]),"\n";
		}
	}
	close(O);
}
close($Det);


#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################

# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------

sub get_fusion_result{
	my ($adetect, $fresult, $ftarget, $fvcf)=@_;
	
	my %fusion;
	my %fusion_id;
	open(T, $ftarget) or die $!;
	while(<T>){
		chomp;
		next if(/^$/ || /^Region/);
		my ($type,$id, $c1, $p1, $c2, $p2)=split;
		$c1 =~s/chr//;
		$c2 =~s/chr//;
		$fusion{$c1}{$p1}{$c2}{$p2}=$id;
		@{$fusion_id{$id}}=($c1, $p1, $c2, $p2);
	}
	close(T);
	
	open(R, $fresult) or die $!;
	<R>;
	while(<R>){
		chomp;
		next if(/^$/ || /\?/);
		my ($br1, $br2, $sup1, $sup2, $total1, $total2)=split;
		my ($c1, $p1)=split /,/, $br1;
		my ($c2, $p2)=split /,/, $br2;
		$c1 =~s/chr//;
		$c2 =~s/chr//;
		if(exists $fusion{$c1}{$p1}{$c2}{$p2}){
			my $id = $fusion{$c1}{$p1}{$c2}{$p2};
			##check
			if($sup1!=$sup2){
				print "Wrong: support1 is not equal with support2, $_\n";
				die;
			}
			my $depa = $sup1;
			my $dept = $total1==$sup1? $total2: ($total2==$sup1? $total1: $total1+$total2); ## 某个断点，若total==sup，说明该断点未设计引物 
			@{$adetect->{$id}}=($depa, $dept, $depa/$dept);
		}
	}
	close(R);

	foreach my $id(keys %fusion_id){
		if(!exists $adetect->{$id}){
			my ($c1, $p1, $c2, $p2) = @{$fusion_id{$id}};
			my $line = `less $fvcf|awk '\$1==\"chr$c1\"'|awk '\$2==$p1'`;
			my ($ad)=$line=~/AD=([\d,]+)/;
			my @ad = defined $ad? split /,/, $ad: ();
			my $line2 = `less $fvcf|awk '\$1==\"chr$c2\"'|awk '\$2==$p2'`;
			my ($ad2)=$line2=~/AD=([\d,]+)/;
			my @ad2 = defined $ad2? split /,/, $ad2:();
			my $dep;
			foreach my $d(@ad, @ad2){
				$dep+=$d;
			}
			if(defined $dep){
				@{$adetect->{$id}}=(0, $dep, 0);
			}
		}
	}
}

sub get_variant_vcf{
	my ($adetect, $fvcf, $ftarget,$fpr)=@_;

	## get target variant
	#11      5247992 5247992 C       A       ID=HBB:CD43;GENE=UNKOWN;STRAND=-;CDS=c.130G>T;AA=UNKOWN
	#11      5247993 5247996 AAAG    -       ID=HBB:CD41-42;GENE=UNKOWN;STRAND=-;CDS=c.126_129delCTTT;AA=UNKOWN
	#11      5247994 5247997 AAGA    -       ID=HBB:CD41-42;GENE=UNKOWN;STRAND=-;CDS=c.126_129delCTTT;AA=UNKOWN
	#11      5248166 5248166 -       G       ID=HBB:CD27/28;GENE=UNKOWN;STRAND=-;CDS=c.84_85insC;AA=UNKOWN
	my %target;
	my %indel_pos;
	my %id;
	if(!-e $ftarget){
		print $ftarget,"\n";
		die;
	}
	open(T, $ftarget) or die $!;
	while(<T>){
		chomp;
		next if(/^$/);
		my ($c, $s, $e, $ref, $alt)=split;
		my ($id)=$_=~/ID=(\S+?);/;
		$id=~s/HBB://;
		$id=~s/HBA2://;
		my $lenr = length $ref;
		my $lena = length $alt;
		my ($pos1, $pos2);
		my $type = "SNP";
		if($alt eq "-"){ #del
			$pos1 = $s;
			$pos2 = $e;
			$type = "Del";
			$s --; ## self
		}elsif($ref eq "-"){
			$pos1 = $s;
			$pos2 = $e+1;
			$type = "Ins";
		}
		if($type ne "SNP"){
			$indel_pos{$c}{$pos1}=1;
			$indel_pos{$c}{$pos2}=1;
			$id{$c}{$pos1}{$id}=1;
			$id{$c}{$pos2}{$id}=1;
		}
		$id{$c}{$s}{$id}=1;
		if(exists $target{$c}{$s}{$id}){
			if($ref ne $target{$c}{$s}{$id}->[3]){
				print "Wrong: $id has more different ref in $ftarget!\n";
				die;
			}
			push @{$target{$c}{$s}{$id}}, $alt;
		}else{
			@{$target{$c}{$s}{$id}}=($type, $pos1, $pos2,$ref, $alt);
		}
	}
	close(T);


	my %primer_spe;
	open(R, $fpr) or die $!;
	while(<R>){
		chomp;
		next if(/^$/ || /^#/);
		my ($c, $s, undef, $primers)=split;
		$c =~s/chr//;
		next if(!exists $target{$c}{$s});
		my @primer = split /,/, $primers;
		foreach my $id(keys %{$target{$c}{$s}}){
			foreach my $p (@primer){
				if($p=~/-U-/){
					$primer_spe{$id}{"U"}{$p}=1;
				}elsif($p=~/-D-/){
					$primer_spe{$id}{"D"}{$p}=1;
				}
			}
		}
	}
	close(R);


	open(I, $fvcf) or die $!;
	my %UD_depth;
	my $adepth;
	while(<I>){
		chomp;
		next if(/^$/ || /^#/);
		my ($c, $p, $ref, $alt, $info)=(split /\t/, $_)[0,1,3,4,7];
		$c=~s/chr//;
		next if(!exists $indel_pos{$c}{$p} && !exists $target{$c}{$p});
#		my $is_next = 1;
#		$is_next = 0 if(exists $indel_pos{$c}{$p} && $alt=~/\<\*\>/);
#		foreach my $var(keys %{$target{$c}{$p}}){
#			my ($id, $tref, $talt, $type) = @{$target{$c}{$p}{$var}};
#			$is_next = 0 if($type eq "SNP" && $alt=~/$talt/);
#			$is_next = 0 if($type ne "SNP" && $alt!~/\<\*\>/);
#		}
#		next if($is_next);
	
		my @id = keys %{$id{$c}{$p}};
		my $id=$id[0];

		##check
		my $primerU = join(",", sort{$a cmp $b} keys %{$primer_spe{$id}{"U"}});
		my $primerD = join(",", sort{$a cmp $b} keys %{$primer_spe{$id}{"D"}});
		foreach my $idt(@id){
			next if($idt eq $id);
			my $U = join(",", sort{$a cmp $b} keys %{$primer_spe{$idt}{"U"}});
			my $D = join(",", sort{$a cmp $b} keys %{$primer_spe{$idt}{"D"}});
			if($U ne $primerU && $D ne $primerD){
				print "Wrong: the pos $c $p has more primers requirement in $fpr, because the depth of the pos is need by nearby indel!\n";
				die;
			}
		}

		if($alt=~/\<\*\>/){
			$adepth = \%{$UD_depth{"SNP"}};
		}else{
			$adepth = \%{$UD_depth{"Indel"}};
		}

		## split
		my ($pr, $adpr)=$info=~/PR=(\S+?);ADPR=(\S+?);/;
		my @pr=split /,/, $pr;
		my @pr_unit =split /:/, $adpr;
		my %adpr;
		for(my $i=0; $i<@pr_unit; $i++){
			my @bunit = split /,/, $pr_unit[$i];
			for(my $j=0; $j<@bunit; $j++){
				my ($b, $d)=split /-/, $bunit[$j];
				$adpr{$pr[$i]}{$b}=$d;
			}
		}
		$adepth->{$c}{$p}{"U"}{"Total"}=0;
		$adepth->{$c}{$p}{"D"}{"Total"}=0;
		$adepth->{$c}{$p}{"Total"}{"Total"}=0;
		foreach my $ud(keys %{$primer_spe{$id}}){
			foreach my $pr(keys %{$primer_spe{$id}{$ud}}){
				foreach my $b (keys %{$adpr{$pr}}){
					my $d = $adpr{$pr}{$b};
					$adepth->{$c}{$p}{$ud}{$b}+=$d;
					$adepth->{$c}{$p}{$ud}{"Total"}+=$d;
					$adepth->{$c}{$p}{"Total"}{$b}+=$d;
					$adepth->{$c}{$p}{"Total"}{"Total"}+=$d;
				}
			}
		}
	}
	close(I);

	foreach my $c (keys %target){
		foreach my $s(keys %{$target{$c}}){
			foreach my $id(keys %{$target{$c}{$s}}){
				next if(!exists $UD_depth{"SNP"}{$c}{$s} && !exists $UD_depth{"Indel"}{$c}{$s});
				my ($type, $pos1, $pos2, $ref, @alt)=@{$target{$c}{$s}{$id}};
				my $uddep = $type eq "SNP"? $UD_depth{"SNP"}{$c}{$s}: $UD_depth{"Indel"}{$c}{$s};
				my ($deptU, $depaU, $deptD, $depaD, $dept, $depa);
				foreach my $alt(@alt){
					$alt = $ref if($type eq "Del");
					$depaU += exists $uddep->{"U"}{$alt}? $uddep->{"U"}{$alt}: 0;
					$depaD += exists $uddep->{"D"}{$alt}? $uddep->{"D"}{$alt}: 0;
				}
				$depa = $depaU+$depaD;
				if($type eq "SNP"){
					$deptU = $uddep->{"U"}{"Total"};
					$deptD = $uddep->{"D"}{"Total"};
					$dept = $uddep->{"Total"}{"Total"};
				}else{
					my ($dpos1U, $dpos1D) = exists $UD_depth{"SNP"}{$c}{$pos1}? ($UD_depth{"SNP"}{$c}{$pos1}{"U"}{"Total"}, $UD_depth{"SNP"}{$c}{$pos1}{"D"}{"Total"}):(0,0);
					my ($dpos2U, $dpos2D) = exists $UD_depth{"SNP"}{$c}{$pos2}? ($UD_depth{"SNP"}{$c}{$pos2}{"U"}{"Total"}, $UD_depth{"SNP"}{$c}{$pos2}{"D"}{"Total"}):(0,0);
					my $deprU = $dpos1U > $dpos2U? $dpos1U: $dpos2U;
					my $deprD = $dpos1D > $dpos2D? $dpos1D: $dpos2D;
					if($type eq "Ins"){
						$dept = $deprU+$deprD;
						$deptU = $deprU;
						$deptD = $deprD;
					}else{
						$deptU = $deprU+$depaU;
						$deptD = $deprD+$depaD;
						$dept = $deptU+$deptD;
					}
				}
				my( $ratio, $ratioU, $ratioD)=($dept>0? $depa/$dept:0, $deptU>0?$depaU/$deptU:0, $deptD>0?$depaD/$deptD:0);
				my $primerU = exists $primer_spe{$id}{"U"}? join(",",sort {$a cmp $b} keys %{$primer_spe{$id}{"U"}}):"-";
				my $primerD = exists $primer_spe{$id}{"D"}? join(",",sort {$a cmp $b} keys %{$primer_spe{$id}{"D"}}):"-";
				next if(exists $adetect->{$id} && $adetect->{$id}->[0] > $depa);
				@{$adetect->{$id}}=($depa,$dept,$ratio, $depaU,$deptU,$ratioU, $depaD,$deptD,$ratioD,$primerU, $primerD);
			}
		}
	}
#	print Dumper %{$adetect};
#	die;
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
  -id  <dir,>   indirs of xxx.vcf and xxx.fusion.final.txt. forced
  -is  <file>   samplesheet,forced
  -ic  <file>   Input config file, [$Bin/config/config.txt]
  -cid <str>	config ID in config file, forced
  -d   <str>    date string like 161206, optinal

  -k  <str>	Key of output file, forced
  -od <dir>	Dir of output file, default ./
  -h		 Help

USAGE
	print $usage;
	exit;
}

