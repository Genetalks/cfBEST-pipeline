
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
my (@indir, $prefix,$fkey,$outdir, $date);
GetOptions(
				"help|?" =>\&USAGE,
				"id:s{,}"=>\@indir,
				"d:s"=>\$date,
				"p:s"=>\$prefix,
				"k:s"=>\$fkey,
				"od:s"=>\$outdir,
				) or &USAGE;
&USAGE unless (@indir and $fkey);

$outdir||="./";
`mkdir $outdir`	unless (-d $outdir);
$outdir=AbsolutePath("dir",$outdir);


if (defined $prefix) {
	`cp "$prefix.stat" "$outdir/$fkey.stat"`;
	open (O,">>$outdir/$fkey.stat") or die $!;
}else{
	open (O,">$outdir/$fkey.stat") or die $!;
	print O join("\t","Sample", "TotalReads", "UniqRatio","BarErrorRatio", "PrimerRatio", "PrimerReads", "SoftClipRatio", "LowQualityRatio", "OffTargetRatio", "PosErrorRatio", "InsertSize", "SEratio", "LowDepth", "ConsensError", "TotalUniq", "TotalUniq/TotalReads", "DupRatio", "RemainRatio", "AverDepth", "VarTotal", "VarPass", "Description"),"\n";
}

foreach my $indir (@indir) {
	my %stat;
#	my ($date)=$indir =~/(\d\d\d\d\d\d)_NS500/;

	## get TotalReads BarErrorRatio 
	get_reads_and_uniq_stat(\%stat, $date, "$indir/statistic/01.Uniq.stat");

	## get_barcode_stat
	get_barcode_stat(\%stat, $date, "$indir/statistic/02.Error_Bar.stat");

	## primer stat
	&get_primer_stat(\%stat, $date, "$indir/statistic/03.Primer.stat");

	## map stat
	&get_map_stat(\%stat, $date, "$indir/statistic/04.Map.stat");

	## group stat
	&get_dedup_stat(\%stat, $date, "$indir/statistic/04.Dedup.stat");

	## Amplicon.stat
	&get_uniq_stat(\%stat, $date, "$indir/statistic/04.Uniq.stat");

	## 04.PrimerDiff.stat
	#&get_remain_stat(\%stat, $date, "$indir/statistic/04.PrimerDiff.stat");

	## 05.Hotspot.stat
	&get_depth_stat(\%stat, $date, "$indir/statistic/05.Hotspot.stat");

	## 00.Sample.des
	&get_sample_description(\%stat, $date, "$indir/statistic/00.Sample.des");

	### output
	foreach my $sample (sort {$a cmp $b} keys %stat) {
		print O join("\t",$sample, @{$stat{$sample}}),"\n";
	}

}
close (O) ;

#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################

# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------

#SampleID        Index   Description
#CT1600012B      CGGACGT ÔçÉ¸50gene£¬µÚ¶þ´Î£¬DSN free
#KYCT058 TTAATAG ÔçÉ¸50gene£¬DSN free
sub get_sample_description {#
	my ($astat, $date, $file) = @_;
	open (I,"$file") or die $!;
	while (<I>) {
		chomp;
		next if(/^$/ || /Sample/);
		my ($SampleID, @Descript)=split /\t/, $_;
		if ($SampleID!~/_/ && defined $date) {
			$SampleID.="_$date";
		}
		push @{$astat->{$SampleID}},@Descript;
	}
	
	close (I) ;
}

##Info   CTA1600008
#Total   1577
#CrowdFilter     3
#DepthFilter     41
#PASS    273
#PolyFilter      33
#VariantFilter   1227



##HotSpot        KYCT102B2       KYCT104B2       KYCT105B2       KYCT107B2  
#Average 513.77  55.32   328.64  1965.57 598.06  404.88  461.89  649.34  143
#chr1.115256529.NRAS.Target      365     8       103     789     454     227
sub get_depth_stat {#
	my ($astat, $date, $file) = @_;
	open (I,"$file") or die $!;
	my $head = <I>;
	my $depth = <I>;
	chomp $head;
	chomp $depth;
	my (undef, @sample)=split /\t/, $head;
	my (undef, @depth)=split /\t/, $depth;
	for (my $i=0; $i<@sample ;$i++) {
		my $SampleID = $sample[$i];
		if ($SampleID!~/_/ && defined $date) {
			$SampleID.="_$date";
		}
		push @{$astat->{$SampleID}},($depth[$i]);
	}
	close (I) ;
}

#SampleID        Total   Filter  Error   Remain  RemainRatio
#KYCT102B2       100161  14171   101     85889   85.75%
#KYCT104B2       9152    569     2       8581    93.76%
sub get_remain_stat {#
	my ($astat, $date, $file) = @_;
	open (I,"$file") or die $!;
	while (<I>) {
		chomp;
		next if(/^$/ || /Sample/);
		my ($SampleID, $RemainRatio)=(split /\t/, $_)[0,5];
		if ($SampleID!~/_/ && defined $date) {
			$SampleID.="_$date";
		}
		push @{$astat->{$SampleID}},($RemainRatio);
	}
	
	close (I) ;
}


#SampleID        TotalUniq       PosUniq BarUniq TotalReads      Duplicate       DupRatio
#KYCT102B2       100161  38455   61706   4019430 3919269 97.51%
#KYCT104B2       9152    5540    3612    72456   63304   87.37%
sub get_uniq_stat {#
	my ($astat, $date, $file) = @_;
	open (I,"$file") or die $!;
	while (<I>) {
		chomp;
		next if(/^$/ || /Sample/);
		my ($SampleID, $TotalUniq, $DupRatio)=(split /\t/, $_)[0,1,6];
		if ($SampleID!~/_/ && defined $date) {
			$SampleID.="_$date";
		}
		push @{$astat->{$SampleID}},($TotalUniq, sprintf("%.2f",$TotalUniq/$astat->{$SampleID}->[0]*100)."%", $DupRatio);
	}
	
	close (I) ;
}

#SampleID        TotalReads      TotalUniq       PosUniq BarUniq Duplicate       DupRatio        SEratio LowDepth        ConsensError    Normal
#PAGBI10BGR2     6958257 17460   6544    10916   6940797 99.75%  0.00%   0.10%   0.34%   99.56%
#PAGBI10BGR3     6740465 13137   6129    7008    6727328 99.81%  0.00%   0.07%   0.67%   99.27%
sub get_dedup_stat {#
	my ($astat, $date, $file) = @_;
	open (I,"$file") or die $!;
	while (<I>) {
		chomp;
		next if(/^$/ || /Sample/);
		my ($SampleID, $SEratio, $LowDepth, $ConsensError)=(split /\t/, $_)[0,7,8,9];
		if ($SampleID!~/_/ && defined $date) {
			$SampleID.="_$date";
		}
		push @{$astat->{$SampleID}},($SEratio, $LowDepth, $ConsensError);
	}
	
	close (I) ;
}


#SampleID        TotalPair       Mapratio        SEratio ProperRatio     SoftClipRatio   DiffChrRatio    LowQualityRatio OffTargetRatio  PosErrorRatio   InsertSize      SD
#PAGBI10BGR2     7152808 100.00% 0.45%   97.96%  12.08%  1.34%   0.24%   2.20%   0.01%   133.07  52.09
#PAGBI10BGR3     6899543 100.00% 0.66%   96.80%  11.66%  2.20%   0.16%   1.59%   0.01%   133.05  50.18
sub get_map_stat {#
	my ($astat, $date, $file) = @_;
	open (I,"$file") or die $!;
	while (<I>) {
		chomp;
		next if(/^$/ || /Sample/);
		my ($SampleID, $SoftClipRatio, $LowQualityRatio, $OffTargetRatio, $PosErrorRatio, $InsertSize)=(split /\t/, $_)[0,5,7,8,9,10];
		if ($SampleID!~/_/ && defined $date) {
			$SampleID.="_$date";
		}
		push @{$astat->{$SampleID}},($SoftClipRatio, $LowQualityRatio, $OffTargetRatio, $PosErrorRatio, $InsertSize);
	}
	
	close (I) ;
}

#Sample  Total   Error   Primer  PrimerRatio
#KYCT102B2       9014093 969137  8044956 89.25%
#KYCT104B2       4816349 1704610 3111739 64.61%
#KYCT105B2       7475489 838684  6636805 88.78%
sub get_primer_stat {#
	my ($astat, $date, $file) = @_;
	open (I,"$file") or die $!;
	while (<I>) {
		chomp;
		next if(/^$/ || /Sample/);
		my ($SampleID,$PrimeRratio, $PrimerReads)=(split /\t/, $_)[0,4,3];
		if ($SampleID!~/_/ && defined $date) {
			$SampleID.="_$date";
		}
		push @{$astat->{$SampleID}},($PrimeRratio, $PrimerReads);
	}
	
	close (I) ;
}

sub get_barcode_stat {#
	my ($astat, $date, $file) = @_;
	open (I,"$file") or die $!;
	while (<I>) {
		chomp;
		next if(/^$/ || /Sample/);
		my ($SampleID,$BarErrorRatio)=(split /\t/, $_)[0,3];
		if ($SampleID!~/_/ && defined $date) {
			$SampleID.="_$date";
		}
		push @{$astat->{$SampleID}},($BarErrorRatio);
	}
	
	close (I) ;
}

sub get_reads_and_uniq_stat {#
	my ($astat, $date, $file) = @_;
	open (I,"$file") or die $!;
	while (<I>) {
		chomp;
		next if(/^$/ || /Sample/);
		my ($SampleID, $TotalReads,$UniqRatio)=(split /\t/, $_)[0,1,3];
		if ($SampleID!~/_/ && defined $date) {
			$SampleID.="_$date";
		}
		push @{$astat->{$SampleID}},($TotalReads,$UniqRatio);
	}
	
	close (I) ;
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
  -id <dir,>   indirs, forced
  -d  <str>    date string like 161206, optinal
  -p  <str>    prefix of pre-existing statistic file, optional

  -k  <str>    Key of output file, forced
  -od <dir>    Dir of output file, default ./
  -h         Help

USAGE
	print $usage;
	exit;
}
