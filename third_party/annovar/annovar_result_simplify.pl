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
my ($fanno, $fvcf, $fkey,$outdir);
GetOptions(
				"help|?" =>\&USAGE,
				"i1:s"=>\$fanno,
				"i2:s"=>\$fvcf,
				"k:s"=>\$fkey,
				"od:s"=>\$outdir,
				) or &USAGE;
&USAGE unless ($fanno and $fvcf and $fkey);

$outdir||="./";
`mkdir $outdir`	unless (-d $outdir);
$outdir=AbsolutePath("dir",$outdir);

my @sample;
open (V, $fvcf) or die $!;
while (<V>) {
	chomp;
	next if(/^$/);
	if ($_=~/^#CHROM/) {
		my @unit = split /\s+/, $_;
		@sample = @unit[9..$#unit];
		last;

	}
}
close (V);


open (O, ">$outdir/$fkey.final.total.txt") or die $!;
open (F, ">$outdir/$fkey.final.filter.txt") or die $!;
open (A, $fanno) or die $!;
while (<A>) {
	chomp;
	next if(/^$/);
	my @info = split /\t/, $_;
	my @retain = @info[0..23,43];

	if ($_=~/Start/) { ##head line
		push @retain, @sample;
	}else{
		for (my $i=56; $i<@info ; $i++) {  ## 57 is the colmun of the first sample 
			if ($info[$i] eq "./.") {
				push @retain,".";
			}else{
				my ($ref_dep, $alt_dep)=$info[$i]=~/:(\d+),(\d+):/;
				my ($total_dep, $freq);
				if (!defined $ref_dep ) {
					($total_dep) = $info[$i]=~/\d\/\d:(\d+)/;
					$ref_dep = "?";
					$alt_dep = "?";
					$freq = "?";
				}else{
					$total_dep = ($ref_dep+$alt_dep);
					$freq = sprintf("%.2f",$alt_dep/$total_dep*100)."%";
				}

				push @retain, $alt_dep.",".$total_dep.",".$freq;
			}
		}
	}
	print O join("\t", @retain),"\n";

	## filter
	next if (($info[5] ne "Func.refGene" && $info[5] ne "exonic") ||  $info[8]=~/nonframeshift/ || $info[8] eq "synonymous SNV"); 
	print F join("\t", @retain),"\n";
}

close (O);
close (F);
close (A);








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
  -i1  <file>   annovar multianno.txt file, forced
  -i2  <file>   vcf file, forced
  -k   <str>    key of Output file, forced

  -od  <dir>    outdir of output file, default ./
  -h         Help

USAGE
	print $usage;
	exit;
}
