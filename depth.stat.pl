#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
require "$Bin/mylib.pl";
my $BEGIN_TIME=time();
my $version="1.0.0";
#######################################################################################

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($bamfile,$outdir,$outfile,$testchr,$testpos);
GetOptions(
				"help|?" =>\&USAGE,
				"i:s"=>\$bamfile,
				"d:s"=>\$outdir,
				"o:s"=>\$outfile,
				) or &USAGE;
&USAGE unless defined($bamfile and $outfile);

$testchr||="chr11";
$testpos||=5247905;
$outdir||="./";
`mkdir $outdir` unless (-d $outdir);
$outdir=&AbsolutePath("dir",$outdir);
open(I, "samtools view $bamfile|") or die $!;
open(OUT,">$outfile");
my %hash;
my %count;
&SHOW_TIME("reading");
while(<I>){
	chomp; my$lines = $_;
	my ($id,$flag,$chr,$pos,$mapQ,$matchInfo,@line) = split(/\t/,$lines);
	#my($rnum, $is_proper_pair, $is_reverse, $is_unmap, $is_munmap, $is_supplementary)=&explain_bam_flag($flag);
	my($rnum, $is_unmap)=&explain_bam_flag($flag);
	if($is_unmap == 1){  # attention: mpileup $7 ne "="
		#print "$id,$flag,$chr,$pos,unmap\n";	
	}else{
		&sep_and_count_array($matchInfo,$chr,$pos,$id,$rnum);
		#&sep_and_count($matchInfo,$chr,$pos,"add");
		#my($start1,$end1) = &sep_and_count_n($matchInfo,$pos,$pos);
		#&check_info_print("checkregion",$chr,$start1,$end1,$id,"add");
#		(exists $hash{$id}{"READ".$rnum})?
#			($hash{$id}{"READ".$rnum}=$hash{$id}{"READ".$rnum}.",$chr:$start1:$end1"):
#			($hash{$id}{"READ".$rnum}="$chr:$start1:$end1");
	}
}
&SHOW_TIME("counting");
foreach my$id(keys %hash){
	foreach my $chr(keys %{$hash{$id}}){
		if( exists $hash{$id}{$chr}{"READ1"} and exists $hash{$id}{$chr}{"READ2"} ){ 
			my @read1 = @{$hash{$id}{$chr}{"READ1"}};
			my @read2 = @{$hash{$id}{$chr}{"READ2"}};
			for(my $i=0; $i<@read1; $i++){	
				my ($start1,$end1) = @{$read1[$i]};
				for (my $j=0; $j<@read2; $j++){
					my ($start2,$end2) = @{$read2[$j]};
					if( $start1 <= $start2 ){
						if( $start2<=$end1 and $end1<=$end2 ){
							&countDel( $chr ,$start2 ,$end1-$start2+1  );
							#&check_info_print("checkregion",$chr1,$start2,$end1,$k1,"del");
						}elsif( $start2<=$end1 and $end1>$end2 ){
							&countDel( $chr ,$start2 ,$end2-$start2+1  );
							#&check_info_print("checkregion",$chr1,$start2,$end2,$k1,"del");
						}
					}else{
						if( $start1<=$end2 and $end2<=$end1 ){
							&countDel( $chr ,$start1 ,$end2-$start1+1  );
							#&check_info_print("checkregion",$chr1,$start1,$end2,$k1,"del");
						}elsif( $start1<=$end2 and $end2>$end1 ){
							&countDel( $chr ,$start1 ,$end1-$start1+1  );
							#&check_info_print("checkregion",$chr1,$start1,$end1,$k1,"del");
						}
					}
				}
			}
		}
	}
}
close I;
&SHOW_TIME("printing:");
foreach my$key(keys %count){
	my $hash2 = $count{$key};
	foreach my $key2 (sort{ $a <=> $b } keys %$hash2){
		print OUT $key,"\t$key2\t$hash2->{$key2}\n";
	}
}
&SHOW_TIME("ending");
#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################

# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------
sub SHOW_TIME {
	#显示当时时间函数，参数内容是时间前得提示信息，为字符串
	my ($str)=@_;
	my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
	my $temp=sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
	print "$str:\t[".$temp."]\n";
}
sub explain_bam_flag{
	my ($flag)=@_;
	my $flag_bin=sprintf("%b", $flag);
	my @flag_bin = split //, $flag_bin;
	my $is_read1 = $flag_bin[-7];
	#my $is_read2 = @flag_bin>=8? $flag_bin[-8]: 0;
	#my $is_supplementary = @flag_bin>=12? $flag_bin[-12]: 0;
	#my $is_proper_pair = $flag_bin[-2];
	#my $is_reverse = $flag_bin[-5];
	my $is_unmap = $flag_bin[-3];
	my $is_munmap = $flag_bin[-4];
	#my $dup = @flag_bin>=11? $flag_bin[-11]: 0;
	my $rnum = $is_read1==1? 1: 2;
	#return($rnum, $is_proper_pair, $is_reverse, $is_unmap, $is_munmap, $is_supplementary);
	return($rnum,$is_unmap);
}
sub sep_and_count_array{
	my($match_info,$chr,$pos,$id,$rnum) = @_;
	#my @match_n = split(/[MDIS]/,$match_info);
	#my @match_str = split(/\d+/,$match_info);
	my @ucigar = split //, $match_info;
	my @match_n;
	my @match_str;
	my $cigar="";
	for(my $i=0; $i<@ucigar; $i++){
		if($ucigar[$i] eq "M" || $ucigar[$i] eq "I" || $ucigar[$i] eq "D" || $ucigar[$i] eq "S"){
			push @match_str, $ucigar[$i];
			push @match_n, $cigar;
			$cigar = "";
		}else{
			$cigar.=$ucigar[$i];
		}
	}
	my $start = $pos;
	for (my $index=0; $index<=$#match_n;$index++){
		if($match_str[$index] eq "S" && $index!=0 && $index!=$#match_n){
			print "Wrong: $match_info 'S' is in middle!\n";
			die;
		}
		my $len = $match_n[$index];
		if($match_str[$index] eq "M" or $match_str[$index] eq "D"){
			&countAdd($chr,$pos,$len );
			$pos = $pos+$len;
		}
	}	
	push @{$hash{$id}{$chr}{"READ".$rnum}},[$start, $pos-1];
}
sub sep_and_count{
	my($match_info,$chr,$pos,$type) = @_;
	my($len,$match_info_new);
	if( $match_info =~ m/^([0-9]+)M(.*)$/ ){
		$len = $1; $match_info_new=$2;
	}elsif( $match_info =~ m/^[0-9]+[IS]([0-9]+)M(.*)$/){
		$len = $1; $match_info_new=$2;
	}
	&count($chr,$pos,$len,$type );
	if( defined $match_info_new ){}else{print "$match_info,$chr,$pos,$type\twhat?\n"; }
	if( length($match_info_new)>1 ){ &sep_and_count_next($match_info_new,$chr,$pos+$len,$type ); }
}
sub sep_and_count_next{
	my($match_info,$chr,$pos,$type) = @_;
	my($len,$match_info_new);
	if( $match_info =~ m/^([0-9]+)M(.*)$/){
		$len = $1; $match_info_new=$2;
		&count($chr,$pos,$len,$type );
		if( length($match_info_new)>1 ){ &sep_and_count_next($match_info_new,$chr,$pos+$len,$type ); }
	}elsif( $match_info =~ m/^([0-9]+)D(.*)$/ ){
		$len = $1; $match_info_new=$2;
		&count($chr,$pos,$len,$type ); ######
		if( length($match_info_new)>1 ){ &sep_and_count_next($match_info_new,$chr,$pos+$len,$type ); }
	}elsif( $match_info =~ m/^([0-9]+)S(.*)$/ ){
		$len = $1; $match_info_new=$2;
		if( length($match_info_new)>1 ){ &sep_and_count_next($match_info_new,$chr,$pos+$len,$type ); }
	}elsif( $match_info =~ m/^([0-9]+)I(.*)$/ ){
		$len = $1; $match_info_new=$2;
		if( length($match_info_new)>1 ){ &sep_and_count_next($match_info_new,$chr,$pos,$type ); }
	}
}
sub sep_and_count_n{
	my($match_info,$n1,$n2) = @_;
	my($len,$match_info_new);
	if( $match_info =~ m/^([0-9]+)M(.*)$/){
		$len = $1; $match_info_new=$2;
	}elsif( $match_info =~ m/^[0-9]+[IS]([0-9]+)M(.*)$/){
		$len = $1; $match_info_new=$2;
	}
	if( length($match_info_new)>1 ){ &sep_and_count_n_next($match_info_new,$n1,$n2+$len-1) }
	else{ return($n1,$n2+$len-1); }
}
sub sep_and_count_n_next{
	my($match_info,$n1,$n2) = @_;
	my($len,$match_info_new);
	if( $match_info =~ m/^([0-9]+)[MD](.*)$/){
		$len = $1; $match_info_new=$2;
		if( length($match_info_new)>1 ){ &sep_and_count_n_next($match_info_new,$n1,$n2+$len) }
		else{ return($n1,$n2+$len); }
	}elsif( $match_info =~ m/^([0-9]+)S(.*)$/ ){
		$len = $1; $match_info_new=$2;
		if( length($match_info_new)>1 ){ &sep_and_count_n_next($match_info_new,$n1+$len,$n2+$len) }
		else{ return($n1,$n2); }
	}elsif( $match_info =~ m/^([0-9]+)I(.*)$/ ){
		$len = $1; $match_info_new=$2;
		if( length($match_info_new)>1 ){ &sep_and_count_n_next($match_info_new,$n1,$n2) }
		else{ return($n1,$n2); }
	}
}

sub check_info_print{
	my($checktype,$chr,$pos,$len,$type,$j)=@_;
	if($checktype eq "checkpos"){
		if( $chr eq $testchr and $j == $testpos ){
			print "$checktype,$chr,$pos,$len,$type,$j\n";
		}
	}elsif($checktype eq "checkregion"){
		if( $chr eq $testchr and $testpos >=$pos and $testpos <= $len ){
			print "$checktype,$chr,$pos,$len,$type,$j\n";
		}
	}
}
sub countAdd{
	my($chr,$pos,$len)=@_;
	foreach my$i(1..$len){
		#&check_info_print("checkpos",$chr,$pos,$len,$type,$j);
		$count{$chr}{$pos + $i -1}++;
	}
}
sub countDel{
	my($chr,$pos,$len)=@_;
	foreach my$i(1..$len){
		#if(not exists$count{$chr}{$j}){print "$chr,$pos,$len,$type,$j\tnotadd\n";}
		$count{$chr}{$pos + $i -1}--;
	}
}

sub USAGE {#
	my $usage=<<"USAGE";
Program: $0 (count base depth without overlap!)
Version: $version
Contact: Wu Guizhi<guizhi.wu\@genetalks.com> 

Usage:
  Options:
  -i  <file>    Input bam file, forced
  -d  <outdir>  Outdir, default ./
  -o  <outfile>	Output file, forced
  -h		 Help

USAGE
	print $usage;
	exit;
}

