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
$|++;
#######################################################################################

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($fdata_sheet, $indir, $outdir, $step, $only);
my ($min_CV, $min_depth, $threads);
my $sample;
my ($SingleEnd, $secBar_len);
my ($mismatch_bar1, $mismatch_bar2);
my $configID;
my ($distance);
my $ref = "/data/bioit/biodata/duyp/bin/hg19/hg19.fasta";
my $fpanel = "$Bin/config/panel/panel.config";
my $fbarcode = "$Bin/config/barcode.config";
GetOptions(
				"help|?" =>\&USAGE,
				"i:s"=>\$fdata_sheet,
				"ip:s"=>\$fpanel,
				"ib:s"=>\$fbarcode,
				"ir:s"=>\$ref,
				"cid:s"=>\$configID,
				"s:s"=>\$sample,
				"indir:s"=>\$indir,
				"SingleEnd:s"=>\$SingleEnd,
				"threads:s"=>\$threads,
				"secBar_len:s"=>\$secBar_len,
				"min_depth:s"=>\$min_depth,
				"min_CV:s"=>\$min_CV,
				"mismatch_bar1:s"=>\$mismatch_bar1,
				"mismatch_bar2:s"=>\$mismatch_bar2,
				"step:s"=>\$step,
				"only:s"=>\$only,
				"od:s"=>\$outdir,
				"distance:i"=>\$distance,
				) or &USAGE;
&USAGE unless ($fdata_sheet and $indir and $sample);
$outdir||="./";
`mkdir $outdir`	unless (-d $outdir);
$outdir=AbsolutePath("dir",$outdir);
$step = defined $step? $step: 1;
$min_CV = defined $min_CV? $min_CV: 0.65;
$min_depth = defined $min_depth? $min_depth: 4;
$threads = defined $threads? $threads: 6;
$secBar_len = defined $secBar_len? $secBar_len: 3;
$mismatch_bar1 = defined $mismatch_bar1? $mismatch_bar1: 1;
$mismatch_bar2 = defined $mismatch_bar2? $mismatch_bar2: 1;
$distance = defined $distance? $distance: 2000;

my %data_config;
my $sh;
my $log_dir = "$outdir/log/";
mkdir $log_dir unless (-d $log_dir);
open ($sh, ">$log_dir/$sample.sh_".time()) or die $!;

#------------------------------------------------------------------
# read config
#------------------------------------------------------------------
&read_data_sheet($fdata_sheet, $fpanel, \%data_config, $configID);

#==================================================================
# fastq uniq
#==================================================================
my $dir_fastq_uniq = "$outdir/01.fastq_uniq";
if ($step == 1) {
	mkdir $dir_fastq_uniq unless (-d $dir_fastq_uniq);
	&fastq_uniq($indir,\%data_config, $dir_fastq_uniq);
	if (!defined $only) {
		$step ++;
	}
}
my $step1_time = time();
print "\nStep1 elapsed time : ",$step1_time-$BEGIN_TIME,"s\n";

#------------------------------------------------------------------
# fastq barcode process
#------------------------------------------------------------------
my $dir_fastq_barcode = "$outdir/02.barcode_process";
if ($step == 2) {
	mkdir $dir_fastq_barcode unless (-d $dir_fastq_barcode);
	&barcode_process($indir, \%data_config, $dir_fastq_barcode, $dir_fastq_uniq);
	if (!defined $only) {
		$step ++;
	}
}
my $step2_time = time();
print "\nStep2 elapsed time : ",$step2_time-$step1_time,"s\n";

#------------------------------------------------------------------
# data primer
#------------------------------------------------------------------
my $dir_data_process = "$outdir/03.data_process";
if ($step == 3) {
	mkdir $dir_data_process unless (-d $dir_data_process);
	&data_process(\%data_config, $dir_data_process);
	if (!defined $only) {
		$step ++;
	}
}
my $step3_time = time();
print "\nStep3 elapsed time : ",$step3_time-$step2_time,"s\n";


#------------------------------------------------------------------
# consensus maker and QC
#------------------------------------------------------------------
my $dir_consensus = "$outdir/04.consensus_make";
if ($step==4) {
	mkdir $dir_consensus unless (-d $dir_consensus);
	&consensus_make(\%data_config, $dir_consensus);
	if (!defined $only) {
		$step ++;
	}
}
my $step4_time = time();
print "\nStep4 elapsed time : ",$step4_time-$step3_time,"s\n";

#------------------------------------------------------------------
# variant detect
#------------------------------------------------------------------
my $dir_variant = "$outdir/05.variant_detect";
if ($step==5) {
	mkdir $dir_variant unless (-d $dir_variant);
	&consensus_variant_detect(\%data_config, $dir_variant, $dir_consensus);
	if (!defined $only) {
		$step ++;
	}
}
my $step5_time = time();
print "\nStep5 elapsed time : ",$step5_time-$step4_time,"s\n";

close ($sh);

#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################

# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------

sub consensus_variant_detect {#
	my ($adc, $dir, $dir_consensus)=@_;
	my $cmd;

	my $fq1 = "$dir_consensus/$sample.consensus.final.R1.fastq";
	my $fq2 = "$dir_consensus/$sample.consensus.final.R2.fastq";

	## bwa and IndelRealigner
	my $RG="\'\@RG\\tID:$sample\\tPL:ILLUMINA\\tSM:$sample\\tDS:ref=hg19,pfx=$ref\'";
	&Run("$Bin/third_party/bwa mem -C -Y -R $RG -t $threads $ref $fq1 $fq2 | $Bin/third_party/samtools view -bS - |$Bin/third_party/samtools sort - -o $dir/$sample.sorted.bam", $sh);
	&Run("$Bin/third_party/samtools index $dir/$sample.sorted.bam", $sh);
	if (!-f "$dir/$sample.target.bed") {
		&Run("ln -s $adc->{$sample}{bedFile} $dir/$sample.target.bed", $sh);
	}
	
	## stat
	$cmd = "perl $Bin/hotspot.stat.pl -i1 $dir/$sample.sorted.bam -i2 $adc->{$sample}{ftarget}";
	if (exists $adc->{$sample}{"fusion_pos"}) {
		$cmd .=" -i3 $adc->{$sample}{fusion_pos}";
	}
	$cmd .=" -k $sample -od $dir/QC/";
	&Run($cmd, $sh); 
	if(-e "$dir/QC/$sample.fusion.bed" && !-z "$dir/QC/$sample.fusion.bed"){
		$adc->{$sample}{fusion_bed}="$dir/QC/$sample.fusion.bed";
	}

	my $dir_variant_sample = "$dir/variants/";
	mkdir $dir_variant_sample unless(-d $dir_variant_sample);

	## detect variants 
	&Run("$Bin/variant/realign_indel -i $dir/$sample.sorted.bam -r $ref -o $dir/$sample.realign.bam", $sh);
	&Run("$Bin/variant/call_snp -i $dir/$sample.realign.bam  -g $ref -rp $adc->{$sample}->{fregion_primer} -o $dir_variant_sample/$sample.vcf", $sh);
		
	if (exists $adc->{$sample}{"ftarget"}) {
		&Run("less $dir_variant_sample/$sample.vcf|awk '\$5!=\"<*>\"'|sed 's/,<\\*>//'|less >$dir_variant_sample/$sample\_variant.vcf", $sh);
		my $dir_annovar_db = "$outdir/annovar_db/";
		`mkdir $dir_annovar_db` unless(-d $dir_annovar_db);
		`ln -s $adc->{$sample}{ftarget} "$dir_annovar_db/hg19_variant_target.txt"` unless(-e "$dir_annovar_db/hg19_variant_target.txt");
		&Run("perl $Bin/third_party/annovar/table_annovar.pl $dir_variant_sample/$sample\_variant.vcf $dir_annovar_db -buildver hg19 -out $dir_variant_sample/$sample -remove -protocol variant_target -operation f -nastring . -vcfinput --otherinfo >$dir_variant_sample/$sample.annovar.log 2>&1", $sh);
		&Run("perl $Bin/variant/annovar_result_simplify_convert.pl --thalassemia -i1 $dir_variant_sample/$sample.hg19_multianno.txt -i2 $adc->{$sample}{ftarget} -s $sample -od $dir_variant_sample/", $sh);
	}

	## fusion
	&Run("$Bin/variant/fusion_detect -i $dir/$sample.sorted.bam -k $sample -od $dir_variant_sample", $sh);
}



sub consensus_make{
	my ($adc, $dir)=@_;
	my $fq1 = (-e "$dir_data_process/$sample.hasPrimer.R1.fastq")? "$dir_data_process/$sample.hasPrimer.R1.fastq" :"$dir_data_process/$sample.hasPrimer.R1.fastq.gz";
	my $fq2 = (-e "$dir_data_process/$sample.hasPrimer.R2.fastq")? "$dir_data_process/$sample.hasPrimer.R2.fastq" : "$dir_data_process/$sample.hasPrimer.R2.fastq.gz";
	my $filter_flag = defined $SingleEnd? "LTP": "SLTP";
	&Run("$Bin/third_party/bwa mem -Y -C -t $threads $ref $fq1 $fq2 | $Bin/cfbest filter -k $dir/$sample -$filter_flag -t $adc->{$sample}->{bedFile} -p $adc->{$sample}->{primerFile} | $Bin/cfbest dedup -s $sample -o $dir -t $threads -d $min_depth -c $min_CV -N 2 -m $mismatch_bar1 -M $mismatch_bar2", $sh);
	if (!defined $SingleEnd) {
		&Run("perl $Bin/consensus_filter.pl -fq1 $dir/$sample.consensus.R1.fastq -fq2 $dir/$sample.consensus.R2.fastq -k $sample -od $dir -distance $distance", $sh);
	}else{
		&Run("ln -s $dir/$sample.consensus.R1.fastq  $dir/$sample.consensus.final.R1.fastq", $sh);
		&Run("ln -s $dir/$sample.consensus.R2.fastq  $dir/$sample.consensus.final.R2.fastq", $sh);
	}
	&Run("$Bin/cfbest cstat -1 $dir/$sample.consensus.final.R1.fastq -p $adc->{$sample}->{primerFile} -o $dir/$sample", $sh);
}

sub data_process{
	my ($adc, $dir)=@_;

	### primer filter
	my	$fq1 = "$dir_fastq_barcode/$sample.clean.R1.fastq";
	my	$fq2 = "$dir_fastq_barcode/$sample.clean.R2.fastq";
	if(!-e $fq1){
		$fq1 = "$dir_fastq_barcode/$sample.corbar.R1.fastq";
		$fq2 = "$dir_fastq_barcode/$sample.corbar.R2.fastq";
	}
	&Run("$Bin/cfbest detprim -1 $fq1 -2 $fq2 -p $adc->{$sample}->{primerFile} -o $dir/$sample ", $sh);
}

sub barcode_process{
	my ($aindir, $adc, $dir, $dir_uniq)=@_;
	my $fq1 = "$dir_uniq/$sample.uniq.R1.fastq";
	my $fq2 = "$dir_uniq/$sample.uniq.R2.fastq";
	if ($adc->{$sample}->{barGroup}=~/NNNNN/){
		my $blen = length($adc->{$sample}->{barGroup});
		&Run("$Bin/cfbest ctrim -1 $fq1 -2 $fq2 -o $dir/$sample", $sh);
	}elsif ($adc->{$sample}->{barGroup}=~/:/){
		my ($bar1, $bar2)=split /:/, $adc->{$sample}->{barGroup};
		&Run("$Bin/cfbest detbar -1 $fq1 -2 $fq2 -t 4 -b $fbarcode -g $bar1 -G $bar2  -o $dir/$sample", $sh);
		&Run("$Bin/cfbest ctrim -1 $dir/$sample.corbar.R1.fastq -2 $dir/$sample.corbar.R2.fastq -o $dir/$sample", $sh);
	}elsif (defined $SingleEnd) {
		&Run("$Bin/cfbest detbar -1 $fq1 -2 $fq2 -B -L $secBar_len -b $fbarcode -g $adc->{$sample}->{barGroup} -o $dir/$sample", $sh);
	}else{
		&Run("$Bin/cfbest detbar -1 $fq1 -2 $fq2 -b $fbarcode -g $adc->{$sample}->{barGroup} -o $dir/$sample", $sh);
		&Run("$Bin/cfbest ctrim -1 $dir/$sample.corbar.R1.fastq -2 $dir/$sample.corbar.R2.fastq -o $dir/$sample", $sh);
	}
}


sub fastq_uniq {#
	my ($indir, $adc, $dir,$dir_barcode) = @_;
	my $fq1;
	my $fq2;
	if($adc->{$sample}->{barGroup}=~/NNNNN/ || $adc->{$sample}->{barGroup}=~/:/){#planC or DoubleBarcode  
		$fq1 = "$outdir/00.call_data/$sample.R1.fastq";
		$fq2 = "$outdir/00.call_data/$sample.R2.fastq";
		&Run("$Bin/cfbest fquniq -1 $fq1 -2 $fq2 -o $dir/$sample", $sh);
	}else{  
		$fq1 = -e "$indir/$sample.R1.fastq.gz"? "$indir/$sample.R1.fastq.gz": "$indir/$sample\_L1*R1.fastq.gz";
		$fq2 = -e "$indir/$sample.R2.fastq.gz"? "$indir/$sample.R2.fastq.gz" : "$indir/$sample\_L1*R2.fastq.gz";
		&Run("$Bin/cfbest fquniq -1 $fq1 -2 $fq2 -o $dir/$sample -C", $sh);
	}
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
  -i        <file>    Input data sheet file(separate by tab), forced
  -ip       <file>    Input panel config file, [$fpanel]
  -ib       <file>    Input barcode config file, [$fbarcode]
  -ir       <file>    Input ref genome file, [$ref]
  -cid      <str>     config ID in config file, forced
  -s        <str>     sample name, forced
  -indir    <dir>     rawdata dir, "Undetermined_*_001.fastq.gz", or "L1/*.R1.fastq.gz"
  --SingleEnd         SingleEnd, read2 as second barcode, optinal

  -secBar_len  <int>     length of r2 seq as barcode when specified --SingleEnd, [3] 
  -min_depth   <int>     min group depth, [4]
  -min_CV     <float>	 min CV, [0.65]
  -mismatch_bar1 <int>   mismatch of the first barcode, [1]
  -mismatch_bar2 <int>   mismatch of the second barcode, [1]

  -threads  <int>    threads num, [6]
  -distance  <int>   maximum cutoff value of read1 and read2 alignment position in same chr(not when --SingleEnd),[2000]
  -step     <int>    step, [1]
	1: data uniq
	2: barcode process
	3: data process
	4: consensus maker
	5: variants detect
  --only            only run the current step, optional

  -od <dir>    outdir of output file, default ./
  -h         Help

USAGE
	print $usage;
	exit;
}
