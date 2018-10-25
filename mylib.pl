use FindBin qw($Bin $Script); 
sub Run{
	my ($cmd, $log_fh, $nodie)=@_;

	if (defined $log_fh && $log_fh ne "") {
		print $log_fh $cmd,"\n";
	}
	print STDERR $cmd, "\n";
	my $ret = system($cmd);
	if (!defined $nodie && $ret) {
		die "Run $cmd failed!\n";
	}
}


sub read_data_sheet{
	my ($fsampsheet, $fconfig, $adata, $cid)=@_;
	open(C, $fconfig) or die $!;
	while (<C>){
		chomp;
		next if(/^$/ || /^#/);
		my ($id, $bed, $primerfile, $fregion_primer, $fusion_pos, $fvariant_target)=split /\t/, $_;
		my $cbin = dirname(&AbsolutePath("file", $fconfig));
		@{$config{$id}}=($cbin."/".$bed, $cbin."/".$primerfile, $cbin."/".$fregion_primer, $cbin."/".$fusion_pos, $cbin."/".$fvariant_target);
	}
	close(C);

	open (S, $fsampsheet) or die $!;
	while (<S>) {
		chomp;
		next if(/^$/ || /^#/ || /^SampleName/);
		my @unit = split /\t/, $_;

		## get configID and variant_before
		my ($configID, $variant_before, $sample, $index, $barInfo, @sample_info);
		if($unit[1]=~/^[ATCGN]+$/){## judge the column is index or not
			($sample, $index, $barInfo, @sample_info)=@unit;
		}elsif($unit[2]=~/^[ATCGN]+$/){
			($sample, $index, $barInfo, @sample_info)=@unit[1..$#unit];
			if(-e $unit[0]){
				$variant_before = $unit[0];
			}else{
				$configID = $unit[0];
			}
		}elsif($unit[3]=~/^[ATCGN]+$/){
			($configID, $variant_before, $sample, $index, $barInfo, @sample_info)=@unit;
		}else{
			print "Wrong samplesheet: $fsampsheet!\n";
			die;
		}
		if(defined $configID){
			$cid =$configID;
		}
		my ($bed, $primerfile, $fregion_primer, $fusion_pos, $ftarget, $fgenotype) = @{$config{$cid}};
		$adata->{$sample}->{"barGroup"}=$barInfo;
		$adata->{$sample}->{"bedFile"}=$bed;
        $adata->{$sample}->{"primerFile"}=$primerfile;
        $adata->{$sample}->{"fregion_primer"}=$fregion_primer;
		$adata->{$sample}->{"ftarget"} = $ftarget;
		$adata->{$sample}->{"fgenotype"} = $fgenotype if(defined $fgenotype);
		$adata->{$sample}->{"index"}=$index;
		@{$adata->{$sample}->{"sample_info"}}=@sample_info;
		if (defined $variant_before) {
			$adata->{$sample}->{"variant_before"}=$variant_before;
		}
		if ($fusion_pos ne "") {
			$adata->{$sample}->{"fusion_pos"}=$fusion_pos;
		}
	}
	close (S);
}

sub AbsolutePath
{		#获取指定目录或文件的决定路径
		my ($type,$input) = @_;

		my $return;

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
1;

