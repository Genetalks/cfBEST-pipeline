## Getting started

	git clone https://github.com/Genetalks/cfBest-pipeline.git
	perl cfBest_pipeline.pl -i demo.config -ir hg19.fasta –ip panel.config -ib barcode.config -cid pID -indir dir/ 
	-od outdir/ -threads 10 -min_depth 4

## Introduction

	CfBest_pipeline is a software package for detecting variants accurately and typing maternal and fetal genotypes 
	for cfBEST which is a barcode-enabled cell-free DNA allelic molecule counting system. It consists of six 
	procedure: data preprocessing(include data merging, barcode trimming, adapter trimming and primer recognizing), 
	Mapping, consensus sequence calling, allele counting and genotyping. The accuracy is solidly guaranteed though 
	adapter trimming precisely by recognizing the overlap between pair ends, realigning indel closed to read end, 
	counting allele depth by specific primers in specific sites, and a precise genotyping model.

## Usage
	-i        	<file>    	Input data sheet file(separate by tab), forced
  	-ir       	<file>    	Input ref genome file(fasta), forced
  	-ip       	<file>    	Input panel config file, [config/panel.config]
  	-ib       	<file>    	Input barcode config file, [config/barcode.config]
  	-cid      	<str>     	config ID in panel config file, optional
  	-indir    	<dir>      	raw fastq dir, "*.R1.fastq.gz"
  	--SingleEnd         		SingleEnd, read2 as second barcode, optinal
  	--Double_Barcode    		Double barcode, optinal
  	--Abnormal_Recognize   		abnormal sample recognize
	
 	-min_depth     	<int>    	depth threshold to filter low depth reads group.[4]
 	-min_CV       	<float>  	min depth CV=(max-sec)/max when call consensus sequence, [0.65]
  	-mdepth_snp    	<int>    	min depth of snp to infer cffDNA, [100]
  	-msize         	<int>    	intersize threshold of consensus sequence,[2000]
  	-secBar_len    	<int>    	length of r2 seq as barcode when specified --SingleEnd, [9]
  	-mismatch_bar1 	<int>    	mismatch of the first barcode, [1]
  	-mismatch_bar2 	<int>    	mismatch of the second barcode, [1]
  	-parallel      	<int>     	jobs num of parallel, [4]
  	-threads       	<int>    	threads num, [6]
  	-step          	<int>    	step, [1]
        				  	0: call data from Undetermined
        					1: data merge
        					2: barcode process
        					3: primer process
        					4: consensus maker
        					5: variants detect
						6: genotyping
        					7: merge statistic and abnormal sample recognize
  	-od 		<dir>    	outdir of output file, default ./
  	-h         			Help

#### Format of data sheet file: 

Col|Name|description
:-:|:---|:---------
1|SampleID|Sample ID
2|Index|Sample Index 
3|BarcodeGroup|Specified barcode group IDs used in the sample, detailed barcode sequence is listed in barcode config file
4|SampleType|Sample type, include gravida-plasma-thalassemia and standards-thalassemia and other. It will infer fetal DNA fraction and genotype for fetal and maternal allele if choose gravida-plasma-thalassemia, and only genotype based on supplied fetal DNA fraction if choose standards-thalassemia, and neither infer fetal fetal DNA fraction nor genotype but only detect mutations when choose other.
5|Description|Sample description
6-8|Genotype|Fetal maternal and paternal genotype, optional, it will judge True or False when supplied. N represent negative, N/N, xx/N, xx/xx represent three heterozygosis state

	Example:
	S1    CTTCGTT B6,B7,B8        gravida-plasma-thalassemia      sample descriptions     -28/N   -28/N   CD41-42/N

#### Format of panel config file:

Col|Name|Description
:-:|:---|:---------
2|TargetBedFile|File, 4 column: chr start end primerID
3|PrimerInfoFile|File, 5 column: chr pos strand primerID primerSeq
4|SpecificPrimerFile|File, 5 column: chr start end primerID distance
5|FusionFile|File, 1) fusion region, 4 column: “Region” chr start end; 2)fusion breakpoints, 6 column: “Break” fusionID chr1 pos1 chr2 pos2
6|TargetDetectFile|File, 6 column: chr start end ref alt info(ID=xxx;GENE=xxx; STRAND=xxx; CDS=xxx;AA=xxx)



#### Format of panel config file:

Col	|Description
:-:	|:----------
1	|Barcode group ID
2	|Barcode sequence

