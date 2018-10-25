#!/share/software/software/R-3.0_install/R-3.0.1/bin/Rscript

library(getopt)
library(dplyr)
#+--------------------
# get options
#+--------------------
spec <- matrix(c(
	'help', 	'h', 	0, "logical", 	"help",
	'input', 	'i', 	1, "character",	"input sample.hotspots.result file, forced.",
	'cffdna',	'c',	1, "character",	"input cffdna of this sample, forced.",
	'outfile', 	'o', 	1, "character",	"outfile of sample.EM.xls with fetalDNA and genotyping info, forced."
), byrow = TRUE, ncol = 5)

opt <- getopt(spec)

#+--------------------
# check options
#+--------------------
if ( !is.null(opt$help) | is.null(opt$input) | is.null(opt$cffdna) | is.null(opt$outfile) ) {
	cat(getopt(spec, usage=TRUE))
	q(status=1)
}

#+--------------------
# some default options
#+--------------------

#+--------------------
# EM - settings
#+--------------------
infer_fetal_genotype_by_cffDNA <- function(Ni, bi, fetal_fraction=0.05, pii=c(0.7, 0.1, 0.1, 0.1)){
	G <- c("AAaa", "AAab", "ABaa", "ABab")
	mG <- c("BBbb", "BBab", "ABbb", "ABab") # the mirror of genotype space

	Ni <- as.numeric(Ni)
	bi <- as.numeric(bi)
	is_mirror <- ifelse(Ni-bi>bi, 0, 1);
	b <- ifelse(Ni-bi>bi, bi, Ni-bi)

	u <- c(0.0001, fetal_fraction/2, 0.5-fetal_fraction/2, 0.5)
	gammak <- sapply(1:length(u), FUN=function(k){ # k = 1..4
		pii[k] * dbinom(as.numeric(b), as.numeric(Ni), u[k])
	})

	i <- match(max(gammak), gammak)
	g <- ifelse(is_mirror, mG[i], G[i])
	p <- gammak/sum(gammak)
	c(g, p[i])
}
format_gt <- function(id, gt){
    formated_gt <- sapply(1:length(id), function(i){
	if (is.na(gt[i])){
	    return(NA)
	}
	switch(gt[i],
	  AAaa = sprintf("%s/%s|%s/%s", "N", "N", "N", "N"),
	  AAab = sprintf("%s/%s|%s/%s", "N","N",id[i],"N"),
	  ABaa = sprintf("%s/%s|%s/%s", id[i], "N", "N", "N"),
	  ABab = sprintf("%s/%s|%s/%s", id[i], "N", id[i], "N"),
	  ABbb = sprintf("%s/%s|%s/%s", id[i], "N", id[i], id[i]),
	  BBab = sprintf("%s/%s|%s/%s", id[i], id[i], id[i], "N"),
	  BBbb = sprintf("%s/%s|%s/%s", id[i], id[i], id[i], id[i])
	)        
    })
    formated_gt
}
format_single_gt <- function(id, gt){
    formated_gt <- sapply(1:length(id), function(i){
	if (is.na(gt[i])){
	    return(NA)
	}
	switch(tolower(gt[i]),
	  aa = sprintf("%s/%s", "N", "N"),
	  ab = sprintf("%s/%s", id[i], "N"),
	  bb = sprintf("%s/%s", id[i], id[i]),
	)        
    })
    formated_gt
}
calculate_aver_depth <- function(df, min_snp_depth = 10, total_col = 6, anno_col = 5){
    hotspot_aver_depth <- NA
    overall_aver_depth <- NA

    snp_idx <- grep("^rs", df[, anno_col], perl=T)
    hotspot <- !((1:nrow(df)) %in% snp_idx)
	if (length(hotspot)){
        hotspot_aver_depth <- mean(df[hotspot, total_col], na.rm = TRUE)
    }
    
    is_snp = !(( 1:nrow(df) ) %in% hotspot);
    depth_larger_than_threshold <- df[, total_col] >= min_snp_depth

    overall_depth_valid <- (is_snp & depth_larger_than_threshold) | (!is_snp)
    overall_aver_depth <- mean(df[overall_depth_valid, total_col], na.rm = T)

    list(hotspot_aver_depth = sprintf("%.0f", hotspot_aver_depth), overall_aver_depth = sprintf("%.0f", overall_aver_depth))
}
# regenotyping  WS/CS/QS  base on SEA's genotyping 

#(1) AAaa (SEA): skip
#(2) AAab (SEA): AAa, AAb(x), ABa, ABb, BBa(x), BBb
#(3) ABaa (SEA): Aaa, Aab, Abb(x), Baa(x), Bab, Bbb
#(4) ABab (SEA): Aa,  Ab,  Ba, Bb 
#(5) ABbb (SEA): A(x)    B(x)
#(6) BBab (SEA): a    b
#(7) BBbb (SEA): 深度为0, skip

get_possible_genotypes <- function(mcpy, fcpy){
   stopifnot(mcpy <= 2 & mcpy >= 0)
   stopifnot(fcpy <= 2 & fcpy >= 0)

   if (mcpy == 2 & fcpy == 2){
      return(c("AAaa", "AAab", "ABaa", "ABab", "ABbb", "BBab", "BBbb"))
   } 

   if (mcpy == 2 & fcpy == 1){
      return(c("AAa", "ABa", "ABb", "BBb"))
   }
   
   if (mcpy == 1 & fcpy == 2){
      return(c("Aaa", "Aab", "Bab", "Bbb"))
   }

   if (mcpy == 1 & fcpy == 1){
      return(c("Aa", "Ab", "Ba", "Bb"))
   }

   if (mcpy == 1 & fcpy == 0){
      return(c("A", "B"))
   }
   if (mcpy == 0 & fcpy == 1){
      return(c("a", "b"))
   }

   return(c())
}

get_maternal_genotype <- function(mcpy, fcpy, mutid){
   stopifnot(mcpy <= 2 & mcpy >= 0)
   stopifnot(fcpy <= 2 & fcpy >= 0)

   if (mcpy == 2 & fcpy == 2){
      return(c(rep("N/N", 2), rep(sprintf("%s/N", mutid), 3), rep(sprintf("%s/%s", mutid, mutid), 2)))
   } 

   if (mcpy == 2 & fcpy == 1){
      return(c("N/N", rep(sprintf("%s/N", mutid), 2), sprintf("%s/%s", mutid, mutid)))
   }
   
   if (mcpy == 1 & fcpy == 2){
      return(rep(c("N", mutid), each = 2))
   }

   if (mcpy == 1 & fcpy == 1){
      return(rep(c("N", mutid), each = 2))
   }

   if (mcpy == 1 & fcpy == 0){
      return(c("N", mutid))
   }
   if (mcpy == 0 & fcpy == 1){
      return(c(NA, NA))
   }

   return(c())
}

get_fetal_genotype <- function(mcpy, fcpy, mutid){
   stopifnot(mcpy <= 2 & mcpy >= 0)
   stopifnot(fcpy <= 2 & fcpy >= 0)

   if (mcpy == 2 & fcpy == 2){
      return(c(rep(c("N/N", sprintf("%s/N", mutid)), 2), sprintf("%s/%s", mutid, murid), sprintf("%s/N", mutid), sprintf("%s/%s", mutid, mutid)))
   } 

   if (mcpy == 2 & fcpy == 1){
      return(rep(c("N", mutid), each=2))
   }
   
   if (mcpy == 1 & fcpy == 2){
      return(c("N/N", rep(sprintf("%s/N", mutid), 2), sprintf("%s/%s", mutid, mutid)))
   }

   if (mcpy == 1 & fcpy == 1){
      return(c("N", mutid, "N", mutid))
   }

   if (mcpy == 1 & fcpy == 0){
      return(c(NA, NA))
   }
   if (mcpy == 0 & fcpy == 1){
      return(c("N", mutid))
   }

   return(c())
}

get_maternal_mutation_allele_count <- function(mcpy, fcpy){
   stopifnot(mcpy <= 2 & mcpy >= 0)
   stopifnot(fcpy <= 2 & fcpy >= 0)

   if (mcpy == 2 & fcpy == 2){
      return(c(0, 0, 1, 1, 1, 2, 2))
   } 

   if (mcpy == 2 & fcpy == 1){
      return(c(0, 1, 1, 2))
   }
   
   if (mcpy == 1 & fcpy == 2){
      return(c(0, 0, 1, 1))
   }

   if (mcpy == 1 & fcpy == 1){
      return(c(0, 0, 1, 1))
   }

   if (mcpy == 1 & fcpy == 0){
      return(c(0, 1))
   }
   if (mcpy == 0 & fcpy == 1){
      return(c(0, 0))
   }
   return(c())
}

get_fetal_mutation_allele_count <- function(mcpy, fcpy){
   stopifnot(mcpy <= 2 & mcpy >= 0)
   stopifnot(fcpy <= 2 & fcpy >= 0)

   if (mcpy == 2 & fcpy == 2){
      return(c(0, 1, 0, 1, 2, 1, 2))
   } 

   if (mcpy == 2 & fcpy == 1){
      return(c(0, 0, 1, 1))
   }
   
   if (mcpy == 1 & fcpy == 2){
      return(c(0, 1, 1, 2))
   }

   if (mcpy == 1 & fcpy == 1){
      return(c(0, 1, 0, 1))
   }

   if (mcpy == 1 & fcpy == 0){
      return(c(0, 0))
   }
   if (mcpy == 0 & fcpy == 1){
      return(c(0, 1))
   }

   return(c())
}

genotyping_given_copy_number <- function(data, maternal_copy_number, fetal_copy_number, fetal_fraction, posterior_prob, annocol, target_loci){
    if (maternal_copy_number == 2 & fetal_copy_number == 2){
        return(data)
    }

    total <- maternal_copy_number * (1 - fetal_fraction) / 2 + fetal_copy_number * fetal_fraction / 2

    gts <- get_possible_genotypes(maternal_copy_number, fetal_copy_number)
    if (is.null(gts)){
        return(data)
    }

    mat_allele_count <- get_maternal_mutation_allele_count(maternal_copy_number, fetal_copy_number)
    fetal_allele_count <- get_fetal_mutation_allele_count(maternal_copy_number, fetal_copy_number)
    expect_mutation_ratio <- (mat_allele_count*(1-fetal_fraction) / 2 + fetal_allele_count * fetal_fraction / 2) / total
    error_rate <- 0.0001
    expect_mutation_ratio <- sapply(expect_mutation_ratio, function(x){
        if (x == 0){
            x <- error_rate
        }
        if (x == 1){
            x <- 1 - error_rate
        }
        x
    })

    print(gts)
    print( expect_mutation_ratio) 

    target_loci_idx <- which(data[, annocol] %in% target_loci)
    if(is.null(target_loci_idx)){
        stop("failed to extract data for target loci!")
    }
    gt_and_pp <- t(sapply(target_loci_idx, FUN = function(i, data, u, G){
        N <- data[i, "DepTotal"]
        if (is.null(N)){
            stop("get DepTotal for target loci failed!")
        }
        b <- data[i, "MutTotal"]
        if (is.null(b)){
            stop("get MutTotal for target loci failed!")
        }
        
        gammak <- sapply(1:length(u), FUN=function(k){
            dbinom(as.numeric(b), as.numeric(N), u[k])
        })
        mutid <- data[i, "MutID"]
        mgts <- get_maternal_genotype(maternal_copy_number, fetal_copy_number, mutid)
        fgts <- get_fetal_genotype(maternal_copy_number, fetal_copy_number, mutid)

        best <- match(max(gammak), gammak)
        matg <- mgts[best]
        fetg <- fgts[best]
        g <- G[best]
        p <- gammak / sum(gammak)
        c(matg, fetg, g, p[best] * posterior_prob)

    }, data = data, u=expect_mutation_ratio, G=gts))
    print(  gt_and_pp )

    # replace old gt and pp 
    #print(data[target_loci_idx, c("MutID", "MutTotal", "DepTotal", "FreqTotal", "Refined.genotype.EM", "p")])
    data[target_loci_idx, c("MaternalGenotype", "FetalGenotype", "Refined.genotype.EM", "p")] <- gt_and_pp
    data
}

regenotype_SNV_mutation_in_CNV_region <- function(data = stop("data is NULL"), fetal_fraction = 0.05, annocol = 5, target_loci=c("WS", "QS", "CS"), given_loci=c("SEA")){
    if (is.null(target_loci) | length(target_loci) == 0){
        return(data)
    } 
    # get sea genotyping 
    given_loci_rowidx <- which(data[, annocol] %in% given_loci)
    if (length(given_loci_rowidx) != 1){
        warning(sprintf("the number of given loci is not 1!\n"))
        return(data)
    }
    given_loci_gt <- data[given_loci_rowidx, "Refined.genotype.EM"]
    if (is.null(given_loci_gt)){
        stop("failed to get given loci's genotype!")
    }
    given_loci_gt_posterior_prob <- data[given_loci_rowidx, "p"]
    if (is.null(given_loci_gt_posterior_prob)){
        stop("failed to get given loci's genotype probability!")
    }

    get_maternal_fetal_copy_number <- function(gt){
        switch(as.character(gt),
            AAaa = c(2, 2),
            AAab = c(2, 1),
            ABaa = c(1, 2),
            ABab = c(1, 1),
            ABbb = c(1, 0),
            BBab = c(0, 1),
            BBbb = c(0, 0),
            stop(sprintf("mapformatted genotype: %s", gt)) 
        )
    }
    mfcpy <- get_maternal_fetal_copy_number(given_loci_gt) 
    genotyping_given_copy_number(data, mfcpy[1], mfcpy[2], fetal_fraction, given_loci_gt_posterior_prob, annocol, target_loci)
}





#+--------------------
# Main
#+--------------------
if (!file.exists(opt$input)) stop(sprintf("input file dose not exist: %s\n", opt$input))
df = read.table(opt$input, header = T, fill = T, sep = "\t", stringsAsFactors = FALSE, comment.char="", quote="")
cffdnaEach = as.numeric(opt$cffdna)
#colnames(df) = c("SampleName", "rsID", "Mutation", "Total", "Ratio", "cffDNA")
#df_index = which(df$Total>1)
df_index = 1:nrow(df)
if( length(df_index) > 0 ){
	HBBid = df$MutID[df_index]
	Total = df$DepTotal[df_index]
	Mutation = df$MutTotal[df_index]
	Ratio = df$FreqTotal[df_index]
	cffdnafile = df$cffDNA[df_index]
	SnpNum = df$SnpNum[df_index]
	AverDepth = df$AverDepth[df_index]
	AverDepthHotspot = df$AverDepthHotspot[df_index]
	Mother = c()
	Fetal = c()
	Refined.genotype.EM = c()
	Pvalue = c()
	for( j in 1:length(Total) ){
		if(length(cffdnafile)!=0){	cffdnaEach=cffdnafile[j]	}
		re = infer_fetal_genotype_by_cffDNA(Total[j], Mutation[j], cffdnaEach)
		#Mother_fetal =  unlist( strsplit( format_gt(HBBid[j],re[1]), "[|]" ) )
		#Mother = c( Mother, Mother_fetal[1] )
		#Fetal = c( Fetal, Mother_fetal[2] )
		Mother = c(Mother, ifelse(substr(HBBid[j], 1, 2) == "rs", substr( re[1], 1, 2), format_single_gt(HBBid[j], substr(re[1], 1, 2))))
		Fetal = c(Fetal, ifelse(substr(HBBid[j], 1, 2) == "rs", substr( re[1], 3, 4), format_single_gt(HBBid[j], substr(re[1], 3, 4))))
		Refined.genotype.EM = c( Refined.genotype.EM, re[1] )
		Pvalue = c( Pvalue, re[2] )
	}
	if(length(cffdnafile)!=0){ cffdnaEach = cffdnafile }
#	print (df$Sample[df_index])
#	print (HBBid)
#	print (Mutation)
#	print (cffdnaEach)
	if(length(SnpNum) == 0){SnpNum = 0}
	if(length(AverDepth) == 0){AverDepth= round( mean(Total),0 )}
	if(length(AverDepthHotspot) ==0){AverDepthHotspot = round( mean(Total,0) )}
	resultEM = data.frame(
		Sample = df$Sample[df_index],
		MutID = HBBid,
		MutTotal = Mutation,
		DepTotal = Total,
		FreqTotal = Ratio,
		fetal.EM = cffdnaEach,
		MaternalGenotype = Mother,
		FetalGenotype = Fetal,
		Refined.genotype.EM = Refined.genotype.EM,
		p = as.numeric(Pvalue),
		stringsAsFactors = FALSE
	)
	colnames(resultEM)=c("Sample","MutID","MutTotal","DepTotal","FreqTotal","fetal.EM","MaternalGenotype","FetalGenotype","Refined.genotype.EM","p")  # ,"SnpNum","AverDepth","AverDepthHotspot")
	resultEM$MutID <- sub("^-","'-",resultEM$MutID)
	resultEM$MaternalGenotype <- sub("^-","'-",resultEM$MaternalGenotype)
	resultEM$FetalGenotype <- sub("^-","'-",resultEM$FetalGenotype)
	write.table(resultEM, file=paste0(opt$outfile,"_3"), quote = F, row.names = F, col.names = T, sep = "\t")
	resultEM <- regenotype_SNV_mutation_in_CNV_region(resultEM, cffdnaEach, 2)
	depth_stat <- calculate_aver_depth(resultEM, 10, 4, 2)
	data <- df[df[, 4] >= 100, c(4, 3, 5)]
	snp_num_used <- nrow(data)
	resultEM <- transform(resultEM,
        FreqTotal = MutTotal / DepTotal,
        SnpNum = snp_num_used,
        AverDepth = depth_stat$overall_aver_depth,
        AverDepthHotspot = depth_stat$hotspot_aver_depth,
        MutID = sub("^-","'-", MutID),
        MaternalGenotype = sub("^-","'-", MaternalGenotype),
        FetalGenotype = sub("^-","'-", FetalGenotype))

	write.table(resultEM, file=opt$outfile, quote = F, row.names = F, col.names = T, sep = "\t")
}


