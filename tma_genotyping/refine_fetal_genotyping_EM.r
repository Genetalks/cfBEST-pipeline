#!/share/software/software/R-3.0_install/R-3.0.1/bin/Rscript

library(getopt)
library(dplyr)

#+--------------------
# get options
#+--------------------
spec <- matrix(c(
	'help', 'h', 0, "logical", "help",
	'verbose', 'v', 2, "integer", "verbose mode, default [1]",
	'input', 'i', 1, "character", "input brief.xls file, forced.",
	'output', 'o', 1, "character", "output brief.EM.xls with fetalDNA and genotyping info, forced.",
	'offset', 's', 2, "integer", "colume offset of total templates in input table, default [6]",
	'annocol', 'a', 2, "integer", "the colume number of dbSNP annotation info in input table, default [5]",
	'mindep', 'd', 2, "integer", "min snp depth, default [100]",
	'infer_by_f', 'e', 0, "logical", "infer fetal genotypes using fetal dna fraction"
), byrow = TRUE, ncol = 5)

opt <- getopt(spec)

#+--------------------
# check options
#+--------------------
if ( !is.null(opt$help) | is.null(opt$input) | is.null(opt$output)) {
	cat(getopt(spec, usage=TRUE))
	q(status=1)
}

#+--------------------
# some default options
#+--------------------

if (is.null(opt$annocol)) opt$annocol <- 5
if (is.null(opt$offset)) opt$offset <- 6
if (is.null(opt$mindep)) opt$mindep <- 100

#+--------------------
# EM - settings
#+--------------------

f <- 0.1  ## the fetal DNA fraction
u <- c(1e-8, f/2, 0.5-f/2, 0.5)
pii <- c(0.25, 0.25, 0.25, 0.25)  

calculate_log_likelihood <- function(data, u = c(1e-8, 0.05, 0.45, 0.5), pii = c(0.25, 0.25, 0.25, 0.25)){
    stopifnot(!is.null(data))
    
    Pi <- sapply(1:nrow(data), FUN=function(i){
        Ni <- data[i, 1]
        bi <- data[i, 2]
        
        Pik <- sapply(1:length(u), FUN=function(k){
            if(!is.na(u[k])){ pii[k] * dbinom(bi, Ni, u[k])}
            else {0}
        })
        log(sum(Pik))
    })
    sum(Pi)
}

grid_search_f <- function(data, u, pii, l0, f_start = 0.05, f_end = 0.5, f_step = 0.01){
    f <- NA
    if (is.na(u[2]) & is.na(u[3])){
        f <- NA
    }else if (is.na(u[2])){
        f <- 2*(u[4]-u[3])
    }else if(is.na(u[3])){
        f <- 2*u[2]
    }else{
        f <- mean(2*u[2], 2*(u[4]-u[3]))
    }
    if(!is.na(f)){
        tmp.u <- u
        silent <- sapply(seq(f_start, f_end, f_step), FUN=function(fi){
            tmp.u[2] <- fi / 2
            tmp.u[3] <- tmp.u[4] - fi / 2
            l1 <- calculate_log_likelihood(data, tmp.u, pii)
            if (l1 > l0){ # log-likelihood improved, update u and f
                u <- tmp.u
                f <- fi
                ## update l1 
                l0 <- l1
            }
        })
    }
    list(u=u, f=f, loglikelihood=l0)
}

fetal_genotyping_EM <- function(data = NULL, u = c(1e-8, 0.05, 0.45, 0.5), pii = c(0.25, 0.25, 0.25, 0.25), f_start = 0.05, f_end = 0.5, f_step = 0.01){
    stopifnot(!is.null(data))
    
    ## init parameters to control EM iteration
    n_iter <- 1  # the number of iterations
    e <- 1e-8    # the stop condition of EM-iteration, threshold for log-likelihood improvement
    max_iter <- 1000  # maximal iteration times
    nG <- length(pii) # the number of genotypes in our model, this is always 4
    T <- nrow(data)  # the number of SNP/INDEL
    Gn <- NULL   # genotype vector for each SNP during the nth iteration, length(Gn) == T
    f <- 0       # the fetal DNA fraction
    best_params <- list(f=f, u=u, pii=pii, Gn=Gn, likelihood=NA, iter=0, nsnp=T)
    # calculate initial log likelihood 
    l0 <- calculate_log_likelihood(data, u, pii)
    while (n_iter <= max_iter){
        
	# the E-step, calculate Gamma using bayes' rule and infer genotypes for each SNP/INDEL
        Gn <- sapply(1:T, FUN=function(i){ # i = 1..T
            Ni <- data[i, 1]
            bi <- data[i, 2]
            gammak <- sapply(1:nG, FUN=function(k){ # k = 1..4
                pii[k] * dbinom(as.numeric(bi), as.numeric(Ni), u[k])
            })
            if (bi / Ni < 0.001){  # too small mutation fraction
		        1	
	        }else{
		        match(max(gammak), gammak)
	        }
        })
     
        # the M-step, update parameters
        ## update pii
        pii <- sapply(1:nG, FUN=function(k){
            length(which(Gn == k)) / T
        })
        ## update u 
        u <- sapply(1:nG, FUN=function(k){
            if (sum(Gn == k) != 0){
                sum(data[, 2] * (Gn == k)) / sum(data[, 1] * (Gn == k))
            }else{
                switch(k, 1e-10, NA, NA, 0.5)
            }   
        })
        
        ## grip search f
        l1 <- calculate_log_likelihood(data, u, pii)
        
        ret <- grid_search_f(data, u, pii, l1, f_start, f_end, f_step)
		u <- ret$u
        f <- ret$f
        l1 <- ret$loglikelihood

		cat(sprintf(">>> The %dth iteration \n", n_iter))
		cat(sprintf("=== params. estimation ===\n"))
		cat(sprintf("Prop. of different genetypes: pi = (%s)\n", paste(round(pii, digits=2), collapse=",")))
		cat(sprintf("Mean estimation of B allele freq.: u = (%s)\n", paste(round(u, digits=2), collapse=",")))
		cat(sprintf("log likelihood: %.2f\n", l1))
		cat(sprintf("estimated fetal DNA fraction: %.6f\n\n", f))

        ## check convergency
        if (l1 - l0 < e){
            break
        }
		best_params <- list(f=f, u=u, pii=pii, Gn=Gn, likelihood=l1, iter=n_iter, nsnp=T)
        # update log likelihood 
        l0 <- l1
        # update iteration times
        n_iter <- n_iter + 1
    }
    #list(f=f, u=u, pii=pii, Gn=Gn, iter=n_iter)
    best_params
}

infer_genotype <- function(Ni, bi, u, pii){
    G <- c("AAaa", "AAab", "ABaa", "ABab")
    mG <- c("BBbb", "BBab", "ABbb", "ABab") # the mirror of genotype space
    
    Ni <- as.numeric(Ni)
    bi <- as.numeric(bi)
    is_mirror <- ifelse(Ni-bi>bi, 0, 1);
    b <- ifelse(Ni-bi>bi, bi, Ni-bi)
    
    gammak <- sapply(1:length(u), FUN=function(k){ # k = 1..4
        pii[k] * dbinom(as.numeric(b), as.numeric(Ni), u[k])
    })
    
    i <- ifelse(b/Ni < 0.001, 1, match(max(gammak), gammak))
    g <- ifelse(is_mirror, mG[i], G[i])
    p <- gammak/sum(gammak)
    c(g, p[i])
}

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

sdbinom <- function(x, N, u, pii){
    pii*dbinom(as.integer(x*N), size=as.integer(N), prob=u)
}

plot_mix_binom <- function(data, model){
    require(ggplot2)
    N <- data[, 1]
    b <- data[, 2]
    
    N_mean <- sapply(1:4, function(k){
	    mean(N[model$Gn==k], na.rm=T)
    })
    colnames(data) <- c("N", "b", "f")

    p <- ggplot(data) + 
    #geom_histogram(aes(x=f,y=..density..),fill="white",color="gray", bins=150) +
    stat_function(fun=sdbinom, args=list(u=model$u[1],N=N_mean[1],pii=model$pii[1]),fill="#DA70D6",geom="polygon", alpha=0.5) +
    stat_function(fun=sdbinom, args=list(u=model$u[2],N=N_mean[2],pii=model$pii[2]),fill="#00FF0080",geom="polygon", alpha=0.5) +
    stat_function(fun=sdbinom, args=list(u=model$u[3],N=N_mean[3],pii=model$pii[3]),fill="#FF000080",geom="polygon", alpha=0.5) +
    stat_function(fun=sdbinom, args=list(u=model$u[4],N=N_mean[4],pii=model$pii[4]),fill="skyblue",geom="polygon", alpha=0.5) +
    xlim(c(-0.5, 1)) + xlab("B allele frequency")
    print(p)
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
    print(expect_mutation_ratio)

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
    print(gt_and_pp)

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

## load input
if (!file.exists(opt$input)) stop(sprintf("input file dose not exist: %s\n", opt$input))

df <- read.table(opt$input, header = T, fill = T, sep = "\t", stringsAsFactors=F, comment.char="", quote="")
flt <- grep("^rs", df[, opt$annocol], perl = T)

df[, opt$offset] <- as.numeric(df[, opt$offset])
df[, opt$offset+1] <- as.numeric(df[, opt$offset+1])
df[, opt$offset+2] <- as.numeric(df[, opt$offset+2])

data <- df[df[, opt$offset+1] >= opt$mindep, c(opt$offset+1, opt$offset, opt$offset+2)]
data[, 2] <- ifelse(data[,3]>0.525, data[,1]-data[,2], data[, 2])
data[, 1] <- as.numeric(data[, 1])
data[, 2] <- as.numeric(data[, 2])

#data <- data[flt, ]
snp_num_used <- nrow(data)

# estimate cffDNA using EM algorithm
theta <- fetal_genotyping_EM(data, u, pii)

# plot mix binorm model 
pdf(sprintf("%s.model.pdf", opt$output))
plot_mix_binom(data, theta)
dev.off()

# result
result <- NULL
if (!is.null(opt$infer_by_f)){
	cat(sprintf("\n#### infer genotypes using deduced fetal dna fration and experimental pii\n"))
	result <- t(unlist(sapply(1:nrow(df), FUN=function(i){
	   infer_fetal_genotype_by_cffDNA(df[i, opt$offset+1], df[i, opt$offset], theta$f) 
	})))

}else{
	result <- t(unlist(sapply(1:nrow(df), FUN=function(i){
	   infer_genotype(df[i, opt$offset+1], df[i, opt$offset], theta$u, theta$pii) 
	})))
}

if (is.null(result)){
    stop("infer fetal fraction failed!")
}

rsid <- gsub("\\w+:([^_]+)_?.*", "\\1", df[, opt$annocol], perl=T)
transform(df, 
    fetal.EM = as.numeric(theta$f),
    MaternalGenotype = ifelse(substr(df[, opt$annocol], 1, 2) == "rs", 
                              substr(result[, 1], 1, 2), 
                              format_single_gt(rsid, substr(result[, 1], 1, 2))),
    FetalGenotype = ifelse(substr(df[, opt$annocol], 1, 2) == "rs", 
                           substr(result[, 1], 3, 4), 
                           format_single_gt(rsid, substr(result[, 1], 3, 4))),
    Refined.genotype.EM = result[, 1],
    p = as.numeric(result[, 2])) %>% 
    mutate_if(is.factor, as.character) -> df  # factor -> character

# regenotype alpha mutation given SEA's genotype 
df <- regenotype_SNV_mutation_in_CNV_region(df, theta$f, opt$annocol)

depth_stat <- calculate_aver_depth(df, 10, opt$offset+1, opt$annocol)
df <- transform(df,
        FreqTotal = MutTotal / DepTotal,
        SnpNum = snp_num_used,
        AverDepth = depth_stat$overall_aver_depth,
        AverDepthHotspot = depth_stat$hotspot_aver_depth,
        MutID = sub("^-","'-", MutID),
        MaternalGenotype = sub("^-","'-", MaternalGenotype),
        FetalGenotype = sub("^-","'-", FetalGenotype))

# save result 
warnings()
write.table(df, file=opt$output, quote = F, row.names = F, sep = "\t")

rdata_file <- sprintf("%s.params.RData", opt$output)
save(theta, file = rdata_file)
