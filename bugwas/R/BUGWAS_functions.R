# File KmerAnalysisMain.R
# Authors: Earle, S. G., Wu, C.-H. and Wilson, D. J.
#
# Copyright (C) 2015 University of Oxford
#
# This file is part of the bugwas R package.
# See the NOTICE file distributed with this work for additional
# information regarding copyright ownership and licensing.
#
# This is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
#  This is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this software package; if not, write to the
# Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
# Boston, MA  02110-1301  USA

################################################################################################
## Compact SNPs down to patterns.
## @gen: SNP data, columns = samples, rows = SNPs. Can contain biallelic or multiallelic SNPs
## @ps: SNP positions for each row in gen
## @prefix: prefix for output files
##
## Outputs:
## XX: unique biallelic SNP patterns (0s and 1s)
## XX.tritetra: unique tri and tetra allelic SNP patterns (0-3s)
## bippat: for each pattern in XX, how many SNPs does it represent
## snppat: for each pattern in XX.tritetra, how many SNPs does it represent
## pattern: for each biallelic SNP, which pattern is it
## pattern.snps: for each tri or tetra allelic SNP, which pattern is it
################################################################################################

compact_SNPs <- function(gen = NULL,
						 ps = NULL,
						 prefix = NULL){
						 	
	# Sanity check
	if(any(gen!="A" & gen!="C" & gen!="G" & gen!="T")){
		stop("\nError: SNP data must be imputed, only containing bases A, C, G and T\n")
	}

	# First pass count the number of A,C,G,T per site
	fa <- gen[, 1]
	m <- matrix(0, 4, length(fa))
	for(i in 1:ncol(gen)) {
		# Read the mapcall file
		fa = gen[, i]
		m[1, fa=="A"] = m[1, fa=="A"]+1
		m[2, fa=="C"] = m[2, fa=="C"]+1
		m[3, fa=="G"] = m[3, fa=="G"]+1
		m[4, fa=="T"] = m[4, fa=="T"]+1
	}

	nalleles <- (m[1, ]>0) + (m[2, ]>0) + (m[3, ]>0) + (m[4, ]>0)
	is.poly <- nalleles>1
	if(length(which(is.poly))==0){
		stop("\nError: no variable sites in SNP data\n")
	}
	if(any(!is.poly)){
		stop("\nError: there are invariant sites in the SNP data\n")
	}
	is.fixed <- !is.poly
	isbiallelic <- nalleles==2
	if(length(which(isbiallelic))==0){
		stop("\nError: no biallelic variable sites in SNP data\n")
	}
	
	
	# Order the alleles at SNPs
	allele.id <- matrix(c("A", "C", "G", "T")[apply(m[1:4, is.poly],
						2, order, decreasing=TRUE)], nrow=4)
	
	# Output filenames
	# biallelic polymorphisms encoded -1 (missing) 0 (allele 0) 1 (allele 1)
	bip_outfile <- paste0(prefix, ".bip.patterns.txt");		
	# positional and allelic information for biallelic polymorphisms
	bipinfo_outfile <- paste0(prefix, ".bipinfo.txt");		
	# Allocate memory for bip and snp, because need to transform so cannot output on the fly
	bip <- matrix(NA,sum(isbiallelic),ncol(gen))
	colnames(bip) <- colnames(gen)
	
	poly <- FALSE
	if(length(which(is.poly))>length(which(isbiallelic))){
		poly <- TRUE
		# Output filenames
		# non-biallelic polymorphisms encoded similarly
		snp_outfile <- paste0(prefix, ".snp.patterns.txt");
		# non-biallelic polymorphism positional information
		snpinfo_outfile <- paste0(prefix, ".snpinfo.txt");
		# Allocate memory for bip and snp, because need to transform so cannot output on the fly
		snp <- matrix(NA, sum(is.poly & !isbiallelic), ncol(gen))
		colnames(snp) <- colnames(gen)
	}
	
	# Second pass over the files: populate bip and snp objects
	for(i in 1:ncol(gen)) {
		# Read the mapcall file
		fa <- gen[, i]
		fa.poly <- fa[is.poly]
	#    This should not be needed (could insert sanity check):
	#    fa[is.poly][fa.poly!="A" & fa.poly!="C" & fa.poly!="G" & fa.poly!="T"] = "N"
		# Populate bip object
		bip[, i] <- -1
		bip[fa[isbiallelic]==allele.id[1, isbiallelic[is.poly]], i] = 0
		bip[fa[isbiallelic]==allele.id[2, isbiallelic[is.poly]], i] = 1
		# Populate snp object
		if(poly){
			snp[,i] <- -1
			snp[fa[is.poly & !isbiallelic]==allele.id[1, !isbiallelic[is.poly]], i] = 0
			snp[fa[is.poly & !isbiallelic]==allele.id[2, !isbiallelic[is.poly]], i] = 1
			snp[fa[is.poly & !isbiallelic]==allele.id[3, !isbiallelic[is.poly]], i] = 2
			snp[fa[is.poly & !isbiallelic]==allele.id[4, !isbiallelic[is.poly]], i] = 3
		}
	}
	
	# Convert BIP and SNP patterns to factors to identify equivalencies
	bip.pat <- factor(apply(bip, 1, paste, collapse=""))
	# Record only unique patterns, and record the pattern equivalence in the bipinfo file
	bip.pat1 <- match(levels(bip.pat), bip.pat)
	
	if(poly){
		# Convert BIP and SNP patterns to factors to identify equivalencies
		snp.pat <- factor(apply(snp, 1, paste, collapse=""))
		# Record only unique patterns, and record the pattern equivalence in the snpinfo file
		snp.pat1 <- match(levels(snp.pat), snp.pat)
	}

	# Output compacted bip and snp objects
	write.table(bip[bip.pat1, ], bip_outfile, row=FALSE, col=TRUE, quote=FALSE, sep="\t")
	# Output info files
	bipinfo <- data.frame("Position"=ps[which(isbiallelic)],
						 "Allele0"=allele.id[1,isbiallelic[is.poly]],
						 "Allele1"=allele.id[2,isbiallelic[is.poly]],
						 "A"=m[1,isbiallelic],
						 "C"=m[2,isbiallelic],
						 "G"=m[3,isbiallelic],
						 "T"=m[4,isbiallelic],
						 "Pattern"=as.numeric(bip.pat)); 
	write.table(bipinfo, bipinfo_outfile, row=FALSE, col=TRUE, quote=FALSE, sep="\t")
	ps.bips <- bipinfo$Position
	bippat <- sapply(1:max(bipinfo$Pattern), function(x)sum(bipinfo$Pattern==x))
	ipat <- bipinfo$Pattern
	rm(bipinfo)
	
	if(poly){
		# Output compacted bip and snp objects
		write.table(snp[snp.pat1, ], snp_outfile, row=FALSE, col=TRUE, quote=FALSE, sep="\t")
		# Output info files
		snpinfo <- data.frame("Position"=ps[which(is.poly & !isbiallelic)],
							 "Allele0"=allele.id[1,!isbiallelic[is.poly]],
							 "Allele1"=allele.id[2,!isbiallelic[is.poly]],
							 "Allele2"=allele.id[3,!isbiallelic[is.poly]],
							 "Allele3"=allele.id[4,!isbiallelic[is.poly]],
							 "A"=m[1,is.poly & !isbiallelic],
							 "C"=m[2,is.poly & !isbiallelic],
							 "G"=m[3,is.poly & !isbiallelic],
							 "T"=m[4,is.poly & !isbiallelic],
							 "Pattern"=as.numeric(snp.pat)); 
		write.table(snpinfo, snpinfo_outfile, row=FALSE, col=TRUE, quote=FALSE, sep="\t")
	}
	
	
	if(poly){
		snppat <- sapply(1:max(snpinfo$Pattern), function(x)sum(snpinfo$Pattern==x))
		ipat.snps <- snpinfo$Pattern
	} else {
		snp <- NULL
		snp.pat1 <- NULL
		snppat <- NULL
		ipat.snps <- NULL
		snpinfo <- NULL
	}

	# Return
	
	return(list("XX" = bip[bip.pat1, ],
				"XX.tritetra" = snp[snp.pat1, ],
				"bippat" = bippat,
				"snppat" = snppat,
				"pattern" = ipat,
				"pattern.snps" = ipat.snps,
				"ps" = ps.bips,
				"ps.snps" = snpinfo$Position,
				"n.triallelic" = length(which(nalleles == 3)),
				"n.tetraallelic" = length(which(nalleles == 4))))

}

################################################################################################
## Compact presence absence matrices down to patterns
## @var: Presence absence matrix, of 0s and 1s
##
## Outputs:
## var: unique patterns of presence and absence
## pattern: for each variant, which pattern it is
## varpat: for each variant pattern, how many variants it represents
################################################################################################

compact_variants <- function(var = NULL){
	
	is.var <- apply(var, 1, function(x) length(unique(x)))
	
	if(any(is.var==1)){
		stop("\nError: there are invariant sites in the presence absence data\n")
	}
	if(any(is.var>2)){
		stop("\nError: presence absence data must be binary\n")
	}
	if(any(var!=0 & var!=1)){
		stop("\nError: presence absence matrix must be coded in 0s and 1s\n")
	}

	var.pat <- factor(apply(var, 1, paste, collapse=""))

	var.pat1 <- match(levels(var.pat), var.pat)
	
	ipat <- match(var.pat, levels(var.pat))
	
	varpat <- sapply(1:max(as.numeric(var.pat)), function(x) sum(as.numeric(var.pat)==x))

	return(list("var" = var[var.pat1, ],
				"pattern" = as.numeric(var.pat),
				"varpat" = varpat))
		
}

################################################################################################
## Rescale variant patterns by how many variants they represent.
## @var: unique variant patterns
## @varpat: for each variant pattern, how many variants it represents
##
## Outputs:
## var: unique patterns of presence and absence rescaled by how many variants each represents
################################################################################################

rescale_variants <- function(var = NULL,
							 varpat = NULL){
	var <- t(var)
	var <- t(t(var) * sqrt(varpat))
	var.f <- colMeans(var)
	var <- t(t(var) - var.f)
	return(var)
}

################################################################################################
## Get tree patterns.
## @tree: phylogeny of all genome sequences in SNP matrix
## @tiporder: vector of sample IDs in the order of the SNP matrix
##
## Outputs:
## pat: branch patterns
## ancestral_edge: tree edges
## ancestral_edge.length: edge lengths
################################################################################################

# Get tree patterns
tree2patterns <- function(tree = NULL,
						  tiporder = NULL) {
						  	
	if(is.null(tiporder)){
		tiporder <- tree$tip.label
	}
	n <- length(tree$tip.label)
	
	
	mtp <- matrix(0,n,n+tree$Nnode)
	
	#Sorting out tree patterns for external nodes.
	mtp[,1:n] <- diag(n)
	
	#Sorting out tree patterns for internal nodes.
	for(i in 1:tree$Nnode) {
		wh <- match(ape::extract.clade(tree, n+i)$tip.label, tree$tip.label)
		mtp[wh, n+i] <- 1
	}
	
	#Reorder the rows to match the order of the individuals
	mtp <- mtp[match(tiporder, tree$tip.label), ]
	
	if(is.null(tree$node.label)){
		tree$node.label = paste("node",1:tree$Nnode+n, sep="")
	}
	
	
	
	mtp.f <- apply(mtp, 2, mean)
	mtp[, mtp.f>0.5] <- 1-mtp[, mtp.f>0.5]
	mtp.f <- apply(mtp, 2, mean)
	
	#Label the rows and columns of tree patterns.
	rownames(mtp) <- tiporder
	colnames(mtp) = c(tiporder, tree$node.label)
	
	
	edge <- match(1:ncol(mtp), tree$edge[, 2])
	edge.length <- tree$edge.length[edge]
	return(list("pat" = mtp, "labelled.tree" = tree,
				"ancestral_edge" = edge,
				"ancestral_edge.length" = edge.length))
}


################################################################################################
## Get correlations between variant patterns and principal components.
## @XX: variant matrix scaled by how many variants each pattern represents
## @pca: principal component analysis
## @npcs: number of principal components (= number of individuals)
##
## Outputs:
## which.pc: for each variant pattern, which PC it is most correlated to
## max.cor.pc: for each variant pattern, the maximum correlation value
################################################################################################

get_correlations <- function (XX = NULL, 
							  pca = NULL,
							  npcs = NULL,
							  id = NULL,
							  all.cor  = FALSE){
	
	
	cor.XX.pca <- cor(XX,pca[, 1:npcs])
	cor.XX.pca[is.na(cor.XX.pca)] = 0
	which.pc <- apply(abs(cor.XX.pca), 1, which.max)
	max.cor.pc <- apply(abs(cor.XX.pca), 1, max)
	if(all.cor){
		return(list("which.pc" = which.pc, "max.cor.pc" = max.cor.pc, "all.cor.pc" = cor.XX.pca))
	}else{
		return(list("which.pc" = which.pc, "max.cor.pc" = max.cor.pc))
	}
	
}


################################################################################################
## Get Bayesian Wald Test inputs.
## @fit.lmm: ridge regression results
## @pca: principal component analysis
## @svd.XX: single value decomposition of biallelic SNPs
##
## Outputs:
## pca.Ebeta
## pca.Vbeta
################################################################################################

get_wald_input <- function(fit.lmm = NULL,
						   pca = NULL,
						   svd.XX = NULL,
						   y = NULL,
						   npcs = NULL,
						   XX = NULL){
	
	# Get full posterior covariance matrix for Bayesian Wald Test
	# Need the full posterior covariance matrix for the Bayesian Wald test,
	# to get the posterior uncertainty for each point estimate
	lambda = fit.lmm$lambda_MLE
	Cstar = diag(lambda * svd.XX$d^2 / (lambda * svd.XX$d^2 + 1))
	# For the null model, the posterior mean and variance (may be slow!)
	astar = t(y)%*%y - t(y)%*%XX%*%fit.lmm$Ebeta
	dstar = length(y)
	tau = as.numeric((dstar-2)/astar)

	rotation = t(pca$rotation[,1:npcs])
	# Based on the PCA rotations of raw genetic diversity
	pca.Ebeta = rotation %*% fit.lmm$Ebeta
 
	rtr = tcrossprod(rotation,rotation); # Should be n (sample size) by n
	rv = rotation %*% svd.XX$v; # Should be n by n
	pca.Vbeta = rv %*% Cstar; # Should be n by n
	pca.Vbeta = tcrossprod(pca.Vbeta,rv); # Should be n by n
	pca.Vbeta = lambda/tau * (rtr - pca.Vbeta); # Should be n by n
	
	return(list("Ebeta" = pca.Ebeta, "Vbeta" = pca.Vbeta))

	
}

################################################################################################
## Get order of principal components by Bayesian Wald Test results.
## @fit.lmm: ridge regression results
## @pca: principal component analysis
## @svd.XX: single value decomposition of biallelic SNPs
##
## Outputs:
## pca.Ebeta
## pca.Vbeta
################################################################################################

get_pc_order <- function(p.pca.bwt = NULL,
						 signif_cutoff = NULL){
	o <- order(p.pca.bwt,decreasing=T)
	# How many are above a significance cut-off
	if(length(which(p.pca.bwt>= signif_cutoff))>20){
		pc.lim <- 1:20
	} else {
		if(length(which(p.pca.bwt>= signif_cutoff))>0){
			pc.lim <- 1:length(which(p.pca.bwt>= signif_cutoff))
		} else {
			pc.lim <- NULL
		}
	}
	
	return(list("pc_order" = o, "pc.lim" = pc.lim))
}

################################################################################################
## Convert tri and tetra allelic SNP patterns into multiple columns for GEMMA.
## @XX: SNP patterns for tri and tetra allelic sites
##
## Outputs:
## List the length of the number of SNPs
## For each SNP, there will be number of variants -1 columns
################################################################################################

convert_tri_tetra <- function(XX = NULL){
	if(length(unique(XX))<3){
		stop("\nError: only tri and tetra allelic SNPs should be run as covariates in GEMMA\n")
	} else if(length(unique(XX))==3){
		dataCOV <- XX
		dataCOV[which(XX==1)] <- 0
		dataCOV[which(XX==2)] <- 1
		XX[which(XX==2)] <- 0
	} else if(length(unique(XX))==4){
		dataCOV <- matrix(0, ncol=2, nrow=length(XX))
		dataCOV[which(XX==2), 1] <- 1
		dataCOV[which(XX==3), 2] <- 1
		XX[which(XX==2 | XX==3)] <- 0
	}
	return(cbind(XX, dataCOV))
	
}


################################################################################################
## LMM functions.

################################################################################################

extract_lambda_lognull <- function(datafiles = NULL){
	a <- scan(datafiles, what=character(0), sep="\n")
	b <- unlist(strsplit(a[13], " "))
	b <- b[length(b)]
	c <- unlist(strsplit(a[17], " "))
	c <- c[length(c)]
	return(list("lambda" = b, "lognull" = c))
}

getLH1 <- function(datafiles = NULL){
	a <- read.table(datafiles, header=T)
	a <- a$logl_H1
	return(a)
}

write_pheno <- function(pheno = NULL, 
						prefix = NULL){
	pheno.file <- paste0(prefix, "_gemma_phenotype.txt")
	cat(pheno, file = pheno.file, sep="\n")
	return(pheno.file)
	
}

get_deviance <- function(LH1 = NULL, lognull = NULL){
	D <- as.numeric(2*(as.numeric(LH1) - as.numeric(lognull)))
	return(D)
}

get_pvals <- function(D = NULL){
	pvals <- pchisq(D["D"], D["num.alleles"]-1, low=F)
	return(pvals)
}


################################################################################################
## Get GEMMA kinship matrix.
## @XX: List of patterns from using compact_SNPs
## @pattern: Vector, for each biallelic SNP, which pattern is it
## @prefix: Output file prefix
## @maf: Minor allele frequency to test using GEMMA (Default: 0, all variants are tested)
## @path: Path to where GEMMA is installed
## @dir: Working directory
## @pheno.file: File containing phenotype to be tested (in same order as columns of XX)
##
## Outputs:
## Path to relatedness matrix built from biallelic SNPs
################################################################################################

get_kinship <- function(XX = NULL,
						pattern = NULL,
						prefix = NULL,
						path = NULL,
						dir = NULL,
						maf = NULL,
						pheno.file = NULL){
	message("Get kinship matrix")
	
	gen.output.file <- paste0(prefix, "_gemma_kinship_genfile.txt")
	snp.output.file <- paste0(prefix, "_gemma_kinship_snpfile.txt")
	
	gen.file <- XX[pattern,]
	gen.file <- cbind(paste0("rs",1:nrow(gen.file)),rep(1,nrow(gen.file)),rep(0,nrow(gen.file)),gen.file)
	snp.file <- cbind(paste0("rs",1:nrow(gen.file)), 1:nrow(gen.file), rep(24,nrow(gen.file)))
	
	write.table(gen.file, file = gen.output.file, row=F, col=F, sep="\t", quote=F)
	write.table(snp.file, file = snp.output.file, row=F, col=F, sep="\t", quote=F)
	
	system(paste0(path, " -g ", gen.output.file, " -p ", pheno.file, " -a ",
				  snp.output.file, " -gk 1 -o ", prefix, "_gemma_relmatrixout", " -maf ",maf))
	
	relmatrix <- paste0(dir, "/output/", prefix, "_gemma_relmatrixout.cXX.txt")
	
	return(relmatrix)
}

################################################################################################
## Run GEMMA software on binary data.
## @XX: Binary variant patterns
## @relmatrix: Path to a relatedness matrix built from biallelic SNP data
## @pattern: For each variant, which pattern (row of XX) is it
## @ps: For each variant, what position is it in a reference genome (if mapped data)
##		Otherwise, an ID number
## @pheno.file: File containing phenotype to be tested (in same order as columns of XX)
## @maf: Minor allele frequency to test using GEMMA (Default: 0, all variants are tested)
## @prefix: Output file prefix
## @path: Path to where GEMMA is installed
## @dir: Working directory
## @process.results: Whether to process the p-value output of LMM or just return lambda and
## 					 the log likelihood under the null
## 					 Defaults to TRUE, which processes LMM p-value results
##
## Outputs:
## Likelihood ratio test results, the log likelihood under the null, lambda
################################################################################################

run_lmm_bi <- function(XX = NULL,
					   relmatrix = NULL,
					   pattern = NULL,
					   ps = NULL,
					   pheno.file = NULL,
					   maf = NULL,
					   prefix = NULL,
					   path = NULL,
					   dir = NULL,
					   process.results = TRUE){
					   	
					   	
	#message("run_lmm_bi")				   	
	
	# Output file names
	gen.output.file <- paste0(prefix, "_gemma_genfile.txt")
	snp.output.file <- paste0(prefix, "_gemma_snpfile.txt")
	
	if(is.null(dim(XX))){
		XX <- matrix(XX,nrow=1)
	}
	
	# Check if the variants are binary
	num.alleles <- apply(XX, 1, function(data) length(unique(data)))
	
	if(length(which(num.alleles>2))!=0){
		stop("\nError: function run_lmm_bi requires binary variants\n")
	}
	gen.file <- cbind(paste0("pattern",1:nrow(XX)),rep(1,nrow(XX)),rep(0,nrow(XX)),XX)
	write.table(gen.file, file = gen.output.file, row=F, col=F, sep="\t", quote=F)
	
	snp.file <- cbind(paste0("pattern",1:nrow(XX)), 1:nrow(XX), rep(24,nrow(XX)))
	write.table(snp.file, file = snp.output.file, row=F, col=F, sep="\t", quote=F)
	
	system(paste0(path, " -g ", gen.output.file, " -p ", pheno.file, " -a ", snp.output.file,
				  " -k ", relmatrix," -lmm 2 -o ", prefix, "_lmmout_patterns"," -maf ", maf))
				  
	#message("LMM calculations completed successfully.")
	lmm.log <- paste0(dir, "/output/", prefix, "_lmmout_patterns.log.txt")
	
	#message(paste(c("GEMMA log file:", lmm.log), collapse=""))
	
	lambda.lognull <- extract_lambda_lognull(lmm.log)
	
	#message("Lambda estimations extracted successfully.")
	#message(paste(c("Extracted lambda: ", lambda.lognull), collapse=""))
	
	if(process.results == FALSE){
		
		return(lambda.lognull)
		
	} else {
		
		assocFile = paste0(dir, "/output/", prefix, "_lmmout_patterns.assoc.txt")
		#message(paste(c("assocFile:", assocFile), collapse=" "))
		lmm <- read.table(assocFile, header=T, sep="\t", as.is=T)
		
		#message(paste(c("Header:", names(lmm)), collapse=" "))
		
		LH1 <- lmm$logl_H1
		
		
		D <- sapply(LH1, get_deviance, lognull=lambda.lognull["lognull"], USE.NAMES=FALSE)
		pvals <- pchisq(as.numeric(D), 1, low=F)
		negLog10 <- -log10(pvals)	
		
		#message("flag1")
	
		m <- match(pattern, lmm$ps)
		pattern <- pattern[!is.na(m)]
		ps <- ps[!is.na(m)]
		
		#message("flag2")
		
		lmm <- cbind(lmm, negLog10)
		lmm <- lmm[match(pattern,lmm$ps), ]
		lmm$ps <- ps
		
		#message("flag3")
		
		write.table(lmm, file = paste0(prefix, "_lmmout_allSNPs.txt"), sep="\t",
					row=F, col=T, quote=F)
		
		return(list("lmm" = lmm, "lognull" = lambda.lognull["lognull"],
					"lambda" = lambda.lognull["lambda"]))
	}
	
}

################################################################################################
## Run GEMMA software on muliallelic data.
## @XX: Variant patterns of multiallelic data
## @relmatrix: Path to a relatedness matrix built from biallelic SNP data
## @pattern: For each variant, which pattern (row of XX) is it
## @ps: For each variant, what position is it in a reference genome (if mapped data)
##		Otherwise, an ID number
## @pheno.file: File containing phenotype to be tested (in same order as columns of XX)
## @maf: Minor allele frequency to test using GEMMA (Default: 0, all variants are tested)
## @prefix: Output file prefix
## @path: Path to where GEMMA is installed
## @dir: Working directory
##
## Outputs:
## Likelihood ratio test results
################################################################################################


run_lmm_multi <- function(XX = NULL,
						  pattern = NULL,
						  relmatrix = NULL,
						  pheno.file = NULL,
						  maf = NULL,
						  ps = NULL,
						  lognull = NULL,
						  prefix = NULL,
						  path = NULL,
						  dir = NULL){
	
	gen.output.file <- paste0(prefix, "_gemma_genfile.txt")
	snp.output.file <- paste0(prefix, "_gemma_snpfile.txt")
	cov.output.file <- paste0(prefix, "_gemma_covfile.txt")
	
	if(is.null(dim(XX))){
		XX <- matrix(XX,nrow=1)
	}
	
	# Check if the variants are multiallelic
	num.alleles <- apply(XX, 1, function(data) length(unique(data)))
	
	if(length(which(num.alleles>2))==0){
		stop("\nError: function run_lmm_multi requires variants of 3 or 4 states\n")
	}
	
	XX <- apply(XX, 1, convert_tri_tetra)
	
	for(i in 1:length(XX)){
		
		gen.info <- cbind(paste0("pattern", i), 1, 0)
		snp.info <- cbind(paste0("pattern", i), i, 24)
		
		gen.file <- matrix(c(gen.info, XX[[i]][, 1]), nrow=1)
		cov.file <- cbind(rep(1,length(XX[[i]][,1])), XX[[i]][, 2:ncol(XX[[i]])])
		
		write.table(gen.file, file = gen.output.file, sep="\t", row=F, col=F, quote=F)
		write.table(snp.info, file = snp.output.file, sep="\t", row=F, col=F, quote=F)
		write.table(cov.file, file = cov.output.file, sep="\t", row=F, col=F, quote=F)
		
		system(paste0(path, " -g ",gen.output.file, " -p ", pheno.file, " -a ", snp.output.file,
					  " -c ", cov.output.file, " -k ", relmatrix, " -lmm 2 -o ", prefix,
					  "_gemma_tritetra_lmmout_pattern", i, " -maf ", maf))
	}
	
	lmm <- paste0(dir, "/output/", prefix, "_gemma_tritetra_lmmout_pattern",
				  1:length(XX), ".assoc.txt")
	lmm.log <- paste0(dir, "/output/", prefix, "_gemma_tritetra_lmmout_pattern",
					  1:length(XX), ".log.txt")
	
	LH1 <- as.numeric(sapply(lmm, getLH1, USE.NAMES=FALSE))
	
	D <- sapply(LH1, get_deviance, lognull=lognull, USE.NAMES=FALSE)
	
	pvals <- apply(data.frame("D" = D, "num.alleles" = num.alleles), 1, get_pvals)
	
	negLog10 <- -log10(pvals)
	
	file.remove(lmm)
	file.remove(lmm.log)
	
	lmm <- data.frame("pattern" = paste0("pattern",1:length(XX)), "ps" = 1:length(XX),
					  "pvals" = pvals, "negLog10" = negLog10)
	write.table(lmm, file = paste0(prefix, "_tritetra_lmmout_patterns.txt"),
				sep="\t", row=F, col=T, quote=F)
						
	
	lmm <- lmm[pattern, ]
	lmm$ps <- ps
	if(any(is.na(lmm$pvals))){
		lmm <- lmm[-which(is.na(lmm$pvals)), ]
	}
	

	write.table(lmm, file = paste0(prefix, "_tritetra_lmmout_allSNPs.txt"),
				sep="\t", row=F, col=T, quote=F)
	
	return(lmm)
	
}

################################################################################################
## Plot manhattan plot ordered on the x-axis by PCs
## @o: PCs in order of significance
## @which.pc: For each pattern, which PC is it most correlated to
## @pattern: For each variant, what pattern it is
## @p.pca.bwt: Bayesian Wald test results for each PC
## @pc.lim: Which PCs are significant by the BWT
## @lmm: LMM results matrix
## @pat.weight: How many variants does each pattern represent
## @prefix: output prefix
## Outputs:
## Manhattan plot with the x-axis ordered by PCs
################################################################################################

plot_pc_manhattan <- function(o = NULL,
							  which.pc = NULL,
							  pattern = NULL,
							  p.pca.bwt = NULL,
							  pc.lim = NULL,
							  negLog10 = NULL,
							  pat.weight = NULL,
							  prefix = NULL,
							  randomCOL = NULL,
							  npcs = NULL){
	
	# Find which PCs have a variant that is most correlated to it
	m <- match(o,unique(which.pc))
	o.pats <- o[which(is.na(m)==FALSE)]
	# Subset the Bayesian Wald Test -log10(p) to those PCs which have kmers most correlated to them
	p.pca.bwt.pats <- p.pca.bwt[o.pats]
	
	
	num.variants <- sapply(o.pats, function(x, pat.weight=NULL, which.pc = NULL)
						   sum(pat.weight[which(which.pc==x)]),pat.weight = pat.weight,
						   which.pc = which.pc, USE.NAMES=FALSE)
	
	cols.pcs <- rep_len(c("#5a5a5a","#c6c6c6"), length.out = length(o.pats))
	if(!is.null(pc.lim)){
		cols.pcs[1:length(pc.lim)] <- randomCOL[1:length(pc.lim)]
	}
	m <- match(which.pc,o.pats)
	cols <- cols.pcs[m[pattern]]
	#cols[which(max.cor.pc[which(which.pc==o.pats[i])]<0.3)]="grey50"
	
	
	pos.gap <- sapply(num.variants, function(x) 10000/x, USE.NAMES=FALSE)
	pos.gap[num.variants > 10000] <- 1
	if(!is.null(pc.lim)){
		pos.gap[is.na(match(o.pats, o[pc.lim]))] <- pos.gap[is.na(match(o.pats, o[pc.lim]))]/ (npcs/10)
	} else {
		pos.gap[21:length(pos.gap)] <- pos.gap[21:length(pos.gap)] / (npcs/10)
	}
	
	pos <- rep(0,length(negLog10))
	max.pos <- 0
	plot.lines <- matrix(rep(0, length(o.pats)*2), ncol=2)
	
	which.pc <- which.pc[pattern]
	
	for(i in 1:length(o.pats)){
	
		s = seq(from = (max.pos+pos.gap[i]), by = pos.gap[i], length.out = num.variants[i])
		if(num.variants[i]>1){
			pos[which(which.pc==o.pats[i])] <- sample(s, num.variants[i], replace=FALSE)
		} else {
			pos[which(which.pc==o.pats[i])] <- s
		}
		max.pos <- max(s)
		plot.lines[i,] <- c((max.pos-pos.gap[i]*num.variants[i]+((num.variants[i]*0.15)*pos.gap[i])),
							(max.pos-((num.variants[i]*0.15)*pos.gap[i])))
		
	}

	pl2(paste0(prefix,"_PC_manhattan"),{
	plot(x = pos, y = negLog10, col = cols, xlab = "", ylab = "LMM -log10(p)", xaxt = "n")
	
	for(i in 1:length(o.pats)){
		lines(x = c(plot.lines[i, ]), y = c(p.pca.bwt.pats[i], p.pca.bwt.pats[i]),
			  type = "l", col = cols.pcs[i], lwd=2)
		if(i<=20){
			text(x = c(plot.lines[i, 1]+((plot.lines[i, 2]-plot.lines[i, 1])/2)),
				 y = c(p.pca.bwt.pats[i]+(max(negLog10)/90)), labels = o.pats[i],
				 col = cols.pcs[i], font = 2, cex = 0.9)
			}
		}
})
	
}


## To run for binary variants only
#plot_pc_manhattan(o = biallelic$pc_order$pc_order, which.pc = gen.var$cor.var$which.pc, pattern = gen.var$pattern, p.pca.bwt = biallelic$p.pca.bwt, pc.lim = biallelic$pc_order$pc.lim, negLog10 = gen.var$lmm$negLog10, pat.weight = gen.var$varpat, prefix=paste0(prefix,"_var"))


## To run for all SNPs (biallelic and tri and tetra allelic)
#new.pat <- triallelic$pattern
#new.pat <- sapply(new.pat, function(x, bippat) x + length(bippat), bippat = biallelic$bippat)
#new.pat <- c(biallelic$pattern, new.pat)

#plot_pc_manhattan(o = biallelic$pc_order$pc_order, 
#                  which.pc = c(biallelic$cor.XX$which.pc, triallelic$cor.tritetra$which.pc), 
#                  pattern = new.pat, 
#                  p.pca.bwt = biallelic$p.pca.bwt, 
#                  pc.lim = biallelic$pc_order$pc.lim, 
#                  negLog10 = c(biallelic$lmm$negLog10, triallelic$negLog10), 
#                  pat.weight = c(biallelic$bippat, triallelic$snppat), 
#                  prefix=paste0(prefix,"_SNPs"))




################################################################################################
## Functions for presence absence matrix
################################################################################################


get_gen_var <- function(var.matrix = NULL,
						XX.ID = NULL,
						logreg.var = NULL,
						lmm.var = NULL,
						relmatrix = NULL,
						pheno.file = NULL,
						maf = NULL,
						prefix = NULL,
						gem.path = NULL,
						output.dir = NULL,
						run.lmm = NULL,
						pca = NULL,
						npcs = NULL){

	genVarList = NULL
	genVarCount = NULL
	
	if(!is.null(var.matrix)){
		# Read in variant patterns (optional)
		# Run check if variant patterns are patterns
		genVarCount = length(var.matrix)
		genVarList = list()
		for(iGenVar in 1:genVarCount){
			genVarList[[iGenVar]] = list()
			gen.var <- read.table(var.matrix[iGenVar], header=T, sep="\t", check.names = F)
			if(any(colnames(gen.var)=="ps")){
				var_ps <- gen.var[,"ps"]
				gen.var <- gen.var[, -which(colnames(gen.var)=="ps")]
			} else {
				var_ps <- 1:nrow(gen.var)
			}
			var.ID <- colnames(gen.var)
			
			if(length(var.ID)!=length(XX.ID)){
				stop(paste0("\nNumber of samples in variant data is not equal to number of samples in phenotype data\n"))
			}
			if(any(is.na(match(XX.ID, var.ID)))){
				stop("\nError: Sample names of SNP data do not match the sample names of variant data\n")
			}
    
    		# If sample names are in a different order, match the ordering
    		gen.var <- gen.var[, match(XX.ID, var.ID)]
    		var.ID <- colnames(var)
    		# Compact to patterns
    		gen.var <- compact_variants(var = gen.var)
    		# Rescale variant patterns by how many variants they represent
    		XX.var <- rescale_variants(var = gen.var$var, varpat = gen.var$varpat)
    		genVarList[[iGenVar]][["gen.var"]] = gen.var
    		genVarList[[iGenVar]][["XX.var"]] = XX.var
    	}
	}
	
	if(!is.null(logreg.var)){
		for(iGenVar in 1:genVarCount){
			logregGenVar <- read.table(logreg.var, header=T, sep="\t", as.is=T)
			if(any(is.na(match(var_ps, logregGenVar$ps)))){
				stop("\nError: variant positions/IDs do not match between var.matrix and logreg.var\n")
			}
			pat1 <- match(unique(logregGenVar$ps), logregGenVar$ps)
			logregGenVar <- logregGenVar[pat1, ]
			o <- order(logregGenVar$ps, decreasing = F)
			logregGenVar <- logregGenVar[o, ]
			genVarList[[iGenVar]][["logreg.var"]] = logregGenVar
		}
	}

	if(!is.null(lmm.var)){
		
		for(iGenVar in 1:genVarCount){
			lmm.gen.var <- read.table(lmm.var[iGenVar], header=T, sep="\t", as.is=T)
			
			if(any(is.na(match(lmm.gen.var$ps,var_ps)))){
				stop("\nError: variant positions/IDs do not match between var.matrix and lmm.var\n")
			}
			pat1 <- match(unique(lmm.gen.var$ps), lmm.gen.var$ps)
			lmm.gen.var <- lmm.gen.var[pat1,]
			o <- order(lmm.gen.var$ps, decreasing = F)
			lmm.gen.var <- lmm.gen.var[o, ]
			
			genVarList[[iGenVar]][["lmm"]][["lmm"]] = lmm.gen.var
		}
	} else if(!is.null(var.matrix) & is.null(lmm.var) & run.lmm){
		for(iGenVar in 1:genVarCount){
			lmm.gen.var <- run_lmm_bi(XX = genVarList[[iGenVar]]$gen.var$var, relmatrix = relmatrix, 
									  pattern = genVarList[[iGenVar]]$gen.var$pattern, 
									  ps = var_ps, pheno.file = pheno.file, 
									  maf = maf, prefix = paste0(prefix, "_var"), 
									  path = gem.path, dir = output.dir)
			genVarList[[iGenVar]][["lmm"]] = lmm.gen.var
		}
	}
	
	if(!is.null(var.matrix)){
  		for(iGenVar in 1:genVarCount){
  			genVarList[[iGenVar]][["cor.var"]] <- get_correlations(XX = XX.var, pca = pca$x, npcs = npcs)
  		}
  		
  		genVars = list()
		for(iGenVar in 1:genVarCount){
			genVar = list("XX" = genVarList[[iGenVar]]$gen.var$var,
            		  	  "pattern" = genVarList[[iGenVar]]$gen.var$pattern,
            		  	  "cor.var" = genVarList[[iGenVar]]$cor.var,
            		  	  "logreg" = genVarList[[iGenVar]]$logreg.var,
            		  	  "lmm" = genVarList[[iGenVar]]$lmm$lmm,
            		  	  "varpat" = genVarList[[iGenVar]]$gen.var$varpat)
           	genVars[[iGenVar]] = genVar
		}
	}
	
	return(genVars)

}





################################################################################################
## Functions for getting tree info
################################################################################################


get_tree <- function(phylo = NULL,
					 prefix = NULL,
					 XX.ID = NULL,
					 pca = NULL,
					 npcs = NULL,
					 allBranchAndPCCor = FALSE){
	# Read in tree
	tree <- ape::read.tree(phylo)
	
	if(any(is.na(match(XX.ID, tree$tip.label)))){
		stop("\nError: Phylogeny sample names do not match the sample names of the SNP data\n")
	}
	
	tree <- phangorn::midpoint(tree)
	ape::write.tree(tree, paste0(prefix, "_midpointrooted_tree.txt"))
	tree <- ape::read.tree(paste0(prefix, "_midpointrooted_tree.txt"))
	treepat <- tree2patterns(tree = tree, tiporder = XX.ID)
	mtp <- treepat$pat
	message("Retrieve all correlations between branches and PCs: ",allBranchAndPCCor)
	cor.tree <- get_correlations(
		XX = mtp, pca = pca$x, npcs = npcs, id = XX.ID, all.cor = allBranchAndPCCor)
		
	if(allBranchAndPCCor){
		branchPCCorFileName = paste(prefix, "allBranchAndPCCor.txt", sep="_" )
		write.table(signif(cor.tree$all.cor.pc, digits = 3), branchPCCorFileName, 
		row.names = T, col.names = T, quote=F, sep="\t")
		ape::write.tree(treepat$labelled.tree, paste0(prefix, "_node_labelled_tree.txt"))
	}

	return(list("tree" = tree, "pattern" = treepat, "cor.tree" = cor.tree))
	
}

################################################################################################
## Functions for biallelic data
################################################################################################

get_biallelic <- function(logreg.bi = NULL,
						  XX.all = NULL,
						  XX = NULL,
						  lmm.bi = NULL,
						  lognull = NULL,
						  lambda = NULL,
						  relmatrix = NULL,
						  pheno.file = NULL,
						  maf = NULL,
						  gem.path = NULL,
						  output.dir = NULL,
						  prefix = NULL,
						  run.lmm = NULL,
						  XX.ID = NULL,
						  pca = NULL,
						  npcs = NULL){

	if(!is.null(logreg.bi)){
		logreg.bi <- read.table(logreg.bi, header=T, sep="\t", as.is=T)
		if(nrow(logreg.bi)!=length(XX.all$pattern)){
			cat(paste0("logreg.bi = ", nrow(logreg.bi), " variants"),"\n")
			cat(paste0("gen = ", length(XX.all$pattern), " biallelic variants"),"\n")
			stop("\nError: number of variants in logreg.bi does not match number of biallelic SNPs in gen")
		}
		if(any(is.na(match(logreg.bi$ps, XX.all$ps)))){
			stop("\nError: variant positions/IDs do not match between logreg.bi and biallelic variants in gen\n")
		}
		if(any(logreg.bi$ps != XX.all$ps)){
			m <- match(XX.all$ps, logreg.bi$ps)
			logreg.bi <- logreg.bi[m,]
		}
	}
	
	
	if(!is.null(lmm.bi)){
		lmm.bi <- read.table(lmm.bi, header=T, sep="\t", as.is=T)
		if(nrow(lmm.bi)!=length(XX.all$pattern)){
			cat(paste0("lmm.bi = ", nrow(lmm.bi), " variants"),"\n")
			cat(paste0("gen = ", length(XX.all$pattern), " biallelic variants"),"\n")
			stop("\nError: number of variants in lmm.bi does not match number of biallelic SNPs in gen")
		}
		if(any(is.na(match(lmm.bi$ps, XX.all$ps)))){
			stop("\nError: variant positions/IDs do not match between lmm.bi and biallelic variants in gen\n")
		}
		if(any(lmm.bi$ps != XX.all$ps)){
			m <- match(XX.all$ps, lmm.bi$ps)
			lmm.bi <- lmm.bi[m,]
		}
		
		if(is.null(lognull) | is.null(lambda)){
			lambda.lognull <- run_lmm_bi(XX = XX.all$XX[1,], relmatrix = relmatrix,
										 pheno.file = pheno.file, maf = maf,
										 prefix = paste0(prefix, "_getlognull"), path = gem.path,
										 dir = output.dir, process.results = FALSE)
			if(is.null(lambda)){
				lambda <- as.numeric(lambda.lognull$lambda)
				cat(paste0("## Lambda = ", lambda), file = paste0(prefix,"_logfile.txt"),
					sep="\n", append = TRUE)
    		}
    		if(is.null(lognull)){
    			lognull <- as.numeric(lambda.lognull$lognull)
    			cat(paste0("## Log likelihood under the null = ", lognull),
    				file = paste0(prefix,"_logfile.txt"), sep="\n", append = TRUE)
    		}
		}
		
	} else if(is.null(lmm.bi) & run.lmm){	
		lmm.bi <- run_lmm_bi(XX = XX.all$XX, relmatrix = relmatrix, pattern = XX.all$pattern,
							 ps = XX.all$ps, pheno.file = pheno.file, maf = maf,
							 prefix = paste0(prefix, "_biallelic"), path = gem.path, dir = output.dir)
							 
		lognull <- as.numeric(lmm.bi$lognull)
		cat(paste0("## Log likelihood under the null = ", lognull),
			file = paste0(prefix,"_logfile.txt"), sep="\n", append = TRUE)
		lambda <- as.numeric(lmm.bi$lambda)
		cat(paste0("## Lambda = ", lambda), file = paste0(prefix,"_logfile.txt"),
			sep="\n", append = TRUE) 
		lmm.bi <- lmm.bi$lmm
	}


	cor.XX <- get_correlations(XX = XX, pca = pca$x, npcs = npcs, id = XX.ID)

	return(list("logreg.bi" = logreg.bi, "lmm.bi" = lmm.bi, "lognull" = lognull,
				"lambda" = lambda, "cor.XX" = cor.XX))

}


################################################################################################
## Functions for tri/tetra allelic data
################################################################################################


get_tritetra <- function(logreg.tri.tetra = NULL,
						 XX.all = NULL,
						 lmm.tri.tetra = NULL,
						 run.lmm = NULL,
						 relmatrix = NULL,
						 pheno.file = NULL,
						 lognull = NULL,
						 prefix = NULL,
						 gem.path = NULL,
						 output.dir = NULL,
						 pca = NULL,
						 npcs = NULL,
						 XX.ID = NULL,
						 maf = NULL){
	
	
	if(!is.null(logreg.tri.tetra)){
		logreg.tri.tetra <- read.table(logreg.tri.tetra, header=T, sep="\t", as.is=T)
		if(nrow(logreg.tri.tetra)!=length(XX.all$pattern.snps)){
			cat(paste0("logreg.tri.tetra = ", nrow(logreg.tri.tetra), " variants"),"\n")
			cat(paste0("gen = ", length(XX.all$pattern.snps), " tri or tetra allelic variants"),"\n")
			stop("\nError: number of variants in logreg.tri.tetra does not match number of tri or tetra allelic SNPs in gen")
		}
		if(any(is.na(match(logreg.tri.tetra$ps, XX.all$ps.snps)))){
			stop("\nError: variant positions/IDs do not match between logreg.tri.tetra and tri and tetra allelic variants in gen\n")
		}
		if(any(logreg.tri.tetra$ps != XX.all$ps.snps)){
			m <- match(XX.all$ps.snps, logreg.tri.tetra$ps)
			logreg.tri.tetra <- logreg.tri.tetra[m,]
		}
	}
		
	if(!is.null(lmm.tri.tetra)){
		lmm.tri.tetra <- read.table(lmm.tri.tetra, header=T, sep="\t")
		if(nrow(lmm.tri.tetra)!=length(XX.all$pattern.snps)){
			cat(paste0("lmm.tri.tetra = ", nrow(lmm.tri.tetra), " variants"),"\n")
			cat(paste0("gen = ", length(XX.all$pattern.snps), " tri or tetra allelic variants"),"\n")
			stop("\nError: number of variants in lmm.tri.tetra does not match number of tri or tetra allelic SNPs in gen")
		}
		if(any(is.na(match(lmm.tri.tetra$ps, XX.all$ps.snps)))){
			stop("\nError: variant positions/IDs do not match between lmm.tri.tetra and tri and tetra allelic variants in gen\n")
		}
		if(any(lmm.tri.tetra$ps != XX.all$ps.snps)){
			m <- match(XX.all$ps.snps, lmm.tri.tetra$ps)
			lmm.tri.tetra <- lmm.tri.tetra[m,]
		}
	} else if(!is.null(XX.all$XX.tritetra) & run.lmm){
		lmm.tri.tetra <- run_lmm_multi(XX = XX.all$XX.tritetra, pattern = XX.all$pattern.snps,
                                 	   relmatrix = relmatrix, pheno.file = pheno.file, maf = maf,
                                 	   ps = XX.all$ps.snps, lognull = lognull,
                                 	   prefix = paste0(prefix, "_tritetra"), path = gem.path,
                                 	   dir = output.dir)
	}

	if(!is.null(XX.all$XX.tritetra)){
	  XX.tritetra <- rescale_variants(var = XX.all$XX.tritetra, varpat = XX.all$snppat)
	  cor.tritetra <- get_correlations(XX = XX.tritetra, pca = pca$x, npcs = npcs,
	  								   id = XX.ID)
	} else {
		cor.tritetra = NULL
	}
	
	return(list("logreg.tri.tetra" = logreg.tri.tetra, "lmm.tri.tetra" = lmm.tri.tetra,
				"cor.tritetra" = cor.tritetra))	
	
}



################################################################################################
## Do PCA
################################################################################################

do_pca <- function(pcs = NULL, XX = NULL, XX.ID = NULL){
	if(is.null(pcs)){
		pca <- prcomp(XX)
	} else {
		## Read in PCs
		pca <- read.table(pcs, header = T, as.is = T)
		if(any(XX.ID != rownames(pca))){
			m <- match(XX.ID, rownames(pca))
			pca <- pca[m, ]
		}
		pca <- list("x" = pca, "rotation" = NULL)
	}
	return(list("pca" = pca))
}




################################################################################################
## Do Wald test
################################################################################################

wald_test <- function(y = NULL,
					  XX = NULL,
					  svd.XX = NULL,
					  lambda = NULL,
					  XX.all = NULL,
					  prefix = NULL,
					  npcs = NULL,
					  pca = NULL){
	
	
	fit.lmm <- ridge_regression(y, XX, svdX=svd.XX,
                            lambda_init=as.numeric(lambda)/sum(XX.all$bippat),
                            maximize=FALSE, skip.var=TRUE)
	
	# Fit the grand null model
	
	fit.0 <- lm(y~1)
	
	# LRT for the LMM null vs grand null
	LRTnullVgrand <- -log10(pchisq(2*(fit.lmm$ML - as.numeric(logLik(fit.0))), 1, low=F)/2)
	cat(paste0("## LRT for the LMM null vs grand null = ", LRTnullVgrand),
    	file = paste0(prefix, "_logfile.txt"), sep="\n", append = TRUE)
	
	# Heritability
	fit.lmm.ypred <- XX %*% fit.lmm$Ebeta
	cat(paste0("## Heritability (R^2) = ", cor(fit.lmm.ypred,y)^2),
    	file=paste0(prefix, "_logfile.txt"), sep="\n", append=TRUE)
	
	# Get full posterior covariance matrix for Bayesian Wald Test
	# Need the full posterior covariance matrix for the Bayesian Wald test,
	# to get the posterior uncertainty for each point estimate
	
	wald_input <- get_wald_input(fit.lmm = fit.lmm, pca = pca, svd.XX = svd.XX,
								 y = y, npcs = npcs, XX = XX)

	# Bayesian Wald Test
	pca.bwt <- wald_input$Ebeta^2/diag(wald_input$Vbeta)
	
	p.pca.bwt <- -log10(exp(1))*pchisq(pca.bwt, 1, low=F, log=T)
	
	cat(paste0("## Bayesian Wald Test for PCs range = ", paste(range(p.pca.bwt), collapse=" ")),
    	file=paste0(prefix, "_logfile.txt"), sep="\n", append=TRUE)
	
	write.table(p.pca.bwt, file = paste0(prefix, "_Bayesian_Wald_Test_negLog10.txt"),
    	        sep="\t", row=T, col = F, quote=F)
	
	# Get order of PCs by Wald test results
	signif_cutoff <- -log10(0.05/npcs)
	
	pc_order <- get_pc_order(p.pca.bwt = p.pca.bwt, signif_cutoff = signif_cutoff)
	# Predict phenotype using effect sizes
	effect <- t(t(XX) * as.vector(fit.lmm$Ebeta))
	pred <- rowSums(effect)
	
	
	return(list("pc_order" = pc_order, "p.pca.bwt" = p.pca.bwt, "pred" = pred,
				"signif_cutoff" = signif_cutoff))
	
}


################################################################################################
## Get SNP data
################################################################################################

get_SNP_data <- function(gen = NULL,
						 pheno = NULL,
						 prefix = NULL){
	
	
	# Read in SNP data

	XX <- read.table(gen,header=T,sep="\t",as.is=T, check.names = F)
	
	if(any(colnames(XX)=="ps")){
		XX_ps <- XX[,"ps"]
		XX <- XX[, -which(colnames(XX)=="ps")]
	} else {
		XX_ps <- 1:nrow(XX)
	}
	XX.ID <- colnames(XX)
	
	# Read in sample IDs and phenotype
	pheno <- read.table(pheno,header=T, as.is = T, sep="\t")
	if(!any(colnames(pheno)=="ID")){
		stop("\nError: phenotype file must have 'ID' column")
	}
	if(!any(colnames(pheno)=="pheno")){
		stop("\nError: phenotype file must have 'pheno' column")
	}
	
	sample_ID <- pheno$ID
	npcs <-length(sample_ID)
	pheno <- pheno$pheno
	y <- pheno[match(XX.ID,sample_ID)]
	
	if(length(y)!=length(XX.ID)){
		stop(paste0("\nNumber of phenotypes is not equal to number of samples in SNP genotypes\n"))
	}
	
	if(any(is.na(match(XX.ID, sample_ID)))){
		stop("\nError: Sample names of SNP data do not match the sample names of the phenotype data\n")
	}
	
	
	# Run compact_SNPs function
	XX.all <- compact_SNPs(gen = XX, ps = XX_ps, prefix = prefix)
	
	return(list("XX.all" = XX.all, "sample_ID" = sample_ID,
				"npcs" = npcs, "y" = y, "XX.ID" = XX.ID))	
}




################################################################################################
## Produce log file
################################################################################################

get_log_file <- function(XX.all = NULL, prefix = NULL){
	# Produce log file
	cat("##", file = paste0(prefix,"_logfile.txt"), sep="\n")
	cat("## BUGWAS", file = paste0(prefix,"_logfile.txt"), sep="\n", append = TRUE)
	cat("##", file = paste0(prefix,"_logfile.txt"), sep="\n", append = TRUE)
	cat("## Command line input: ", file = paste0(prefix,"_logfile.txt"), sep="\n", append = TRUE)
	cat("##", file = paste0(prefix,"_logfile.txt"), sep="\n", append = TRUE)
	cat(paste0("## # Biallelic sites = ", length(XX.all$pattern)),
		file = paste0(prefix,"_logfile.txt"), sep="\n", append = TRUE)
	cat(paste0("## # Tri allelic sites = ", XX.all$n.triallelic),
		file = paste0(prefix,"_logfile.txt"), sep="\n", append = TRUE)
	cat(paste0("## # Tetra allelic sites = ", XX.all$n.tetraallelic),
		file = paste0(prefix,"_logfile.txt"), sep="\n", append = TRUE)
	cat("##", file = paste0(prefix,"_logfile.txt"), sep="\n", append = TRUE)

}

##################################################################################
## Extract input arguments.
## @argName: Name of the argument to be extracted
## @commandLineInputs: A n x 2 matrix, where n is the number of command line inputs.
## @default: default value of the input argument if not provided
## The first column of the matrix contains the names of the arguments, and the second
## column contains the arguments values.
##################################################################################
extractInputArgument = function(arg = NULL, name = NULL, default = NULL, canBeNULL = FALSE, checkExist = FALSE){
  
  if(is.null(arg)){
    arg = default
  }
  
  if(!canBeNULL & is.null(arg)){
    stop(paste0(name, " variable cannot be null!"))
    
  }
  
  if(!is.null(arg) & checkExist){
    checkExistence(arg)
  }
  
  return(arg)
}

checkExistence = function(filePath = NULL){
  
  doesNotExist = which(!file.exists(filePath))

  if(length(doesNotExist) > 0){
    stop(paste("The file", filePath[doesNotExist]," does not exist!", collapse="\n", sep=""))
  }
  
}















