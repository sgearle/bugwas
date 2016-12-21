# File BUGWAS_modular.R
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
#### 									BUGWAS				   							   #####
################################################################################################


################################################################################################
## Library dependencies
################################################################################################
# library(ape)
# library(phangorn)

# source("/dipro/mmm/ana/Saur/AbxGWAS/bugwas/BUGWAS_functions.R")
# source("/home/wilson/R/ridge regression.R")
# source("/dipro/mmm/ana/Saur/AbxGWAS/bugwas/all_plots_new.R")
# source("/dipro/mmm/ana/Saur/AbxGWAS/bugwas/jessie_plots_code.R")
                

#' Lineage and locus tests for bacterial GWAS
#'
#' This function tests for locus effects using GEMMA and lineage effects using a bayesian wald test for haploid data
#' @param gen A file name specified by either a variable of mode character, or a double-quoted string, containing imputed haploid SNP data. Rows are SNPs, and columns are samples, with the first column being SNP positions. Column headers must contain 'ps' for the SNP positions with the others being the sample names. This must contain biallelic SNPs, but can also contain tri and tetra-allelic SNPs. Required argument.
#' @param pheno A file name specified by either a variable of mode character, or a double-quoted string, containing a column of sample names with header 'ID' and a column of the binary phenotype (coded by 0s and 1s) with column header 'pheno'. Required argument.
#' @param phylo A file name specified by either a variable of mode character, or a double-quoted string, containing a phylogeny of the samples, with the same names matching with arguments gen and pheno. Required argument.
#' @param prefix Output file prefix. Required argument.
#' @param gem.path A file path specified by either a variable of mode character, or a double-quoted string. gem.path is the file path to the software GEMMA (version >= ?). Required argument.
#' @param pcs A file name specified by either a variable of mode character, or a double-quoted string, containing the principle components of the data. Column names should be 'PC1' to 'PCn' and row names should be the sample names.
#' @param lmm.bi A file name specified by either a variable of mode character, or a double-quoted string, containing GEMMA results (ending '.assoc.txt') for the biallelic SNPs in argument 'gen'.
#' @param lmm.tri.tetra A file name specified by either a variable of mode character, or a double-quoted string, containing GEMMA results for the tri and tetra allelic SNPs in argument 'gen'. This must contain column headers 'ps' for SNP positions/IDs, 'pvals' for p-values and 'negLog10' for -log10(p).
#' @param logreg.bi A file name specified by either a variable of mode character, or a double-quoted string, containing logistic regression -log10(p) for the biallelic SNPs with column names 'ps' for SNP positions/IDs and 'negLog10' for -log10(p).
#' @param logreg.tri.tetra A file name specified by either a variable of mode character, or a double-quoted string, containing logistic regression -log10(p) for the tri and tetra allelic SNPs with column names 'ps' for SNP positions/IDs and 'negLog10' for -log10(p).
#' @param var.matrix A vector of file names specified by double-quoted strings. The files should contain presence absence matrices, with rows being variants (of 0s and 1s), and columns being samples. Column headers must contain 'ps' for variant positions/IDs with the others being the sample names.
#' @param logreg.var A vector of file names specified by double-quoted strings, of files containing logistic regression -log10(p-value) results for the presence absence matrices. Column names must contain 'ps' for variant positions/IDs and 'negLog10' for the -log10(p).
#' @param lmm.var A vector of file names specified by double-quoted strings, of files containing GEMMA results for the presence absence matrices.
#' @param cutOffCor Correlation cut-off for assigning and colouring variants by Principal Components (Default = 0, variants are coloured by the PC they are most correlated with).
#' @param run.lmm Whether to run GEMMA (Default = TRUE).
#' @param maf Minor allele frequency for GEMMA (Default = 0, all varaints are tested).
#' @param relmatrix A file name specified by either a variable of mode character, or a double-quoted string of a file containing the GEMMA relatedness matrix of the samples created from biallelic SNPs. The individual ordering must be in the same order as the column names in argument 'gen'.
#' @param lognull The log likelihood under the null from GEMMA.
#' @param lambda Lambda from GEMMA.
#' @param output.dir Output file directory.
#' @param creatingAllPlots Whether to create all bugwas plots. Default = TRUE.
#' @param allBranchAndPCCor Whether or not to retreive correlation matrix between branches and PCs. Default = FALSE.
#' @keywords bacteria GWAS locus lineage wald GEMMA
#' @export
#' @examples
#' lin_loc()
#' ## An example of running lin_loc with the minimum required inputs
#' ## Assuming gemma is installed in the present working directory
#' gen <- system.file("extdata", "gen.txt", package = "bugwas")
#' pheno <- system.file("extdata", "pheno.txt", package = "bugwas")
#' phylo <- system.file("extdata", "tree.txt", package = "bugwas")
#' prefix <- "test_bugwas"
#' gem.path <- "./gemma"
#' data <- lin_loc(gen = gen, pheno = pheno, phylo = phylo, prefix = prefix, gem.path = gem.path)


lin_loc <- function(gen = NULL,
                   pheno = NULL,
                   phylo = NULL,
                   prefix = NULL,
                   gem.path = NULL,
                   pcs = NULL,
                   lmm.bi = NULL,
                   lmm.tri.tetra = NULL,
                   logreg.bi = NULL,
                   logreg.tri.tetra = NULL,
                   var.matrix = NULL,
                   logreg.var = NULL,
                   lmm.var = NULL,
                   cutOffCor = 0,
                   run.lmm = TRUE,
                   maf = 0,
                   relmatrix = NULL,
                   lognull = NULL,
                   lambda = NULL,
                   output.dir = getwd(),
                   creatingAllPlots = TRUE,
                   allBranchAndPCCor = FALSE,
                   runTriTetrallelic = TRUE){
    
  gen = extractInputArgument(arg = gen, checkExist = TRUE)
  pheno = extractInputArgument(arg = pheno, checkExist = TRUE)
  phylo = extractInputArgument(arg = phylo, checkExist = TRUE)
  prefix = extractInputArgument(arg = prefix)
  gem.path = extractInputArgument(arg = gem.path, checkExist = TRUE)
  pcs = extractInputArgument(arg = pcs, canBeNULL = TRUE)
  lmm.bi = extractInputArgument(arg = lmm.bi, canBeNULL = TRUE, checkExist = TRUE)
  lmm.tri.tetra = extractInputArgument(arg = lmm.tri.tetra, canBeNULL = TRUE, checkExist = TRUE)
  logreg.bi = extractInputArgument(arg = logreg.bi, canBeNULL = TRUE, checkExist = TRUE)
  logreg.tri.tetra = extractInputArgument(arg = logreg.tri.tetra, canBeNULL = TRUE, checkExist = TRUE)
  var.matrix = extractInputArgument(arg = var.matrix, canBeNULL = TRUE, checkExist = TRUE)
  logreg.var = extractInputArgument(arg = logreg.var, canBeNULL = TRUE, checkExist = TRUE)
  lmm.var = extractInputArgument(arg = lmm.var, canBeNULL = TRUE, checkExist = TRUE)
  cutOffCor = extractInputArgument(arg = cutOffCor, default = 0)
  run.lmm = extractInputArgument(arg = run.lmm, default = TRUE)
  maf = extractInputArgument(arg = maf, default = 0)
  relmatrix = extractInputArgument(arg = relmatrix, canBeNULL = TRUE, checkExist = TRUE)
  lognull = extractInputArgument(arg = lognull, canBeNULL = TRUE)
  lambda = extractInputArgument(arg = lambda, canBeNULL = TRUE)
  output.dir = extractInputArgument(arg = output.dir, default = getwd())
  creatingAllPlots = extractInputArgument(arg = creatingAllPlots, default = TRUE)
  allBranchAndPCCor = extractInputArgument(arg = allBranchAndPCCor, default = FALSE)
  runTriTetrallelic = extractInputArgument(arg = runTriTetrallelic, default = TRUE)

	SNPdata <- get_SNP_data(gen = gen, pheno = pheno, prefix = prefix)
	XX.all <- SNPdata$XX.all
	sample_ID <- SNPdata$sample_ID
	npcs <- SNPdata$npcs
	y <- SNPdata$y
	XX.ID <- SNPdata$XX.ID
	rm(SNPdata)
	
	get_log_file(XX.all = XX.all, prefix = prefix)
	
	# Read in LMM results, if no LMM results to read in, run LMM
	if(run.lmm){
		pheno.file <- write_pheno(pheno = y, prefix = prefix)
		if(is.null(relmatrix)){
			message("Calculating kinship matrix.")
			relmatrix <- get_kinship(XX = XX.all$XX,
							 	 	 pattern = XX.all$pattern,
							 	 	 prefix = prefix,
								 	 path = gem.path,
								 	 dir = output.dir,
							 	 	 maf = maf,
                           	 	 	 pheno.file = pheno.file)
            message("Kinship matrix calculated successfully.")
		}
	}
	
	
	XX <- rescale_variants(var = XX.all$XX, varpat = XX.all$bippat)
	message("Rescaled variants.")
	
	svd.XX <- svd(XX)
	message("Single value decomposition complete.")

	# PCA on the bips
	pca <- do_pca(pcs = pcs, XX = XX, XX.ID = XX.ID)
	message("Principle component analysis complete.")
	
	biallelic <- get_biallelic(logreg.bi = logreg.bi,
		 					   XX.all = XX.all,
						  	   XX = XX,
						  	   lmm.bi = lmm.bi,
						  	   lognull = lognull,
						  	   lambda = lambda,
						  	   relmatrix = relmatrix,
						  	   pheno.file = pheno.file,
						  	   maf = maf,
						  	   gem.path = gem.path,
						  	   output.dir = output.dir,
						  	   prefix = prefix,
						  	   run.lmm = run.lmm,
						  	   XX.ID = XX.ID,
						  	   pca = pca$pca,
						  	   npcs = npcs)
	message("Biallelic data processed successfully.")

	if(!is.null(XX.all$XX.tritetra) & runTriTetrallelic){						  
		tritetra <- get_tritetra(logreg.tri.tetra = logreg.tri.tetra,
						 	 	XX.all = XX.all,
						 	 	lmm.tri.tetra = lmm.tri.tetra,
						 	 	run.lmm = run.lmm,
						 	 	relmatrix = relmatrix,
						 	 	pheno.file = pheno.file,
						 	 	lognull = biallelic$lognull,
						 	 	prefix = prefix,
						 	 	gem.path = gem.path,
						 	 	output.dir = output.dir,
						 	 	pca = pca$pca,
						 	 	npcs = npcs,
						 	 	XX.ID = XX.ID,
						 	 	maf = maf)
		message("tri/tetrallelic data processed successfully.")
	} else {
		tritetra = NULL
	}

	# Get list containing all presence/absence matrix info
	if(!is.null(var.matrix)){
		genVars <- get_gen_var(var.matrix = var.matrix,
					  	  XX.ID = XX.ID,
					  	  logreg.var = logreg.var,
					  	  lmm.var = lmm.var,
					  	  relmatrix = relmatrix,
					  	  pheno.file = pheno.file,
					  	  maf = maf,
					  	  prefix = prefix,
					  	  gem.path = gem.path,
					  	  output.dir = output.dir,
					  	  run.lmm = run.lmm,
					  	  pca = pca$pca,
					  	  npcs = npcs)
		message("General binary data processed successfully.")
	} else {
		genVars = NULL
	}

	# Get list of all tree info
	treeInfo <- get_tree(phylo = phylo,
						prefix = prefix,
						XX.ID = XX.ID,
						pca = pca$pca,
						npcs = npcs,
						allBranchAndPCCor = allBranchAndPCCor)
	message("Tree data processed successfully.")


	# Ridge regression

	wald <- wald_test(y = y,
					  XX = XX,
					  svd.XX = svd.XX,
					  lambda = biallelic$lambda,
					  XX.all = XX.all,
					  prefix = prefix,
					  npcs = npcs,
					  pca = pca$pca)
					  
	rm(list=c("XX", "svd.XX"))

	# Put everything into lists
	biallelic <- list("pattern" = XX.all$pattern,
                 	 "cor.XX" = biallelic$cor.XX,
                 	 "npcs" = npcs,
                 	 "pheno" = y,
                 	 "logreg" = biallelic$logreg.bi,
                 	 "lmm" = biallelic$lmm.bi,
                 	 "ps" = XX.all$ps,
                 	 "pred" = wald$pred,
                 	 "pca" = pca$pca,
                 	 "pc_order" = wald$pc_order,
                 	 "p.pca.bwt" = wald$p.pca.bwt,
                 	 "bippat" = XX.all$bippat,
                 	 "id" = XX.ID) # New  - Jessie need to add

	triallelic <- list("pattern" = XX.all$pattern.snps,
                  	  "cor.tritetra" = tritetra$cor.tritetra,
                  	  "logreg" = tritetra$logreg.tri.tetra,
                  	  "lmm" = tritetra$lmm.tri.tetra,
                  	  "pattern" = XX.all$pattern.snps,
                  	  "snppat" = XX.all$snppat) # New  - Jessie need to add

	config <- list("prefix" = prefix,
			  	  "signif_cutoff" = wald$signif_cutoff,
			  	  "cutoffCor" = cutOffCor)

	if(creatingAllPlots){
		all_plots(biallelic = biallelic, triallelic = triallelic,
		  	  genVars = genVars, treeInfo = treeInfo, config = config)
	}
		  	  
	return(list("biallelic" = biallelic, "triallelic" = triallelic, "config" = config, "treeInfo" = treeInfo, "genVars" = genVars))

}


# all_plots(biallelic = biallelic, triallelic = triallelic, config = config, treeInfo = treeInfo, genVars = genVars)
