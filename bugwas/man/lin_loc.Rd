% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BUGWAS_modular.R
\name{lin_loc}
\alias{lin_loc}
\title{Lineage and locus tests for bacterial GWAS}
\usage{
lin_loc(gen = NULL, pheno = NULL, phylo = NULL, prefix = NULL,
  gem.path = NULL, pcs = NULL, lmm.bi = NULL, lmm.tri.tetra = NULL,
  logreg.bi = NULL, logreg.tri.tetra = NULL, var.matrix = NULL,
  logreg.var = NULL, lmm.var = NULL, cutOffCor = 0, run.lmm = TRUE,
  maf = 0, relmatrix = NULL, lognull = NULL, lambda = NULL,
  output.dir = getwd(), creatingAllPlots = TRUE,
  allBranchAndPCCor = FALSE, runTriTetrallelic = TRUE)
}
\arguments{
\item{gen}{A file name specified by either a variable of mode character, or a double-quoted string, containing imputed haploid SNP data. Rows are SNPs, and columns are samples, with the first column being SNP positions. Column headers must contain 'ps' for the SNP positions with the others being the sample names. This must contain biallelic SNPs, but can also contain tri and tetra-allelic SNPs. Required argument.}

\item{pheno}{A file name specified by either a variable of mode character, or a double-quoted string, containing a column of sample names with header 'ID' and a column of the binary phenotype (coded by 0s and 1s) with column header 'pheno'. Required argument.}

\item{phylo}{A file name specified by either a variable of mode character, or a double-quoted string, containing a phylogeny of the samples, with the same names matching with arguments gen and pheno. Required argument.}

\item{prefix}{Output file prefix. Required argument.}

\item{gem.path}{A file path specified by either a variable of mode character, or a double-quoted string. gem.path is the file path to the software GEMMA (version >= ?). Required argument.}

\item{pcs}{A file name specified by either a variable of mode character, or a double-quoted string, containing the principle components of the data. Column names should be 'PC1' to 'PCn' and row names should be the sample names.}

\item{lmm.bi}{A file name specified by either a variable of mode character, or a double-quoted string, containing GEMMA results (ending '.assoc.txt') for the biallelic SNPs in argument 'gen'.}

\item{lmm.tri.tetra}{A file name specified by either a variable of mode character, or a double-quoted string, containing GEMMA results for the tri and tetra allelic SNPs in argument 'gen'. This must contain column headers 'ps' for SNP positions/IDs, 'pvals' for p-values and 'negLog10' for -log10(p).}

\item{logreg.bi}{A file name specified by either a variable of mode character, or a double-quoted string, containing logistic regression -log10(p) for the biallelic SNPs with column names 'ps' for SNP positions/IDs and 'negLog10' for -log10(p).}

\item{logreg.tri.tetra}{A file name specified by either a variable of mode character, or a double-quoted string, containing logistic regression -log10(p) for the tri and tetra allelic SNPs with column names 'ps' for SNP positions/IDs and 'negLog10' for -log10(p).}

\item{var.matrix}{A vector of file names specified by double-quoted strings. The files should contain presence absence matrices, with rows being variants (of 0s and 1s), and columns being samples. Column headers must contain 'ps' for variant positions/IDs with the others being the sample names.}

\item{logreg.var}{A vector of file names specified by double-quoted strings, of files containing logistic regression -log10(p-value) results for the presence absence matrices. Column names must contain 'ps' for variant positions/IDs and 'negLog10' for the -log10(p).}

\item{lmm.var}{A vector of file names specified by double-quoted strings, of files containing GEMMA results for the presence absence matrices.}

\item{cutOffCor}{Correlation cut-off for assigning and colouring variants by Principal Components (Default = 0, variants are coloured by the PC they are most correlated with).}

\item{run.lmm}{Whether to run GEMMA (Default = TRUE).}

\item{maf}{Minor allele frequency for GEMMA (Default = 0, all varaints are tested).}

\item{relmatrix}{A file name specified by either a variable of mode character, or a double-quoted string of a file containing the GEMMA relatedness matrix of the samples created from biallelic SNPs. The individual ordering must be in the same order as the column names in argument 'gen'.}

\item{lognull}{The log likelihood under the null from GEMMA.}

\item{lambda}{Lambda from GEMMA.}

\item{output.dir}{Output file directory.}

\item{creatingAllPlots}{Whether to create all bugwas plots. Default = TRUE.}

\item{allBranchAndPCCor}{Whether or not to retreive correlation matrix between branches and PCs. Default = FALSE.}
}
\description{
This function tests for locus effects using GEMMA and lineage effects using a bayesian wald test for haploid data
}
\examples{
lin_loc()
## An example of running lin_loc with the minimum required inputs
## Assuming gemma is installed in the present working directory
gen <- system.file("extdata", "gen.txt", package = "bugwas")
pheno <- system.file("extdata", "pheno.txt", package = "bugwas")
phylo <- system.file("extdata", "tree.txt", package = "bugwas")
prefix <- "test_bugwas"
gem.path <- "./gemma"
data <- lin_loc(gen = gen, pheno = pheno, phylo = phylo, prefix = prefix, gem.path = gem.path)
}
\keyword{GEMMA}
\keyword{GWAS}
\keyword{bacteria}
\keyword{lineage}
\keyword{locus}
\keyword{wald}

