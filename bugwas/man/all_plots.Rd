% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/all_plots.R
\name{all_plots}
\alias{all_plots}
\title{Generates all plots.}
\usage{
all_plots(biallelic = NULL, triallelic = NULL, genVars = NULL,
  treeInfo = NULL, config = NULL)
}
\arguments{
\item{biallelic}{A list called 'biallelic' created from the lin_loc function}

\item{triallelic}{A list called 'triallelic' created from the lin_loc function}

\item{genVars}{A list called 'genVars' created from the lin_loc function}

\item{treeInfo}{A list called 'treeInfo' created from the lin_loc function}

\item{config}{A list called 'config' created from the lin_loc function}
}
\description{
This function generates all the plots
}
\examples{
data <- lin_loc(gen = gen, pheno = pheno, phylo = phylo, 
 prefix = prefix, gem.path = gem.path)
all_plots(biallelic = data$biallelic, triallelic = data$triallelic, 
	genVars = data$genVars, treeInfo = data$treeInfo, config = data$config)
}
\keyword{plot}

