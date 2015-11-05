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

#' Plot the a Manhattan plot organised by the signficance of PCs.
#' 
#' This function generates a Manhattan plot organised by the signficance of the principal components.
#' @param biallelic A list called 'biallelic' created from the lin_loc function. It is a required input.
#' @param triallelic A list called 'triallelic' created from the lin_loc function. It is a required input.
#' @param config A list called 'config' created from the lin_loc function. It is a required input.
#' @param colourPalette A vector of colours colour the significant principal components identified by the Bayesian Wald test (see testGenomeWidePCs). If this is NULL then colours are chosen from a default colour palette.
#' @param p.genomewidepc A matrix of the significant principal component and their correlation with lineages. This is Bayesian Wald test results produced by the function testGenomeWidePCs. If this is NULL then testGenomeWidePCs is called to generate the required test results. 
#' @keywords SNP
#' @keywords PCA
#' @keywords Manhattan-plot
#' @export
#' @examples
#' data <- lin_loc(gen = gen, pheno = pheno, phylo = phylo, 
#'  prefix = prefix, gem.path = gem.path)
#' snpPCManhattanPlot(biallelic = data$biallelic,  triallelic = data$triallelic,
#'  config = data$config)
snpPCManhattanPlot = function(config, biallelic, triallelic, colourPalette = NULL){
							 	
	if(is.null(colourPalette)){
		colourPalette = getColourPalette(p.pca.bwt = biallelic$p.pca.bwt, 
										signifCutOff = config$signif_cutoff, 
										pc.lim = biallelic$pc_order$pc.lim)	
	}
	## To run for all SNPs (biallelic and tri and tetra allelic)
	
	ipat.snps = triallelic$pattern
	bippat = biallelic$bippat
	ipat = biallelic$pattern
	
	new.pat <- ipat.snps
	new.pat <- sapply(new.pat, function(x, pat){ x + length(pat)}, pat = bippat)
	new.pat <- c(ipat, unlist(new.pat))
	
  
  .plot_pc_manhattan(o = biallelic$pc_order$pc_order, 
                   which.pc = c(biallelic$cor.XX$which.pc, triallelic$cor.tritetra$which.pc), 
                   pattern = new.pat, 
                   p.pca.bwt = biallelic$p.pca.bwt, 
                   pc.lim = biallelic$pc_order$pc.lim, 
                   negLog10 = c(biallelic$lmm$negLog10, triallelic$lmm$negLog10), 
                   pat.weight = c(biallelic$bippat, triallelic$snppat), 
                   prefix=paste0(config$prefix,"_SNPs"),
                   colourPalette = colourPalette,
                   npcs = biallelic$npcs)

}