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

#' Generates Manhattan plots for a SNP GWAS
#' 
#' This function generates the Manhattan plot(s) for a SNP GWAS.
#' @param biallelic A list called 'biallelic' created from the lin_loc function. It is a required input.
#' @param triallelic A list called 'triallelic' created from the lin_loc function. It is a required input.
#' @param config A list called 'config' created from the lin_loc function. It is a required input.
#' @param colourPalette A vector of colours colour the significant principal components identified by the Bayesian Wald test (see testGenomeWidePCs). If this is NULL then colours are chosen from a default colour palette.
#' @keywords SNP
#' @keywords Manhattan-plot
#' @export
#' @examples
#' data <- lin_loc(gen = gen, pheno = pheno, phylo = phylo, 
#'  prefix = prefix, gem.path = gem.path)
#' snpManhattanPlot(biallelic = data$biallelic, triallelic = data$triallelic, 
#'  config = data$config)
snpManhattanPlot = function(biallelic, triallelic, config, colourPalette = NULL){

	sampleCount = length(biallelic$pheno)
	
	if(is.null(colourPalette)){
		colourPalette = getColourPalette(p.pca.bwt = biallelic$p.pca.bwt, signifCutOff = config$signif_cutoff, pc.lim = biallelic$pc_order$pc.lim)	
	}
	
	prefix = config$prefix
	o = biallelic$pc_order$pc_order
	which.pc = biallelic$cor.XX$which.pc
	ipat = biallelic$pattern
	max.cor.pc = biallelic$cor.XX$max.cor.pc
	cutoffCor = config$cutoffCor
	ipat.snps = triallelic$pattern
	which.pc.tritetra = triallelic$cor.tritetra$which.pc
	max.cor.pc.tritetra = triallelic$cor.tritetra$max.cor.pc
	cor.XX = biallelic$cor.XX
	fit.lm = biallelic$logreg
	fit.lm.tritetra = triallelic$logreg
	fit.lmm = biallelic$lmm
	fit.lmm.tritetra = triallelic$lmm
	
	snpColours = .getSNPColours(sampleCount = sampleCount,
                           colourPalette = colourPalette,
                           colouredPCs =  o[1:20],
                           which.pc = which.pc,
                           ipat = ipat,
                           max.cor.pc = max.cor.pc,
                           cutoffCor = cutoffCor,
                           which.pc.tritetra = which.pc.tritetra,
                           ipat.snps = ipat.snps,
                           max.cor.pc.tritetra = max.cor.pc.tritetra,
                           cor.XX = cor.XX)
  
 
  
  bipCount = 0
  ttpCount = 0
  if(!is.null(fit.lmm)){
  	bipCount = nrow(fit.lmm)
  }
  if(!is.null(fit.lmm.tritetra)){
  	ttpCount = nrow(fit.lmm.tritetra)
  }
  snpType = c(rep(1,bipCount),rep(2,ttpCount))
  
  ## If available, plot logistic regression P-values on a Manhattan plot.
  logregPvalues = NULL
  if(!is.null(fit.lm)){
    logregPos = c(fit.lm$ps,fit.lm.tritetra$ps)
    logregPvalues = c(fit.lm$negLog10, fit.lm.tritetra$negLog10)
    .manhattanPlot(prefix = paste0(prefix,"_ManhattanRawPvalues"),
                  snpPos = logregPos,
                  pValues = logregPvalues,
                  snpType = snpType,
                  main = "Logistic Regression SNPs Manhattan Plot",
                  col = c(snpColours$bip, snpColours$ttp))
  }
  
  
  ## Plot LMM P-values on a Manhattan plot.
  lmmPos = c(fit.lmm$ps,fit.lmm.tritetra$ps)
  lmmPvalues = c(-log10(fit.lmm$p_lrt),-log10(as.numeric(fit.lmm.tritetra$pvals)))
  .manhattanPlot(prefix = paste0(prefix,"_ManhattanLMMPvalues"),
                snpPos = lmmPos,
                pValues = lmmPvalues,
                snpType = snpType,
                main = "LMM SNPs Manhattan Plot",
                col = c(snpColours$bip, snpColours$ttp))
}

## Get the colouring of the SNPs for manhattanp plot
.getSNPColours = function(sampleCount = NULL,
                         colourPalette = NULL,#colourPalette
                         colouredPCs = NULL, # o[1:20]
                         which.pc = NULL,
                         ipat = NULL,
                         max.cor.pc = NULL,
                         cutoffCor = NULL,
                         which.pc.tritetra = NULL,
                         ipat.snps = NULL,
                         max.cor.pc.tritetra = NULL,
                         cor.XX = NULL){
 
  COL = rep("grey50", sampleCount)
  
  COL[colouredPCs] = colourPalette
  COL = COL[which.pc][ipat]
  COL[max.cor.pc[ipat] < cutoffCor] = "grey50"
  COL2 = rep("grey50", sampleCount)
  COL2[colouredPCs] = colourPalette
  COL2 = COL2[which.pc.tritetra][ipat.snps]
  COL2[max.cor.pc.tritetra[ipat.snps] < cutoffCor] = "grey50"
  
  return(list(bip = COL ,ttp = COL2))
  
}