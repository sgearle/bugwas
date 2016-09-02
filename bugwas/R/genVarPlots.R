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

#' Plots for GWAS on general variants
#' 
#' This function generates the various Manhattan plots for general variants.
#' @param biallelic A list called 'biallelic' created from the lin_loc function. It is a required input.
#' @param config A list called 'config' created from the lin_loc function. It is a required input.
#' @param genVars A list called 'genVars' created from the lin_loc function. It is a required input.
#' @param colourPalette A vector of colours colour the significant principal components identified by the Bayesian Wald test (see testGenomeWidePCs). If this is NULL then colours are chosen from a default colour palette.
#' @keywords Manhattan-plot
#' @export
#' @examples
#' data <- lin_loc(gen = gen, pheno = pheno, phylo = phylo, 
#'  prefix = prefix, gem.path = gem.path)
#' genVarPlots(genVars = data$genVars, biallelic = data$biallelic, config = data$config)
genVarPlots = function(genVars, biallelic, config, colourPalette = NULL){
	
	o = biallelic$pc_order$pc_order
	p.pca.bwt = biallelic$p.pca.bwt
	pc.lim = biallelic$pc_order$pc.lim
	
	if(is.null(colourPalette)){
		colourPalette = getColourPalette(p.pca.bwt = p.pca.bwt, 
									signifCutOff = config$signif_cutoff, 
									pc.lim = pc.lim)	
	}
	
	sampleCount = length(biallelic$pheno)
	
	.genVarPlots(genVars = genVars, 
				o = o, 
				p.pca.bwt = p.pca.bwt, 
				pc.lim = pc.lim, 
                colourPalette = colourPalette, 
                prefix = config$prefix, 
                npcs = biallelic$npcs, 
                sampleCount = sampleCount,
                cutoffCor = config$cutoffCor)
	
	
}


.genVarPlots = function(genVars = NULL, o = NULL, p.pca.bwt =NULL, pc.lim = NULL, 
                       npcs = NULL, colourPalette = NULL, prefix = NULL,
                       sampleCount = NULL, cutoffCor = NULL){
  
  variantCount = length(genVars)
  variantNames = names(genVars)
  
  
  if(is.null(variantNames)){
    variantNames = paste0("genVar",c(1:variantCount))
  }
  
  
  
  
  
  for(iGenVar in 1:variantCount){
    count = length(genVars[[iGenVar]]$lmm$negLog10)
    
    
    genVarCol = .getGenVarColours(sampleCount = sampleCount,
                     colourPalette = colourPalette,
                     which.pc =  genVars[[iGenVar]]$cor.var$which.pc,
                     varpat = genVars[[iGenVar]]$varpat,
                     max.cor.pc = genVars[[iGenVar]]$cor.var$max.cor.pc,
                     cutoffCor = cutoffCor,
                     o = o)
    
    if(!is.null(genVars[iGenVar[iGenVar]]$logreg)){
      .manhattanPlot(prefix = paste0(prefix,"_", variantNames[iGenVar], "_ManhattanRawPvalues"),
                    snpPos = c(1:count),
                    pValues = genVars[[iGenVar]]$logreg$negLog10,
                    snpType = rep(1,count),
                    main = "Logistic Regression SNPs Manhattan Plot",
                    col = genVarCol[genVars[[iGenVar]]$pattern],
                    xlab = "Variant Index")
    }
    
    
    .manhattanPlot(prefix = paste0(prefix,"_", variantNames[iGenVar], "_ManhattanLMMPvalues"),
                  snpPos = c(1:length(genVars[[iGenVar]]$lmm$negLog10)),
                  pValues = genVars[[iGenVar]]$lmm$negLog10,
                  snpType = rep(1,count),
                  main = "LMM SNPs Manhattan Plot",
                  col = genVarCol[genVars[[iGenVar]]$pattern],
                  xlab = "Variant Index")
    
    .plot_pc_manhattan(o = o, 
                      which.pc = genVars[[iGenVar]]$cor.var$which.pc, 
                      pattern = genVars[[iGenVar]]$pattern, 
                      p.pca.bwt = p.pca.bwt, 
                      pc.lim = pc.lim, 
                      negLog10 = genVars[[iGenVar]]$lmm$negLog10, 
                      pat.weight = genVars[[iGenVar]]$varpat, 
                      prefix=paste0(prefix,"_", variantNames[iGenVar]),
                      colourPalette = colourPalette,
                      npcs = npcs)
    
    
    
  
  }
  
}


.getGenVarColours = function(sampleCount = NULL,
                         colourPalette = NULL,#colourPalette
                         which.pc = NULL,
                         varpat = NULL,
                         max.cor.pc = NULL,
                         cutoffCor = NULL,
                         o = NULL){
                         	  
  COL = rep("grey50", sampleCount)
  
  COL[o[1:20]] = colourPalette
  COL = COL[which.pc]
  COL[max.cor.pc < cutoffCor] = "grey50"  
  
  return(COL)
  
}