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

#' Barplot of Bayesian Wald Test on principal components
#' 
#' This function generates the barplot of Bayesian Wald Test on principal components.
#' @param biallelic A list called 'biallelic' created from the lin_loc function. It is a required input.
#' @param config A list called 'config' created from the lin_loc function. It is a required input.
#' @param treeInfo A list called 'treeInfo' created from the lin_loc function. It is a required input.
#' @param colourPalette A vector of colours colour the significant principal components identified by the Bayesian Wald test (see testGenomeWidePCs). If this is NULL then colours are chosen from a default colour palette.
#' @param p.genomewidepc A matrix of the significant principal component and their correlation with lineages. This is Bayesian Wald test results produced by the function testGenomeWidePCs. If this is NULL then testGenomeWidePCs is called to generate the required test results. 
#' @keywords Bayesian-Wald-test
#' @keywords PCA
#' @keywords Barplot
#' @export
#' @examples
#' data <- lin_loc(gen = gen, pheno = pheno, phylo = phylo, 
#'  prefix = prefix, gem.path = gem.path)
#' testGenomeWidePCs(config = data$config, biallelic = data$biallelic)
bayesianWaldTestPCsBarplot = function(config, biallelic, treeInfo, colourPalette = NULL, p.genomewidepc = NULL){
	
	if(is.null(colourPalette)){
		colourPalette = getColourPalette(p.pca.bwt = biallelic$p.pca.bwt, signifCutOff = config$signif_cutoff, pc.lim = biallelic$pc_order$pc.lim)	
	}
	
	if(is.null(p.genomewidepc)){
		p.genomewidepc = testGenomeWidePCs(config = config, biallelic = biallelic)
	}
	
	o = biallelic$pc_order$pc_order
	which.mtp.pc = treeInfo$cor.tree$which.pc
	pc.lim = biallelic$pc_order$pc.lim
	m = match(o[pc.lim], which.mtp.pc)
	
	.BayesianWaldTestPCsBarplot(prefix = config$prefix,
                             p.pca.bwt = biallelic$p.pca.bwt,
                             colourPalette = colourPalette,
                             o = o,
                             m = m,
                             p.genomewidepc = p.genomewidepc,
                             pc.lim = biallelic$pc_order$pc.lim)

	
}

.BayesianWaldTestPCsBarplot = function(prefix = NULL,
                                      p.pca.bwt = NULL,
                                      colourPalette = NULL,
                                      o = NULL,
                                      m = NULL,
                                      p.genomewidepc = NULL,
                                      pc.lim = NULL){
  
  barDensity = rep(NA, length(colourPalette))
  barDensity[which(is.na(m) == TRUE)] = 30
  
  # Barplot of Bayesian Wald Test for PCs
  .pl(paste0(prefix,"_barplot_BayesianWald_PCs"),{
    b=barplot(as.vector(p.pca.bwt[o][1:20]),
              col = colourPalette,
              xlab = "Principal Component", ylab = "-log10(p-value)",
              main="Bayesian Wald Test", names=o[1:20],
              cex.lab=1.5, cex.axis=1.5,cex.main=2)
  
    colourPalette2 = colourPalette
    colourPalette2[which(is.na(m)==TRUE)]="black"
    b = barplot(as.vector(p.pca.bwt[o][1:20]),
              col=colourPalette2, xlab="",
              main="Bayesian Wald Test",
              cex.lab=1.5,cex.axis=1.5, 
              names=o[1:20],cex.main=2,
              density=barDensity, add=TRUE, axes=FALSE)
  
    if(length(pc.lim)>0){
      text(b,as.vector(p.pca.bwt[o][1:20]), 
           ifelse(as.numeric(p.genomewidepc[,2])>=2 & as.numeric(p.genomewidepc[,2])<5,"*",""),
           pos=3, cex=2, xpd=NA)
      text(b,as.vector(p.pca.bwt[o][1:20]),
           ifelse(as.numeric(p.genomewidepc[,2])>=5,"**",""),
           pos=3, cex=2, xpd=NA)
      legend("topright",
             c("Non genome-wide PCs","",
               expression(paste("*    ",italic('p')," ≤ 0.01"),collapse=""),
               expression(paste("**  ",italic('p')," ≤ 1e-5",collapse=""))), 
             bty="n",cex=2)
    }
    
  })
}