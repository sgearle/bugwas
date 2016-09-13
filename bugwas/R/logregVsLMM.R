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

#' Plot of the P-values of logistic regression versus those of LMM
#' 
#' This function generates the plot of the P-values of logistic regression versus those of SNP GWAS.
#' @param biallelic A list called 'biallelic' created from the lin_loc function. It is a required input.
#' @param triallelic A list called 'triallelic' created from the lin_loc function. It is a required input.
#' @param config A list called 'config' created from the lin_loc function. It is a required input.
#' @param colourPalette A vector of colours colour the significant principal components identified by the Bayesian Wald test (see testGenomeWidePCs). If this is NULL then colours are chosen from a default colour palette.
#' @keywords SNP
#' @keywords Scatter-plot
#' @keywords P-values
#' @export
#' @examples
#' data <- lin_loc(gen = gen, pheno = pheno, phylo = phylo, 
#'  prefix = prefix, gem.path = gem.path)
#' logregVsLMM(biallelic = data$biallelic, triallelic = data$triallelic, config = data$config)
logregVsLMM = function(config, biallelic, triallelic, colourPalette = NULL){
	sampleCount = length(biallelic$pheno)
	
	if(is.null(colourPalette)){
		colourPalette = getColourPalette(p.pca.bwt = biallelic$p.pca.bwt, 
									signifCutOff = config$signif_cutoff, 
									pc.lim = biallelic$pc_order$pc.lim)	
	}
	
	o = biallelic$pc_order$pc_order
	which.pc = biallelic$cor.XX$which.pc
	ipat = biallelic$pattern
	max.cor.pc = biallelic$cor.XX$max.cor.pc
	cutoffCor = config$cutoffCor
	which.pc.tritetra = triallelic$cor.tritetra$which.pc
	ipat.snps = triallelic$pattern
	max.cor.pc.tritetra = triallelic$cor.tritetra$max.cor.pc
	cor.XX = biallelic$cor.XX
	
	fit.lm = biallelic$logreg
	fit.lm.tritetra = triallelic$logreg
	fit.lmm = biallelic$lmm
	fit.lmm.tritetra = triallelic$lmm
	
	logregPvalues = c(fit.lm$negLog10, fit.lm.tritetra$negLog10)
	lmmPvalues = c(-log10(fit.lmm$p_lrt),-log10(as.numeric(fit.lmm.tritetra$pvals)))
	snpType = c(rep(1,nrow(fit.lmm)),rep(2,nrow(fit.lmm.tritetra)))
	logregPos = c(fit.lm$ps, fit.lm.tritetra$ps)
	lmmPos = c(fit.lmm$ps, fit.lmm.tritetra$ps)
	
	if(!identical(sort(logregPos),sort(lmmPos))){
		stop("The SNPs from the logistic regression do not match those from LMM.")
	}
	
	logregPvalues = logregPvalues[order(logregPos)]
	lmmPvalues = lmmPvalues[order(lmmPos)]
	snpType = snpType[order(lmmPos)]

	
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
                           
	.logregVsLMM(prefix = config$prefix,
				col = c(snpColours$bip, snpColours$ttp)[order(lmmPos)],
                         logregPvalues = logregPvalues,
                         lmmPvalues =  lmmPvalues,
                         snpType = snpType,
                         pcOrder = o,
                         pc.lim = biallelic$pc_order$pc.lim,
                         colourPalette = colourPalette)

}

## Plots the logreg results against LMM results
## Coloured by PC with which it has the highest correlation
## If correlation below a cutoff (e.g. 0.5) then coloured grey
## which.pc is in the order of patterns, so need to change it 
## so that which.pc is is order of SNPs
.logregVsLMM = function(prefix = NULL, #prefix,"_",name
         col = NULL, #=c(COL,COL2)
         logregPvalues = NULL, # c(-log10(fit.lm),fit.lm.tritetra[,ncol(fit.lm.tritetra)])
         lmmPvalues = NULL, # c(-log10(gemma$p_lrt),-log10(as.numeric(gemma.tritetra$pvals)))
         snpType = NULL, # c(rep(1,length(fit.lm)),rep(2,nrow(fit.lm.tritetra)))
         pcOrder = NULL, # o
         pc.lim = NULL,
         colourPalette = NULL){ #colourPalette
  .pl2(paste0(prefix,"_SNPs_logreg_v_LMM_PCs_greyLowCor"),{
    
    par(oma = c(1, 1, 1, 10))
    plot(logregPvalues, lmmPvalues,
        col = col, pch= snpType, ,cex.lab=2,
        xlab="Raw -log10(p-values)", 
        ylab="LMM -log10(p-values)",
        main="All SNPs logistic regression v LMM")
  
    abline(0,1,lty=2)
    par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
    plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
    if(length(pc.lim)>0){
      legend("right",
           c(paste("PC", pcOrder[pc.lim]), "Other", "Biallelic", "Multiallelic"),
           pch = c(rep(22,length(pc.lim)+1), 21, 24), 
           pt.bg = c(colourPalette[pc.lim],"grey50","white","white"),
           bty="n", cex=1.8, pt.cex=3, xpd = TRUE, inset = c(0,0))
    } else {
      legend("right",c("Biallelic","Multiallelic"),
           pch=c(21,24),pt.bg=c("white","white"),
           bty="n", cex=1.8, pt.cex=3, xpd = TRUE, inset = c(0,0))
    }	
  })
}