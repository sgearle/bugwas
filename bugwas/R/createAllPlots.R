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

.createAllPlots = function(prefix = NULL, 
genVars = NULL, cutoffCor = NULL,  npcs = NULL, phenotype = NULL,
                          pca = NULL, fit.lm = NULL, fit.lmm = NULL, fit.lm.tritetra = NULL,
                          fit.lmm.tritetra = NULL, gemma = NULL, gemma.tritetra = NULL, 
                          ipat = NULL, ipat.snps = NULL, pred2 = NULL, cor.XX = NULL, cor.tritetra = NULL,
                          which.pc = NULL, max.cor.pc = NULL, 
                          which.pc.tritetra = NULL, max.cor.pc.tritetra = NULL,
                          which.mtp.pc = NULL, max.mtp.cor.pc = NULL,
                          XX.comid = NULL, pc_order = NULL, o = NULL, pc.lim = NULL, p.pca.bwt = NULL,
                          signifCutOff = NULL, bippat = NULL, snppat = NULL,
                          tree = NULL, treepat = NULL){
              
  
  # Get a random sample of colours to use for all plots, equal to the number of significant PCs
  colourPalette = getColourPalette(p.pca.bwt = p.pca.bwt, signifCutOff = signifCutOff, pc.lim = pc.lim)
  
  sampleCount = length(phenotype)
  m = match(o[pc.lim], which.mtp.pc)
  
  #Bayesian Wald test for genome-wide PCs
  p.genomewidepc = .testGenomeWidePCs(prefix = prefix,
                                    pc.lim = pc.lim,
                                     pca = pca,
                                     bippat = bippat,
                                     ipat = ipat,
                                     o = o)
  message("Bayesian Wald test for genome-wide PCs has been completed successfully.")
  
  #The barplot for the Bayesian wald test for genome-wide PCs
  .BayesianWaldTestPCsBarplot(prefix = prefix,
                              p.pca.bwt = p.pca.bwt,
                              colourPalette = colourPalette,
                              o = o,
                              m = m,
                              p.genomewidepc = p.genomewidepc,
                              pc.lim = pc.lim)
   message("The barplot for the Bayesian wald test for genome-wide PCs has been completed successfully.")
  
 
  
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
  
  #The Manhattan plot for SNP GWAS using logistic regerssion
  logregPvalues = NULL
  logregPos = NULL
  if(!is.null(fit.lm)){
  	
  	logregPos = c(fit.lm$ps,fit.lm.tritetra$ps)
    logregPvalues = c(fit.lm$negLog10,fit.lm.tritetra$negLog10)
    .manhattanPlot(prefix = paste0(prefix,"_ManhattanRawPvalues"),
                  snpPos = logregPos,
                  pValues = logregPvalues,
                  snpType = snpType,
                  main = "Logistic Regression SNPs Manhattan Plot",
                  col = c(snpColours$bip, snpColours$ttp))
                  
    message("The Manhattan plot for SNP GWAS using logistic regerssion has been completed successfully.")
                  
  }
  
  
  #The Manhattan plot for SNP GWAS using LMM
  lmmPos = c(fit.lmm$ps,fit.lmm.tritetra$ps)
  lmmPvalues = c(-log10(fit.lmm$p_lrt),-log10(as.numeric(fit.lmm.tritetra$pvals)))
  .manhattanPlot(prefix = paste0(prefix,"_ManhattanLMMPvalues"),
                snpPos = lmmPos,
                pValues = lmmPvalues,
                snpType = snpType,
                main = "LMM SNPs Manhattan Plot",
                col = c(snpColours$bip, snpColours$ttp))
                
  message("The Manhattan plot for SNP GWAS using LMM has been completed successfully.")
  
  
  #The plot of logistic regression P-values vs. LMM P-values for SNP GWAS
  if(!is.null(logregPvalues)){
    .logregVsLMM(prefix = prefix,
                         col = c(snpColours$bip, snpColours$ttp),
                         logregPvalues = logregPvalues,
                         lmmPvalues =  lmmPvalues,
                         snpType = snpType,
                         pcOrder = o,
                         pc.lim = pc.lim,
                         colourPalette = colourPalette)
                         
     message("The plot of logistic regression P-values vs. LMM P-values for SNP GWAS has been completed successfully.")                    
  }
  
  
  
  
  # Plot the individuals by their top two significant additive PCs
  .plotIndividualBy2PCs(pc1 = o[1], pc2 = o[2],
                       pc1.scores = pca$x[,o[1]], pc2.scores = pca$x[,o[2]],
                       prefix= prefix, phenotype)
                       
   message("The reduced space plot of the sample on the top two significant additive PCs been completed successfully.")                    
                       
  
  #The plot with true and predicted phenotype on the tree
  .trueAndPredPhenoOnTreePlot(prefix = prefix, tree = tree, which.mtp.pc = unlist(which.mtp.pc), #Check with SGE
                             max.mtp.cor.pc = max.mtp.cor.pc, cutoffCor = cutoffCor, treepat = treepat,
                             pcOrder = o, p.genomewidepc = p.genomewidepc, phenotype = phenotype, 
                             XX.comid = XX.comid, colourPalette = colourPalette,
                             pc.lim = pc.lim, pred2 = pred2)
                             
   message("The plot with true and predicted phenotype on the tree has been completed successfully.")                         
  
  
  #The plots of PCs loadings
  .pcLoadingsPlot(prefix = prefix, pca = pca, pc.lim = pc.lim, 
                 pcOrder = o, ipat = ipat, bippat = bippat, 
                 bipPos = fit.lmm$ps)
  
  message("The plots of PCs loadings have been completed successfully.")
  

  ## To run for all SNPs (biallelic and tri and tetra allelic)
  new.pat <- ipat.snps
  new.pat <- sapply(new.pat, function(x, pat){ x + length(pat)}, pat = bippat)
  new.pat <- c(ipat, unlist(new.pat))
  
  #The Manhattan plot organised by PCs for SNP GWAS
  .plot_pc_manhattan(o = o, 
                   which.pc = c(cor.XX$which.pc, cor.tritetra$which.pc), 
                   pattern = new.pat, 
                   p.pca.bwt = p.pca.bwt, 
                   pc.lim = pc_order$pc.lim, 
                   negLog10 = c(fit.lmm$negLog10, fit.lmm.tritetra$negLog10), 
                   pat.weight = c(bippat, snppat), 
                   prefix=paste0(prefix,"_SNPs"),
                   colourPalette = colourPalette,
                   npcs = npcs)
   message("The Manhattan plot organised by PCs for SNP GWAS has been completed successfully.")                
  
  
  #The plots for general genetic variants
  if(!is.null(genVars)){
    .genVarPlots(genVars = genVars, o = o, p.pca.bwt = p.pca.bwt, pc.lim = pc.lim, 
                colourPalette = colourPalette, prefix = prefix, npcs = npcs, 
                sampleCount = sampleCount,cutoffCor = cutoffCor)
     message("The plots for general genetic variants have been completed successfully.")            
    
   
  }
  
}

