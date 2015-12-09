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

#' Generates all plots.
#' 
#' This function that plots the true and predicted phenotypes on a tree.
#' @param biallelic A list called 'biallelic' created from the lin_loc function. It is a required input.
#' @param config A list called 'config' created from the lin_loc function. It is a required input.
#' @param treeInfo A list called 'treeInfo' created from the lin_loc function. It is a required input.
#' @param p.genomewidepc A matrix of the significant principal component and their correlation with lineages. This is Bayesian Wald test results produced by the function testGenomeWidePCs. If this is NULL then testGenomeWidePCs is called to generate the required test results. 
#' @param colourPalette A vector of colours colour the significant principal components identified by the Bayesian Wald test (see testGenomeWidePCs). If this is NULL then colours are chosen from a default colour palette.
#' @keywords Phenotype
#' @keywords Phylogram
#' @export
#' @examples
#' data <- lin_loc(gen = gen, pheno = pheno, phylo = phylo, 
#'  prefix = prefix, gem.path = gem.path)
#' trueAndPredPhenoOnTreePlot(biallelic = data$biallelic, treeInfo = data$treeInfo,
#'  config = data$config)
trueAndPredPhenoOnTreePlot = function(config, biallelic, treeInfo, 
										p.genomewidepc = NULL, colourPalette = NULL){
	if(is.null(p.genomewidepc)){
		p.genomewidepc = testGenomeWidePCs(config = config, biallelic = biallelic)
	}
	
	if(is.null(colourPalette)){
		colourPalette = getColourPalette(p.pca.bwt = biallelic$p.pca.bwt, signifCutOff = config$signif_cutoff, pc.lim = biallelic$pc_order$pc.lim)	
	}
	
	.trueAndPredPhenoOnTreePlot(prefix = config$prefix, 
								tree = treeInfo$tree, 
								which.mtp.pc = unlist(treeInfo$cor.tree$which.pc), #Check with SGE
								max.mtp.cor.pc = treeInfo$cor.tree$max.cor.pc, 
								cutoffCor = config$cutoffCor, 
								treepat = treeInfo$pattern,
								pcOrder = biallelic$pc_order$pc_order, 
								p.genomewidepc = p.genomewidepc, 
								phenotype =  biallelic$pheno, 
								XX.comid = biallelic$id, 
								colourPalette = colourPalette,
								pc.lim = biallelic$pc_order$pc.lim, 
								pred2 = biallelic$pred)
}

.trueAndPredPhenoOnTreePlot = function(prefix = NULL, #prefix,"_",name
         tree = NULL,
         which.mtp.pc = NULL,
         max.mtp.cor.pc = NULL,
         cutoffCor = NULL,
         treepat = NULL,
         pcOrder = NULL, # o
         p.genomewidepc = NULL,
         phenotype = NULL, #y
         XX.comid = NULL,
         colourPalette = NULL, #colourPalette
         pc.lim = NULL,
         pred2 = NULL){ #m
         	
    
	m = match(pcOrder[pc.lim], which.mtp.pc)
  
  
  
  branch.col.pc = .getBranchColors(treepat = treepat, tree = tree, phenotype = phenotype, 
                          pcOrder = pcOrder, colourPalette = colourPalette, which.mtp.pc = which.mtp.pc, 
                          max.mtp.cor.pc = max.mtp.cor.pc, cutoffCor = cutoffCor)
  
  tree.eq = .get.tree.eq(tree)

  pred3 = (pred2-min(phenotype-mean(phenotype)))/diff(range(phenotype-mean(phenotype)))
  
  edgeLabels = rep("",length(tree$edge))
  edgeLabelsCOL = rep("white",length(tree$edge))
  add.signif = NULL
  
  
  if(length(pc.lim)>0){
    sigInd = .getSigIndicators(pc.lim = pc.lim, pcOrder = pcOrder, p.genomewidepc = p.genomewidepc, 
                     tree.eq = tree.eq, branch.col.pc = branch.col.pc, colourPalette = colourPalette)
    edgeLabels = sigInd$edgeLabels
    edgeLabelsCOL = sigInd$edgeLabelsCOL
    add.signif = sigInd$add.signif
    
  }
  
  
  # Added st() around pred3 to get between 0 and 1
  predCOLS=colorRamp(c("grey","black"))
  predCOLS=predCOLS(.st(pred3))
  predCOLS =rgb(predCOLS, maxColorValue=256)
  
  
  .trueAndPredPhenoOnTree(prefix = prefix, phenotype = phenotype, tree.eq = tree.eq,
                         branch.col.pc = branch.col.pc, tree = tree, XX.comid = XX.comid,
                         predCOLS = predCOLS, 
                         edgeLabels = edgeLabels, edgeLabelsCOL = edgeLabelsCOL,
                         add.signif = add.signif,
                         pc.lim = pc.lim, colourPalette = colourPalette, matchedPCs = m)
  
}

.getBranchColors = function(tree = NULL,
                           treepat = NULL,
                           phenotype = NULL,
                           pcOrder = NULL,
                           colourPalette = NULL,
                           which.mtp.pc = NULL,
                           max.mtp.cor.pc = NULL,
                           cutoffCor = NULL){
  
  branch.col.pc = rep("grey50",nrow(tree$edge))
  bcols = rep("grey50",length(phenotype))
  bcols[pcOrder[1:20]] = colourPalette
  bcols = bcols[which.mtp.pc]
  bcols[max.mtp.cor.pc<cutoffCor] = "grey50"
  branch.col.pc = bcols[match(1:length(branch.col.pc), treepat$ancestral_edge)]
  
  return(branch.col.pc)
}

.get.tree.eq = function(tree = NULL){
  tree.eq = tree
  tree.eq$edge.length = sqrt(tree$edge.length)
  tree.eq$tip.label = rep("__———__",length(tree$tip.label))
  return(tree.eq)
}


.getSigIndicators = function(pc.lim = NULL,
                            pcOrder = NULL,
                            tree.eq = NULL,
                            branch.col.pc = NULL,
                            colourPalette = NULL,
                            p.genomewidepc = NULL){
  
  edgeCount = length(tree.eq$edge)
  edgeLabels = rep("", edgeCount)
  edgeLabelsCOL = rep("white", edgeCount)
  
  
  
  for(i in 1:length(pc.lim)){
    
    a = which.max(tree.eq$edge.length[which(branch.col.pc==colourPalette[i])])
    if(colourPalette[i]!="grey50"){
      edgeLabels[which(branch.col.pc==colourPalette[i])[a]]=paste0("PC",pcOrder[pc.lim][i])
      edgeLabelsCOL[which(branch.col.pc==colourPalette[i])[a]]=colourPalette[i]
    }
  }
  
  
  add.signif=paste("PC",pcOrder[pc.lim])
  for(i in 1:length(add.signif)){
    if(as.numeric(p.genomewidepc[i,2])>=2 & as.numeric(p.genomewidepc[i,2])<5){
      add.signif[i] = paste("PC", pcOrder[pc.lim[i]],"*")
    } else if(as.numeric(p.genomewidepc[i,2])>=5){
      add.signif[i] = paste("PC", pcOrder[pc.lim[i]],"**")
    }
  }
  
  return(list(edgeLabels = edgeLabels, edgeLabelsCOL = edgeLabelsCOL, add.signif = add.signif))
}

.trueAndPredPhenoOnTree = function(prefix = NULL,
         phenotype = NULL,
         tree.eq = NULL,
         branch.col.pc = NULL,
         tree = NULL,
         XX.comid = NULL,
         predCOLS = NULL,
         edgeLabels = NULL,
         edgeLabelsCOL = NULL,
         add.signif = NULL,
         pc.lim = NULL,
         colourPalette = NULL,
         matchedPCs = NULL){
  # Plot tree with true and predicted phenotypes on tips, with internal branches coloured by PCs
  .pl(paste0(prefix, "_tree_branchescolouredbyPC"),{
    n = length(phenotype)
    par(oma = c(1, 1, 1, 5))
    ape::plot.phylo(tree.eq, 
         type = "fan", show.tip.label = TRUE, edge.col = branch.col.pc, 
         tip.col=c("gray","black")[1+phenotype-min(phenotype)][match(tree$tip.label, XX.comid)],
         lwd=3, adj=0.5, lab4ut = "axial")
    xx = get("last_plot.phylo", envir = ape::.PlotPhyloEnv)$xx[1:n][match(XX.comid,tree$tip.label)]
    yy = get("last_plot.phylo", envir = ape::.PlotPhyloEnv)$yy[1:n][match(XX.comid,tree$tip.label)]
    
    points(1.08*xx, 1.08*yy, col=predCOLS, xpd=TRUE, lwd=.5, cex=1.5)
    ape::edgelabels(edgeLabels,frame="none",bg="none",col=edgeLabelsCOL)
    par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
    plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
    if(length(pc.lim)>0){
      # Turns the colour of the PCs which don't correlate with the tree to grey
      colourPalette[which(is.na(matchedPCs)==TRUE)]="grey50"
      legend("right",
             c(add.signif,"Other"),
             fill = c(colourPalette[pc.lim], "grey50"),
             bty="n", cex=1.8, xpd = TRUE, inset=c(0,0))
    }
  })
}

.st = function(x){ 
  (x-min(x,na.rm=TRUE))/diff(range(x,na.rm=TRUE))
}
