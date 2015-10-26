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
#' This function generates all the plots
#' @param biallelic A list called 'biallelic' created from the lin_loc function
#' @param triallelic A list called 'triallelic' created from the lin_loc function
#' @param genVars A list called 'genVars' created from the lin_loc function
#' @param treeInfo A list called 'treeInfo' created from the lin_loc function
#' @param config A list called 'config' created from the lin_loc function
#' @keywords plot
#' @export
#' @examples
#' data <- lin_loc(gen = gen, pheno = pheno, phylo = phylo, 
#'  prefix = prefix, gem.path = gem.path)
#' all_plots(biallelic = data$biallelic, triallelic = data$triallelic, 
#' 	genVars = data$genVars, treeInfo = data$treeInfo, config = data$config)
all_plots = function(biallelic = NULL,
                     triallelic = NULL,
                     genVars = NULL,
                     treeInfo = NULL,
                     config = NULL){
  
  
  .createAllPlots(prefix = config$prefix,
                 genVars = genVars,
                 cutoffCor = config$cutoffCor,
                 npcs = biallelic$npcs,
                 phenotype = biallelic$pheno,
                 pca = biallelic$pca,
                 fit.lm = biallelic$logreg,
                 fit.lmm = biallelic$lmm,
                 fit.lm.tritetra = triallelic$logreg,
                 fit.lmm.tritetra = triallelic$lmm,
                 gemma = biallelic$lmm,
                 gemma.tritetra = triallelic$lmm,
                 ipat = biallelic$pattern,
                 ipat.snps = triallelic$pattern,
                 pred2 = biallelic$pred,
                 cor.XX = biallelic$cor.XX,
                 cor.tritetra = triallelic$cor.tritetra,
                 which.pc =  biallelic$cor.XX$which.pc,
                 max.cor.pc = biallelic$cor.XX$max.cor.pc,
                 which.pc.tritetra = triallelic$cor.tritetra$which.pc,
                 max.cor.pc.tritetra = triallelic$cor.tritetra$max.cor.pc,
                 which.mtp.pc = treeInfo$cor.tree$which.pc,
                 max.mtp.cor.pc = treeInfo$cor.tree$max.cor.pc,
                 XX.comid = biallelic$id,
                 pc_order = biallelic$pc_order,
                 o = biallelic$pc_order$pc_order,
                 pc.lim = biallelic$pc_order$pc.lim,
                 p.pca.bwt = biallelic$p.pca.bwt,
                 signifCutOff = config$signif_cutoff,
                 bippat = biallelic$bippat,
                 snppat = triallelic$snppat,
                 tree = treeInfo$tree,
                 treepat = treeInfo$pattern);
  
}





