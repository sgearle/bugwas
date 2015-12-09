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

#' Plot of the sample on the first two principal components
#' 
#' This function generates a plot of the sample on the first two principal components.
#' @param biallelic A list called 'biallelic' created from the lin_loc function. It is a required input.
#' @param config A list called 'config' created from the lin_loc function. It is a required input.
#' @keywords Reduced-space-plot
#' @export
#' @examples
#' data <- lin_loc(gen = gen, pheno = pheno, phylo = phylo, 
#'  prefix = prefix, gem.path = gem.path)
#' plotIndividualBy2PCs(biallelic = data$biallelic, config = data$config)
plotIndividualBy2PCs = function(biallelic = NULL, config = NULL){
	o = biallelic$pc_order$pc_order
	pca = biallelic$pca
	# Plot the individuals by their top two significant additive PCs
  .plotIndividualBy2PCs(pc1 = o[1], pc2 = o[2],
                       pc1.scores = pca$x[,o[1]], pc2.scores = pca$x[,o[2]],
                       prefix= config$prefix, phenotype = biallelic$pheno)
	
}

# Plot the individuals by their top two significant additive PCs
.plotIndividualBy2PCs = function(pc1 = NULL, #o[1]
                                pc2 = NULL, #o[2]
                                pc1.scores = NULL,#pca$x[,o[1]]
                                pc2.scores = NULL, #pca$x[,o[2]]
                                prefix= NULL,
                                phenotype = NULL){ # prefix,"_",name,"
  
  
  .pl(paste0(prefix,"_indiv_first2signifPCs"),{
    # Plot the individuals on the first two additive value PCs
    plot(pc1.scores+rnorm(length(pc1.scores), 0, sd(pc1.scores)/10),
         pc2.scores+rnorm(length(pc2.scores), 0, sd(pc2.scores)/10),
         col=c("blue","red")[1+phenotype], lwd=.5, cex=1.5 ,
         xlab = paste0("PC ",pc1), ylab=paste0("PC ",pc2))
  })
  
}