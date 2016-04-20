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

#' Plot PC loadings of all SNPs.
#' 
#' This function plots the loadings of all SNPs for each significant principle component identified by the Bayesian Wald test.
#' @param biallelic A list called 'biallelic' created from the lin_loc function. It is a required input.
#' @param config A list called 'config' created from the lin_loc function. It is a required input.
#' @keywords Scatter-plot
#' @export
#' @examples
#' data <- lin_loc(gen = gen, pheno = pheno, phylo = phylo, 
#'  prefix = prefix, gem.path = gem.path)
#' pcLoadingPlot(config = data$config, biallelic = data$biallelic)
pcLoadingPlot = function(config, biallelic){
	.pcLoadingsPlot(prefix = config$prefix, 
					pca = biallelic$pca, 
					pc.lim = biallelic$pc_order$pc.lim, 
					pcOrder = biallelic$pc_order$pc_order,
					ipat = biallelic$pattern,
					bippat = biallelic$bippat, 
					bipPos = biallelic$ps)

}

.pcLoadingsPlot = function(prefix = NULL, #prefix,"_",name,"
         pca = NULL,
         pc.lim = NULL,
         pcOrder = NULL,
         ipat = NULL,
         bippat = NULL,
         bipPos = NULL){
  # Output plots of the PC loadings
  dirName = paste0(prefix, "_PCloadings")
  dir.create(dirName)

  if(length(pc.lim)>0){
    for(i in 1:length(pc.lim)){
      a=pca$rotation[,pcOrder[pc.lim[i]]]
      a=a/(bippat)
      
      # Keep nmut if want to colour the loadings by homoplasy
      # If colouring my homoplasy, need to put back colours into plot()
      
    
      png(paste0(dirName,"/",prefix,"_PC",pcOrder[pc.lim[i]],"_weighted_loadings.png"),
          width=1000,height=700)
      plot(x = as.numeric(bipPos), y = as.numeric(a[ipat]), 
           xlab = "Position in genome", 
           ylab = paste0("PC",pcOrder[pc.lim[i]]," weighted loadings"),
           cex=0.6)
      dev.off()	
    }
  
  }
}