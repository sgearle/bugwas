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

#' Testing the genome-wide principal components
#' 
#' This function performs the genome-wide principal components.
#' @param biallelic A list called 'biallelic' created from the bugwas function.
#' @param config A list called 'config' created from the bugwas function.
#' @keywords Bayesian-Wald-test
#' @keywords PCA
#' @return The p-value of the Bayesian Wald test for the genome-wide effect of principal components.
#' @export
#' @examples
#' data <- lin_loc(gen = gen, pheno = pheno, phylo = phylo, 
#'  prefix = prefix, gem.path = gem.path)
#' testGenomeWidePCs(config = data$config, biallelic = data$biallelic)
testGenomeWidePCs = function(config = NULL, biallelic = NULL){
	p.genomewidepc = .testGenomeWidePCs(prefix = config$prefix,
						pc.lim = biallelic$pc_order$pc.lim,
						pca = biallelic$pca,
						bippat = biallelic$bippat,
						ipat = biallelic$pattern,
						o = biallelic$pc_order$pc_order)
	return(p.genomewidepc)
}

.testGenomeWidePCs = function(prefix = NULL,
                             pc.lim = NULL,
                             pca = NULL,
                             bippat = NULL,
                             ipat = NULL,
                             o = NULL){
  
  p.genomewidepc = rep(NA,length(pc.lim))
  
  for(i in 1:length(p.genomewidepc)) {
    X = abs(pca$rotation[ipat,o[pc.lim[i]]]/bippat[ipat])
    # pl({plot(X)})
    ngp = 20
    gp = rep(1:20,each=ceiling(length(X)/20))[1:length(X)]
    E = mean(X)
    O = sapply(1:ngp,function(GP)mean(X[gp==GP]))
    N = 1.*as.vector(table(gp))
    ssq  = sapply(1:ngp,function(GP)var(X[gp==GP]))/sqrt(N)
    CHISQ = sum((O-E)^2/ssq)
    p.genomewidepc[i] = pchisq(CHISQ,ngp-1,low=F)
    #cat("Done",i,"\n")
  }
  
  p.genomewidepc = cbind(paste0("PC",o[pc.lim]),-log10(p.genomewidepc))
  write.table(p.genomewidepc,file=paste0(prefix,"_genomewidePCtest.txt"),row=F,col=F,sep="\t",quote=F)
  
  return(p.genomewidepc)
}


