#' Testing the genome-wide principal components for all PCs.
#' 
#' This function identifies the non-genome-wide principal components.
#' @param biallelic A list called 'biallelic' created from the bugwas function
#' @param config A list called 'config' created from the bugwas function
#' @keywords Bayesian-Wald-test
#' @keywords PCA
#' @keywords All
#' @export
#' @return The p-value of the Bayesian Wald test for the genome-wide effect of all principal components.
#' @examples
#' data <- lin_loc(gen = gen, pheno = pheno, phylo = phylo, 
#'  prefix = prefix, gem.path = gem.path)
#' testGenomeWidePCsAll(config = data$config, biallelic = data$biallelic)
testGenomeWidePCsAll = function(config = NULL, biallelic = NULL){

	prefix = config$prefix
	pca = biallelic$pca
	bippat = biallelic$bippat
	ipat = biallelic$pattern

	pcCount = ncol(pca$rotation)
	p.genomewidepc = rep(NA,pcCount)
  
  for(i in 1:pcCount) {
    X = abs(pca$rotation[ipat,i]/bippat[ipat])
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
  
  p.genomewidepc = cbind(paste0("PC", c(1:pcCount)),-log10(p.genomewidepc))
  write.table(p.genomewidepc,file=paste0(prefix,"_genomewidePCtestAll.txt"),row=F,col=F,sep="\t",quote=F)
  
  return(p.genomewidepc)
}
