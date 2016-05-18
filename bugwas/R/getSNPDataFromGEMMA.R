#' This function retrieves the binary data from GEMMA input files and phenotype files.
#' 
#' This function generates the Manhattan plot(s) for a SNP GWAS.
#' @param gemmaGenFile The path to the GEMMA gen input file. It is a required input.
#' @param gemmaSnpFile The path to the GEMMA snp input file. It is a required input.
#' @param pheno The path of the phenotype file with contains sample ids in a column and phenotypes in another. The columns containing the ids and phenotypes have the repective headings id and phenotype.
#' @param prefix Prefix of the output files.
#' @keywords GEMMA
#' @keywords Data
#' @export
#' @examples
#' data <- getSNPDataFromGEMMA(gemmaGenFile = gemmaGenFile, gemmaSnpFile = gemmaSnpFile, 
#'								pheno = pheno, prefix = prefix)
#' @return A list of information on the biallelic snps and phenotypes.
getSNPDataFromGEMMA <- function(gemmaGenFile = NULL, 
								gemmaSnpFile = NULL, 
								pheno = NULL,
								prefix = NULL){
	
	
	# Read in sample IDs and phenotype
	pheno <- read.table(pheno,header=T, as.is = T, sep="\t")
	if(!any(colnames(pheno)=="id")){
		stop("\nError: phenotype file must have 'id' column")
	}
	if(!any(colnames(pheno)=="phenotype")){
		stop("\nError: phenotype file must have 'phenotype' column")
	}
	
	sample_ID <- pheno$id
	npcs <-length(sample_ID)
	pheno <- pheno$phenotype
	
	warning("This method assumes that the phenotypes are in the same order as the gemma input files.")
	
	
	
	# Run compact_SNPs function
	XX.all <- getDataFromGemma(prefix = prefix, gemmaGenFile = gemmaGenFile, gemmaSnpFile = gemmaSnpFile, id = sample_ID)
	
	return(list("XX.all" = XX.all, "sample_ID" = sample_ID,
				"npcs" = npcs, "y" = pheno, "XX.ID" = sample_ID))	
}