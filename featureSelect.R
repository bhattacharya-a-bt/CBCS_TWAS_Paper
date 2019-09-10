###################################################################################
## FUNCTION: featureSelect                                                       ##
## Inputs: pheno = full vector of expression for gene of interest                ##
##         onlyCis = T if only cisSNPs considered, F for full cis-trans eQTL     ##
##         snp = genotype dosages as data.table                                  ##
##         geneInt = gene of interest                                            ##
##         filesnpspos = file name for SNP locations                             ##
##         filegenepos = file name for gene locations                            ##
##         FDRcut = FDR-adjusted p-value cutoff for trans eSNPs                  ##
##                                                                               ##
## Output: list of y.train, W, list of SNPs in W.train                           ##
##                                                                               ##
##                                                                               ##
##                                                                               ##
## Note: For cross-validation with FDRcut > 0, eQTLs were constructed on         ##
##       5-fold partition of data to prevent over-fitting from                   ##
##       overall trans-eQTLs                                                     ##
###################################################################################


featureSelect <- function(pheno,
                          onlyCis,
                          snp,
                          geneInt,
                          filesnpspos,
                          filegenepos,
                          FDRcut=0,
                          trans_results){
  
  
  ## GATHER ALL CIS SNPS TO GENE
  require(data.table)
  snpspos = fread(filesnpspos)
  snpX = subset(snpspos,chr == 'chrX')
  genepos = fread(filegenepos)
  thisGene = subset(genepos,geneid==geneInt)
  theseSNPs = subset(snpspos,chr==thisGene$chr)
  theseSNPs = subset(theseSNPs,pos >= thisGene$left - 5e5)
  theseSNPs = subset(theseSNPs,pos <= thisGene$right + 5e5)
  snpListSmall = theseSNPs$snpid
  
  # COLLECT TRANS eSNPs
  trans <- fread(trans_results)
  trans <- subset(trans,gene == geneInt)
  eSNPs = trans$SNP[trans$FDR <= FDRcut]
  
  if (!onlyCis) { snpListSmall <- union(snpListSmall,eSNPs)}
  
  # DESIGN MATRIX
  colnames(snp)[1] = 'SNP'
  snp = snp[snp$SNP %in% snpListSmall,]
  snp = snp[!duplicated(snp$SNP)]
  W <- t(as.matrix(snp[snp$SNP %in% snpListSmall,-1]))
  snpListSmall <- snpListSmall[snpListSmall %in% snp$SNP]
  
  return(list(y = pheno,
              W = W,
              snps = snpListSmall))
  
}