###################################################################################
## FUNCTION: modelFit                                                            ##
## Inputs: listOutput  = a single list of y, W, list of SNPs in W                ##
##         snpspos = SNP locations data.table                                    ##
##         genfile = file name for .gen file                                     ##
##         samplefile = file name for .sample file                               ##
##         bedfile = file name for .bed file, NO FILE EXTENSION NEEDED           ##
##         windowSize = window size in SNPs                                      ##
##         numSNPShift = number of SNPs to shift after considering a window      ##
##         ldThresh = LD threshold                                               ##
##         outfile = output file for pruned SNPs, NO EXTENSION                   ##
##         method = method with greatest CV R2                                   ##
##         covFile = covariates to regress out from y                            ##
##                                                                               ##
## Output: returns model information for a given gene                            ##
###################################################################################

modelFit <- function(listOutput,
                       snpspos,
                       genfile,
                       samplefile,
                       bedfile,
                       windowSize,
                       numSNPShift,
                       ldThresh,
                       outFile,
                       method,
                       sex = F,
                       covFile){
  
  ## TAKE IN OUTPUT FROM featureSelect
  W = as.matrix(listOutput$W)
  y = as.vector(listOutput$y.train)
  snpList = as.vector(listOutput$snps)
  
  ## CREATE PLINK FORMAT GENOTYPE FILES
  onlyThese <- subset(snpspos,snpid %in% snpList)
  W = W[,snpList %in% onlyThese$snpid]
  snpList = snpList[snpList %in% onlyThese$snpid]
  df = as.data.frame(matrix(nrow = 1,ncol =2))
  ids = row.names(W)
  geno = as.data.frame(matrix(ncol=nrow(W)+4,nrow = ncol(W)))
  colnames(geno) <- c('SNP','Pos','A1','A2',ids)
  geno[,5:ncol(geno)] = t(W)
  geno$SNP = snpList
  onlyThese <- snpspos[snpspos$snpid %in% geno$SNP,]
  geno <- geno[geno$SNP %in% snpspos$snpid,]
  onlyThese <- onlyThese[match(onlyThese$snpid,geno$SNP),]
  chr <- unlist(lapply(strsplit(onlyThese$chr,'r'),function(x) as.character(x[2])))
  geno$Pos <- onlyThese$pos
  geno$A1 <- unlist(lapply(strsplit(geno$SNP,':'),function(x) as.character(x[3])))
  geno$A2 <- unlist(lapply(strsplit(geno$SNP,':'),function(x) as.character(x[4])))
  chr_dosage <- cbind(chr,geno)
  rm(geno)
  
  
  new.levels <- c('1 0 0','0 1 0','0 0 1')
  matrix.alleles <- as.matrix(chr_dosage[,6:ncol(chr_dosage)] + 1)
  impute2.format <- matrix(new.levels[matrix.alleles],ncol=ncol(matrix.alleles))
  gen <- cbind(chr_dosage[,1:5],impute2.format)
  gen[is.na(gen)] <- '<NA>' 
  
  require(data.table)
  fwrite(gen,genfile,row.names=FALSE, 
         col.names = FALSE, quote = FALSE, sep='\t')
  
  sample <- as.data.frame(matrix(nrow=(ncol(chr_dosage)-4),ncol=5))
  colnames(sample) <- c('ID_1','ID_2','missing','gender','pheno')
  sample$ID_1[2:nrow(sample)] <- paste('1A',2:nrow(sample),sep='')
  sample$ID_2[2:nrow(sample)] <- colnames(chr_dosage)[6:ncol(chr_dosage)]
  sample$missing[2:nrow(sample)] <- 0
  sample$gender[2:nrow(sample)] <- 2
  sample$pheno[2:nrow(sample)] <- y
  sample[1,] <- c(0,0,0,'D','P')
  write.table(sample,samplefile,row.names=FALSE, 
              col.names = TRUE, quote = FALSE)
  
  system(paste('/nas/longleaf/home/abhattac/plink','--gen',genfile,'--sample',samplefile,'--make-bed','--out',bedfile))
  
  a = fread(paste0(bedfile,'.fam'))
  a$V5 = 2
  fwrite(a,paste0(bedfile,'.fam'),col.names=F,row.names=F,quote=F,sep='\t')
  
  ### h2 calculation is done with GCTA-LDMS. Follow code at
  ### https://cnsgenomics.com/software/gcta/#GREMLinWGSorimputeddata
  
  
  ### LD PRUNING WITH PLINK
  system(paste('plink','--bfile',bedfile,
               '--indep-pairwise',windowSize,numSNPShift,ldThresh,
               '--out',outFile))  
  
  file.remove(genfile,samplefile)
  
  s <- as.character(fread(paste0(outFile,'.prune.in'),header=F)$V1)
  
  W <- W[,which(snpList %in% s)]
  q <- snpList[snpList %in% s]
  snpList <- q[match(s,q)]
  onlyThese <- subset(onlyThese,snpid %in% snpList)
  snpList = snpList[snpList %in% onlyThese$snpid]
  
  system(paste('plink','--bfile',bedfile,
               '--exclude',paste0(outFile,'.prune.out'),'--make-bed',
               '--out',bedfile))
  
  ### FIT WITH ELASTIC NET 
  if (method == 'Elastic net'){
    require(glmnet)
    set.seed(1218)
    rr <- cv.glmnet(y = y,x=W,alpha=0.5,nfolds=5,standardize=F,intercept=T)
    df = data.frame(Feature=as.character(c('Intercept',snpList)),
                    Chromosome = as.character(c('-',onlyThese$chr)),
                    Position = as.numeric(c('-',onlyThese$pos)),
                    Beta=as.vector(coef(rr,s='lambda.min')))
    colnames(df) = c('Feature','Chromosome','Position','Beta')
    df = subset(df,df$Beta != 0)
    meth = 'Model trained with Elastic-net'
  }
  
  ### IF ELASTIC NET IS FULLY SPARSE, MODEL IS NON-INFORMATIVE
  if (nrow(df) == 1){method$Method[1] = 'BLUP'}
  
  ### FIT WITH rrBLUP
  if (method$Method[1] == 'BLUP'){  
    set.seed(1218)
    require(rrBLUP)
    if (ncol(W) != 0){
      rr <- mixed.solve(y,Z=W)
      df = data.frame(Feature=as.character(c('Intercept',snpList)),
                      Chromosome = as.character(c('-',onlyThese$chr)),
                      Position = as.numeric(c('-',onlyThese$pos)),
                      Beta=as.numeric(c(rr$beta,rr$u)))
      colnames(df) = c('Feature','Chromosome','Position','Beta')}
    if (ncol(W) == 0){
      df = data.frame(Feature=as.character(c('Intercept',snpList)),
                      Chromosome = as.character(c('-',onlyThese$chr)),
                      Position = as.numeric(c('-',onlyThese$pos)),
                      Beta=as.numeric(c(0,rep(0,length(snpList)))))
      colnames(df) = c('Feature','Beta')}
    meth = 'Model trained with rrBLUP-LMM'
    df = subset(df,df$Beta != 0)
  }
  
  ### FIT WITH GEMMA
  if (method$Method[1] == 'LMM'){
    gemmaName = bedfile
    system(paste0('gemma/bin/gemma -miss 1 -maf 0 -r2 1 -rpace 1000 -bfile ',
                  bedfile,' -bslmm 2 -o ',paste0(gemmaName,'test.gemma.lmm')))
    par.lmm = read.table(paste0(gemmaName,'test.gemma.lmm.param.txt'),
                         head=T,as.is=T)
    lmm.wt =  rep(NA,length(snpList))
    m = match( snpList, par.lmm$rs)
    m.keep = !is.na(m)
    m = m[m.keep]
    lmm.wt[m.keep] = (par.lmm$alpha + par.lmm$beta * par.lmm$gamma)[m]
    df = data.frame(Feature=as.character(snpList),
                    Chromosome = as.character(onlyThese$chr),
                    Position = as.numeric(onlyThese$pos),
                    Beta=as.numeric(lmm.wt))
    meth = 'Model trained with GEMMA-LMM'
  }
  
  
  return(list(BetaMatrix = df,
              Method = meth))
  
}