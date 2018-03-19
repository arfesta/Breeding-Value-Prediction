##Filter Methods####
library(parallel); library(genefilter)
load("/mnt/ARF/step2.normalization/192.expt.data.Rdata")
load("/mnt/ARF/step2.normalization/192.norm.log.count.lig.Rdata")

##Create Family mean counts and phenos####
fam.phenos <- data.frame(unique(expt.dat.192[,c("vol","ID")]))
rownames(fam.phenos) <- 1:57

salmon.fam.norm.counts <- mclapply(1:57,function(each.fam){
  these.rows <- which(expt.dat.192$ID %in% fam.phenos$ID[each.fam])
  if(length(these.rows) == 1){
    2^salmon.bio.norm.counts[these.rows,]
  } else {
    apply(2^salmon.bio.norm.counts[these.rows,],2,mean)
  }
},mc.cores =30)
salmon.fam.norm.counts <- do.call(rbind,salmon.fam.norm.counts)
rownames(salmon.fam.norm.counts) <- fam.phenos$ID
colnames(salmon.fam.norm.counts) <- colnames(salmon.bio.norm.counts)

##Generate general summary stats on txpts####
fam.stats <- vector("list")
fam.stats$txpt.sum <- apply(salmon.fam.norm.counts,2,sum)
fam.stats$txpt.sd <- apply(salmon.fam.norm.counts,2,sd)
fam.stats$txpt.mean <- apply(salmon.fam.norm.counts,2,mean)
fam.stats$txpt.dispersion <- fam.stats$txpt.sd^2/fam.stats$txpt.mean 
fam.stats$txpt.coeff_var <- fam.stats$txpt.sd/fam.stats$txpt.mean 

these.txpts <- (which(fam.stats$txpt.dispersion > 1))
these.txpts2 <- which(fam.stats$txpt.coeff_var > 1)
length(which(these.txpts2 %in% these.txpts))
rm(these.txpts,these.txpts2)
##First filter on dispersion index####
fam.filt.exp <- salmon.fam.norm.counts[, which(fam.stats$txpt.dispersion > 1)]
dim(fam.filt.exp)

##Then filter on correlation greater than .99####
bigcor <- function(x, nblocks = 3, verbose = TRUE, ...){
  library(ff, quietly = TRUE)
  NCOL <- ncol(x)
  
  ## test if ncol(x) %% nblocks gives remainder 0
  if (NCOL %% nblocks != 0) stop("Choose different 'nblocks' so that ncol(x) %% nblocks = 0!")
  
  ## preallocate square matrix of dimension
  ## ncol(x) in 'ff' single format
  #corMAT <- ff(NA,vmode="single",dim = c(NCOL, NCOL))
  
  ## split column numbers into 'nblocks' groups
  SPLIT <- split(1:NCOL, rep(1:nblocks, each = NCOL/nblocks))
  
  ## create all unique combinations of blocks
  COMBS <- expand.grid(1:length(SPLIT), 1:length(SPLIT))
  COMBS <- t(apply(COMBS, 1, sort))
  COMBS <- unique(COMBS)
  
  ## iterate through each block combination, calculate correlation matrix
  ## between blocks and store them in the preallocated matrix on both
  ## symmetric sides of the diagonal
  # for (i in 1:nrow(COMBS)) {
  corMAT <-  lapply(1:nrow(COMBS),function(i){
    COMB <- COMBS[i, ]
    G1 <- SPLIT[[COMB[1]]]
    G2 <- SPLIT[[COMB[2]]]
    if (verbose) cat("Block", COMB[1], "with Block", COMB[2], "\n")
    flush.console()
    COR <- cor(x[, G1], x[, G2], ...)
    #corMAT[G1, G2] <- COR
    #corMAT[G2, G1] <- t(COR)
    #COR <- NULL
    return(COR)
  })
  
  gc()
  return(corMAT)}

ex.cor <- bigcor(x=(fam.filt.exp))

final.cor <- ex.cor[[1]]
final.cor2 <- cbind(final.cor,ex.cor[[2]],ex.cor[[3]])
final.cor <- NULL; gc()
final.cor3 <- cbind(t(ex.cor[[2]]),t(ex.cor[[4]]),t(ex.cor[[5]]))
final.cor4 <- cbind(t(ex.cor[[3]]),t(ex.cor[[5]]),t(ex.cor[[6]]))
ex.cor <- NULL ; rm(ex.cor); gc()

final.cor <- rbind(final.cor2,final.cor3,final.cor4)
final.cor2 <- NULL; final.cor3 <- NULL; final.cor4 <- NULL; gc()

#save(final.cor,file="/mnt/txpt.cor.mat.RData",compress=T)
txpt.index <- 1:ncol(fam.filt.exp)
keep <- c()
unique.txpts <- colnames(fam.filt.exp)
for(this.txpt in 1:length(unique.txpts)) {
rm.these <- which(abs(final.cor[txpt.index[this.txpt],]) >= .99 )[-which(txpt.index[[this.txpt]] %in% txpt.index[this.txpt])]
if(length(rm.these) == 0 ){
  txpt.index <- txpt.index
} else {
  txpt.index <- txpt.index[-c(rm.these)]
}}

fam.filt.exp2 <- fam.filt.exp[,c(txpt.index)]


#genefilter####
#Anova
#colFtests
#cv
#genefinder - can be used on smaller subset to identify groups of expression
#var.cutoff

#CV filter####
f1 <- cv(a=1,b=Inf,na.rm=T)
keep.these <- which(genefilter(t(fam.filt.exp2),f1) == T)

#Anova filter with ID (how well does a model with only family ID explain each txpt?)######
bio.filt.exp2 <- salmon.bio.norm.counts[,colnames(fam.filt.exp2)]
f4 <- as.factor(expt.dat.192$ID)
r4 = rowFtests(t(bio.filt.exp2), f4)
hist(p.adjust(r4[,2],method="bonferroni"))
length(which(p.adjust(r4[,2],method="bonferroni") < .01))

keep.these2 <- which(genefilter(test.dat,f2) == T)


ffun_combined <- filterfun(f1, f2)
keep.these3 <- which(genefilter(test.dat,ffun_combined) == T)




###GLMNET#####
tmp.dat.fam2 <- salmon.fam.norm.counts[,names(keep.these3)]
library(glmnet)
tmp.dat.fam2 <- fam.filt.exp2

ll <- (which(apply(fam.filt.exp2,2,sd )/apply(fam.filt.exp2,2,mean ) >= 1)) 
ll.1 <- ll[which(ll %in% which(p.adjust(r4[,2],method = "bonferroni") < .01))]
gl.pred <- unlist(mclapply(1:57,function(each.fam){
  gl.mod <- glmnet(x=log2(tmp.dat.fam2[-each.fam,ll.1]),y=(fam.phenos$vol[-each.fam]),alpha=1)
  
  the.p <- predict(gl.mod,log2(tmp.dat.fam2[,ll.1]),type="response")
  
  the.p[each.fam,ncol(the.p)] - (mean(fam.phenos$vol[-each.fam]))
},mc.cores=30))

plot(gl.pred,fam.phenos$vol,ylab="YTrue",xlab="YPred",pch=19)

cor(gl.pred,fam.phenos$vol)^2
