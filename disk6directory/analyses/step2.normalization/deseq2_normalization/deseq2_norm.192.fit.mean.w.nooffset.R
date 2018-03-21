#Set WD to top level RNA seq analysis####
setwd("/repos/Breeding-Value-Prediction/")
packages <- c("readr","dplyr","parallel","reshape2","readxl","tibble", "DESeq2")
sapply(packages,library,character.only=T)
rm(packages);gc()

#Load salmon count data and experimental data####
load("/mnt/ARF/step1/salmon_count_data.RData")
load("/mnt/ARF/step1/salmon_expt_data.RData") 
expt.dat.720 <- as_tibble(exp.dat); rm(exp.dat)

#Load techincal rep data frame with volume information for all 720 samples####
bio.exp.dat <- as_tibble(read_excel("./disk6directory/resources/phenos/All.xlsx")) %>%
  mutate(vol = round(vol,digits = 1)) %>%
  filter(ht > 0) %>%
  dplyr::select(animal,ID,vol,ht) %>%
  dplyr::rename(animal.id=animal) %>%
  distinct()
head(bio.exp.dat)

#Merge experimental data from salmon with volume and family ID's####
expt.dat.576 <- inner_join(x = bio.exp.dat,y = expt.dat.720,by="animal.id")
expt.dat.576$lane.id[1:144] <- paste("ew_",expt.dat.576$lane.id[1:144],sep="")
expt.dat.576$lane.id[140:145]

#Merge raw counts with volume and family ID's (to ensure correct order as above join)####
salmon.raw.counts <- data.frame("animal.id" = expt.dat.720$animal.id,t(salmon.counts$counts),stringsAsFactors = F) %>%
  inner_join(x = bio.exp.dat,y = ., by=c("animal.id")) %>%
  subset(x=.,select=-c(animal.id,ID,vol,ht))
salmon.raw.counts <- data.frame(salmon.raw.counts)

#Normalize all of the data together####
salmon.raw.counts <- data.frame(apply(salmon.raw.counts,2,as.integer))

expt.dat.tmp <- data.frame(expt.dat.576[,c("animal.id","true.index","lane.id")],stringsAsFactors = T)

bio.raw.counts <- mclapply(1:192,function(each.bio){
  these.rows <- which(expt.dat.tmp$animal.id %in% unique(expt.dat.tmp$animal.id)[each.bio])
  apply(salmon.raw.counts[these.rows,],2,sum)
},mc.cores=30)
bio.raw.counts <- do.call(rbind,bio.raw.counts)
dim(bio.raw.counts)

salmon.counts <- NULL; salmon.raw.counts <- NULL;
rm(salmon.counts,salmon.raw.counts); gc()

og.txi <- vector("list"); og.txi$counts = t(bio.raw.counts); og.txi$countsFromAbundance="no"

expt.dat.192 <- unique(expt.dat.576[,c("animal.id","ID","vol")])

expt.dat.192 <- within(expt.dat.192,expr = {
  animal.id <- factor(animal.id); ID <- factor(ID); vol <- as.numeric(vol)
})

dds <- DESeqDataSetFromMatrix(countData = og.txi$counts,colData = expt.dat.192,design = ~ 0 + ID)
gc()
dds.all <- DESeq(dds,parallel = T,fitType = "mean")

save.image(file="/mnt/deseq2_norm.192.filt.mean.w.nooffset.RData")


######

library(DESeq2)
library(OmicKriging)
new.vst.counts <- vst(object = dds.all,blind = F)
new.vst.countsblind <- vst(object = dds.all,blind = T)
new.vst.counts <- assay(new.vst.counts)
new.vst.countsblind <- assay(new.vst.countsblind)

test.cts <- new.vst.countsblind
cor.mat <- cor(test.cts)
rownames(cor.mat) <- 1:192; colnames(cor.mat) <- 1:192

tmp.pheno <- data.frame(expt.dat.192); rownames(tmp.pheno) <- 1:192

pred <- c(); truth <- c()
for(each.fam in 1:57) {
  these.rows <- which(expt.dat.192$ID %in% unique(expt.dat.192$ID)[each.fam])
  tmp <- okriging(idtest = these.rows,idtrain = c(1:192)[-these.rows],corlist = list(cor.mat),H2vec = c(1),pheno = tmp.pheno,phenoname = "vol")
  pred <- c(pred,mean(tmp[,2]))
  truth <- c(truth,mean(tmp[,3]))
}
plot(pred,truth)
cor(pred,truth)^2
#new.vst.counts with blind=F is better


new.cv <- new.sd/new.mean
hist(new.cv)
rm.txpts <- which(new.cv < .13)

test.cts <- (new.vst.counts[-rm.txpts,])
cor.mat <- cor(test.cts)
rownames(cor.mat) <- 1:192; colnames(cor.mat) <- 1:192

tmp.pheno <- data.frame(expt.dat.192); rownames(tmp.pheno) <- 1:192

pred <- c(); truth <- c()
for(each.fam in 1:57) {
  these.rows <- which(expt.dat.192$ID %in% unique(expt.dat.192$ID)[each.fam])
  tmp <- okriging(idtest = these.rows,idtrain = c(1:192)[-these.rows],corlist = list(cor.mat),H2vec = c(1),pheno = tmp.pheno,phenoname = "vol")
  pred <- c(pred,mean(tmp[,2]))
  truth <- c(truth,mean(tmp[,3]))
}
plot(pred,truth)
cor(pred,truth)^2

########
new.fam.counts <- mclapply(1:57,function(each.fam){
  these.rows <- which(expt.dat.192$ID %in% unique(expt.dat.192$ID)[each.fam])
  if(length(these.rows) == 1){ new.counts[,these.rows]} else {
    apply(new.counts[,these.rows],1,mean)
  }
},mc.cores=25)
new.fam.counts <- do.call(rbind,new.fam.counts)
dim(new.fam.counts)

fam.phenos <- unique(expt.dat.192[,c("ID","vol")])
fam.phenos <- data.frame(fam.phenos); rownames(fam.phenos) <- 1:57


new.sum <- apply(new.fam.counts,2,sum)
rm.txpts <- (which(new.sum <= 576))
new.sd <- apply(new.fam.counts,2,sd)
rm.txpts <- (which(new.sd <= 1))
new.mean <- apply(new.fam.counts,2,mean)
new.var <- apply(new.fam.counts,2,var)
new.cv <- new.sd/new.mean
rm.txpts <- (which(new.sum < 1))

cor.mat <- cor(t(log2(new.fam.counts[,-rm.txpts] +1)))
rownames(cor.mat) <- 1:57; colnames(cor.mat) <- 1:57

fam.pred <- unlist(lapply(1:57,function(each.fam){
  okriging(idtest = each.fam,idtrain = c(1:57)[-each.fam],corlist = list(cor.mat),H2vec = c(1),pheno = fam.phenos,phenoname = "vol")[,2]
}))

plot(fam.pred,fam.phenos$vol)
cor(fam.pred,fam.phenos$vol)^2


rld <- rlog(dds.all)

new.assay <- assay(vsd)
dim(new.assay)


library(OmicKriging)
new.sum <- apply(new.assay,1,sum)
rm.txpts <- (which(new.sum <= 576))
new.sd <- apply(new.assay,1,sd)
rm.txpts <- (which(new.sd <= 1))
new.mean <- apply(new.assay,1,mean)
new.var <- apply(new.assay,1,var)
new.cv <- new.sd/new.mean
rm.txpts <- (which(new.var < 1))


cor.mat <- cor(t(scale(t(2^new.assay[-rm.txpts,]),center=T,scale=T)))
rownames(cor.mat) <- 1:192; colnames(cor.mat) <- 1:192

tmp.pheno <- data.frame(expt.dat.192); rownames(tmp.pheno) <- 1:192

pred <- c(); truth <- c()
for(each.fam in 1:57) {
  these.rows <- which(expt.dat.192$ID %in% unique(expt.dat.192$ID)[each.fam])
  
  tmp <- okriging(idtest = these.rows,idtrain = c(1:192)[-these.rows],corlist = list(cor.mat),H2vec = c(1),pheno = tmp.pheno,phenoname = "vol")
  pred <- c(pred,mean(tmp[,2]))
  truth <- c(truth,mean(tmp[,3]))
}
plot(pred,truth)
cor(pred,truth)^2




new.fam.vst.counts <- mclapply(1:57,function(each.fam){
  these.rows <- which(expt.dat.192$ID %in% unique(expt.dat.192$ID)[each.fam])
  if(length(these.rows) == 1){ new.assay[,these.rows]} else {
    apply(new.assay[,these.rows],1,mean)
  }
},mc.cores=25)
new.fam.vst.counts <- do.call(rbind,new.fam.vst.counts)



fam.phenos <- unique(expt.dat.192[,c("ID","vol")])
fam.phenos <- data.frame(fam.phenos); rownames(fam.phenos) <- 1:57


new.sum <- apply(new.fam.vst.counts,2,sum)
rm.txpts <- (which(new.sum <= 576))
new.sd <- apply(new.fam.vst.counts,2,sd)
rm.txpts <- (which(new.sd <= 1))
new.mean <- apply(new.fam.vst.counts,2,mean)
new.var <- apply(new.fam.vst.counts,2,var)
new.cv <- new.sd/new.mean
rm.txpts <- (which(new.var <= 1))

cor.mat <- cor(t(scale(new.fam.vst.counts[,-rm.txpts],center=T,scale=T)))
rownames(cor.mat) <- 1:57; colnames(cor.mat) <- 1:57

fam.pred <- unlist(lapply(1:57,function(each.fam){
  okriging(idtest = each.fam,idtrain = c(1:57)[-each.fam],corlist = list(cor.mat),H2vec = c(1),pheno = fam.phenos,phenoname = "vol")[,2]
}))

plot(fam.pred,fam.phenos$vol)
cor(fam.pred,fam.phenos$vol)^2


