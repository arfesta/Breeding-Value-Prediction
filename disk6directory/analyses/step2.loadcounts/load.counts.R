# This script is responsible for loading the raw count files produced by salmon
setwd("/media/disk6/ARF/RNASEQ/")
load("/media/disk6/ARF/RNASEQ/Breeding-Value-Prediction/disk6directory/resources/expt.dat.720.RData")
rm.cols <- c(6,8,9,10,11,34,35,36)
expt.dat.240 <- unique(expt.dat.720[,-c(rm.cols)])
expt.dat.240 <- droplevels.data.frame(expt.dat.240)
salmon_bio_path <- c(paste0("./counts/86kSalmon/bio_EW/Sample_",expt.dat.240$animal_id[1:96]),paste0("./counts/86kSalmon/bio_LGEP/Sample_",expt.dat.240$animal_id[97:240]))
files <- file.path(salmon_bio_path, "quant.sf")
names(files) <- paste0(expt.dat.240$animal_id,".")
all(file.exists(files))

library(tximport);library(readr)
txi <- tximport(files, type = "salmon", txOut = T,varReduce = T)
save(txi,file = "/media/disk6/ARF/RNASEQ/counts/86kSalmon/bio_rep_txi_count_data_object.RData",compress=T)
names(txi)
colnames(txi$counts)

# summary of counts
bio.rep.sum <- apply(txi$counts,2,sum)
summary(bio.rep.sum)
hist(bio.rep.sum)

library(dplyr)
expt.dat.240$bio.rep.sum = bio.rep.sum
sum.by.index <- expt.dat.240 %>%
  group_by(index_seq) %>%
  summarise(Frequency = sum(bio.rep.sum))
hist(sum.by.index$Frequency)

# Check out EW exclusive batch due to bad samples:

ew.bio.rep.sum <- apply(txi$counts[,1:96],2,sum)
ew.txpt.sum <- apply(txi$counts[,1:96],1,sum)
#hist(ew.tech.rep.sum)
#which.min(ew.tech.rep.sum)
#hist(ew.txpt.sum)
#plot(ew.txpt.sum)
#which.max(ew.txpt.sum)
summary(ew.txpt.sum)
rm.zero <- (which(ew.txpt.sum < 10))

mypc <- prcomp(x = t(txi$counts[-rm.zero,1:96]),scale. = T,center = T)
plot(mypc$x[,1],mypc$x[,2])
#Large chunk of outliers..which samps are they?
rm.cols <- sort(names(which(mypc$x[,1] > 0 & mypc$x[,2] < 0)))
rm.fams <- which(colnames(txi$counts) %in% rm.cols)

dim(txi$abundance); dim(txi$counts); dim(txi$variance); dim(txi$length)
txi$abundance = txi$abundance[,-rm.fams]
txi$counts = txi$counts[,-rm.fams]
txi$variance = txi$variance[,-rm.fams]
txi$length = txi$length[,-rm.fams]

expt.dat.227 <- droplevels.data.frame(expt.dat.240[-rm.fams,])


# Do same for LGEP
lgep.samps <- which(expt.dat.227$batch=="LGEP")
lgep.tech.rep.sum <- apply(txi$counts[,lgep.samps],2,sum)
lgep.txpt.sum <- apply(txi$counts[,lgep.samps],1,sum)
#hist(ew.tech.rep.sum)
#which.min(ew.tech.rep.sum)
#hist(ew.txpt.sum)
#plot(ew.txpt.sum)
#which.max(ew.txpt.sum)
summary(lgep.txpt.sum)
rm.zero <- (which(lgep.txpt.sum < 10))

mypc <- prcomp(x = t(txi$counts[-rm.zero,lgep.samps]),scale. = T,center = T)
plot(mypc$x[,1],mypc$x[,2])

rm.fams <- names(which(table(as.character(expt.dat.227$fam_id)) < 3))
rm.rows <- which(as.character(expt.dat.227$fam_id) %in% rm.fams)

dim(txi$abundance); dim(txi$counts); dim(txi$variance); dim(txi$length)
txi$abundance = txi$abundance[,-rm.rows]
txi$counts = txi$counts[,-rm.rows]
txi$variance = txi$variance[,-rm.rows]
txi$length = txi$length[,-rm.rows]

expt.dat.206 <- droplevels.data.frame(expt.dat.227[-rm.rows,])

library(DESeq2)
str(expt.dat.206)
rownames(expt.dat.206) <- colnames(txi$counts)
test.fam <- which(expt.dat.206$fam_id == c("N08001"))
test.dat <- droplevels.data.frame(expt.dat.206[test.fam,])

test.txi <- txi
test.txi$abundance = test.txi$abundance[,test.fam]
test.txi$counts = test.txi$counts[,test.fam]
test.txi$variance = test.txi$variance[,test.fam]
test.txi$length = test.txi$length[,test.fam]

test.dat$rep = as.factor(paste0(test.dat$fam_id,"_", test.dat$batch))

dds.txi.test <- DESeqDataSetFromTximport(test.txi,
                                   colData = test.dat,
                                   design = ~ index_seq + rep)

N101027


test.fam <- which(expt.dat.206$fam_id == c("N101027"))
test.dat <- droplevels.data.frame(expt.dat.206[test.fam,])

test.txi <- txi
test.txi$abundance = test.txi$abundance[,test.fam]
test.txi$counts = test.txi$counts[,test.fam]
test.txi$variance = test.txi$variance[,test.fam]
test.txi$length = test.txi$length[,test.fam]

test.dat$rep = as.factor(paste0(test.dat$fam_id,"_", test.dat$batch))

dds.txi.test <- DESeqDataSetFromTximport(test.txi,
                                   colData = test.dat,
                                   design = ~ rep)

test.fam <- which(expt.dat.206$fam_id == c("N111095"))
test.dat <- droplevels.data.frame(expt.dat.206[test.fam,])

test.txi <- txi
test.txi$abundance = test.txi$abundance[,test.fam]
test.txi$counts = test.txi$counts[,test.fam]
test.txi$variance = test.txi$variance[,test.fam]
test.txi$length = test.txi$length[,test.fam]

test.dat$rep = as.factor(paste0(test.dat$fam_id,"_", test.dat$batch))
test.dat$an = as.factor(c(1,2,3,1,2,3))

dds.txi.test <- DESeqDataSetFromTximport(test.txi,
                                   colData = test.dat,
                                   design = ~  rep)


library(BiocParallel)
register(MulticoreParam(50))
keep <- rowSums(counts(dds.txi.test)) >= 10
dds.txi.test <- dds.txi.test[keep,]
dds.txi.test <- DESeq(object = dds.txi.test,minReplicatesForReplace = Inf,parallel = T)
resultsNames(dds.txi.test)
res <- results(dds.txi.test,name = "rep_N111095_LGEP_vs_N111095_EW",
               lfcThreshold = 1.2,pAdjustMethod = "bonferroni",alpha = .05)
summary(res)
mcols(res$description)