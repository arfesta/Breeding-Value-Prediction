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

salmon.raw.lengths <- data.frame("animal.id" = expt.dat.720$animal.id,t(salmon.counts$length),stringsAsFactors = F) %>%
  inner_join(x = bio.exp.dat,y = ., by=c("animal.id")) %>%
  subset(x=.,select=-c(animal.id,ID,vol,ht))
salmon.raw.lengths <- data.frame(salmon.raw.lengths)

#Normalize all of the data together####
salmon.raw.counts <- data.frame(apply(salmon.raw.counts,2,as.integer))
raw.sum <- apply(salmon.raw.counts,2,sum)
rm.txpt <- which(raw.sum <= 1)
salmon.raw.counts <- salmon.raw.counts[,-rm.txpt]
salmon.raw.lengths <- salmon.raw.lengths[,-rm.txpt]
expt.dat.tmp <- data.frame(expt.dat.576[,c("animal.id","true.index","lane.id")],stringsAsFactors = T)


og.txi <- vector("list"); og.txi$counts = t(salmon.raw.counts); og.txi$countsFromAbundance="no"; og.txi$length = t(salmon.raw.lengths)

expt.dat.tmp <- within(expt.dat.tmp,expr = {
  animal.id <- factor(animal.id); true.index <- factor(true.index); lane.id <- factor(lane.id)
})


dds <- DESeqDataSetFromTximport(txi=og.txi,colData = expt.dat.tmp,design = ~ animal.id)
salmon.counts <- NULL; rm(salmon.counts); gc()
dds.all <- DESeq(dds,parallel = T,fitType = "mean")
