# Set WD to top level RNA seq analysis####
setwd("/repos/Breeding-Value-Prediction/")
packages <- c("readr","dplyr","parallel","reshape2","readxl","tibble", "DESeq2")
sapply(packages,library,character.only=T)
rm(packages);gc()

# Load salmon count data and experimental data####
load("/mnt/ARF/step1/salmon_count_data.RData")
load("/mnt/ARF/step1/salmon_expt_data.RData") 
expt.dat.720 <- as_tibble(exp.dat); rm(exp.dat)

# Load techincal rep data frame with volume information for all 720 samples####
bio.exp.dat <- as_tibble(read_excel("./disk6directory/resources/phenos/All.xlsx")) %>%
  mutate(vol = round(vol,digits = 1)) %>%
  filter(ht > 0) %>%
  dplyr::select(animal,ID,vol,ht) %>%
  dplyr::rename(animal.id=animal) %>%
  distinct()
head(bio.exp.dat)

# Merge experimental data from salmon with volume and family ID's####
expt.dat.576 <- inner_join(x = bio.exp.dat,y = expt.dat.720,by="animal.id")
expt.dat.576$lane.id[1:144] <- paste("ew_",expt.dat.576$lane.id[1:144],sep="")
expt.dat.576$lane.id[140:145]


# Merge raw counts with volume and family ID's (to ensure correct order as above join)####
salmon.counts$counts <- data.frame("animal.id" = expt.dat.720$animal.id,t(salmon.counts$counts),stringsAsFactors = F) %>%
  inner_join(x = bio.exp.dat,y = ., by=c("animal.id")) %>%
  subset(x=.,select=-c(animal.id,ID,vol,ht))
salmon.counts$counts <- data.frame(salmon.counts$counts)

salmon.counts$length <- data.frame("animal.id" = expt.dat.720$animal.id,t(salmon.counts$length),stringsAsFactors = F) %>%
  inner_join(x = bio.exp.dat,y = ., by=c("animal.id")) %>%
  subset(x=.,select=-c(animal.id,ID,vol,ht))
salmon.counts$length <- data.frame(salmon.counts$length)

# Prep counts and col dat####
salmon.counts$counts <- data.frame(apply(salmon.counts$counts,2,as.integer))
expt.dat.tmp <- data.frame(expt.dat.576[,c("animal.id","ID","true.index","lane.id")],stringsAsFactors = T)
expt.dat.tmp$tech.id <- rep(1:3,times=192)
expt.dat.tmp <- within(expt.dat.tmp, {
  animal.id <- factor(animal.id); ID <- factor(ID); true.index <- factor(true.index); lane.id <- factor(lane.id); tech.id <- factor(tech.id)
})
salmon.counts$abundance = NULL; salmon.counts$infReps = NULL;

#toy example: #####
toy.cts <- salmon.counts$counts[,3232:3240]
toy.lng <- salmon.counts$length[,3232:3240]

toy.tx <- vector("list"); toy.tx$counts <- t(toy.cts); toy.tx$length <- t(toy.lng); toy.tx$countsFromAbundance = "no"

test.design <- model.matrix(~ 0+ true.index + lane.id + ID,expt.dat.tmp)
rm.cols <- which(apply(test.design,2,sum) == 0)
if(length(rm.cols) > 0) {test.design <- test.design[,-rm.cols]}
ncol(test.design)
head(test.design)


dds.toy <- DESeqDataSetFromTximport(txi=toy.tx,colData = expt.dat.tmp,design = ~ 0 + true.index + lane.id + ID)

gc()
dds.toy <- DESeq(dds.toy,parallel = T)

resultsNames(dds.toy)
res <- results(dds.toy,contrast = c("ID","N07002","N11061"))
res

norm.cts <- counts(dds.toy,normalized=T)

# Create DESeq object ####
salmon.counts$counts <- t(salmon.counts$counts); salmon.counts$length = t(salmon.counts$length)
dds <- DESeqDataSetFromTximport(txi=salmon.counts,colData = expt.dat.tmp,design = ~ 0 + true.index + lane.id + ID)
gc()


nrow(dds)
#86008
length(which(apply(counts(dds),1,min) > 1))
dds.test <- dds[apply(counts(dds),1,min) > 1,]
nrow(dds.test)
#21561

vsd <- vst(dds.test,blind=FALSE)

dds.test <- estimateSizeFactors(dds.test)
nl2 <- log2(counts(dds.test, normalized=TRUE))

library(ggplot2)
df <- as_data_frame( log2(counts(dds.test, normalized=TRUE))[,c(1,432)]) %>% mutate(transformation = "log2(x + 1)")
df2 <- as_data_frame(assay(vsd)[,c(1,432)]) %>% mutate(transformation = "vst")
colnames(df) <- colnames(df2)
df <- rbind(df,df2)
colnames(df)[1:2] <- c("x", "y")  

ggplot(df, aes(x = x, y = y)) + geom_hex(bins = 80) +
  coord_fixed() + facet_grid( . ~ transformation)  


library("pheatmap")
library("RColorBrewer")

sampleDists.vsd <- as.matrix(dist(t(assay(vsd))))
sampleDists.nl2 <- as.matrix(dist(t(nl2)))

sampleDistMatrix <- sampleDists.nl2
  #rownames(sampleDistMatrix) <- paste( vsd$ID, vsd$animal.id, sep = " - " )
rownames(sampleDistMatrix) <- paste( vsd$ID)
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)

plotPCA(nl2, intgroup = c("dex", "cell"))

#dds <- DESeq(dds,parallel = T)

#save.image(file = "/mnt/deseq2_normalization.576.RData",compress=T)

##Remove offset and use only counts#####
rm(nl2,df,salmon.counts)

dds <- DESeqDataSetFromMatrix(countData =salmon.counts$counts,colData = expt.dat.tmp,design = ~ 0 + true.index + lane.id + ID)
gc()
dds <- DESeq(dds,parallel = T)
