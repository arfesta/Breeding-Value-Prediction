---
title: "Prediction of breeding values using LMM"
author: "Adam Festa"
date: "5/8/2018"
output: html_document
editor_options: 
  chunk_output_type: console
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

## Script objective


```{r cars}
summary(cars)
```

## Load data

```{r pressure, echo=FALSE}
# Load counts
load("/mnt/media/disk6/ARF/RNASEQ/counts/86kSalmon/Step3_LMM_Counts.RData")

# Load phenos
load("/mnt/media/disk6/ARF/RNASEQ/Breeding-Value-Prediction/disk6directory/resources/expt.dat.720.RData")
# Change factor to character and subset 3 columns
expt.dat.720$fam_id <- as.character(expt.dat.720$fam_id)
select.columns <- c("fam_id","Volume","Height")

# Remove missing phenotypes and subset for phenos which have counts
fam.phenos <- unique(expt.dat.720[-which(is.na(expt.dat.720$Volume) == T),select.columns])
fam.phenos <- fam.phenos[which(fam.phenos$fam_id %in% rownames(asreml.counts)),]
rownames(fam.phenos) <- fam.phenos$fam_id

# Subset each count data to contain only rows which have phenos
select.rows <- which(rownames(asreml.counts) %in% fam.phenos$fam_id)

asreml.counts <- asreml.counts[select.rows,]
```


## Summary stats

```{r}
sum.of.counts.per.gene <- apply(asreml.counts,2,sum)
sum.of.counts.per.fam <- apply(asreml.counts,1,sum)
hist(sum.of.counts.per.fam)
```

## Predictions

### Omic Kriging

```{r}
suppressPackageStartupMessages(library(OmicKriging))
true.val <- sapply(1:57,function(each.fam){
 fam.phenos$Volume[each.fam] - mean(fam.phenos$Volume[-each.fam])})
```

#### All transcripts

```{r}
cor.mat <- cor(t(log2(asreml.counts + 1)))

the.predictions <- sapply(1:57,function(each.fam){
  test.fam <- rownames(fam.phenos)[each.fam]
  train.fam <- rownames(fam.phenos)[-each.fam]
  okriging(idtest = test.fam,idtrain = train.fam,corlist = list(cor.mat),H2vec = c(.99),pheno = fam.phenos,phenoname = "Volume")[,2] - mean(fam.phenos$Volume[-each.fam])
})

cor(the.predictions,true.val)^2
```
