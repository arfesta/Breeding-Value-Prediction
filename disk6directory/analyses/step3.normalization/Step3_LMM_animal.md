---
title: 'Alt Approaches Step 3: LMM Normalization'
author: "Adam Festa"
date: "April 22, 2018"
output: 
  html_document: 
    toc: yes
    theme: paper
    highlight: kate
editor_options: 
  chunk_output_type: console
---

## Load the txi counts and phenos 

* The technical replicate counts and phenotype data were previously prepared for normalization.  For more details on the preparation see Step 2 readme.

```
load("/media/titan/disk6/ARF/RNASEQ/counts/86kSalmon/Step2_TECH_Load_Counts.RData")
```

## Obtain pedigree information

* Generate a csv file of unique family names to obtain pedigree information by subsetting the unique ID's from the phenotype data frame 

```
unique.fams <- unique(unlist(strsplit(x = unique(as.character(load.counts_tech$phenos$fam_id)) , split = "x")))
unique.fams <- data.frame("ID"=unique.fams)
#Write out a CSV file which can be uploaded to TIPRoot to access pedigree
    #write_csv(unique.fams,col_names = T,path ="/media/disk6/ARF/RNASEQ/Breeding-Value-Prediction/disk6directory/resources/phenos/ufams.v2.csv")
```

* Steps to Generate Ped:
    * Accessed tip root
    * pedigree infomration --> pedigree and relatives
    *  upload csv written up and generate ancestors
    *  download file and put in /resources/pedigree/
    *  Open file and delete the comment column

## Load pedigree 

  Subset the unique families from the phenotype data fram
  
```
    unique.fams <- unique(unlist(strsplit(x = unique(as.character(load.counts_tech$phenos$fam_id)),split = "x")))
    unique.fams <- data.frame("ID"=unique.fams)
  pedigree <-read_xlsx("/media/titan/disk6/ARF/RNASEQ/Breeding-Value-Prediction/disk6directory/resources/pedigree/pedigree_2018_04_19_20_36_15.xlsx")
  as.character(unique.fams[,1])[which(!as.character(unique.fams[,1]) %in% as.character(pedigree$Id))]
  #Note *** The pedigree returned from tiproot does not contain any WG families
```

## Add missing OP families to the pedigree (WG fams)

```
    add.fams <- as.character(unique.fams[,1])[which(!as.character(unique.fams[,1]) %in% as.character(pedigree$Id))]
  #These are the families which are missing
    head(add.fams)
  #Create data frame for families which are missing from the pedigree
    missing.ped <- data.frame("Id"=add.fams,"Parent 1"=NA,"Parent 2"=NA,"County"=NA,"State"=rep("WG",length(add.fams)))
    colnames(missing.ped) <- c("Id","Parent 1","Parent 2","County","State")
  #Rbind the missing pedigree to the one obtained from tiproot
    pedigree <- rbind(pedigree,missing.ped)
    rm(unique.fams,add.fams,missing.ped)
  #Note**Some families are missing a State while most are WG one is N11011..what's up with that?
    pedigree$Id[which(is.na(pedigree$State) == T)]
```
## Create progeny pedigree

```
    prog.pedigree <- unique(load.counts_tech$phenos[,c("animal_id","p1","p2")])
    colnames(prog.pedigree) <- c("animal","sire","dam")
```

## Complete pedigree 

  Create A inverse from pedigree
  
```
    head(pedigree)
    head(prog.pedigree)
    asreml.ped <-data.frame(pedigree[,1:3])
    colnames(asreml.ped) <- c("animal","sire","dam")
    asreml.ped <- rbind(asreml.ped,prog.pedigree)
    asreml.ped$dam[which(asreml.ped$dam == "NA")] <- NA
    asreml.ped <- within(asreml.ped,{animal <- factor(animal); sire <- factor(sire); dam <- factor(dam)})
    rownames(asreml.ped) <- asreml.ped[,1]
library(asreml)
    ainv2 <- asreml.Ainverse(pedigree = asreml.ped)$ginv
```

## Normalize all of the data together

```
    pre.norm.counts <- data.frame(apply(load.counts_tech$txi_object$counts,2,as.integer))
      colnames(pre.norm.counts) <- load.counts_tech$phenos$animal_id
    pre.norm.length <- load.counts_tech$txi_object$length
      colnames(pre.norm.length) <- load.counts_tech$phenos$animal_id
    pre.norm.var <- load.counts_tech$txi_object$variance
      colnames(pre.norm.var) <- load.counts_tech$phenos$animal_id
    expt.dat.tmp <- data.frame(load.counts_tech$phenos[,c("animal_id","fam_id","index_seq","batch","lane")],stringsAsFactors = T)
```

```
pt <- proc.time()
    asreml.counts <- mclapply(1:nrow(pre.norm.counts),function(each.txpt){
      #each.txpt=323
    tip.subset <- data.frame(cbind("txpt"=unlist(pre.norm.counts[each.txpt,]),expt.dat.tmp,"length" = pre.norm.length[each.txpt,]))
      #hist(tip.subset$txpt)

    tip.subset <- within(tip.subset,{txpt <- log2(as.integer(txpt)+1); 
      animal_id <- factor(animal_id); 
      batch <- factor(batch);
      fam_id <- factor(fam_id); 
      true.index <- factor(index_seq); 
      length <- as.numeric(length);
      lane <- factor(lane)})

    tmp1 <- asreml(fixed= txpt ~1 + batch + lane + index_seq, random = ~ ped(animal_id,var=T,init=1),
               ginverse=list(animal_id=ainv2),data=tip.subset,maxiter=50,trace=F,
               family=asreml.gaussian(link = "identity"))
   
     tmp1$coefficients$random + tmp1$coefficients$fixed['(Intercept)']
      },mc.cores=23)
proc.time() - pt
```

153 seconds for 20K txpts (approx 2.5 minutes * 4 sets of 20k = ~10 minutes for normalization)


## Combine count output and save 

```
asreml.counts <- do.call(cbind,asreml.counts)
mypc <- prcomp(2^asreml.counts,scale. = T,center = T)
mypc <- mypc$x
plot(mypc[,1:2])

 #edit.names <- gsub(pattern = "fam_id, var = T)_",replacement = "",x = rownames(asreml.counts))
#      edit.names <- gsub("ped\\(", "", edit.names)
    
#rownames(asreml.counts) <- edit.names
colnames(asreml.counts) <- rownames(pre.norm.counts)
asreml.counts <- 2^asreml.counts
save(asreml.counts,file = "../Step3_LMM_animal.RData",compress=T)
```