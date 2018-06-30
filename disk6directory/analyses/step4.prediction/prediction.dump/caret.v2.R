# Caret
library(caretEnsemble);library(caret); library(Metrics);library(parallel)
#setwd("~/Documents/Grad_Projects/Breeding-Value-Prediction")
setwd("~/Documents/Grad_Projects/BV_Prediction/mnt/media/disk6/ARF/shared")

### Try 012 ####
out.012 <- read.table("./snps/012/Q30.snps.nomiss.maf05.012",header = F)
rownames(out.012) <- as.character(read.table("./snps/012/Q30.snps.nomiss.maf05.012.indv")[,1])
colnames(out.012) <- paste0(as.character(read.table("./snps/012/Q30.snps.nomiss.maf05.012.pos")[,1]),"_",read.table("./snps/012/Q30.snps.nomiss.maf05.012.pos")[,2])


rel.colnames <- gsub(rownames(out.012),pattern = "_",replacement = "")
# Identify those matching the EW samples and replace the following statement
rel.colnames.ew <- gsub(rel.colnames,pattern = ".ew.bio.rep.merge.bam",replacement = "")
# We need to add "1000" to these names in order to match the phenos, so here just subset the ones which when the string is forced to be numeric, it has numbers
ew.bio.reps <- which(!is.na(as.numeric(rel.colnames.ew)))
# The remaining 192 samples are LGEP
lgep.bio.reps <- c(1:192)[-ew.bio.reps]
# Now add 1000 to the EW sample names
mc.ew <- as.numeric(rel.colnames.ew[ew.bio.reps]) + 1000


# Subset the lgep bio reps, replace the matching pattern, and remove the index character attached to the sample number
mc.lgep <- rel.colnames[lgep.bio.reps]
mc.lgep <- gsub(pattern = "*.lgep.bio.rep.merge.bam",replacement = "",x = mc.lgep)
mc.lgep <- substr(mc.lgep,1,nchar(mc.lgep)-1)

# Finally replace the original column names with the adjusted colnames and remove vars
rel.colnames[ew.bio.reps] <- as.character(mc.ew)
rel.colnames[lgep.bio.reps] <- mc.lgep
rel.colnames
rownames(out.012) <- rel.colnames
rm(rel.colnames,rel.colnames.ew,ew.bio.reps,lgep.bio.reps,mc.ew,mc.lgep)

## Load phenotype data ####
load("~/Documents/Grad_Projects/Breeding-Value-Prediction/disk6directory/resources/expt.dat.720.RData")
str(expt.dat.720)
expt.dat.720$animal_id <- as.character(expt.dat.720$animal_id)
expt.dat.720$fam_id <- as.character(expt.dat.720$fam_id)
expt.dat.720$batch <- as.character(expt.dat.720$batch)

# Create bio phenos from individuals with no missing Volume and who have unique, fam*animal*batch*vol combo
bio.phenos <- unique(expt.dat.720[which(!is.na(expt.dat.720$Volume) ==T),c("fam_id","animal_id","batch","Volume")])

## Match bio phenos to relmat ####
# Subset bio.phenos to only include individuals within the snp mat
keep.index <-  which(rownames(out.012) %in% bio.phenos$animal_id)
out.012.subset <- out.012[keep.index,]

bio.pheno.subset <- bio.phenos[which(bio.phenos$animal_id %in% rownames(out.012.subset)),]

identical(bio.pheno.subset$animal_id[match(rownames(out.012.subset),bio.pheno.subset$animal_id)],rownames(out.012.subset))

bio.pheno.subset <- bio.pheno.subset[match(rownames(out.012.subset),bio.pheno.subset$animal_id),]
rownames(bio.pheno.subset) <- bio.pheno.subset$animal_id

## Pred setup #####
fam.phenos <- unique(bio.pheno.subset[,c("fam_id","Volume")])
rownames(fam.phenos) <- fam.phenos$fam_id
u.fams <- unique(fam.phenos$fam_id)

# Filter BCV ####
# This near zero var will check for things which have less than 10 or less unique values and 95/5 most common/least common (19:1 ratio)
check.vars <- nearZeroVar(x = out.012.subset)
#check.vars2 <- nearZeroVar(x = out.012.subset,freqCut = 90/5 ,uniqueCut =10)
out.012.subset.0var <- out.012.subset[,check.vars2]

#test.find.cor <- findCorrelation(out.012.subset,verbose = T)
var.snps <- apply(out.012.subset.0var,2,var)
sd.snps <- apply(out.012.subset.0var,2,sd)
mean.snps <- apply(out.012.subset.0var,2,mean)
bcv.snps <- sd.snps/mean.snps

hist(bcv.snps[which(sd.snps > .3)])
hist(sd.snps[which(bcv.snps > 3)])
select.snps <- which(bcv.snps > 3 & sd.snps > .3)
train.dat.bio <- out.012.subset.0var[,select.snps]


fam.012 <- matrix(,nrow=nrow(fam.phenos),ncol=ncol(out.012.subset))
for(each.fam in 1:length(u.fams)) {
  these.rows <- which(bio.pheno.subset$fam_id %in% u.fams[each.fam])
  if(length(these.rows) < 2) { fam.012[each.fam,] <- unlist(out.012.subset[these.rows,])
  } else {
    fam.012[each.fam,] <- apply(out.012.subset[these.rows,],2,mean)
  }
}

fam.012.cormat <- matrix(,nrow=nrow(fam.phenos),ncol=ncol(train.dat.bio))
for(each.fam in 1:length(u.fams)) {
  these.rows <- which(bio.pheno.subset$fam_id %in% u.fams[each.fam])
  if(length(these.rows) < 2) { fam.012.cormat[each.fam,] <- unlist(train.dat.bio[these.rows,])
  } else {
    fam.012.cormat[each.fam,] <- apply(train.dat.bio[these.rows,],2,mean)
  }
}
rownames(fam.phenos) <- u.fam

### Caret ####
# Stacking Algorithms - Run multiple algos in one call.((


# Create groups for training
#grp.list <- vector("list")
#new.l <- 1:56
#for(each in 1:7){
#grp.list[[each]] <- sample(new.l,replace = F,size = 8)
#new.l <- new.l[-unlist(grp.list[[each]])]
#}

# Set up the train control
trainControl <- trainControl(method="repeatedcv",repeats=100,allowParallel = T, savePredictions="all",verboseIter = TRUE,search="grid")
#algorithmList <- c('earth','bagEarth','glmnet','xgbDART','xgbLinear','glmboost','treebag','rf','BstLm')
#algorithmList <- c('bayesglm','glm')'svmLinear2'superpc'rrlda, rlm ,foba,ridge,RRF,rfRules,glm.nb

algorithmList <- c('glmnet','glm')

models <- caretList(y=train.p,x=train.dat,  trControl=trainControl, methodList=algorithmList,continue_on_fail = T)
results <- resamples(models)
summary(results)

# Box plots to compare models
scales <- list(x=list(relation="free"), y=list(relation="free"))
bwplot(results, scales=scales)

pred <- vector()
for(check.samp in 1:56){
  set.seed(100)
train.p <- fam.phenos$Volume[-check.samp]
train.dat <- fam.012.cormat[-check.samp,-c(25,  26,  30, 187, 203,  94)]
colnames(train.dat) <- 1:ncol(train.dat)
#mod1 <- train(y=train.p,x = train.dat,trControl=trainControl,method="ridge")
mod1 <- train(y=train.p,x = train.dat,trControl=trainControl,method="glmnet")
#mod1 <- train(y=train.p,x = train.dat,trControl=trainControl,method="BstLm")
mod1
test.dat <- fam.012.cormat[,-c(25,  26,  30, 187, 203,  94)]; colnames(test.dat) <- colnames(train.dat)
pred <- c(pred,predict(mod1,test.dat)[check.samp])
plot(pred,fam.phenos$Volume[1:check.samp])
}
abline(h=mean(fam.phenos$Volume),v=mean(pred))
cor(pred,fam.phenos$Volume)^2
plot(pred-fam.phenos$Volume,fam.phenos$Volume)
plot(pred[-c(2,15,32,48)],fam.phenos$Volume[-c(2,15,32,48)])





## Caret pt1. 
```{r setup.pred}
crossover.fams <- unique(bio.pheno.subset[,c("fam_id","Volume","batch")])
crossover.fams <- names(which(table(crossover.fams$fam_id) ==2))

# Split test and train set by batch
ew.fams <- unique(bio.pheno.subset$fam_id[which(bio.pheno.subset$batch == "EW")])
lgep.fams <- unique(bio.pheno.subset$fam_id[which(bio.pheno.subset$batch == "LGEP")])

# Remove crossover fams from test set
ew.fams <- ew.fams[-which(ew.fams %in% crossover.fams)]

train.p <- fam.phenos$Volume[which(fam.phenos$fam_id %in% lgep.fams)]
train.dat <- fam.012.cormat[which(fam.phenos$fam_id %in% lgep.fams),]
colnames(train.dat) <- 1:ncol(train.dat)
```
trainControl <- trainControl(method="LOOCV",allowParallel = T, savePredictions="all",verboseIter = TRUE,search="grid")
pred <- vector()
for(check.samp in 1:56){
  set.seed(100)
  #train.p <- fam.phenos$Volume[which(fam.phenos$fam_id %in% lgep.fams)]
  #train.dat <- fam.012.cormat[which(fam.phenos$fam_id %in% lgep.fams),]
  train.p <- fam.phenos$Volume[-check.samp]
  train.dat <- fam.012.cormat[-check.samp,]
  colnames(train.dat) <- 1:ncol(train.dat)
  #mod1 <- train(y=train.p,x = train.dat,trControl=trainControl,method="ridge")
  mod1 <- train(y=train.p,x = train.dat,trControl=trainControl,method="glmnet")
  #mod1 <- train(y=train.p,x = train.dat,trControl=trainControl,method="BstLm")
  #mod1
  test.dat <- fam.012.cormat; colnames(test.dat) <- colnames(train.dat)
  pred <- c(pred,predict(mod1,test.dat)[check.samp])
  #pred <- c(pred,predict(mod1,test.dat)[which(fam.phenos$fam_id %in% ew.fams)])
  plot(pred,fam.phenos$Volume[1:check.samp])
}
abline(h=mean(fam.phenos$Volume),v=mean(pred))

