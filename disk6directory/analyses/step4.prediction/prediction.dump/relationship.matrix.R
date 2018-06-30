## Project Setup ####
library(reshape2); library(parallel); library(OmicKriging); library(glmnet)
#setwd("~/Documents/Grad_Projects/Breeding-Value-Prediction")
setwd("~/Documents/Grad_Projects/BV_Prediction/mnt/media/disk6/ARF/shared")

# Load relatedness matrix ####
relatedness2 <- acast(read.table("./snps/relmats/Q30.snps.nomiss.maf05.relatedness2",header = T), INDV1 ~ INDV2 , value.var="RELATEDNESS_PHI")
# Collect sample names
rel.colnames <- gsub(colnames(relatedness2),pattern = "_",replacement = "")
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
colnames(relatedness2) <- rel.colnames
rownames(relatedness2) <- rel.colnames
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
keep.index <-  which(colnames(relatedness2) %in% bio.phenos$animal_id)
relatedness2.subset <- relatedness2[keep.index,keep.index]

bio.pheno.subset <- bio.phenos[which(bio.phenos$animal_id %in% colnames(relatedness2.subset)),]

identical(bio.pheno.subset$animal_id[match(colnames(relatedness2.subset),bio.pheno.subset$animal_id)],colnames(relatedness2.subset))

bio.pheno.subset <- bio.pheno.subset[match(colnames(relatedness2.subset),bio.pheno.subset$animal_id),]

rownames(bio.pheno.subset) <- bio.pheno.subset$animal_id

##### Now test with OK ####
ok.phenos <- bio.pheno.subset; rownames(ok.phenos) <- 1:nrow(ok.phenos)
fam.phenos <- unique(ok.phenos[,c("fam_id","Volume")])
rownames(fam.phenos) <- fam.phenos$fam_id
u.fams <- unique(fam.phenos$fam_id)

# Predict LOO with bio reps ####
cor.mat <- relatedness2.subset
cor.mat <- cor(cor.mat)
rownames(cor.mat) <- 1:nrow(bio.pheno.subset); colnames(cor.mat) <- 1:nrow(bio.pheno.subset)

out <- c()
for(each.bio in 1:nrow(bio.pheno.subset)){
  tester <- as.character(each.bio)
  train <- as.character(c(1:nrow(ok.phenos))[-each.bio])
  out <- c(out,okriging(idtest = tester,idtrain = train,corlist = list(cor.mat),H2vec = c(.9),pheno = ok.phenos,phenoname = "Volume")[,2])
}
plot(out,bio.pheno.subset$Volume)

# Predict LOO families with bio reps ####
cor.mat <- relatedness2.subset
cor.mat <- cor(cor.mat)
rownames(cor.mat) <- 1:nrow(bio.pheno.subset); colnames(cor.mat) <- 1:nrow(bio.pheno.subset)

out <- c()
for(each.bio in 1:length(u.fams)){
  tester <- as.character(which(ok.phenos$fam_id %in% u.fams[each.bio]))
  train <- as.character(c(1:nrow(ok.phenos))[-which(ok.phenos$fam_id %in% u.fams[each.bio])])
  out <- c(out,mean(okriging(idtest = tester,idtrain = train,corlist = list(cor.mat),H2vec = c(.9),pheno = ok.phenos,phenoname = "Volume")[,2]))
}
plot(out,fam.phenos$Volume)
abline(h=mean(fam.phenos$Volume),v=mean(out))
cor(out,fam.phenos$Volume)^2



# Create fam mean relationship matrix ####
out <- vector("list")
for(family in 1:56){
#subset the family
  these.rows <- which(bio.pheno.subset$fam_id == u.fams[family])
  this.fam1 <- relatedness2.subset[these.rows,]
# Now for that single family, create a empty vector
  fam.vector <- rep(0,length(u.fams))
# Replace the indivdiual with 1
  fam.vector[family] <- 1
# Find out which are still 0
  still.empty <- which(fam.vector == 0)
  if(length(this.fam1) < 193){
    for(each in 1:length(still.empty)){
    other.ufam <- u.fams[still.empty[each]]
    these.cols <- which(bio.pheno.subset$fam_id == other.ufam)
    fam.vector[each] <- mean((this.fam1[these.cols]))
  }} else {
  for(each in 1:length(still.empty)){
  other.ufam <- u.fams[still.empty[each]]
  these.cols <- which(bio.pheno.subset$fam_id == other.ufam)
  fam.vector[still.empty[each]] <- mean(unlist(this.fam1[,these.cols]))
  }}
  
  out[[family]] <- fam.vector
  print(family)}
fam.snp.mat <- do.call("rbind",out)

# Predict LOO with family mean relationships ####
rm.these <- c(12,27,33,52)
cc <- cor(fam.snp.mat[-rm.these,-rm.these])
colnames(cc) <- u.fams[-rm.these]
rownames(cc) <- u.fams[-rm.these]
test.phenos <- fam.phenos[-rm.these,]
r.fams <- u.fams[-rm.these]
out <- c()
for(each.fam in 1:length(r.fams)){
  tester <- as.character(r.fams[each.fam])
  train <- r.fams[-c(which(r.fams %in% tester))] 
  out <- c(out,(okriging(idtest = tester,idtrain = train,corlist = list(cc),H2vec = c(1),pheno = test.phenos,phenoname = "Volume")[,2]))
}
plot(out,test.phenos$Volume)
abline(h=mean(test.phenos$Volume),v=mean(out))
cor(out,test.phenos$Volume)^2







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

##### Now test with OK ####
ok.phenos <- bio.pheno.subset; rownames(ok.phenos) <- 1:nrow(ok.phenos)
fam.phenos <- unique(ok.phenos[,c("fam_id","Volume")])
rownames(fam.phenos) <- fam.phenos$fam_id
u.fams <- unique(fam.phenos$fam_id)

# Predict LOO with bio reps ####
cor.mat <- train.dat.bio
cor.mat <- cor(t(cor.mat))
rownames(cor.mat) <- 1:nrow(bio.pheno.subset); colnames(cor.mat) <- 1:nrow(bio.pheno.subset)

out <- c()
for(each.bio in 1:nrow(bio.pheno.subset)){
  tester <- as.character(each.bio)
  train <- as.character(c(1:nrow(ok.phenos))[-each.bio])
  out <- c(out,okriging(idtest = tester,idtrain = train,corlist = list(cor.mat),H2vec = c(.9),pheno = ok.phenos,phenoname = "Volume")[,2])
}
plot(out,bio.pheno.subset$Volume)

# Predict LOO families with bio reps ####
cor.mat <- out.012.subset
cor.mat <- cor(t(cor.mat))
rownames(cor.mat) <- 1:nrow(bio.pheno.subset); colnames(cor.mat) <- 1:nrow(bio.pheno.subset)

out <- c()
for(each.bio in 1:length(u.fams)){
  tester <- as.character(which(ok.phenos$fam_id %in% u.fams[each.bio]))
  train <- as.character(c(1:nrow(ok.phenos))[-which(ok.phenos$fam_id %in% u.fams[each.bio])])
  out <- c(out,mean(okriging(idtest = tester,idtrain = train,corlist = list(cor.mat),H2vec = c(1),pheno = ok.phenos,phenoname = "Volume")[,2]))
}
plot(out,fam.phenos$Volume)
abline(h=mean(fam.phenos$Volume),v=mean(out))
cor(out,fam.phenos$Volume)^2

# Check variance of snps ####
var.snps <- apply(out.012.subset,2,var)
summary(var.snps)
sd.snps <- apply(out.012.subset,2,sd)
mean.snps <- apply(out.012.subset,2,mean)

bcv.snps <- sd.snps/mean.snps

hist(bcv.snps)

## Pred with filters ####
rm.var2 <- which(var.snps > .2 )
rm.var <- which(bcv.snps < 3 )
cor.mat <- out.012.subset[,-c(unique(c(rm.var,rm.var2)))]
cor.mat <- cor(t(cor.mat))
rownames(cor.mat) <- 1:nrow(bio.pheno.subset); colnames(cor.mat) <- 1:nrow(bio.pheno.subset)
out <- c()
for(each.bio in 1:length(u.fams)){
  tester <- as.character(which(ok.phenos$fam_id %in% u.fams[each.bio]))
  train <- as.character(c(1:nrow(ok.phenos))[-which(ok.phenos$fam_id %in% u.fams[each.bio])])
  out <- c(out,mean(okriging(idtest = tester,idtrain = train,corlist = list(cor.mat),H2vec = c(.9),pheno = ok.phenos,phenoname = "Volume")[,2]))
}
plot(out,fam.phenos$Volume)
abline(h=mean(fam.phenos$Volume),v=mean(out))
cor(out,fam.phenos$Volume)^2



# Test glmnet #####
cor.mat <- as.matrix(out.012.subset[,-c(unique(c(rm.var,rm.var2)))])

gl.mod <- glmnet(x = cor.mat,y = bio.pheno.subset$Volume,alpha = .5)
gl.pred <- predict.glmnet(gl.mod,newx = cor.mat)
plot(gl.pred[,ncol(gl.pred)],bio.pheno.subset$Volume)

# Create fam relationships ####
fam.012.cormat <- matrix(,nrow=nrow(fam.phenos),ncol=ncol(cor.mat))
for(each.fam in 1:length(u.fams)) {
  these.rows <- which(bio.pheno.subset$fam_id %in% u.fams[each.fam])
  if(length(these.rows) < 2) { fam.012.cormat[each.fam,] <- unlist(cor.mat[these.rows,])
  } else {
    fam.012.cormat[each.fam,] <- apply(cor.mat[these.rows,],2,mean)
  }
}
rownames(fam.phenos) <- u.fams

## Get variance estimates####
fam.var.snps <- apply(fam.012,2,var)
fam.sd.snps <- apply(fam.012,2,sd)
fam.mean.snps <- apply(fam.012,2,mean)
fam.bcv.snps <- fam.sd.snps/fam.mean.snps
hist(fam.bcv.snps)

hist(fam.var.snps[which(fam.bcv.snps > 3)])

### Predict with these ####

#cc <- cor(t(fam.012[,which(fam.bcv.snps > 3)]))
cc <- cor(t(fam.012.cormat[,-c(25,  26,  30, 187, 203,  94)]))
colnames(cc) <- u.fams; rownames(cc) <- u.fams
test.phenos <- fam.phenos
out <- c()
for(each.fam in 1:length(u.fams)){
  tester <- as.character(u.fams[each.fam])
  train <- u.fams[-c(which(u.fams %in% tester))] 
  out <- c(out,(okriging(idtest = tester,idtrain = train,corlist = list(cc),H2vec = c(1),pheno = test.phenos,phenoname = "Volume")[,2]))
}
plot(out,test.phenos$Volume)
abline(h=mean(test.phenos$Volume),v=mean(out))
cor(out,test.phenos$Volume)^2


#### Pred with glment family #####
filter <- which(bcv.snps > 3 & var.snps < .2)
pred.out <- c()

for(each.bio in 1:length(u.fams)){
  tester <- which(fam.phenos$fam_id %in% u.fams[each.bio])
  this.this <- fam.012[-tester,filter]
  gl.mod <- glmnet(x = this.this,y = fam.phenos$Volume[-c(tester)],alpha = .5,standardize = F)
  gl.pred <- predict.glmnet(gl.mod,fam.012[,filter])
  pred.out <- c(pred.out,gl.pred[each.bio,ncol(gl.pred)])
}
print(cor(pred.out,fam.phenos$Volume)^2)
plot(pred.out,fam.phenos$Volume)
abline(h=mean(fam.phenos$Volume),v=mean(pred.out))



alpha.seq <- seq(0,1,.01)
out.cor <- vector()
for(each.seq in 1:length(alpha.seq)) {
  pred.out <- c()
  thisa <- alpha.seq[each.seq]
      for(each.bio in 1:length(u.fams)){
        tester <- which(fam.phenos$fam_id %in% u.fams[each.bio])
        this.this <- fam.012.cormat[-tester,]
        gl.mod <- glmnet(x = this.this,y = fam.phenos$Volume[-c(tester)],alpha = .5,standardize = F)
        gl.pred <- predict.glmnet(gl.mod,fam.012.cormat)
        pred.out <- c(pred.out,gl.pred[each.bio,ncol(gl.pred)])
        }
  print(cor(pred.out,fam.phenos$Volume)^2)
  out.cor <- c(out.cor,cor(pred.out,fam.phenos$Volume)^2)
}
plot(alpha.seq,out.cor)


library(caret)
pred.out <- c()
for(each in 1:56){
out <- knnreg(x = fam.012.cormat[-each,],y=fam.phenos$Volume[-each],k = 13)
pred.out <- c(pred.out,predict(out,fam.012.cormat)[each])
}
plot(pred.out,fam.phenos$Volume)
colnames(fam.012.cormat) <- 1:ncol(fam.012.cormat)
pcn <- avNNet(x = fam.012.cormat[-1,],y = fam.phenos$Volume[-1],size=1,repeats = 1,bag = T)
predict(pcn,fam.012.cormat,type = "raw")

library(caretEnsemble)

# Stacking Algorithms - Run multiple algos in one call.
trainControl <- trainCotrol(method="boot",verboseIter = T,allowParallel = T,number=1, savePredictions="all")
                             

algorithmList <- c('earth','bagEarth','glmnet','xgbDART','xgbLinear','glmboost','treebag','rf','logicBag')
algorithmList <- c('bayesglm','BstLm','gamboost')
algorithmList <- c('bayesglm','glm')'svmLinear2'superpc'rrlda, rlm ,foba,ridge,RRF,rfRules,glm.nb
set.seed(100)
this.dat <- data.frame(cbind("Volume"=fam.phenos$Volume,fam.012.cormat),stringsAsFactors=F)
models <- caretList(Volume ~ ., data=this.dat, trControl=trainControl, methodList=algorithmList) 
results <- resamples(models)
summary(results)
# Box plots to compare models
scales <- list(x=list(relation="free"), y=list(relation="free"))
bwplot(results, scales=scales)

library(randomForest)

pred <- c()
for(each in 1:56){
 these.rows <- which(bio.pheno.subset$fam_id %in% u.fams[each])
dat <- cor.mat[-these.rows,]
ph <- bio.pheno.subset$Volume[-these.rows]
rf.mod <- randomForest(x = dat,y=ph,do.trace = T,ntree = 300,nPerm = 2,)
pred <- c(pred,predict(rf.mod,cor.mat)[these.rows])
}
plot(pred,bio.pheno.subset$Volume)


plot(predict(rf.mod,fam.012.cormat),fam.phenos$Volume)
