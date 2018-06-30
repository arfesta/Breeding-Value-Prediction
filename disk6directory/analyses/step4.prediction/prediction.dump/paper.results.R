# RUN ONCE ####
#Loads raw snp data and saves with rownames and column names (ONLY RUN ONCE)
initial.loadsnps <- function(){
  bio.012.df <- read.table("./snps/012/Q30.snps.nomiss.maf05.012",header = F,row.names = 1)
  rownames(bio.012.df) <- as.character(read.table("./snps/012/Q30.snps.nomiss.maf05.012.indv")[,1])
  colnames(bio.012.df) <- paste0(as.character(read.table("./snps/012/Q30.snps.nomiss.maf05.012.pos")[,1]),
                                 "_",
                                 read.table("./snps/012/Q30.snps.nomiss.maf05.012.pos")[,2])
  
  # Remove any extra "_" characters in the sample names
  rel.colnames <- gsub(rownames(bio.012.df),pattern = "_",replacement = "")
  # Identify those matching the EW samples and replace the following statement
  rel.colnames.ew <- gsub(rel.colnames,pattern = ".ew.bio.rep.merge.bam",replacement = "")
  
  # We need to add "1000" to these names in order to match the phenos, 
  # so here just subset the ones which when the string is forced to be numeric, it has numbers
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
  head(rel.colnames)
  rownames(bio.012.df) <- rel.colnames
  rm(rel.colnames,rel.colnames.ew,ew.bio.reps,lgep.bio.reps,mc.ew,mc.lgep)
  save(bio.012.df,file="/mnt/bio.012.df.RData",compress=T)
}
## This function loads the SNP data created above as well as the Txpt data and gets them in matching order
load.bio.data <- function(){
  # Load counts
  load("/mnt/media/disk6/ARF/RNASEQ/counts/86kSalmon/Step3_LMM_animal.RData")
  head(rownames(asreml.counts))
  asreml.counts <- asreml.counts[-c(1:94),]
  head(rownames(asreml.counts))
  
  # Load phenos
  load("/mnt/media/disk6/ARF/RNASEQ/Breeding-Value-Prediction/disk6directory/resources/expt.dat.720.RData")
  # Change factor to character and subset 3 columns
  expt.dat.720$animal_id <- as.character(expt.dat.720$animal_id)
  expt.dat.720$fam_id <- as.character(expt.dat.720$fam_id)
  expt.dat.720$batch <- as.character(expt.dat.720$batch)
  select.columns <- c("fam_id","animal_id","Volume","Volume Accuracy","Height","index_seq","batch")
  
  # Remove missing phenotypes and subset for phenos which have counts
  bio.phenos <- unique(expt.dat.720[-which(is.na(expt.dat.720$Volume) == T),select.columns])
  
  # Edit rownames of asreml.counts
  sub.names <- gsub(pattern = "animal_id, var = T)_",replacement = "",x = rownames(asreml.counts))
  sub.names <- gsub("ped\\(", "", sub.names)
  
  # Subset each count data to contain only rows which have phenos
  select.rows <- which(sub.names %in% as.character(bio.phenos$animal_id))
  
  asreml.counts <- asreml.counts[select.rows,]
  
  select.rows <- which(as.character(bio.phenos$animal_id) %in% sub.names)
  bio.phenos <- bio.phenos[select.rows,]
  
  final.r.names <- gsub(pattern = "animal_id, var = T)_",replacement = "",x = rownames(asreml.counts))
  final.r.names <- gsub("ped\\(", "", final.r.names)
  head(as.character(bio.phenos$animal_id))
  head(final.r.names)
  
  rownames(asreml.counts) <- paste0(final.r.names,".",bio.phenos$index_seq)
  rownames(bio.phenos) <- paste0(bio.phenos$animal_id,".",bio.phenos$index_seq)
  
  rm(final.r.names,select.columns,select.rows,sub.names)
  
  ## Load SNPs
  ### Load SNPs
  load("/mnt/bio.012.df.RData")
  ## Load phenotype data ####
  #load("~/Documents/Grad_Projects/Breeding-Value-Prediction/disk6directory/resources/expt.dat.720.RData")
  load("/mnt/Breeding-Value-Prediction/disk6directory/resources/expt.dat.720.RData")
  expt.dat.720$animal_id <- as.character(expt.dat.720$animal_id)
  expt.dat.720$fam_id <- as.character(expt.dat.720$fam_id)
  expt.dat.720$batch <- as.character(expt.dat.720$batch)
  
  # Create bio phenos from individuals with no missing Volume and who have unique, fam*animal*batch*vol combo
  bio.phenos.snps <- unique(expt.dat.720[which(!is.na(expt.dat.720$Volume) ==T),c("fam_id","animal_id","batch","Volume")])
  
  
  # Index the snp data frame to only include individuals which have phenotypes
  keep.index <-  which(rownames(bio.012.df) %in% bio.phenos.snps$animal_id)
  
  # Keep only those rows
  bio.012.df <- bio.012.df[keep.index,]
  
  bio.phenos.snps <- bio.phenos.snps[which(bio.phenos.snps$animal_id %in% rownames(bio.012.df)),]
  identical(bio.phenos.snps$animal_id[match(rownames(bio.012.df),bio.phenos.snps$animal_id)],rownames(bio.012.df))
  
  bio.phenos.snps <- bio.phenos.snps[match(rownames(bio.012.df),bio.phenos.snps$animal_id),]
  rownames(bio.phenos.snps) <- bio.phenos.snps$animal_id
  rm(keep.index)
  
  crossover <- which(bio.phenos$animal_id %in% bio.phenos.snps$animal_id)
  
  bio.phenos <- bio.phenos[crossover,]
  asreml.counts <- asreml.counts[crossover,]
  
  
  bio.phenos$animal_id
  bio.phenos.snps$animal_id
  
  identical(bio.phenos.snps$animal_id,bio.phenos$animal_id[match(bio.phenos.snps$animal_id,bio.phenos$animal_id)])
  asreml.counts <- asreml.counts[match(bio.phenos.snps$animal_id,bio.phenos$animal_id),]
  
  rownames(asreml.counts) <- rownames(bio.012.df)
  bio.phenos <- bio.phenos[match(bio.phenos.snps$animal_id,bio.phenos$animal_id),]
  rownames(bio.phenos) <- rownames(bio.phenos.snps)
  
  save(bio.phenos,file="/mnt/bio.phenos.RData",compress=T)
  save(asreml.counts,file="/mnt/asreml.counts.RData",compress=T)
  save(bio.012.df,file="/mnt/bio.012.df.RData",compress=T)
  
  counts <- list("asreml.counts"=asreml.counts,"bio.phenos"=bio.phenos,"bio.012.df"=bio.012.df)
  return(counts)
}
bio.data <- load.bio.data()

# Using the bio.data, now create the family mean data
create.family.data <- function(bio.data.list=NULL){
  bio.phenos <- bio.data.list$bio.phenos
  asreml.counts <- bio.data.list$asreml.counts
  bio.012.df <- bio.data.list$bio.012.df
  
  fam.phenos <- unique(bio.phenos[,c("fam_id","Volume")])
  rownames(fam.phenos) <- fam.phenos$fam_id
  u.fams <- unique(fam.phenos$fam_id)
  save(fam.phenos,file="/mnt/fam.phenos.RData",compress=T)
  
  fam.012.df <- mclapply(1:56,function(each.fam){
    merge.bio.reps <-  which(bio.phenos$fam_id %in% u.fams[each.fam])
    if(length(merge.bio.reps) == 1){ bio.012.df[merge.bio.reps,]} else {
      apply(bio.012.df[merge.bio.reps,],2,mean)}
  },mc.cores=30)
  fam.012.df <- do.call(rbind,fam.012.df)
  rownames(fam.012.df) <- u.fams
  colnames(fam.012.df) <- colnames(bio.012.df)
  save(fam.012.df,file="/mnt/fam.012.df.RData",compress=T)
  
  fam.asreml.counts <- matrix(,nrow = 56,ncol=ncol(asreml.counts))
  for(each.fam in 1:56){
    merge.bio.reps <-  which(bio.phenos$fam_id %in% u.fams[each.fam])
    if(length(merge.bio.reps) == 1){ fam.asreml.counts[each.fam,] <- asreml.counts[merge.bio.reps,]} else {
      fam.asreml.counts[each.fam,] <- apply(asreml.counts[merge.bio.reps,],2,mean)}
  }
  rownames(fam.asreml.counts) <- u.fams
  save(fam.asreml.counts,file="/mnt/fam.asreml.counts.RData",compress=T)
  
  fam.data <- list("fam.012.df"=fam.012.df,"fam.phenos"=fam.phenos,"fam.asreml.counts"=fam.asreml.counts)
  return(fam.data)
}
fam.data <- create.family.data(bio.data.list = bio.data)

rm(bio.data,fam.data,initial.loadsnps,load.bio.data,create.family.data);gc()

## All files are saved as outputs and can be loaded below ####
load("/mnt/bio.012.df.RData")
load("/mnt/fam.012.df.RData")
load("/mnt/asreml.counts.RData")
load("/mnt/fam.asreml.counts.RData")
load("/mnt/bio.phenos.RData")
load("/mnt/fam.phenos.RData")
load("/mnt/fam.hm.012.df.RData")
load("/mnt/bio.hm.012.df.RData")

# With all starting files loaded we can now generate the first plot ####

# Create even test groups ####
# Index the phenotypes so they are ind decreasing order
pheno.order <- order(fam.phenos$Volume,decreasing = T)
first <- seq(1,56,7)
last <- seq(7,56,7)

se <- c(3223,42435,23534,5235,3423,1345,21541,35432)
pheno.groups <- lapply(1:8,function(each.group){
set.seed(se[each.group])
    sample(pheno.order[first[each.group]:last[each.group]])
})
pheno.groups <- do.call(rbind,pheno.groups)

# Each of the 7 groups contain 8 individuals spread across the population of phenotypes

# ALL TRANSCRIPTS OR SNPS ####
library(OmicKriging); library(Metrics)
txpt.cor <- cor(t(scale(log2(fam.asreml.counts+1),center=F)))
rownames(txpt.cor) <- rownames(fam.phenos); colnames(txpt.cor) <- rownames(fam.phenos)

snp.cor <- cor(t(fam.hm.012.df))
#snp.cor <- cor(t(fam.012.df))
rownames(snp.cor) <- rownames(fam.phenos); colnames(snp.cor) <- rownames(fam.phenos)

# Identify crossover families between batches
crossover.fams <- unique(bio.phenos[,c("fam_id","Volume","batch")])
crossover.fams <- names(which(table(crossover.fams$fam_id) ==2))
# Subset the unique list of families that are in the LGEP or EW
ew.fams <- unique(bio.phenos$fam_id[which(bio.phenos$batch == "EW")])
lgep.fams <- unique(bio.phenos$fam_id[which(bio.phenos$batch == "LGEP")])
# Remove LGEP families from the EW families (i.e. remove crossover families)
ew.fams <- ew.fams[-which(ew.fams %in% crossover.fams)]


set.seed(424)
random.txpt.df <- fam.asreml.counts[sample(1:56,56),]
rownames(random.txpt.df) <- fam.phenos$fam_id
txpt.random.cor <- cor(t(scale(log2(random.txpt.df + 1),center=F)))
colnames(txpt.random.cor) <- fam.phenos$fam_id; rownames(txpt.random.cor) <- fam.phenos$fam_id


set.seed(424)
random.snp.df <- fam.012.df[sample(1:56,56),]
rownames(random.snp.df) <- fam.phenos$fam_id
snp.random.cor <- cor(t(random.snp.df))
colnames(snp.random.cor) <- fam.phenos$fam_id; rownames(snp.random.cor) <- fam.phenos$fam_id


# IMAGE-7FOLD: All txpts 7-Fold with OmicKriging ####
all.predictions <- c();all.truth <- c(); group.cor <- c(); group.rmse <- c()
for(each.group in 1:7){
  test.group <- rownames(fam.phenos)[pheno.groups[,each.group]]
  train.group <- rownames(fam.phenos)[-pheno.groups[,each.group]]
  
 pred <-  okriging(idtest = test.group,idtrain = train.group,corlist = list(txpt.cor),H2vec = c(1),pheno = fam.phenos,phenoname = "Volume")
 all.predictions <- c(all.predictions,pred[,2]); all.truth <- c(all.truth,pred[,3])
 group.cor <- c(group.cor,cor(pred[,2],pred[,3],method="spearman"))
 group.rmse <- c(group.rmse,rmse(pred[,2],pred[,3]))
}
#plot(all.predictions,all.truth)
#plot(group.cor); plot(group.rmse)
#cor(all.predictions,all.truth)^2
gp1.rmse <- cbind("model"=rep("transcripts",7),"rmse"=group.rmse)
gp1.cor <- cbind("model"=rep("transcripts",7),"cor"=group.cor,"rmse"=group.rmse) 
all.predictions1 <- all.predictions
# IMAGE-7FOLD: All snps 7-fold with OmicKriging ####
all.predictions <- c();all.truth <- c(); group.cor <- c(); group.rmse <- c()
for(each.group in 1:7){
  test.group <- rownames(fam.phenos)[pheno.groups[,each.group]]
  train.group <- rownames(fam.phenos)[-pheno.groups[,each.group]]
  
  pred <-  okriging(idtest = test.group,idtrain = train.group,corlist = list(snp.cor),H2vec = c(1),pheno = fam.phenos,phenoname = "Volume")
  all.predictions <- c(all.predictions,pred[,2]); all.truth <- c(all.truth,pred[,3])
  group.cor <- c(group.cor,cor(pred[,2],pred[,3],method="spearman"))
  group.rmse <- c(group.rmse,rmse(pred[,3],pred[,2]))
  plot(pred[,2],pred[,3])
}
plot(all.predictions,all.truth)
plot(group.cor); plot(group.rmse)
cor(all.predictions,all.truth)^2
gp2.rmse <- cbind("model"=rep("snps",7),"rmse"=group.rmse)
gp2.cor <- cbind("model"=rep("snps",7),"cor"=group.cor,"rmse"=group.rmse)   
all.predictions2 <- all.predictions
#IMAGE-7fOLD: All snps + All txpts with OmicKriging ####
all.predictions <- c();all.truth <- c(); group.cor <- c(); group.rmse <- c()
for(each.group in 1:7){
  test.group <- rownames(fam.phenos)[pheno.groups[,each.group]]
  train.group <- rownames(fam.phenos)[-pheno.groups[,each.group]]
  
  pred <-  okriging(idtest = test.group,idtrain = train.group,corlist = list(snp.cor,txpt.cor),H2vec = c(.98,.02),pheno = fam.phenos,phenoname = "Volume")
  all.predictions <- c(all.predictions,pred[,2]); all.truth <- c(all.truth,pred[,3])
  group.cor <- c(group.cor,cor(pred[,2],pred[,3],method="spearman"))
  group.rmse <- c(group.rmse,rmse(pred[,2],pred[,3]))
}

plot(all.predictions,all.truth)
gp3.rmse <- cbind("model"=rep("transcripts*snps",7),"rmse"=group.rmse)
gp3.cor <- cbind("model"=rep("transcripts*snps",7),"cor"=group.cor,"rmse"=group.rmse) 
all.predictions3 <- all.predictions

## 7-FOLD Image: ####
all.dat.7fold <- data.frame(rbind(gp1.cor,gp2.cor,gp3.cor))
all.dat.7fold$group <- as.factor(as.character(rep(x = 1:7,times=3)))
all.dat.7fold$rmse <- as.numeric(as.character(all.dat.7fold$rmse))
all.dat.7fold$cor <- as.numeric(as.character(all.dat.7fold$cor))

ggboxplot(all.dat.7fold, x = "group", y = "cor",
               color = "model", palette = "jco",
               add = "jitter",short.panel.labs = F,shape="model",size=3)
ggboxplot(all.dat.7fold, x = "group", y = "cor",
          color = "group", palette = "jco",short.panel.labs = F,show.legend=F) + geom_point(aes(shape=model),size=2) 

p <- ggplot(data=all.dat.7fold,aes(x=group, y=cor,colour=model)) + geom_point(size=3) + 
  geom_path(data = all.dat.7fold, aes(x=as.numeric(as.character(group)),y=cor)) 

p + ggplot(data=all.dat.7fold,aes(x=group, y=rmse,colour=model)) + geom_point(size=3) + 
  geom_path(data = all.dat.7fold, aes(x=as.numeric(as.character(group)),y=rmse))


qplot(group, cor, data = all.dat.7fold,color=model,shape=model,size=1)

qplot(model, rmse, data = all.dat.7fold, 
      geom= "violin", fill = model)

final.dat <- data.frame(cbind("true.val"=rep(all.truth,times=3), 
"pred.val"=c(all.predictions1,all.predictions2,all.predictions3),
"model"=c(rep("transcripts",56),rep("snps",56),rep("transcripts*snps",56))))
final.dat$true.val <- as.numeric(as.character(final.dat$true.val))
final.dat$pred.val <- as.numeric(as.character(final.dat$pred.val))

plot(final.dat$true.val,final.dat$pred.val)
(plot5 <- ggplot(aes(x = true.val, y = pred.val, group = model, colour = model),
                 data = final.dat) +
    geom_line() + geom_point())
t.test(as.numeric(gp1.cor[,2]),as.numeric(gp2.cor[,2]),paired = T)
plot(all.dat.7fold$cor,all.dat.7fold$rmse)

# IMAGE-LGEP: All txpts lgep vs. ew with OmicKriging ####
test.group <- ew.fams
train.group <- lgep.fams
pred <-  okriging(idtest = test.group,idtrain = train.group,corlist = list(txpt.cor),H2vec = c(1),pheno = fam.phenos,phenoname = "Volume")
plot(pred[,2],pred[,3])
pred1 <- pred[,2]
lgep.rmse1 <- rmse(pred[,2],pred[,3])
lgep.cor1 <- cor(pred[,2],pred[,3],method = "spearman")
# IMAGE-LGEP: All snps lgep vs. ew with OmicKriging ####
test.group <- ew.fams
train.group <- lgep.fams
pred <-  okriging(idtest = test.group,idtrain = train.group,corlist = list(snp.cor),H2vec = c(.3),pheno = fam.phenos,phenoname = "Volume")
pred2 <- pred[,2]
lgep.rmse2 <- rmse(pred[,2],pred[,3])
lgep.cor2 <- cor(pred[,2],pred[,3],method = "spearman")
# IMAGE-LGEP: All snps +txpts vs. ew with OmicKriging ####
pred <-  okriging(idtest = test.group,idtrain = train.group,corlist = list(snp.cor,txpt.cor),H2vec = c(.97,.03),pheno = fam.phenos,phenoname = "Volume")
pred3 <- pred[,2]
lgep.rmse3 <- rmse(pred[,2],pred[,3])
lgep.cor3 <- cor(pred[,2],pred[,3],method = "spearman")

##LGEP IMAGE #####
final.dat <- data.frame(cbind("true.val"=rep(pred[,3],times=3), 
                              "pred.val"=c(pred1,pred2,pred3),
                              "model"=c(rep("transcripts",11),rep("snps",11),rep("transcripts*snps",11))))
final.dat$true.val <- as.numeric(as.character(final.dat$true.val))
final.dat$pred.val <- as.numeric(as.character(final.dat$pred.val))

plot(final.dat$true.val,final.dat$pred.val)
ggplot(aes(x = true.val, y = pred.val, group = model, colour = model),
                 data = final.dat) +  geom_point() + geom_line()

# across h2vecs ####
weights <- seq(.01,1,.01); rmse.list <- c(); cor.list <- c() ; spearman.list <- c()
for(this.w in 1:length(weights)){
  weight <- weights[this.w]
  all.predictions <- c();all.truth <- c(); group.cor <- c(); group.rmse <- c()
  for(each.group in 1:7){
    test.group <- rownames(fam.phenos)[pheno.groups[,each.group]]
    train.group <- rownames(fam.phenos)[-pheno.groups[,each.group]]
    
    pred <-  okriging(idtest = test.group,idtrain = train.group,corlist = list(txpt.cor),H2vec = c(weight),pheno = fam.phenos,phenoname = "Volume")
    all.predictions <- c(all.predictions,pred[,2]); all.truth <- c(all.truth,pred[,3])
    group.cor <- c(group.cor,cor(pred[,2],pred[,3],method = "spearman"))
    group.rmse <- c(group.rmse,rmse(pred[,2],pred[,3]))
  }
  
  cor.list <- c(cor.list, cor(all.predictions,all.truth))
  rmse.list <- c(rmse.list, rmse(all.predictions,all.truth))
  spearman.list <- c(spearman.list, cor(all.predictions,all.truth,method="spearman"))
}
plot(weights,cor.list)
plot(weights,spearman.list)
plot(weights,rmse.list)

# across h2vecs ####
weights <- seq(.01,1,.01); rmse.list <- c(); cor.list <- c() ; spearman.list <- c()
for(this.w in 1:length(weights)){
  weight <- weights[this.w]
  all.predictions <- c();all.truth <- c(); group.cor <- c(); group.rmse <- c()
  for(each.group in 1:7){
    test.group <- rownames(fam.phenos)[pheno.groups[,each.group]]
    train.group <- rownames(fam.phenos)[-pheno.groups[,each.group]]
    
    pred <-  okriging(idtest = test.group,idtrain = train.group,corlist = list(snp.cor),H2vec = c(weight),pheno = fam.phenos,phenoname = "Volume")
    all.predictions <- c(all.predictions,pred[,2]); all.truth <- c(all.truth,pred[,3])
    group.cor <- c(group.cor,cor(pred[,2],pred[,3],method = "spearman"))
    group.rmse <- c(group.rmse,rmse(pred[,2],pred[,3]))
  }
  
  cor.list <- c(cor.list, cor(all.predictions,all.truth))
  rmse.list <- c(rmse.list, rmse(all.predictions,all.truth))
  spearman.list <- c(spearman.list, cor(all.predictions,all.truth,method="spearman"))
}
plot(weights,cor.list)
plot(weights,spearman.list)
plot(weights,rmse.list)


# Vary weights across spectrum All snps + All txpts with OmicKriging  ####
snp.weights <- seq(.01,.99,.01); rmse.list <- c(); cor.list <- c() ; spearman.list <- c()
for(each.weight in 1:length(snp.weights)){
  s.weight <- snp.weights[each.weight]
  t.weight <- 1- s.weight
  
  all.predictions <- c();all.truth <- c(); 
  for(each.group in 1:7){
    test.group <- rownames(fam.phenos)[pheno.groups[,each.group]]
    train.group <- rownames(fam.phenos)[-pheno.groups[,each.group]]
    
    pred <-  okriging(idtest = test.group,idtrain = train.group,corlist = list(snp.cor,txpt.cor),H2vec = c(s.weight,t.weight),pheno = fam.phenos,phenoname = "Volume")
    all.predictions <- c(all.predictions,pred[,2]); all.truth <- c(all.truth,pred[,3])
  }
  cor.list <- c(cor.list, cor(all.predictions,all.truth))
  rmse.list <- c(rmse.list, rmse(all.predictions,all.truth))
  spearman.list <- c(spearman.list, cor(all.predictions,all.truth,method="spearman"))
}
plot(snp.weights,cor.list)
which.max(cor.list)
plot(snp.weights,spearman.list)
plot(snp.weights,rmse.list)

# across h2vecs ####
weights <- seq(.01,1,.01); rmse.list <- c(); cor.list <- c() ; spearman.list <- c()
for(this.w in 1:length(weights)){
  weight <- weights[this.w]
  all.predictions <- c();all.truth <- c(); group.cor <- c(); group.rmse <- c()
  test.group <- ew.fams
  train.group <- lgep.fams
  
  pred <-  okriging(idtest = test.group,idtrain = train.group,corlist = list(txpt.cor),H2vec = c(weight),pheno = fam.phenos,phenoname = "Volume")
  
  cor.list <- c(cor.list, cor(pred[,2],pred[,3]))
  rmse.list <- c(rmse.list, rmse(pred[,2],pred[,3]))
  spearman.list <- c(spearman.list, cor(pred[,2],pred[,3],method="spearman"))
}
plot(weights,cor.list)
plot(weights,spearman.list)
plot(weights,rmse.list)

# across h2vecs####
weights <- seq(.01,1,.01); rmse.list <- c(); cor.list <- c() ; spearman.list <- c()
for(this.w in 1:length(weights)){
  weight <- weights[this.w]
  all.predictions <- c();all.truth <- c(); group.cor <- c(); group.rmse <- c()
  test.group <- ew.fams
  train.group <- lgep.fams
  
  pred <-  okriging(idtest = test.group,idtrain = train.group,corlist = list(snp.cor),H2vec = c(weight),pheno = fam.phenos,phenoname = "Volume")
  
  cor.list <- c(cor.list, cor(pred[,2],pred[,3]))
  rmse.list <- c(rmse.list, rmse(pred[,2],pred[,3]))
  spearman.list <- c(spearman.list, cor(pred[,2],pred[,3],method="spearman"))
}
plot(weights,cor.list)
plot(weights,spearman.list)
plot(weights,rmse.list)





# Vary weights across spectrum All snps + All txpts LGEP OK #####
test.group <- ew.fams
train.group <- lgep.fams
# Vary weights across spectrum
snp.weights <- seq(.01,.99,.01); rmse.list <- c(); cor.list <- c() ; spearman.list <- c()
for(each.weight in 1:length(snp.weights)){
  s.weight <- snp.weights[each.weight]
  t.weight <- 1- s.weight
  pred <-  okriging(idtest = test.group,idtrain = train.group,corlist = list(snp.cor,txpt.cor),H2vec = c(s.weight,t.weight),pheno = fam.phenos,phenoname = "Volume")
  cor.list <- c(cor.list, cor(pred[,2],pred[,3]))
  rmse.list <- c(rmse.list, rmse(pred[,2],pred[,3]))
  spearman.list <- c(spearman.list, cor(pred[,2],pred[,3],method="spearman"))
}
#plot(snp.weights,cor.list)
#plot(snp.weights,spearman.list)
#plot(snp.weights,rmse.list)
# IMAGE-10K: 10K 8 random sample true data #####
# TRANSCRIPTS
all.predictions <- c();all.truth <- c(); group.cor <- c(); group.rmse <- c();group.spearman <- c()
samp.seed <- sample(14234:1053453,size = 10000)
for(each.group in 1:10000){
  set.seed(seed = samp.seed[each.group])
  this.samp <- sample(1:56,8)
  test.group <- rownames(fam.phenos)[this.samp]
  train.group <- rownames(fam.phenos)[-this.samp]
  
  pred <-  okriging(idtest = test.group,idtrain = train.group,corlist = list(txpt.cor),H2vec = c(1),pheno = fam.phenos,phenoname = "Volume")
  all.predictions <- c(all.predictions,pred[,2]); all.truth <- c(all.truth,pred[,3])
  group.spearman <- c(group.spearman,cor(pred[,2],pred[,3],method = "spearman"))
  group.cor <- c(group.cor,cor(pred[,2],pred[,3]))
  group.rmse <- c(group.rmse,rmse(pred[,2],pred[,3]))
}
#plot(all.predictions,all.truth)
hist(group.cor); 
hist(group.spearman)
plot(group.rmse)
cor(all.predictions,all.truth)^2

# SNPS
all.predictions <- c();all.truth <- c(); group.cor <- c(); group.rmse <- c();group.spearman <- c()
samp.seed <- sample(14234:1053453,size = 10000)
for(each.group in 1:10000){
  set.seed(seed = samp.seed[each.group])
  this.samp <- sample(1:56,8)
  test.group <- rownames(fam.phenos)[this.samp]
  train.group <- rownames(fam.phenos)[-this.samp]
  
  pred <-  okriging(idtest = test.group,idtrain = train.group,corlist = list(snp.cor),H2vec = c(1),pheno = fam.phenos,phenoname = "Volume")
  all.predictions <- c(all.predictions,pred[,2]); all.truth <- c(all.truth,pred[,3])
  group.spearman <- c(group.spearman,cor(pred[,2],pred[,3],method = "spearman"))
  group.cor <- c(group.cor,cor(pred[,2],pred[,3]))
  group.rmse <- c(group.rmse,rmse(pred[,2],pred[,3]))
}
#plot(all.predictions,all.truth)
hist(group.cor); 
hist(group.spearman)
#plot(group.rmse)
cor(all.predictions,all.truth)^2

# IMAGE-10K: 10K 8 randomized sample fake data #####
# TRANSCRIPTS
all.predictions <- c();all.truth <- c(); group.cor <- c(); group.rmse <- c();group.spearman <- c()
samp.seed <- sample(14234:1053453,size = 10000)
for(each.group in 1:10000){
  set.seed(seed = samp.seed[each.group])
  this.samp <- sample(1:56,8)
  test.group <- rownames(fam.phenos)[this.samp]
  train.group <- rownames(fam.phenos)[-this.samp]
  
  pred <-  okriging(idtest = test.group,idtrain = train.group,corlist = list(txpt.random.cor),H2vec = c(1),pheno = fam.phenos,phenoname = "Volume")
  all.predictions <- c(all.predictions,pred[,2]); all.truth <- c(all.truth,pred[,3])
  group.spearman <- c(group.spearman,cor(pred[,2],pred[,3],method = "spearman"))
  group.cor <- c(group.cor,cor(pred[,2],pred[,3]))
  group.rmse <- c(group.rmse,rmse(pred[,2],pred[,3]))
}
#plot(all.predictions,all.truth)
hist(group.cor); 
hist(group.spearman)
hist(group.rmse)
cor(all.predictions,all.truth)^2

#SNPS
all.predictions <- c();all.truth <- c(); group.cor <- c(); group.rmse <- c();group.spearman <- c()
samp.seed <- sample(14234:1053453,size = 10000)
for(each.group in 1:10000){
  set.seed(seed = samp.seed[each.group])
  this.samp <- sample(1:56,8)
  test.group <- rownames(fam.phenos)[this.samp]
  train.group <- rownames(fam.phenos)[-this.samp]
  
  pred <-  okriging(idtest = test.group,idtrain = train.group,corlist = list(snp.random.cor),H2vec = c(1),pheno = fam.phenos,phenoname = "Volume")
  all.predictions <- c(all.predictions,pred[,2]); all.truth <- c(all.truth,pred[,3])
  group.spearman <- c(group.spearman,cor(pred[,2],pred[,3],method = "spearman"))
  group.cor <- c(group.cor,cor(pred[,2],pred[,3]))
  group.rmse <- c(group.rmse,rmse(pred[,2],pred[,3]))
}
#plot(all.predictions,all.truth)
hist(group.cor); 
hist(group.spearman)
hist(group.rmse)
cor(all.predictions,all.truth)^2


## Prediction of Breeding Values 

### Predictions of Breeding Values across batch
#### Pre filtering predictions
      * Txpts(txpt.matrix), SNPs (rel.matrix), Both
      * 10k fold 7-fold random draw data structure
#### Filtering Predictions


### Prediction of Breeding Values among batch
#### Pre filtering predictions
       * Txpts(txpt.matrix), SNPs (rel.matrix), Both
#### Filtering Predictions


## All files are saved as outputs and can be loaded below ####
load("/mnt/bio.012.df.RData")
load("/mnt/fam.012.df.RData")
load("/mnt/asreml.counts.RData")
load("/mnt/fam.asreml.counts.RData")
load("/mnt/bio.phenos.RData")
load("/mnt/fam.phenos.RData")
load("/mnt/fam.hm.012.df.RData")
load("/mnt/bio.hm.012.df.RData")
u.fams <- unique(fam.phenos$fam_id)

# With all starting files loaded we can now generate the first plot ####

# Create even test groups ####
# Index the phenotypes so they are ind decreasing order
pheno.order <- order(fam.phenos$Volume,decreasing = T)
first <- seq(1,56,7)
last <- seq(7,56,7)

se <- c(3223,42435,23534,5235,3423,1345,21541,35432)
pheno.groups <- lapply(1:8,function(each.group){
  set.seed(se[each.group])
  sample(pheno.order[first[each.group]:last[each.group]])
})
pheno.groups <- do.call(rbind,pheno.groups)
# Each of the 7 groups contain 8 individuals spread across the population of phenotypes


# For each of the test groups, run an anova on the train.test ####

# Generate anova scores for groups 
group.snps.anova.scores <- mclapply(1:7,function(each.group){
  test.set <- which(bio.phenos$fam_id %in% u.fams[c(pheno.groups[,each.group])])
  train.set <- bio.hm.012.df[-test.set,]
  train.y <- bio.phenos$Volume[-test.set]
  apply(train.set,2,function(each.snp) anovaScores(each.snp,train.y))
},mc.cores=10)

group.txpt.anova.scores <- mclapply(1:7,function(each.group){
  test.set <- which(bio.phenos$fam_id %in% u.fams[c(pheno.groups[,each.group])])
  train.set <- asreml.counts[-test.set,]
  train.y <- bio.phenos$Volume[-test.set]
  apply(train.set,2,function(each.snp) anovaScores(each.snp,train.y))
},mc.cores=10)

save(group.snps.anova.scores,file="/mnt/group.snps.anova.scores.RData",compress=T)
save(group.txpt.anova.scores,file="/mnt/group.txpt.anova.scores.RData",compress=T)

# IMAGE-7FOLD: All txpts 7-Fold with OmicKriging ####
group.snps.anova.adj <- lapply(group.snps.anova.scores,function(x) p.adjust(x,method = "fdr"))
group.txpt.anova.adj <- lapply(group.txpt.anova.scores,function(x) p.adjust(x,method="fdr"))


all.predictions <- c();all.truth <- c(); group.cor <- c(); group.rmse <- c()
for(each.group in 1:7){
  txpt.cor <- cor(t(scale(center=F,log2(fam.asreml.counts[,which(group.txpt.anova.adj[[each.group]] < .05)]))))
  test.group <- rownames(fam.phenos)[pheno.groups[,each.group]]
  train.group <- rownames(fam.phenos)[-pheno.groups[,each.group]]
  
  pred <-  okriging(idtest = test.group,idtrain = train.group,corlist = list(txpt.cor),H2vec = c(.8),pheno = fam.phenos,phenoname = "Volume")
  all.predictions <- c(all.predictions,pred[,2]); all.truth <- c(all.truth,pred[,3])
  group.cor <- c(group.cor,cor(pred[,2],pred[,3],method="spearman"))
  group.rmse <- c(group.rmse,rmse(pred[,2],pred[,3]))
}
plot(all.predictions,all.truth)
plot(group.cor); plot(group.rmse)
cor(all.predictions,all.truth)^2


gp1.rmse <- cbind("model"=rep("transcripts",7),"rmse"=group.rmse)
gp1.cor <- cbind("model"=rep("transcripts",7),"cor"=group.cor,"rmse"=group.rmse) 
all.predictions1 <- all.predictions
# IMAGE-7FOLD: All snps 7-fold with OmicKriging ####
all.predictions <- c();all.truth <- c(); group.cor <- c(); group.rmse <- c()
for(each.group in 1:7){
  txpt.cor <- cor(t(fam.hm.012.df[,which(group.snps.anova.adj[[each.group]] < .3)]))
  test.group <- rownames(fam.phenos)[pheno.groups[,each.group]]
  train.group <- rownames(fam.phenos)[-pheno.groups[,each.group]]
  
  pred <-  okriging(idtest = test.group,idtrain = train.group,corlist = list(txpt.cor),H2vec = c(1),pheno = fam.phenos,phenoname = "Volume")
  all.predictions <- c(all.predictions,pred[,2]); all.truth <- c(all.truth,pred[,3])
  group.cor <- c(group.cor,cor(pred[,2],pred[,3],method="spearman"))
  group.rmse <- c(group.rmse,rmse(pred[,2],pred[,3]))
}
plot(all.predictions,all.truth)
plot(group.cor); plot(group.rmse)
cor(all.predictions,all.truth)^2
#gp1.rmse <- cbind("model"=rep("transcripts",7),"rmse"=group.rmse)
#gp1.cor <- cbind("model"=rep("transcripts",7),"cor"=group.cor,"rmse"=group.rmse) 
#all.predictions1 <- all.predictions
#IMAGE-7fOLD: All snps + All txpts with OmicKriging ####
all.predictions <- c();all.truth <- c(); group.cor <- c(); group.rmse <- c()
for(each.group in 1:7){
  txpt.cor <- cor(t(scale(center=F,log2(fam.asreml.counts[,which(group.txpt.anova.adj[[each.group]] < .8)]))))
  snp.cor <- cor(t(fam.hm.012.df[,which(group.snps.anova.adj[[each.group]] < .8)]))
  test.group <- rownames(fam.phenos)[pheno.groups[,each.group]]
  train.group <- rownames(fam.phenos)[-pheno.groups[,each.group]]

  pred <-  okriging(idtest = test.group,idtrain = train.group,corlist = list(snp.cor,txpt.cor),H2vec = c(.98,.02),pheno = fam.phenos,phenoname = "Volume")
  all.predictions <- c(all.predictions,pred[,2]); all.truth <- c(all.truth,pred[,3])
  group.cor <- c(group.cor,cor(pred[,2],pred[,3],method="spearman"))
  group.rmse <- c(group.rmse,rmse(pred[,2],pred[,3]))
}

plot(all.predictions,all.truth)
cor(all.predictions,all.truth)^2

gp3.rmse <- cbind("model"=rep("transcripts*snps",7),"rmse"=group.rmse)
gp3.cor <- cbind("model"=rep("transcripts*snps",7),"cor"=group.cor,"rmse"=group.rmse) 
all.predictions3 <- all.predictions



tC <- trainControl(method="LGOCV",allowParallel = T,p=.75,savePredictions="all",verboseIter = F,search="grid")

algorithmList <- c("glmnet","BstLm","ranger","svmLinear2","simpls")

the.predictions <- lapply(1:7,function(each.group){
#txpt.cor <- scale(center=F,log2(fam.asreml.counts[,which(group.txpt.anova.adj[[each.group]] < .05)]))
snp.cor <- fam.hm.012.df[,which(group.snps.anova.adj[[each.group]] < .1)]
test.group <- c(pheno.groups[,each.group])
train.group <- c(1:56)[-pheno.groups[,each.group]]
t3 <- train(y=fam.phenos$Volume[train.group], x=snp.cor[train.group,],algorithmList[each.group],metric = "RMSE", trControl = tC)
cor(predict(object = t3,newdata = snp.cor[test.group,]),fam.phenos$Volume[test.group])

     


t2 <- mclapply(1:length(algorithmList),function(i){ 
  t3 <- train(y=LGEP.train.p, x=LGEP.train.dat,algorithmList[i], trControl = tC)
},mc.cores=20)


