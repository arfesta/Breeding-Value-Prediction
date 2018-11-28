# Next look at 70 fold 
library(caret)
load("~/Desktop/Alt Approaches/Breeding-Value-Prediction/shared/ordered.data.RData")
names(ordered.data)
set.seed(24826)
all.folds <- createMultiFolds(y=ordered.data$fam.phenos$Volume,times = 10,k = 7)
# Each Rep corresponds to a single 7-fold CV (Rep 1: folds 1-7 rep1; Rep2: folds 1-7 rep2)

train.index.matrix <- do.call(rbind,all.folds)
test.index.matrix <- do.call(rbind,lapply(1:70,function(x) {
  c(1:56)[-c(train.index.matrix[x,])]
}))

train.fam.matrix <- do.call(rbind,lapply(1:70,function(x) {
  as.character(ordered.data$fam.phenos$fam_id)[c(train.index.matrix[x,])]
}))
test.fam.matrix <- do.call(rbind,lapply(1:70,function(x) {
  as.character(ordered.data$fam.phenos$fam_id)[-c(train.index.matrix[x,])]
}))


log.bio.counts <- scale(log2(ordered.data$bio.counts + 1),center = F)

pt <- proc.time()
txpt.anova.scores <- mclapply(1:70,function(each.group){
  select <- which(as.character(ordered.data$bio.phenos$fam_id) %in% train.fam.matrix[each.group,])
  train.dat <- log.bio.counts[select,]
  train.p <- ordered.data$bio.phenos$Volume[select]
  pval.trains <- apply(train.dat,2,function(z) anovaScores(x = z,y = train.p))
  p.adjust(pval.trains,method = "fdr")
},mc.cores = 8)
proc.time() - pt
#user    system   elapsed 
#12659.345   154.509  2029.395 

pt <- proc.time()
snp.anova.scores <- mclapply(1:70,function(each.group){
  select <- which(as.character(ordered.data$bio.phenos$fam_id) %in% train.fam.matrix[each.group,])
  train.dat <- ordered.data$bio.snps[select,]
  train.p <- ordered.data$bio.phenos$Volume[select]
  pval.trains <- apply(train.dat,2,function(z) anovaScores(x = z,y = train.p))
  p.adjust(pval.trains,method = "fdr")
},mc.cores = 8)
proc.time() - pt
#user   system  elapsed 
#8411.899  100.505 1337.628 

save.image(file = "../Breeding-Value-Prediction/shared/step4_create.70.RData",compress=T)
