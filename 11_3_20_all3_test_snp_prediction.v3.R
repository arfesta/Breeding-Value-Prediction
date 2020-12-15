{
setwd("/media/blackbox/pine.families/")
  ALL3 <- read.table("./out.012")[,-1]
  pos <- read.table("./out.012.pos")
  indv <- read.table("out.012.indv")
#q.all <- read.table("./bam.wrg/hq.05_192_maf05.95_nomiss_dp10_q30.012.012")[,-1]
#pos <- read.table("./bam.wrg/hq.05_192_maf05.95_nomiss_dp10_q30.012.012.pos")
#keep <- which(pos[,1] %in% as.character(cc$V1))
pos <- paste0(pos[,1],"_",pos[,2])
colnames(ALL3) <- pos
#indv <- read.table("./bam.wrg/hq.05_192_maf05.95_nomiss_dp10_q30.012.012.indv")
rownames(ALL3) <- indv[,1]

#q.all <- q.all[,keep]

#old.snps.pos <- read.table("/media/disk6/ARF/RNASEQ/shared/snps/86k/192.biorep.Q30.snps.no.miss.maf05.ann.high.mod.012.pos")
#old.snps.pos <- paste0(old.snps.pos[,1],"_",old.snps.pos[,2])
rm(pos,indv)
q.all <- ALL3
r.fams <- gsub(".....fq.gz","",x=rownames(q.all)[-c(1:192)])
u.fams <- unique(r.fams)
batch3 <- q.all[-c(1:192),]
new.dat <- do.call(rbind,parallel::mclapply(1:length(u.fams),function(x){
  apply(batch3[which(r.fams %in% u.fams[x]),],2,mean)
},mc.cores=22))
dim(new.dat)
colnames(new.dat) <- colnames(q.all)
rownames(new.dat) <- u.fams
rm(r.fams,u.fams)


r.names <- rownames(q.all)[c(1:192)]
ew.names <- gsub(pattern = ".ew.bio.rep.merge.bam",replacement = "",x = r.names[1:48])
ew.names <- gsub(pattern="_",replacement = "",ew.names)
ew.names <- paste0(as.numeric(ew.names) + 1000,".")
lgep.names <- gsub(pattern = "..lgep.bio.rep.merge.bam",replacement = "",x = r.names[-c(1:48)])
lgep.names <- paste0(lgep.names,".")
old.snps.names <- c(ew.names,lgep.names)
rownames(q.all)[c(1:192)] <- old.snps.names

rm(ew.names,lgep.names,old.snps.names,r.names)


load("/media/disk6/ARF/RNASEQ/BV-Prediction/Breeding-Value-Prediction/disk6directory/resources/expt.dat.720.RData")
colnames(expt.dat.720)
bio.phenos <- unique(expt.dat.720[,c("animal_id","fam_id","Volume")])
keep.cols <- which(paste0(bio.phenos$animal_id,".") %in% rownames(q.all))

q.all <- q.all[c(1:192),]
bio.phenos <- bio.phenos[keep.cols,]
rownames(q.all);bio.phenos$animal_id
identical(paste0(as.character(bio.phenos$animal_id),".")[match(rownames(q.all),paste0(as.character(bio.phenos$animal_id),"."))],
          rownames(q.all))

bio.phenos <- bio.phenos[match(rownames(q.all),paste0(as.character(bio.phenos$animal_id),".")),]
rownames(bio.phenos) <- paste0(as.character(bio.phenos$animal_id),".")
q.all <- q.all[-c(which(is.na(bio.phenos$Volume))),]
bio.phenos <- bio.phenos[-c(which(is.na(bio.phenos$Volume))),]
rm(keep.cols,expt.dat.720)

u.fams <- unique(as.character(bio.phenos$fam_id))
fam.snp.mat <- do.call(rbind,parallel::mclapply(1:length(u.fams),function(x) {
  these.rows <- which(as.character(bio.phenos$fam_id) %in% u.fams[x])
  if(length(these.rows) > 1) {
    apply(q.all[these.rows,],2,mean) } else {
      q.all[these.rows,]
    }
},mc.cores=56))
rownames(fam.snp.mat) <- u.fams
fam.snp.mat <- as.matrix(fam.snp.mat)

fam.phenos <- unique(bio.phenos[,c(2,3)])
rm(bio.phenos,q.all,u.fams)
}

## Identify MAF of snps in beach 1,2 and then 3
a0=apply(batch3,2,function(x) length(which(x == 0))) 
a1=apply(batch3,2,function(x) length(which(x == 1))) 
a2=apply(batch3,2,function(x) length(which(x == 2))) 
batch3_maf <- rbind(a0,a1,a2)
a22= a2*2
minor_count = sapply(1:19896,function(x)sum(a1[x],a2[x]))
maf_batch3 <- minor_count/(66*2)

a0=apply(q.all,2,function(x) length(which(x == 0))) 
a1=apply(q.all,2,function(x) length(which(x == 1))) 
a2=apply(q.all,2,function(x) length(which(x == 2))) 

minor_count = sapply(1:19896,function(x)sum(a1[x],a2[x]))
maf_batch12 <- minor_count/(192*2)

plot(maf_batch3,maf_batch12)
maf_residual = maf_batch3-maf_batch12
plot(maf_residual)
drop.these <- which(abs(maf_residual) <.1) 
library(ggfortify)
pca_res_aLl <- prcomp(test.dat2[,-drop.these])
autoplot(pca_res_aLl)

### PRedict ####
library(glmnet)
identical(colnames(fam.snp.mat),colnames(new.dat))
test.dat <- as.matrix(rbind(fam.snp.mat,new.dat))
rm(new.dat)

library(readr)
Batch3_phenos <- read_csv("~/Batch3_phenos.csv")
test.dat2 <- test.dat[c(1:56,match(as.character(Batch3_phenos$RNA_ID),rownames(test.dat))),]

colnames(Batch3_phenos)[c(1,4)] <- c("fam_id","Volume")

all_phenos <- rbind(fam.phenos,Batch3_phenos[,c(1,4)])

a=0; std=F
loo_1_f <- unlist(parallel::mclapply(1:78,function(x){
   train.dat <- test.dat2[-c(x),-drop.these]
  train.y <- all_phenos$Volume[-c(x)]
  gl.mod <- glmnet::glmnet(x = train.dat,y = train.y,alpha = a,standardize=std,nlambda = 100)
  cc <- which.min(apply(predict(gl.mod,train.dat)[,-1],2,function(x) caret::RMSE(x,obs = train.y)))
predict(gl.mod,test.dat2[c(x,57:78),-drop.these],s = gl.mod$lambda[cc])[1,1]
},mc.cores=56))
plot(loo_1_f[57:78],all_phenos$Volume[57:78])

loo_1_f <- (parallel::mclapply(1:56,function(x){
   train.dat <- test.dat2[-c(x,57:78),-drop.these]
  train.y <- all_phenos$Volume[-c(x,57:78)]
  gl.mod <- glmnet::glmnet(x = train.dat,y = train.y,alpha = a,standardize=std,nlambda = 100)
  cc <- which.min(apply(predict(gl.mod,train.dat)[,-1],2,function(x) caret::RMSE(x,obs = train.y)))
predict(gl.mod,test.dat2[c(x,57:78),-drop.these],s = gl.mod$lambda[cc])
},mc.cores=56))
mean_ans <- apply(do.call(cbind,loo_1_f),1,mean)[-1]
plot(mean_ans,all_phenos$Volume[57:78])

new.dat <- test.dat2[c(57:78),]
a=1; std=F
loo_22_1_f <- unlist(parallel::mclapply(1:22,function(x){
   train.dat <- new.dat[-c(x),]
  train.y <- Batch3_phenos$Volume[-x]
  gl.mod <- glmnet::glmnet(x = train.dat,y = train.y,alpha = a,standardize=std,nlambda = 100)
  cc <- which.min(apply(predict(gl.mod,train.dat)[,-1],2,function(x) caret::RMSE(x,obs = train.y)))
predict(gl.mod,new.dat,s = gl.mod$lambda[cc])[x,1]
},mc.cores=22))
plot(loo_22_1_f,Batch3_phenos$Volume)
cor(loo_22_1_f,Batch3_phenos$Volume)

a=1; std=T
loo_56_1_f <- unlist(parallel::mclapply(1:56,function(x){
   train.dat <- fam.snp.mat[-c(x),]
  train.y <- fam.phenos$Volume[-x]
  gl.mod <- glmnet::glmnet(x = train.dat,y = train.y,alpha = a,standardize=std,nlambda = 100)
  cc <- which.min(apply(predict(gl.mod,train.dat)[,-1],2,function(x) caret::RMSE(x,obs = train.y)))
predict(gl.mod,fam.snp.mat,s = gl.mod$lambda[cc])[x,1]
},mc.cores=56))
cor(loo_56_1_f,fam.phenos$Volume)
plot(loo_

alphas <- c(0,.5,1)
non_stand_preds <- list("vector")
stand_preds <- list("vector")

for(each in 1:length(alphas)){
  a= alphas[each]
  std=F
non_stand_preds[[each]] <- do.call(cbind,parallel::mclapply(1:56,function(x){
  train.dat <- test.dat[-c(x,57:78),]
  train.y <- fam.phenos$Volume[-x]
  gl.mod <- glmnet(x = train.dat,y = train.y,alpha = a,standardize=std,nlambda = 100)
  cc <- which.min(apply(predict(gl.mod,train.dat)[,-1],2,function(x) RMSE(x,obs = train.y)))
predict(gl.mod,test.dat[c(x,57:78),],s = gl.mod$lambda[cc])
},mc.cores=56))

std=T
 stand_preds[[each]] <- do.call(cbind,parallel::mclapply(1:56,function(x){
  train.dat <- test.dat[-c(x,57:78),]
  train.y <- fam.phenos$Volume[-x]
  gl.mod <- glmnet(x = train.dat,y = train.y,alpha = a,standardize=std,nlambda = 100)
  cc <- which.min(apply(predict(gl.mod,train.dat)[,-1],2,function(x) caret::RMSE(x,obs = train.y)))
predict(gl.mod,test.dat[c(x,57:78),],s = gl.mod$lambda[cc])
},mc.cores=56))
print(each)
}


rm(a,std,each,alphas)

unlist(lapply(stand_preds,function(x) cor(x[1,],fam.phenos$Volume)))
unlist(lapply(non_stand_preds,function(x) cor(x[1,],fam.phenos$Volume)))

stand_preds_new <- do.call(cbind,lapply(stand_preds,function(x) apply(x[-1,],1,mean)))

non_stand_preds_new <- do.call(cbind,lapply(non_stand_preds,function(x) apply(x[-1,],1,mean)))


stand_preds_new <- stand_preds_new[match(check[,1],rownames(stand_preds_new)),]



non_stand_preds_new <- non_stand_preds_new[match(check[,1],rownames(non_stand_preds_new)),]



predictions <- do.call(cbind,predictions)
cor(predictions[1,],fam.phenos$Volume)
fam= rownames(test.dat)[57:78]
bv = apply(predictions[-1,],1,mean)
preds <- data.frame(fam,bv)


cor(predictions[1,],fam.phenos$Volume)

select.txpts <- as.numeric(as.character(names(table(unlist(predictions)))))

 train.dat <- fam.snp.mat[,select.txpts]
  train.y <- fam.phenos$Volume
gl.mod <- glmnet(x = train.dat,y = train.y,alpha = 1,standardize=F)
gl.mod <- coef(gl.mod)
gl.mod <- apply(gl.mod,1,function(x) length(unique(x)))[-1]
  these.txpts <-  which(gl.mod > 1)
gl.mod <- cv.glmnet(train.dat[,these.txpts],fam.phenos$Volume,nfolds = 10,alpha=1,type.measure = "mae",standardize=F,grouped=T)

test.dat <- test.dat[,select.txpts]
the.preds <- predict(gl.mod,test.dat[c(57:78),these.txpts])

cor(the.preds,check$`16`)

predictions <- do.call(cbind,predictions)
fam= names(apply(predictions[-1,],1,mean))
bv = apply(predictions[-1,],1,mean)
preds <- data.frame(fam,bv)

new.dat.pred.mean <- apply(predictions[-1,],1,mean)
new.dat.pred.sd <- apply(predictions[-1,],1,sd)
n <- 56
error <- (qnorm(0.975)*new.dat.pred.sd)/sqrt(n)
l <- new.dat.pred.mean -error
r <- new.dat.pred.mean +error

new.dat.pred.min <- apply(predictions[-1,],1,min)
new.dat.pred.max <- apply(predictions[-1,],1,max)
new.preds <- cbind("mean"=new.dat.pred.mean,"left"=l,"right"=r,"min"=new.dat.pred.min,"max"=new.dat.pred.max)
new.dat.pred.range <- new.dat.pred.max - new.dat.pred.min

#plot(preds2,fam.phenos$Volume)
preds2 <- do.call(cbind,preds2)
cor(preds2[1,],fam.phenos$Volume)

apply(preds2[-1,],1,mean)
fam.snp.mat2 <- fam.snp.mat[,which(colnames(fam.snp.mat) %in% paste0(new.snp.pos[,1],"_",new.snp.pos[,2]))]

test.dat <- as.matrix(fam.snp.mat)
preds <- unlist(mclapply(1:56,function(x){
  train.dat <- as.matrix(fam.snp.mat[-x,])
  train.y <- fam.phenos$Volume[-x]
gl.mod <- glmnet(x = train.dat,y = train.y,alpha = 1,standardize=F)
gl.mod <- coef(gl.mod)
gl.mod <- apply(gl.mod,1,function(x) length(unique(x)))[-1]
  these.txpts <-  which(gl.mod > 1)
gl.mod <- cv.glmnet(train.dat[,these.txpts],train.y,nfolds = nrow(train.dat),alpha=1,type.measure = "mae",standardize=F,grouped=F)

predict(gl.mod,test.dat[,these.txpts])[x,1]
},mc.cores=56))

plot(preds,fam.phenos$Volume);cor(preds,fam.phenos$Volume)
pred <- data.frame("family"=rownames(preds2)[-1],"pred"=pred[,1])
cc <- cor(t(q.all))
plot(cc[lower.tri(cc)])

rownames(new.dat)
new.dat <- new.dat[match(check[,1],rownames(new.dat)),]
preds <- unlist(mclapply(1:22,function(x){
  train.dat <- new.dat[-x,]
  train.y <- check$`16`[-x]
gl.mod <- glmnet(x = train.dat,y = train.y,alpha = 1,standardize=F)
gl.mod <- coef(gl.mod)
gl.mod <- apply(gl.mod,1,function(x) length(unique(x)))[-1]
  these.txpts <-  which(gl.mod > 1)
gl.mod <- cv.glmnet(train.dat[,these.txpts],train.y,nfolds = nrow(train.dat),alpha=1,type.measure = "mae",standardize=F,grouped=F)

predict(gl.mod,new.dat[,these.txpts])[x,1]
},mc.cores=22))


new.fam.cor <- cor(t(new.preds))
