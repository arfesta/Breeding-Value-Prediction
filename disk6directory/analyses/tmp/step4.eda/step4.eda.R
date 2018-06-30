  load("/media/disk6/ARF/RNASEQ/counts/86kSalmon/Step3_LMM_Counts.RData")
  load("/media/disk6/ARF/RNASEQ/Breeding-Value-Prediction/disk6directory/resources/expt.dat.720.RData")

  library(OmicKriging);library(glmnet)
  
  rownames(asreml.counts)
fam.phenos <- unique(expt.dat.720[,c("fam_id","Volume","Volume Accuracy",
                                     "Height","Height Accuracy",
                                     "Strt", "Strt Accuracy")])
fam.phenos <- fam.phenos[which(as.character(fam.phenos$fam_id) %in% rownames(asreml.counts)),]
fam.phenos <- fam.phenos[which(is.na(fam.phenos$Volume) == F),]
asreml.counts <- asreml.counts[which(rownames(asreml.counts) %in% as.character(fam.phenos$fam_id)),]

fam.phenos <- fam.phenos[match(rownames(asreml.counts),as.character(fam.phenos$fam_id)),]
identical(as.character(fam.phenos$fam_id),rownames(asreml.counts))
#TRUE
rownames(fam.phenos) <- as.character(fam.phenos$fam_id)
log.asreml <- scale(log2(asreml.counts+1),center=F)
the.pred <- unlist(mclapply(1:57,function(x){
gl.test <- glmnet(x = log.asreml[-x,],y = fam.phenos$Volume[-x],alpha = 0.4,standardize = F)
p.vals <- predict(gl.test,log.asreml)
p.vals[x,ncol(p.vals)]
},mc.cores=57))
cor(the.pred,fam.phenos$Volume)^2



