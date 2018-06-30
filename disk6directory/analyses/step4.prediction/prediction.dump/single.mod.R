the.p <- seq(.01,.15,.01)

for(each in 1:length(the.p)){
LGEP.train.dat <- fam.012[which(fam.phenos$fam_id %in% lgep.fams),which(p.adjust(anova_scores) < the.p[each])]
colnames(LGEP.train.dat) <- colnames(fam.012[,which(p.adjust(anova_scores) < the.p[each])])
LGEP.train.p <- fam.phenos$Volume[which(fam.phenos$fam_id %in% lgep.fams)]
t2 <- train(y=LGEP.train.p, x=LGEP.train.dat,method = "glmnet", trControl = trainControl(method = "LOOCV"),standardize=F)
t2
the.ew.predictions <- predict(t2,fam.012[which(fam.phenos$fam_id %in% ew.fams),which(p.adjust(anova_scores) < the.p[each])])
plot(the.ew.predictions,fam.phenos$Volume[which(fam.phenos$fam_id %in% ew.fams)])
abline(h=mean(fam.phenos$Volume[which(fam.phenos$fam_id %in% ew.fams)]),v=mean(the.ew.predictions))
}

bio.pheno.subset$v



fam.phenos2 <- unique(bio.pheno.subset[,c("fam_id","Volume","binary_vol")])
rownames(fam.phenos2) <- fam.phenos2$fam_id
u.fams <- unique(fam.phenos2$fam_id)

ctrl <- rfeControl(functions = rfFuncs,
                   method = "repeatedcv",
                   repeats = 1,
                   verbose = FALSE)

the.pred <- c()

the.pred <- mclapply(1:56,function(each.fam){
td <- fam.012[-each.fam,which(p.adjust(family_anova_scores[[each.fam]]) < .5)]
tp <- as.factor(fam.phenos2$binary_vol[-each.fam])
test.d <- fam.012[,which(p.adjust(family_anova_scores[[each.fam]]) < .5)]
t2 <- train(y=tp, x=td,method = "spls", trControl = trainControl(method = "boot"))
predict(t2,test.d)[each.fam]
},mc.cores=20)
confusionMatrix(reference = as.factor(fam.phenos2$binary_vol), data = unlist(the.pred), mode='everything')

gl.rfe <- rfe(x = td,y = tp,rfeControl = ctrl,method="glmnet")
pls.rfe <- rfe(x = td,y = tp,rfeControl = ctrl,method="pls")
rrg.rfe <- rfe(x = td,y = tp,rfeControl = ctrl,method="RRFglobal")
pen.rfe <- rfe(x = td,y = tp,rfeControl = ctrl,method="penalized")



t2 <- train(y=LGEP.train.p, x=LGEP.train.dat,method = "glmnet", trControl = trainControl(method = "LOOCV"),standardize=F)
t2