#load part 1####
load("~/Desktop/corrected.identifying.top.DE.RDA.pt1.RData")
library(OmicKriging)
##Calculate for fdr is correct subset those less than .01 ####
comparisons.fdr.01 <- vector("list",55)
for(i in 1:55){
  a <- length(corrected.full.comparisons.fdr[[i]])
  for(each in 1:a){
    comparisons.fdr.01[[i]][[each]] <-  rownames(corrected.full.comparisons[[i]][[each]])[which(corrected.full.comparisons.fdr[[i]][[each]] < .01)] 
  }
  print(i)}

total.sig.count.fdr.01 <- vector("list",55)
for(i in 1:55){
  a <- length(corrected.full.comparisons.fdr[[i]])
  for(each in 1:a){
    total.sig.count.fdr.01[[i]][[each]] <-  length(comparisons.fdr.01[[i]][[each]])
  }
  print(i)}

summary(unlist(total.sig.count.fdr.01))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#1.0   134.8   329.5   561.1   751.2  5727.0 
which.min(unlist(total.sig.count.fdr.01))
#comparison number 36

total.sig.mat.fdr.01 <- matrix(NA,55,56)
beg <- 1; end <- 55; match <- 2; all <- 55
for(i in 1:55){
  c <- beg; d <- end
  total.sig.mat.fdr.01[i,match:56] <- unlist(total.sig.count.fdr.01)[c:d]
  print(i); beg <- d + 1; end <- d+(all-i); match <- match+1}

##Calculate for fdr is correct subset those less than .05 ####
comparisons.fdr.05 <- vector("list",55)
for(i in 1:55){
  a <- length(corrected.full.comparisons.fdr[[i]])
  for(each in 1:a){
    comparisons.fdr.05[[i]][[each]] <-  rownames(corrected.full.comparisons[[i]][[each]])[which(corrected.full.comparisons.fdr[[i]][[each]] < .05)] 
  }
  print(i)}

total.sig.count.fdr.05 <- vector("list",55)
for(i in 1:55){
  a <- length(corrected.full.comparisons.fdr[[i]])
  for(each in 1:a){
    total.sig.count.fdr.05[[i]][[each]] <-  length(comparisons.fdr.05[[i]][[each]])
  }
  print(i)}

summary(unlist(total.sig.count.fdr.05))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   7.0   324.5   746.5  1088.0  1522.0  8259.0 

total.sig.mat.fdr.05 <- matrix(NA,55,56)
beg <- 1; end <- 55; match <- 2; all <- 55
for(i in 1:55){
  c <- beg; d <- end
  total.sig.mat.fdr.05[i,match:56] <- unlist(total.sig.count.fdr.05)[c:d]
  print(i); beg <- d + 1; end <- d+(all-i); match <- match+1
}

##Calculate for fdr is correct subset those less than .1 ####
comparisons.fdr.1 <- vector("list",55)
for(i in 1:55){
  a <- length(corrected.full.comparisons.fdr[[i]])
  for(each in 1:a){
    comparisons.fdr.1[[i]][[each]] <-  rownames(corrected.full.comparisons[[i]][[each]])[which(corrected.full.comparisons.fdr[[i]][[each]] < .1)] 
  }
  print(i)}

total.sig.count.fdr.1 <- vector("list",55)
for(i in 1:55){
  a <- length(corrected.full.comparisons.fdr[[i]])
  for(each in 1:a){
    total.sig.count.fdr.1[[i]][[each]] <-  length(comparisons.fdr.1[[i]][[each]])
  }
  print(i)
}

summary(unlist(total.sig.count.fdr.1))
#     Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#     17.0   505.2  1099.0  1512.0  2127.0  9780.0 

total.sig.mat.fdr.1 <- matrix(NA,55,56)
beg <- 1; end <- 55; match <- 2; all <- 55
for(i in 1:55){
  c <- beg; d <- end
  total.sig.mat.fdr.1[i,match:56] <- unlist(total.sig.count.fdr.1)[c:d]
  print(i); beg <- d + 1; end <- d+(all-i); match <- match+1}

#Create comparison number 1-1540 by DE pair being compared####
comp.number.mat <- matrix(NA,55,56)
a <- 1; b <- 55; c <- 2; g <- 53
for(i in 1:55){
  comp.number.mat[i,c:56] <- a:b
  a <- b+1; b <- a+g; c <- c+1; g <- g-1; print(i)
}

#get family mean for counts and phenos####
phenos.fam <- aggregate(x=phenos$vol,by=list(phenos$group),mean)
rownames(phenos.fam) <- 1:56
fam.logcpm2 <- aggregate(t(logcpm2),by=list((phenos$group)),mean)
fam.logcpm2 <- fam.logcpm2[,-1]
scaled.fam.logcpm2 <- scale(fam.logcpm2)
identical(as.numeric(colnames(logcpm2)),phenos$animal)

#Create unrelated groups for CV####
cv.folds <- vector("list")
cv.folds[[1]] <- c(1,11,26,35,36,46,54)
cv.folds[[2]] <- c(2,9,12,20,27,37,47)
cv.folds[[3]] <- c(3,13,14,24,28,38,48)
cv.folds[[4]] <- c(4,15,29,39,43,45,49)
cv.folds[[5]] <- c(5,8,16,30,40,50,56)
cv.folds[[6]] <- c(6,17,22,31,33,41,51)
cv.folds[[7]] <- c(7,18,23,32,42,52,55)
cv.folds[[8]] <- c(10,19,21,25,34,44,53)
sort(unlist(cv.folds))
identical(as.numeric(colnames(logcpm2)),phenos$animal)


##Use different FDR thresholds to test prediction ranges

#use .1 threshold for prediction####
fam.cor.1 <- vector()
for(fam in 1:55){
  num.tests <- length(corrected.full.comparisons[[fam]])
  for(test.group in 1:num.tests){
    test <- unique(c((comparisons.fdr.1[[fam]][[test.group]])))
    log.DE <- scaled.fam.logcpm2[,match(test,colnames(scaled.fam.logcpm2))]
    cor.fam <- cor(t(log.DE))
    colnames(cor.fam) <- 1:56; rownames(cor.fam) <- 1:56
    all.fam <- 1:56
    pred <- vector();true <- vector(); all.cors <- vector()
    for(i in 1:8){
      set <- cv.folds[[i]]
      tmp <- okriging(idtest=set,idtrain = all.fam[-set],
                      corlist=list(cor.fam),H2vec=c(.99),pheno=phenos.fam,phenoname="x")
      pred <- c(pred,tmp[,2])
      true <- c(true,tmp[,3])}
    fam.cor.1 <- c(fam.cor.1,cor(pred,true))}
  print(fam)}

#use .05 threshold for prediction####
fam.cor.05 <- vector()
for(fam in 1:55){
  num.tests <- length(corrected.full.comparisons[[fam]])
  for(test.group in 1:num.tests){
    test <- unique(c((comparisons.fdr.05[[fam]][[test.group]])))
    log.DE <- scaled.fam.logcpm2[,match(test,colnames(scaled.fam.logcpm2))]
    cor.fam <- cor(t(log.DE))
    colnames(cor.fam) <- 1:56; rownames(cor.fam) <- 1:56
    all.fam <- 1:56
    pred <- vector();true <- vector(); all.cors <- vector()
    for(i in 1:8){
      set <- cv.folds[[i]]
      tmp <- okriging(idtest=set,idtrain = all.fam[-set],
                      corlist=list(cor.fam),H2vec=c(.99),pheno=phenos.fam,phenoname="x")
      pred <- c(pred,tmp[,2])
      true <- c(true,tmp[,3])}
    fam.cor.05 <- c(fam.cor.05,cor(pred,true))}
  print(fam)}

#use .01 threshold for prediction####
fam.cor.01.ht <- vector()
for(fam in 1:55){
  num.tests <- length(corrected.full.comparisons[[fam]])
  for(test.group in 1:num.tests){
    test <- unique(c((comparisons.fdr.01[[fam]][[test.group]])))
    if(length(test)==1){fam.cor.01.ht <- c(fam.cor.01.ht,NA)} else {
    log.DE <- scaled.fam.logcpm2[,match(test,colnames(scaled.fam.logcpm2))]
    cor.fam <- cor(t(log.DE))
    colnames(cor.fam) <- 1:56; rownames(cor.fam) <- 1:56
    all.fam <- 1:56
    pred <- vector();true <- vector(); all.cors <- vector()
    for(i in 1:8){
      set <- cv.folds[[i]]
      tmp <- okriging(idtest=set,idtrain = all.fam[-set],
                      corlist=list(cor.fam),H2vec=c(.99),pheno=phenos.fam,phenoname="ht")
      pred <- c(pred,tmp[,2])
      true <- c(true,tmp[,3])}
    fam.cor.01.ht <- c(fam.cor.01.ht,cor(pred,true))}}
  print(fam)}

#compare the three different thresholds####
hist(fam.cor.01)
hist(fam.cor.05)
hist(fam.cor.1)
boxplot(fam.cor.01)
boxplot(fam.cor.05)
boxplot(fam.cor.1)

summary(fam.cor.01);summary(fam.cor.05); summary(fam.cor.1)
 #Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#-0.3940  0.1522  0.2409  0.2299  0.3210  0.7595       1 
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#-0.3465  0.1856  0.2508  0.2525  0.3282  0.7088 
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#-0.2585  0.1907  0.2524  0.2543  0.3200  0.6220 
boxplot(fam.cor.1[365:412])
which(fam.cor.01 %in% sort(fam.cor.01,decreasing = T)[1:10])
which(fam.cor.01 > .6)
which(fam.cor.05 %in% sort(fam.cor.05,decreasing = T)[1:25])
which(fam.cor.05 > .6)
which(fam.cor.1 %in% sort(fam.cor.1,decreasing = T)[1:10])
which(fam.cor.1 > .6)


#Look at NCxFL x AL 8x17 with fdr####
test <- unique(c((comparisons.fdr.1[[8]][[9]])))
log.DE <- scaled.fam.logcpm2[,match(test,colnames(scaled.fam.logcpm2))]
cor.fam <- cor(t(log.DE))
colnames(cor.fam) <- 1:56; rownames(cor.fam) <- 1:56
all.fam <- 1:56
pred <- vector();true <- vector(); all.cors <- vector()
for(i in 1:8){
  set <- cv.folds[[i]]
  tmp <- okriging(idtest=set,idtrain = all.fam[-set],
                  corlist=list(cor.fam),H2vec=c(.99),pheno=phenos.fam,phenoname="x")
  print(cor(tmp[,2],tmp[,3]))
  pred <- c(pred,tmp[,2])
  true <- c(true,tmp[,3])}

plot(pred,true)

##Lets see how many times transcripts shows up for a family####
all.transcripts <- unlist(comparisons.fdr.05)
all.transcripts <- table(all.transcripts)
str(all.transcripts)
all.txpts <- matrix(NA,27949,1)
rownames(all.txpts) <- names(all.transcripts)
all.txpts[,1] <- all.transcripts
hist(all.txpts[,1])
summary(all.txpts[,1])
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#1.00   15.00   42.00   59.97   85.00  600.00 

#Try different number of transcripts to see how they do for prediction####
test <- rownames(all.txpts)[which(all.txpts[,1] > 1)]

phenos.fam <- aggregate(phenos[,4:7],by=list(phenos$group),mean)
log.DE <- scaled.fam.logcpm2[,match(comparisons.fdr.01[[15]][[34]],colnames(scaled.fam.logcpm2))]
top.txpt.de.fdr.01.ht <- colnames(log.DE)
#save(top.txpt.de.fdr.01.ht,file="~/Desktop/top.txpt.defdr.01.ht.RDA")
cor.fam <- cor(t(log.DE))
colnames(cor.fam) <- 1:56; rownames(cor.fam) <- 1:56
all.fam <- 1:56
pred <- vector();true <- vector()
for(i in 1:56){
  #set <- cv.folds[[i]]
  set <- i
  tmp <- okriging(idtest=set,idtrain = all.fam[-set],
                  corlist=list(cor.fam),H2vec=c(.99),pheno=phenos.fam,phenoname="ht")
  print(cor(tmp[,2],tmp[,3]))
  pred <- c(pred,tmp[,2])
  true <- c(true,tmp[,3])}
cor(pred,true)^2
summary(lm(pred~true))[8]
plot(pred,true,ylab = "true value",xlab="prediction")
legend(x='topleft', legend=paste("R^2 =",round((cor(pred,true)^2),2),sep=""),bty = "n",cex = .75)
pdf(file = "~/Desktop/10vs.49.vol.pdf",width = 4,height = 4)
plot(pred,true,xlab = "", ylab="", cex.axis=.8, cex.lab=1.1, cex=.5,
     ylim = c(90,145), xlim=c(95,145))
title(xlab = "prediction", ylab="true value", mgp=c(2,1,0),cex.lab=1)
legend(x='topleft', legend=paste("R^2 =",round((cor(true,pred)^2),2),sep=""),bty = "n",cex = .75)
dev.off()

cor(pred,true)^2

###subset all transcripts which pass .05 threshold####
all.fam.de <- vector("list",56)
for(fam in 1:55){
  d=fam
  fam.de <- vector()
  for(i in 1:fam){
    d <- d-1
    if(d > 0){
      fam.de <- c(fam.de,unlist(comparisons.fdr.05[[i]][[d]]))
    } else {
      fam.de <- c(fam.de,unlist(comparisons.fdr.05[[fam]]))
    }
  }
  all.fam.de[[fam]] <- fam.de}

fam=56
d=56
fam.de <- vector()
for(i in 1:55){
  d <- d-1
  if(d > 0){
    fam.de <- c(fam.de,unlist(comparisons.fdr.05[[i]][[d]]))}}

all.fam.de[[56]] <- fam.de

#all.fam.de contains the list of DE transcripts called for each family

#lets see how many DE per family
fam.de.count <- vector()
for(i in 1:56){
  fam.de.count <- c(fam.de.count,length(all.fam.de[[i]]))
}

summary(fam.de.count)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#3747   37390   54480   59860   77120  185500 
boxplot(fam.de.count)
which(fam.de.count > 185000)
#family 33 has most DE transcripts across 55 comparisons
#family 1 has least DE transcripts across 55 comparisons



unique.all.fam.de <- vector("list",56)
for(fam in 1:55){
  d=fam
  fam.de <- vector()
  for(i in 1:fam){
    d <- d-1
    if(d > 0){
      fam.de <- unique(c(fam.de,unlist(comparisons.fdr.05[[i]][[d]])))
    } else {
      fam.de <- unique(c(fam.de,unlist(comparisons.fdr.05[[fam]])))
    }
  }
  unique.all.fam.de[[fam]] <- fam.de}

fam=56
d=56
fam.de <- vector()
for(i in 1:55){
  d <- d-1
  if(d > 0){
    fam.de <- unique(c(fam.de,unlist(comparisons.fdr.05[[i]][[d]])))}}

unique.all.fam.de[[56]] <- fam.de

#lets see how many DE per family
unique.fam.de.count <- vector()
for(i in 1:56){
  unique.fam.de.count <- c(unique.fam.de.count,length(unique.all.fam.de[[i]]))
}

summary(unique.fam.de.count)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#1600   11610   13480   13620   16560   22220 
plot(unique.fam.de.count,fam.de.count)#roughly exponential?
boxplot(unique.fam.de.count)
which.max(unique.fam.de.count)
#family 7 has most unique DE transcripts across 55 comparisons
#family 1 has least unique DE transcripts across 55 comparisons

##
for(each in 1:56){
  test <- unique.all.fam.de[[each]]
  log.DE <- fam.logcpm2[,match(test,colnames(fam.logcpm2))]
  cor.fam <- cor(t(log.DE))
  colnames(cor.fam) <- 1:56; rownames(cor.fam) <- 1:56
  all.fam <- 1:56
  pred <- vector();true <- vector(); all.cors <- vector()
  for(i in 1:8){
    set <- cv.folds[[i]]
    tmp <- okriging(idtest=set,idtrain = all.fam[-set],
                    corlist=list(cor.fam),H2vec=c(.99),pheno=phenos.fam,phenoname="x")
    #print(cor(tmp[,2],tmp[,3]))
    pred <- c(pred,tmp[,2])
    true <- c(true,tmp[,3])}
  print(cor(pred,true))
  plot(pred,true)}
#all roughly r=.3


##identify transcripts from bad or good families and test those####
bad.transcripts <- unique(c(comparisons.fdr.05[[1]][[1]],comparisons.fdr.05[[1]][[2]],comparisons.fdr.05[[1]][[3]],
                            comparisons.fdr.05[[2]][[1]],comparisons.fdr.05[[2]][[2]],
                            comparisons.fdr.05[[3]][[1]]))


good.transcripts <- unique(c(comparisons.fdr.05[[55]][[1]],
                             comparisons.fdr.05[[54]][[2]],comparisons.fdr.05[[54]][[1]],
                             comparisons.fdr.05[[53]][[3]],comparisons.fdr.05[[53]][[2]],
                             comparisons.fdr.05[[53]][[1]]))


good.v.bad.txpts <- unique(c(comparisons.fdr.05[[1]][[55]],comparisons.fdr.05[[1]][[54]],comparisons.fdr.05[[1]][[53]],
                             comparisons.fdr.05[[2]][[54]],comparisons.fdr.05[[2]][[53]],comparisons.fdr.05[[2]][[52]],
                             comparisons.fdr.05[[3]][[53]],comparisons.fdr.05[[3]][[52]],comparisons.fdr.05[[3]][[51]]))


length(good.transcripts)#5399
length(bad.transcripts) #5355
length(good.v.bad.txpts) #4779

length(which(good.transcripts %in% bad.transcripts))#1146



test <- good.v.bad.txpts[!(good.v.bad.txpts %in% unique(c(good.transcripts,bad.transcripts)))]
log.DE <- fam.logcpm2[,match(test,colnames(scaled.fam.logcpm2))]
cor.fam <- cor(t(log.DE))
colnames(cor.fam) <- 1:56; rownames(cor.fam) <- 1:56
all.fam <- 1:56
pred <- vector();true <- vector(); all.cors <- vector()
for(i in 1:8){
  set <- cv.folds[[i]]
  tmp <- okriging(idtest=set,idtrain = all.fam[-set],
                  corlist=list(cor.fam),H2vec=c(.99),pheno=phenos.fam,phenoname="x")
  #print(cor(tmp[,2],tmp[,3]))
  pred <- c(pred,tmp[,2])
  true <- c(true,tmp[,3])}
print(cor(pred,true))
plot(pred,true)
#.45 for scaled
#.4 for unscaled



#AL ;17
#SCxTX; 21
#FLxTX; 13,38,43
#GAxTX; 35
#VA; 1
#GA; 23,33,34,36,38,48
#FL; 14,16,32,39,44,45,46,47,49,50,51,52,53,54,55
#NC; 2,3,4,5,7,15,20
#SC; 6,10,11,12,18,22,24,25,26,2728,29,30,31,37,40,41,42,56
#SCxGA; 9
#SCxSC; 19
#FLxFL; 38
#SCxTX; 21
#VA; 1

AL <- 17
SC.TX <- 21
FL.TX <- c(13,38,43)
VA <- 1
GA <-  c(23,33,34,36,38,48)
FL<- c( 14,16,32,39,44,45,46,47,49,50,51,52,53,54,55)
NC <- c(2,3,4,5,7,15,20)
SC <- c(6,10,11,12,18,22,24,25,26,2728,29,30,31,37,40,41,42,56)
SC.GA <- 9
SC.SC <-  19
FL.FL <-  38
SC.TX <-  21

###save
save.image("~/Desktop/corrected.identifying.top.DE.RDA.pt2.RData")
#AL vs. Diff Regions####

#vs. SC.TX #.139
all.fam[[AL]][[SC.TX]] #.139
ALa <- all.fam[[17]]
mean(ALa[c(13,37,42)]) #.18 vs. FL.TX
all.fam[[AL]][[1]]  # .0945 vs. VA

mean(ALa[c(23,33,34,36,38,48)]) #. 151 vs. GA
mean(ALa[c( 14,16,31,38,43,44,45,46,48,49,50,51,52,53,54)]) #. 151 vs. GA
mean(ALa[c(23,33,34,36,38,48)]) #. 151 vs. GA
mean(ALa[c(23,33,34,36,38,48)]) #. 151 vs. GA


mean(ALa[c(23,33,34,36,38,48)]) #. 151 vs. GA
mean(ALa[c(23,33,34,36,38,48)]) #. 151 vs. GA
mean(ALa[c(23,33,34,36,38,48)]) #. 151 vs. GA
mean(ALa[c(23,33,34,36,38,48)]) #. 151 vs. GAall.fam[[AL]][[SC.TX]] 
all.fam[[AL]][[SC.TX]] 
all.fam[[AL]][[SC.TX]] 
all.fam[[AL]][[SC.TX]] 
all.fam[[AL]][[SC.TX]] 
all.fam[[AL]][[SC.TX]] 
all.fam[[AL]][[SC.TX]] 
all.fam[[AL]][[SC.TX]] 

#est_pita_67707934
#AT3G62550.1  | Symbols:  | Adenine nucleotide alpha hydr...    82   2e-15  1 
#est_Pita_61779901
#AT4G11360.1  | Symbols: RHA1B | RING-H2 finger A1B | chr...    34   0.26   1
#est_pita_51498052
#AT5G14780.1  | Symbols: FDH | formate dehydrogenase | ch...   248   1e-65  2 
#isotig09939
#AT5G39040.1  | Symbols: ALS1, ATTAP2, TAP2 | transporter...   264   2e-70  3 
#est_Pita_66977002
#AT1G44740.1  | Symbols:  | unknown protein; FUNCTIONS IN...    34   0.57   1
#est_pita_66981047
#AT5G07990.1  | Symbols: TT7, CYP75B1, D501 | Cytochrome ...   150   2e-49  3 
#est_pita_67551644
#AT4G22020.1  | Symbols:  | pseudogene, hypothetical prot...    30   0.064  2 


##TRY PCA OR SVA####
#first identify top 10 correlations
fam.cor.05[which(fam.cor.05 %in% sort(fam.cor.05,decreasing=T)[1:10])]
which(fam.cor.05 %in% sort(fam.cor.05,decreasing=T)[1:10])
#47  318  373  386  393  483  794  836 1261 1352
total <- c(comparisons.fdr.05[[1]][[47]],comparisons.fdr.05[[7]][[3]],comparisons.fdr.05[[8]][[9]],
  comparisons.fdr.05[[8]][[22]],comparisons.fdr.05[[8]][[29]],comparisons.fdr.05[[10]][[24]],
  comparisons.fdr.05[[17]][[34]],comparisons.fdr.05[[18]][[37]],comparisons.fdr.05[[32]][[21]],
  comparisons.fdr.05[[37]][[2]])
tab.total <- table(total)
test <- names(which(tab.total > 2))
test <- unique(c(comparisons.fdr.05[[1]][[47]],comparisons.fdr.05[[7]][[3]],comparisons.fdr.05[[8]][[9]],
       comparisons.fdr.05[[8]][[22]],comparisons.fdr.05[[8]][[29]],comparisons.fdr.05[[10]][[24]],
       comparisons.fdr.05[[17]][[34]],comparisons.fdr.05[[18]][[37]],comparisons.fdr.05[[32]][[21]],
       comparisons.fdr.05[[37]][[2]]))
dim(scaled.fam.logcpm2)
fam.logcpm2.7.10 <- scaled.fam.logcpm2[,which(colnames(fam.logcpm2) %in% test)]
#fam.logcpm2.7.10 <- fam.logcpm2[,which(colnames(fam.logcpm2) %in% c(corrected.full.comparisons.fdr.05[[32]][[20]]))]

p.7.10 <- prcomp(fam.logcpm2.7.10)
str(p.7.10$rotation)
print(head(p.7.10))
summary(p.7.10)
plot(p.7.10,type="l")
biplot(p.7.10)
plot(p.7.10$x[-10,1],p.7.10$x[-10,2])
s.7.10 <- svd(x = fam.logcpm2.7.10)
plot(s.7.10$v[,1],s.7.10$v[,2])
plot(s.7.10$u[,1])
plot(s.7.10$d^2/sum(s.7.10$d^2), xlim = c(0, 56), type = "b", pch = 16, xlab = "principal components", 
     ylab = "variance explained")
plot(s.7.10$v[,1],s.7.10$v[,2])
library(randomForest)
randomForest()

out <- apply(fam.logcpm2,MARGIN = 1,FUN = function(x) sum(p.7.10$rotation[,1]*x)) 

pc.fam.cpm <- prcomp(fam.logcpm2.7.10, retx = TRUE, center = TRUE, scale. = TRUE)
uniq.fam <- as.character(unique(phenos$group))
fam.phenos <- data.frame(fam=uniq.fam,ht=phenos$ht[match(uniq.fam,phenos$group)],
                         vol=phenos$vol[match(uniq.fam,phenos$group)])
summary(fam.phenos$ht)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 96.8   105.2   108.6   108.0   111.2   119.2 
summary(fam.phenos$vol)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#    91.5   109.0   119.6   118.8   128.2   145.7 
vol.grps <- list(q1=which(fam.phenos$vol<109.0),
                 q2=intersect(which(fam.phenos$vol>=109.0),which(fam.phenos$vol<119.6)),
                 q3=intersect(which(fam.phenos$vol>=119.6),which(fam.phenos$vol<128.2)),
                 q4=which(fam.phenos$vol>=128.2))
max(pc.fam.cpm$x[,1]) # 33
min(pc.fam.cpm$x[,1]) # -20
max(pc.fam.cpm$x[,2]) # 25
min(pc.fam.cpm$x[,2]) # -14
plot(pc.fam.cpm$x[vol.grps$q1,1],pc.fam.cpm$x[vol.grps$q1,2],col="red",xlim=c(-21,35),
     ylim=c(-15,26))
points(pc.fam.cpm$x[vol.grps$q2,1],pc.fam.cpm$x[vol.grps$q2,2],col="orange")     
points(pc.fam.cpm$x[vol.grps$q3,1],pc.fam.cpm$x[vol.grps$q3,2],col="purple")     
points(pc.fam.cpm$x[vol.grps$q4,1],pc.fam.cpm$x[vol.grps$q4,2],col="blue")     
# No visible pattern of differences. Try calculating median value of PCs for each group,
# plot those to look for differences.  
med.q1 <- apply(pc.fam.cpm$x[vol.grps$q1,],2,median)
med.q2 <- apply(pc.fam.cpm$x[vol.grps$q2,],2,median)
med.q3 <- apply(pc.fam.cpm$x[vol.grps$q3,],2,median)
med.q4 <- apply(pc.fam.cpm$x[vol.grps$q4,],2,median)
max(c(med.q1,med.q2,med.q3,med.q4)) # 30
min(c(med.q1,med.q2,med.q3,med.q4)) # -50
plot(1:56,med.q1,col="red",pch=19,ylim=c(-3,3))
points(1:56,med.q2,col="orange",pch=20)
points(1:56,med.q3,col="purple",pch=16)
points(1:56,med.q4,col="blue",pch=20)
# PC2 and PC4 show scatter of median values - plot those against each other
plot(pc.fam.cpm$x[vol.grps$q1,3],pc.fam.cpm$x[vol.grps$q1,5],col="red",xlim=c(-10,10),
     ylim=c(-10,10),pch=19,xlab="PC2",ylab="PC4")
points(pc.fam.cpm$x[vol.grps$q2,3],pc.fam.cpm$x[vol.grps$q2,5],col="magenta",pch=20)     
points(pc.fam.cpm$x[vol.grps$q3,3],pc.fam.cpm$x[vol.grps$q3,5],col="purple",pch=16)     
points(pc.fam.cpm$x[vol.grps$q4,3],pc.fam.cpm$x[vol.grps$q4,5],col="blue",pch=20)  
# Some separation of groups is visible that accounts for differences in medians,
# but still lots of overlap. Patterns (if any) seem to be too complex to see in
# pairwise analysis of PCs.

