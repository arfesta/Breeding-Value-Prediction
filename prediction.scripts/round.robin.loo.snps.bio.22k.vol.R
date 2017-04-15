#Set WD####
setwd("/media/fairfax/sg3/Adam.RNA.SNPs/192.relatedness2")
library(parallel); library(qvalue); library(OmicKriging)

phenos <- read.table(file="/media/fairfax/seagate/Adam.RNA/analyses/combined/phenos/all.phenos.192.csv", header=T, sep= ",", row.names=1)
phenos192 <- phenos[order(phenos$animal, decreasing=F),]
phenos192$animal[1:48] <- phenos192$animal[1:48]-900; phenos192$animal[49:192] <- phenos192$animal[49:192]-10000
phenos192$group[which(phenos192$group==25)] <- 24
phenos192$group[which(phenos192$group %in% c(26:57))] <- phenos192$group[which(phenos192$group %in% c(26:57))] -1
phenos.192 <- phenos192 ; rownames(phenos.192) <- 1:192
phenos.fam <- aggregate(phenos192[,4:6],by=list(phenos192$group),mean)
rm(phenos)

df <- read.table("with.snp.eff/with.no.missing.snps/out.012", header=F,stringsAsFactors = F,row.names=1,na.strings = "-1")
k22.counts <- df
k22.counts <- as.matrix(sapply(k22.counts,as.numeric))
snp.fam <- aggregate(k22.counts,by=list(phenos192$group),mean,na.rm=T)
k22.snp.fam.df <- snp.fam[,-1]

rm(df,snp.fam)



lm.bio.rep.snp.22k.vol <- mclapply(1:56,FUN=function(x){
  select <- as.numeric(which(phenos192$group== x))
  test.logcpm <- k22.counts[-select,]
  test.phenos <- phenos192[-select,]
  rest <- c(1:56)[-x]
  test.apply <- lapply(1:55,function(x){
    r <- rest[x]
    se <- as.numeric(which(test.phenos$group == r))
    tester <- test.logcpm[-se,]
    tester.p <- test.phenos$vol[-se]
    te.2 <- apply(tester,2, function(x){
      lml <- lm(tester.p ~ x)
      tr <- unlist(summary(lml)[4])
      tr[8]})
    te.2})
  test.apply},mc.cores = 56)


all.55.22.qvalues <- mclapply(1:56,function(x){
  vars <- lm.bio.rep.snp.22k.vol[[x]]
  list.qvalues <- lapply(1:55,function(x){
    var2 <- vars[[x]]
    var.q <- qvalue(var2)
    var.qvalues <- var.q$qvalues
    var.qvalues
  })
  l.qvalues <- unlist(list.qvalues)
  l.qvalues
  #c <- table(names(l.qvalues)[which(l.qvalues <= .015)])                       
  #names(c)[which(c == 55)]
},mc.cores=56) 

se <- seq(.05,.06,.001)
for(each in 1:length(se)){
  all.vars <- vector(); l.vars <- vector()
  pred <- vector() ; true <- vector()
  for(i in 1:56){
    c <- table(names(all.55.22.qvalues[[i]])[which(all.55.22.qvalues[[i]] <= .005)])                       
    vars <- names(c)[which(c >= 55)]
    all.vars <- c(all.vars,vars)
    l.vars <- c(l.vars,length(vars))
    cor.matrix <- cor(t(k22.snp.fam.df[,match(try,colnames(k22.snp.fam.df))])); colnames(cor.matrix) <- 1:56; rownames(cor.matrix) <- 1:56
    test <-  i
    tmp <- okriging(idtest=test,idtrain=seq(1,56,1)[-c(test)],corlist = list(cor.matrix),
                    H2vec = c(.9),pheno = phenos.fam,phenoname = "vol")
    pred <- c(pred,tmp[,2])
    true <- c(true,tmp[,3])
    print(cor(tmp[,2],tmp[,3]))}
  plot(pred,true)
  print(cor(pred,true))^2
  print(se[[each]])
}
abline(h=mean(true),col="red")
abline(v=mean(pred),col="red")
#best result: qvalue .05; weight .8; .47
pdf(file = "~/Desktop/snp.lmoo.vol.pdf",width = 4,height = 4)
plot(pred,true,xlab = "", ylab="", cex.axis=.8, cex.lab=1.1, cex=.5,
     ylim = c(95,120), xlim=c(103,113))
title(xlab = "prediction", ylab="true value", mgp=c(2,1,0),cex.lab=1)
legend(x='topleft', legend=paste("R^2 =",round((cor(true,pred)^2),2),sep=""),bty = "n",cex = .75)
dev.off()

summary(l.vars)
length(all.vars)
try <- unique(all.vars)

snp.vol.list.vars <- vector("list")
for(i in 1:56){
  c <- table(names(all.55.22.qvalues[[i]])[which(all.55.22.qvalues[[i]] <= .005)])                       
  vars <- names(c)[which(c >= 55)]
  snp.vol.list.vars[[i]] <- vars}
