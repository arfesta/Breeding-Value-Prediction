#This codes bring in countfiles from sailfish and then combines them into all good Tip samples
##Bring in countfiles from .sf files after sailfish and put into data set####
#creates a vector of file names within the working direcotry that have the associated pattern
files <- as.list(dir("/media/seagate2/Adam.RNASEQ/countfiles/LGEP/countfiles",recursive=T, full.names=T, pattern="quant_bias_corrected.sf"))
#function which processes each file by reading in the table and then subsetting for 7 because that contains
#read counts
processfile<- function(f) { a<- read.table(f)[c("V7")]}

#using apply(list style) for all of the files and running the function process file
#then converting to data frame
#gene name is in every other column which is unneccesary, so it is removed by first creating a vector of numbers
#which correspond to the column that we want to keep in the data set (read counts) 
library(parallel)
lgep.countdata <- mclapply(files,processfile,mc.cores=22 ); lgep.countdata <- do.call(cbind,lgep.countdata)
rownames(lgep.countdata) <- read.table(files[[1]])[c("V1")][1:86008,1]
files <- as.list(dir("/media/seagate2/Adam.RNASEQ/countfiles/LGEP/countfiles",recursive=T, full.names=F, pattern="quant_bias_corrected.sf"))
a <- regmatches(unlist(files), regexpr("[0-9].*[0-9]", unlist(files))); a <- gsub( "2_", "_", a)
colnames(lgep.countdata) <- a ; rm(a,files)
lgep.countdata <- lgep.countdata[,sort(colnames(lgep.countdata))]


files <- as.list(dir("/media/seagate2/Adam.RNASEQ/countfiles/EastWest/raw.countfiles",recursive=T, full.names=T, pattern="quant_bias_corrected.sf"))
ew.countdata <- mclapply(files,processfile,mc.cores=22 ); ew.countdata <- do.call(cbind,ew.countdata)
files <- as.list(dir("/media/seagate2/Adam.RNASEQ/countfiles/EastWest/raw.countfiles",recursive=T, full.names=F, pattern="quant_bias_corrected.sf"))
a <- regmatches(unlist(files), regexpr("[0-9].*[0-9]", unlist(files)))
B <- gsub('_S[[:digit:]][[:digit:]]_','_',a); colnames(ew.countdata) <- B
B <- gsub('_S[[:digit:]]_','_',colnames(ew.countdata)); colnames(ew.countdata) <- B
ew.countdata <- ew.countdata[ ,sort(colnames(ew.countdata))]
rm(a,B,files,processfile)

#sorting countdata for nolumn names
full.countdata <- cbind(lgep.countdata,ew.countdata)
full.countdata <- full.countdata[,sort(colnames(full.countdata))]
#rm(ew.countdata,lgep.countdata)


all.data <- read.csv(file="/media/seagate/Adam.RNA/analyses/combined/phenos/All.tech.reps.rfsrc.csv", header=T, sep=",")
expt.data <- all.data[order(all.data$animal),]
Bad.samples <- which(expt.data$ht==0)
expt.data <- expt.data[-Bad.samples,]
ord <- match(expt.data$row.name,colnames(full.countdata))
tip.countdata<- as.matrix(full.countdata[,ord])
identical(colnames(tip.countdata),as.character(expt.data$row.name))
#TRUE
dim(tip.countdata) #86008x576
#filter to include transcripts with at least greater than 1 count per technical rep
tip.countdata <- tip.countdata[rowSums(tip.countdata) > 576,]
dim(tip.countdata) #46782x576

###Load pedigree and create a inverse####
library(asreml)
pedigree <- read.csv(file="/media/seagate/Adam.RNA/analyses/combined/phenos/Full.pedigree.csv",header=T)
pedigree <- within(pedigree, {animal <- factor(animal); dam <- factor(dam)
sire <- factor(sire)})
yo <- unique(match(expt.data$animal,pedigree$animal))
pedigree <- pedigree[c(1:97,yo),]
ainv <- asreml.Ainverse(pedigree)$ginv
##########
practice <- cbind(expt.data[,c(2:4,9)],t(tip.countdata1))

# First taking log2 of counts and binding with expt data;  convert animal to factor####
log.counts.samp <- data.frame(cbind(practice[,1:4], log2(practice[,5:46786] + 1)))
log.counts.samp <- within(log.counts.samp, {animal <- factor(animal)})

###Run ASREML normalization on all families ####
tip.full.sample.est <- matrix(NA,nrow=192,ncol=46782)
rownames(tip.full.sample.est)<-as.character(pedigree[98:289,1])
colnames(tip.full.sample.est)<-names(log.counts.samp)[5:46786]

log.lik.lane0.5 <- list("vector")
#Run ASREML 46782
for(i in 1:46782){ y <- log.counts.samp[,i+4]
tmp <- asreml(fixed = y ~ lane + index, random = ~ ide(animal,var=T,init=1),ginverse=list(animal=ainv),
              data=log.counts.samp,maxiter=30)
tip.full.sample.est[,i] <- tmp$coefficients$random[98:289]+tmp$coefficients$fixed["(Intercept)"]
log.lik.lane0.5[i] <- tmp$loglik; rm(y);print(i)}

#save(tip.full.sample.est,file="/media/seagate/Adam.RNA/analyses/combined/tip.full.sample.est.rda")
