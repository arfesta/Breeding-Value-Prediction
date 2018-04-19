#Set WD to top level RNA seq analysis####
setwd("/media/disk6/ARF/RNASEQ/analyses")
packages <- c("readr","sommer","dplyr","parallel","pedigree","reshape2","readxl","Metrics","caret")
sapply(packages,library,character.only=T)
rm(packages);gc()

#Load salmon count data and experimental data####
load("./step1/salmon_count_data.RData")
load("./step1/salmon_expt_data.RData") 
expt.dat.720 <- as_tibble(exp.dat.720); rm(exp.dat.720)

#Load techincal rep data frame with volume information for all 720 samples####
bio.exp.dat <- as_tibble(read_excel("../resources/phenos/All.xlsx")) %>%
          mutate(vol = round(vol,digits = 1)) %>%
          filter(ht > 0) %>%
          dplyr::select(animal,ID,vol,ht) %>%
          dplyr::rename(animal.id=animal) %>%
          distinct()
head(bio.exp.dat)

ew.fams <- unique(unlist(bio.exp.dat[which(bio.exp.dat$animal.id < 1000),"ID"]))
lgep.fams <- unique(unlist(bio.exp.dat[which(bio.exp.dat$animal.id > 1000),"ID"]))
ew.lgep.fams <- ew.fams[which(ew.fams %in% lgep.fams)]
ew.fams <- ew.fams[-c(which(ew.fams %in% ew.lgep.fams))]
lgep.fams <- lgep.fams[-c(which(lgep.fams %in% ew.lgep.fams))]

#Merge experimental data from salmon with volume and family ID's####
expt.dat.576 <- inner_join(x = bio.exp.dat,y = expt.dat.720,by="animal.id")
expt.dat.576$lane.id[1:144] <- paste("ew_",expt.dat.576$lane.id[1:144],sep="")
expt.dat.576$lane.id[140:145]


#Merge raw counts with volume and family ID's (to ensure correct order as above join)####
salmon.raw.counts <- data.frame("animal.id" = expt.dat.720$animal.id,t(salmon.counts$counts),stringsAsFactors = F) %>%
  inner_join(x = bio.exp.dat,y = ., by=c("animal.id")) %>%
  subset(x=.,select=-c(animal.id,ID,vol,ht))
rm(expt.dat.720,salmon.counts,bio.exp.dat); gc()
#Load pedigree:####
pedigree <- read_excel("../resources/pedigree/pedigree_2018_02_18_17_28_36.xlsx") %>%
  mutate(animal.id=rep(1:length(unique(Id)))) %>%
  mutate(p1.animal.id=rep(NA,length(Id))) %>%
  mutate(p2.animal.id=rep(NA,length(Id)))

pa.num <- list(p1.T = which(is.na(pedigree$`Parent 1`) == F),
               p2.T = which(is.na(pedigree$`Parent 2`) == F))

for(each.p in 1:length(pa.num$p1.T)){pedigree$p1.animal.id[pa.num$p1.T[each.p]] <- which(pedigree$Id == pedigree$`Parent 1`[pa.num$p1.T[each.p]])}
for(each.p in 1: length(pa.num$p2.T)){pedigree$p2.animal.id[pa.num$p2.T[each.p]] <- which(pedigree$Id == pedigree$`Parent 2`[pa.num$p2.T[each.p]])}
rm(pa.num,each.p)
head(pedigree)

#Create progeny pedigree####
colnames(expt.dat.576)
expt.dat.192 <- unique(subset.data.frame(expt.dat.576,select = c("animal.id","ID","vol","ht")))

progeny.pedigree <- as_data_frame(cbind("Id" = expt.dat.192$ID,
                                        "Parent 1" = rep(NA,192),
                                        "Parent 2" = rep(NA,192),
                                        "County" = rep(NA,192),
                                        "State" = rep(NA,192),
                                        "animal.id"=expt.dat.192$animal.id,
                                        "p1.animal.id"= rep(NA,192),
                                        "p2.animal.id" = rep(NA,192)))

#Families which are crosses have "x" in the middle, let's split that:
prog.ped.fam <- sapply(expt.dat.192$ID,function(x) strsplit(x = x,split = "x"))
head(prog.ped.fam)
#Each biological rep has at least 1 parent because it is at a minimum OP.  Subset that first parent
progeny.pedigree$`Parent 1` <- unlist(lapply(prog.ped.fam,function(x) x[[1]]))

#Only FS crosses have a parent 2.  Subset those FS crosses and isolate parent 2.
prog.ped.p2.list <- which(sapply(1:192,function(x) length(prog.ped.fam[[x]]) > 1) == T)
progeny.pedigree$`Parent 2`[prog.ped.p2.list] <- sapply(prog.ped.fam[prog.ped.p2.list],function(x) x[2])
rm(prog.ped.p2.list,prog.ped.fam)
head(progeny.pedigree)

#create numeric par id
progeny.pedigree$p1.animal.id <- sapply(1:192,function(x) which(pedigree$Id == progeny.pedigree$`Parent 1`[x]))
progeny.pedigree$p2.animal.id  <- sapply(1:192,function(x) {
  this.par <- which(pedigree$Id == progeny.pedigree$`Parent 2`[x])
  if(length(this.par) < 1){ NA} else {this.par}
})

progeny.pedigree$County  <- sapply(1:192,function(x) {pedigree$County[which(pedigree$Id == progeny.pedigree$`Parent 1`[x])]})
progeny.pedigree$State  <- sapply(1:192,function(x) {pedigree$State[which(pedigree$Id == progeny.pedigree$`Parent 1`[x])]})

#Join parent and progeny pedigree, create A matrix####
full.pedigree <- rbind(pedigree,progeny.pedigree)
A.ped <- data.frame(cbind("animal.id"=(as.character(full.pedigree$animal.id)),
                          "p1.animal.id"=(as.character(full.pedigree$p1.animal.id)),
                          "p2.animal.id"=(as.character(full.pedigree$p2.animal.id))))
str(A.ped)
makeA(A.ped,which = rep(T,nrow(A.ped)))
A <- read.table("A.txt") %>% acast(.,V1~V2, value.var="V3"); 
dimnames(A) <- list(A.ped$animal.id,A.ped$animal.id); 

order.Amat <- match(as.character(expt.dat.192$animal.id),rownames(A))
identical(rownames(A)[order.Amat],as.character(expt.dat.192$animal.id))
#TRUE
A <- A[order.Amat,order.Amat]
A[upper.tri(A)] <- t(A)[upper.tri(A)]
expt.dat.192$animal.id <- as.character(expt.dat.192$animal.id)
rm(order.Amat,pedigree,progeny.pedigree)

#Normalize all of the data together####
#Now remove txpts which have mean 0 (No counts)
fam.phenos <- unique(expt.dat.192[,c("ID","vol")])

remove.zeros <- apply(salmon.raw.counts,2,sd)
rm.txpt <-which(remove.zeros <= 1)
removal.stats <- apply(salmon.raw.counts[,rm.txpt],2,sum)
hist(removal.stats)
summary(removal.stats)

unique.animals <- unique(expt.dat.576$animal.id)
tmp.dat <- salmon.raw.counts[,-rm.txpt]; 

tip.dat <- data.frame("index.id" = as.factor(expt.dat.576$true.index),
                      "animal.id" = as.factor(expt.dat.576$animal.id),
                      "lane.id"=as.factor(expt.dat.576$lane.id),
                      "batch" = as.factor(c(rep(1,144),rep(2,432)))
)

salmon.bio.norm.counts <- mclapply(1:ncol(tmp.dat),function(each.txpt) {
  tip.dat$txpt <- unlist(log2(tmp.dat[,each.txpt]+1))

  mm.mod<- mmer2(fixed = txpt ~ index.id + lane.id ,random = ~ g(animal.id),G = list(animal.id=A),
                 data = tip.dat,iters=30,silent = T)
  mm.mod$u.hat$`g(animal.id)`[,1] +  mm.mod$beta.hat[1,1]
},mc.cores=60)

r.names <- names(salmon.bio.norm.counts[[1]])
salmon.bio.norm.counts <- do.call(cbind,salmon.bio.norm.counts)
rownames(salmon.bio.norm.counts) <- r.names; rm(r.names)
colnames(salmon.bio.norm.counts) <- colnames(tmp.dat)
save(salmon.bio.norm.counts,file = "./step2.normalization/dplyr.tech.bio.norm/192.norm.log.count.lig.Rdata",compress=T)
save(expt.dat.192,file = "./step2.normalization/dplyr.tech.bio.norm/192.expt.data.Rdata",compress=T)
save.image(file = "./step2.normalization/dplyr.tech.bio.norm/dplyr.normalization.all.bio.RData",compress=T)
