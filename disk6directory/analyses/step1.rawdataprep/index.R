#####
setwd("/media/disk6/ARF/RNASEQ/trimmedfiltreads/86k")
 the.call <-  paste('cd EW; rm all_indicies.txt; while read file; do gzip -cd $file |head -100 | grep "^@D00555" | cut -d : -f 10 | sort | uniq -c | sort -nr  >> all_indicies.txt; done < all.ew.reps',sep="")
system(the.call,ignore.stderr = T)
ew.index.table <- read.table("./EW/all_indicies.txt")
ew.index.table <- table(ew.index.table[,2])
true.EW.index <- names(which(ew.index.table == 12))


the.call <-  paste0('cd LGEP; rm all_indicies.txt; while read file; do gzip -cd $file |head -200 | grep "^@HISEQ" | cut -d : -f 10 | sort | uniq -c | sort -nr  >> all_indicies.txt; done < all.lgep.reps')
system(the.call,ignore.stderr = T)

lgep.index.table <- read.table("./LGEP/all_indicies.txt")
lgep.index.table <- table(lgep.index.table[,2])
true.LGEP.index <- names(which(lgep.index.table == 18))

table(c(true.EW.index,true.LGEP.index))
#Same indicies were used
#ACAGTG ACTGAT ACTTGA AGTCAA AGTTCC ATCACG ATGTCA ATTCCT CAGATC CCGTCC CGATGT CGTACG CTTGTA GAGTGG GATCAG GCCAAT GGCTAC GTCCGC 
#     2      2      2      2      2      2      2      2      2      2      2      2      2      2      2      2      2      2 
#GTGAAA GTGGCC GTTTCG TAGCTT TGACCA TTAGGC 
#     2      2      2      2      2      2 
#     
identical(true.EW.index,true.LGEP.index)
#TRUE
#
all.ew.reps <- read.table("./EW/all.ew.reps")
#We have the true barcodes for each batch now lets assign to samples
#ew.index.table <- rep(NA,288)
ew.indicies <- vector()
for(each.file in 1:nrow(all.ew.reps) ) {
  file <- as.character(all.ew.reps$V1[each.file])

  the.call <- paste0('gzip -cd ./EW/', file,' | head -10000 | grep "^@D00555" | cut -d : -f 10 > ../out.txt')
  system(the.call)
  index.count <- table(read.table(file = "/media/disk6/ARF/RNASEQ/trimmedfiltreads/out.txt")[,1])
index.name <- names(index.count)[which(names(index.count) %in% true.EW.index)]
ew.indicies <- c(ew.indicies,index.name)
}
 
table(ew.indicies)
#ACAGTG ACTGAT ACTTGA AGTCAA AGTTCC ATCACG ATGTCA ATTCCT CAGATC CCGTCC CGATGT CGTACG CTTGTA GAGTGG GATCAG GCCAAT GGCTAC GTCCGC 
#    12     12     12     12     12     12     12     12     12     12     12     12     12     12     12     12     12     12 
#GTGAAA GTGGCC GTTTCG TAGCTT TGACCA TTAGGC 
#    12     12     12     12     12     12 
#    
all.lgep.reps <- read.table("./LGEP/all.lgep.reps")
lgep.indicies <- vector()
for(each.file in 1:nrow(all.lgep.reps) ) {
  file <- as.character(all.lgep.reps$V1[each.file])

  the.call <- paste0('gzip -cd ./LGEP/', file,' | head -10000 | grep "^@HISEQ" | cut -d : -f 10 > ../out.txt')
  system(the.call)
  index.count <- table(read.table(file = "/media/disk6/ARF/RNASEQ/trimmedfiltreads/out.txt")[,1])
index.name <- names(index.count)[which(names(index.count) %in% true.LGEP.index)]
lgep.indicies <- c(lgep.indicies,index.name)
}
 
table(lgep.indicies)
#ACAGTG ACTGAT ACTTGA AGTCAA AGTTCC ATCACG ATGTCA ATTCCT CAGATC CCGTCC CGATGT CGTACG CTTGTA GAGTGG GATCAG GCCAAT GGCTAC GTCCGC 
#    18     18     18     18     18     18     18     18     18     18     18     18     18     18     18     18     18     18 
#GTGAAA GTGGCC GTTTCG TAGCTT TGACCA TTAGGC 
#    18     18     18     18     18     18 

#Output index and save file#####
index=c(ew.indicies,lgep.indicies)
file_name=c(as.character(all.ew.reps[,1]),as.character(all.lgep.reps[,1]))
tech.rep.index <- cbind("index"=index,"file_name"=file_name)
save(tech.rep.index,file = "/media/disk6/ARF/RNASEQ/Breeding-Value-Prediction/disk6directory/resources/exptdesign/sequencing/index.RData",compress=T)

