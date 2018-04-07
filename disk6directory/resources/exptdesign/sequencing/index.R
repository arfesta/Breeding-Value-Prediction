#####
setwd("/media/disk6/ARF/RNASEQ/trimmedfiltreads")
 the.call <-  paste('cd EW; rm all_indicies.txt; while read file; do gzip -cd $file |head -100 | grep "^@D00555" | cut -d : -f 10 | sort | uniq -c | sort -nr  >> all_indicies.txt; done < all.ew.reps',sep="")
system(the.call,ignore.stderr = T)
ew.index.table <- read.table("./EW/all_indicies.txt")
ew.index.table <- table(ew.index.table[,2])
true.EW.index <- names(which(ew.index.table == 12))


the.call <-  paste('cd LGEP; rm all_indicies.txt; while read file; do gzip -cd $file |head -100 | grep "^@HISEQ" | cut -d : -f 10 | sort | uniq -c | sort -nr  >> all_indicies.txt; done < all.lgep.reps',sep="")
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

#We have the true barcodes for each batch now lets assign to samples
i = 1
ew.index.table <- rep(NA,288)
while(anyNA(ew.index.table)){
 the.call <-  paste('cd EW; rm all_indicies.txt; while read file; do gzip -cd $file |head -1000 | grep "^@D00555" | sed "',as.character(i +1),'q;d" | cut -d : -f 10 | sort | uniq -c | sort -nr  >> all_indicies.txt; done < all.ew.reps',sep="")
system(the.call,ignore.stderr = T)
ew.index <- read.table("./EW/all_indicies.txt")
correct <- which(ew.index[,2] %in% true.EW.index)
ew.index.table[correct] <- as.character(ew.index[correct,2])
i <- i + 1
}
table(ew.index.table)
#ACAGTG ACTGAT ACTTGA AGTCAA AGTTCC ATCACG ATGTCA ATTCCT CAGATC CCGTCC CGATGT CGTACG CTTGTA GAGTGG GATCAG GCCAAT GGCTAC GTCCGC 
#    12     12     12     12     12     12     12     12     12     12     12     12     12     12     12     12     12     12 
#GTGAAA GTGGCC GTTTCG TAGCTT TGACCA TTAGGC 
#    12     12     12     12     12     12 
#    

i=1
lgep.index.table <-  rep(NA,432)
while(anyNA(lgep.index.table) ){
  the.call <-  paste('cd ./LGEP; rm all_indicies.txt; while read file; do gzip -cd $file |head -1000 | grep "^@HISEQ" | sed "',as.character(i +1),'q;d" |  cut -d : -f 10 | sort | uniq -c | sort -nr  >> all_indicies.txt; done < all.lgep.reps',sep="")
system(the.call,ignore.stderr = T)

lgep.index <- read.table("./LGEP/all_indicies.txt")
correct <- which(lgep.index[,2] %in% true.LGEP.index)
lgep.index.table[correct] <- as.character(lgep.index[correct,2])
i <- i+1
}
table(lgep.index.table)
#ACAGTG ACTGAT ACTTGA AGTCAA AGTTCC ATCACG ATGTCA ATTCCT CAGATC CCGTCC CGATGT CGTACG CTTGTA GAGTGG GATCAG GCCAAT GGCTAC GTCCGC 
#    18     18     18     18     18     18     18     18     18     18     18     18     18     18     18     18     18     18 
#GTGAAA GTGGCC GTTTCG TAGCTT TGACCA TTAGGC 
#    18     18     18     18     18     18 

#Output index and save file#####
ew.filenames <- read.table("./EW/all.ew.reps")
lgep.filenames <- read.table("./LGEP/all.lgep.reps")
index=c(ew.index.table,lgep.index.table)
file_name=c(as.character(ew.filenames[,1]),as.character(lgep.filenames[,1]))
tech.rep.index <- cbind("index"=index,"file_name"=file_name)
save(tech.rep.index,file = "../resources/exptdesign/index.RData",compress=T)

