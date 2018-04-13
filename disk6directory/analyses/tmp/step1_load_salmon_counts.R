#Set WD to top level RNA seq analysis
setwd(dir = "/media/disk6/ARF/RNASEQ/")

#All count file names are "quant.sf", we can get the sample names from each respective file path:
countfiles = list.files(path = "./counts/86kSalmon",pattern = "quant.sf",recursive = T,full.names = T)
head(countfiles)

#Extract meta data from file path for all 720 technical replicates:####
#Extract the lane by splitting each file name on "/" and selecting the 5th element, which corresponds to lane id.
ew.exp.dat <- lapply(1:288,function(each.file){
  split.file.name <- strsplit(x=countfiles[each.file],split = "/")[[1]]
  split.samp.name <- strsplit(split.file.name[6],split = "_")[[1]]
  
  lane.name <- split.file.name[5]
  samp.name <- split.samp.name[1]
  
  
  index <- split.samp.name[2]
  index <- as.numeric(gsub(pattern = "S",replacement = "",x = index))
  if(index < 25){index <- index} else if(index > 24 & index < 49){index <- index - 24}
  if(index > 48 & index < 73){index <- index - 24*2}
  if(index > 72) {index <- index - 24*3}
  index <- LETTERS[index]
  
  return(list("sample.name" = samp.name,
              "lane.id" = lane.name,
              "index.id" = index))
})
head(ew.exp.dat)
ew.exp.dat <- data.frame(stringsAsFactors = F,
      "animal.id" = as.numeric(sapply(ew.exp.dat,function(each.file) smp <- each.file$sample.name,simplify =T)) + 100,
      "lane.id" = sapply(ew.exp.dat,function(each.file) smp <- each.file$lane.id,simplify =T),
      "index.id" = sapply(ew.exp.dat,function(each.file) smp <- each.file$index.id,simplify =T))

lgep.exp.dat <-  lapply(289:720,function(each.file){
  split.file.name <- strsplit(x=countfiles[each.file],split = "/")[[1]]
  
  lane.name <- split.file.name[5]
  samp.name <- (gsub("([0-9]+).*$", "\\1", strsplit(split.file.name[6],split = "_")[[1]][2]))
  index <- regmatches(split.file.name[6], regexpr("[^.](?=\\.)", split.file.name[6], perl = TRUE))
  
  return(list("sample.name" = samp.name,
              "lane.id" = lane.name,
              "index.id" = index))
})

lgep.exp.dat <- data.frame(stringsAsFactors = F,
                           "animal.id" = as.numeric(sapply(lgep.exp.dat,function(each.file) smp <- each.file$sample.name,simplify =T)) + 1000,
                           "lane.id" = sapply(lgep.exp.dat,function(each.file) smp <- each.file$lane.id,simplify =T),
                           "index.id" = sapply(lgep.exp.dat,function(each.file) smp <- each.file$index.id,simplify =T))
exp.dat <- rbind(ew.exp.dat,lgep.exp.dat)
str(exp.dat)
rm(ew.exp.dat,lgep.exp.dat)

#Load counts using tximport####
library(tximport)
salmon.counts <- tximport(files=countfiles,type = "salmon",txIn = T,countsFromAbundance = "no",txOut = T)

#The length matrix contains the average transcript length for each gene which can be used as an offset for gene-level analysis.
#Load in the correct index data frame; match  file names first to make sure correct assignment
load("./resources/exptdesign/index.RData")

countfiles = list.files(path = "./counts/86kSalmon",pattern = "quant.sf",recursive = T,full.names = T)
head(countfiles)
split.file.name <- strsplit(x=countfiles[1],split = "/")[[1]]
split.file.name.i <- strsplit(x=as.character(tech.rep.index[,2])[1],split = "/")[[1]]

#We want to collect 5&6 from the countfiles and 2&3 from the index data frame
count.file.split <- unlist(lapply(1:720,function(each.file){
  split.file.name <- strsplit(x=countfiles[each.file],split = "/")[[1]]
  split.file.name <- split.file.name[c(5,6)]
  split.file.name <- paste(split.file.name[1],split.file.name[2],sep="_")}))
  
index.file.split <- unlist(lapply(1:720,function(each.file){
split.file.name.i <- strsplit(x=as.character(tech.rep.index[each.file,2]),split = "/")[[1]]
  split.file.name.i <- split.file.name.i[c(2,3)]
  split.file.name.i <- paste(split.file.name.i[1],split.file.name.i[2],sep="_")}))

#match file names and replace index in exp dat
new.order <- match(count.file.split,index.file.split)
identical(count.file.split,index.file.split[new.order]) 
exp.dat$true.index <- tech.rep.index[new.order,1]
#Save loaded counts and expt data:
save(exp.dat,file="./analyses/step1/salmon_expt_data.RData",compress=T)
save(salmon.counts,file="./analyses/step1/salmon_count_data.RData",compress=T)
