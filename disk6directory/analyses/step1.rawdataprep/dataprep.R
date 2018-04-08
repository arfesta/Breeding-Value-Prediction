# Purpose of this script is to combine files in experimental folder for input to salmon and modeling ####

# Index.R was ran to confirm/identify the indicies used in experimental design
 load("/media/disk6/ARF/RNASEQ/Breeding-Value-Prediction/disk6directory/resources/exptdesign/sequencing/index.RData")

# Using the file.name path and index, we can generate an experimental design table
# First split the file path for each technical replicate read fastq file
file.path.split <- strsplit(x = tech.rep.index[,2],split = "/")
(file.path.split[[1]])

# file.path.split is a list equal to the length of tech rep files. The second element in each of the lists is the lane.  For EW we will paste "ew_" in front of the lane to make the lane names unique from LGEP
lane.names <- unlist(lapply(file.path.split,function(x) x[2]))
lane.names[1:288] <- paste0("ew_",lane.names[1:288])
lane.names[c(1,4,290,500)]

# The third element of each list is the complete sample name with file extension. 
comp.sample.names <- unlist(lapply(file.path.split,function(x) x[3]))

# Split the complete name and retrieve sample id # for EW
split.samp.names <- strsplit(comp.sample.names,split = "_")
split.samp.names[[1]]
ew.sample.ids <- unlist(lapply(split.samp.names[1:288],function(x) x[1]))  
head(ew.sample.ids)

# LGEP have index and extension still attached which must be removed
split.samp.names[[720]]
lgep.sample.ids <- (gsub("([0-9]+).*$", "\\1", unlist(lapply(split.samp.names[-c(1:288)],function(x) x[2]))))

# We can obtain the letter by using regmatch to place as index letter
index.let <- regmatches(unlist(lapply(split.samp.names[-c(1:288)],function(x) x[2])), regexpr("[^.](?=\\.)", unlist(lapply(split.samp.names[-c(1:288)],function(x) x[2])), perl = TRUE))

# First 288 samples are EW and the next 432 are LGEP
batch.name <- c(rep("EW",288),rep("LGEP",432))

# Add 1000 to numeric EW sample ids so that they are unique from LGEP
animal.id <- as.character(c(as.numeric(ew.sample.ids) + 1000, as.numeric(lgep.sample.ids)))

# Store original sample ids
sample.id <- c(ew.sample.ids,lgep.sample.ids)

# Temporarily join and subset the unique combo of index_sequence and index_letter from LGEP to see which sequences correspond to which letters
tmp.dat <- unique(cbind("index_seq"=tech.rep.index[-c(1:288),1],index.let))
head(tmp.dat)

# Assign each tech rep an index letter based on matching to the unique subset defined above
index.let <- tmp.dat[match(tech.rep.index[,1],tmp.dat[,1]),2]

# Create unique copy of numeric sample id + new letter index
samp.index <- paste0(sample.id,index.let)

# Load library sizes, estimate FL's and assign to tech reps
library(readr)
EW_lib_sizes <- rbind(read_csv("/media/disk6/ARF/RNASEQ/Breeding-Value-Prediction/disk6directory/resources/exptdesign/sequencing/EW_lib_sizes_1.csv"),read_csv("/media/disk6/ARF/RNASEQ/Breeding-Value-Prediction/disk6directory/resources/exptdesign/sequencing/EW_lib_sizes_2.csv"))

LGEP_pool_info <- read_csv("/media/disk6/ARF/RNASEQ/Breeding-Value-Prediction/disk6directory/resources/exptdesign/sequencing/LGEP_pool_info.csv")

# The LGEP pool currently only has pool id in columns, we need the lane name to join with sample id so that we can assign FL lengths
# Identify which lanes need to have a 0 added in front (i.e. lane01 instead of lane1)
add.zero <- which(LGEP_pool_info$lane_number < 10)
LGEP_pool_info$lane_id= c(paste0("lane",0,LGEP_pool_info$lane_number[add.zero]),paste0("lane",LGEP_pool_info$lane_number[-add.zero]))

# Estimate FL's from library sizes:
EW_lib_sizes$FL = EW_lib_sizes$`Average Size [bp]` - 122
LGEP_pool_info$FL = LGEP_pool_info$lib_size - 122
sd.LGEP_fl <- sd(LGEP_pool_info$FL)

# Check that assignment will result in correct order
identical(samp.index[-c(1:288)],LGEP_pool_info$smp.name[match(paste0(lane.names[-c(1:288)],samp.index[-c(1:288)]),paste0(LGEP_pool_info$lane_id,LGEP_pool_info$smp.name))])
lgep.fl.lengths <- LGEP_pool_info$FL[match(paste0(lane.names[-c(1:288)],samp.index[-c(1:288)]),paste0(LGEP_pool_info$lane_id,LGEP_pool_info$smp.name))]

mean.EW_fl <- mean(EW_lib_sizes$FL); sd.EW_fl <- sd(EW_lib_sizes$FL)

ew.fl.lengths <- EW_lib_sizes$FL[match(sample.id[1:288],as.character(EW_lib_sizes$`Sample Description`))]
ew.fl.lengths[is.na(ew.fl.lengths)] <- mean.EW_fl
ew.fl.lengths <- as.integer(ew.fl.lengths)

fl.lengths <- c(ew.fl.lengths,lgep.fl.lengths)

# Add Family ID's
library(readxl)
All <- read_excel("/media/disk6/ARF/RNASEQ/Breeding-Value-Prediction/disk6directory/resources/phenos/All.xlsx")

# Subset EW fam ID's and get numeric sample id
ew.samp.fam.id <- unique(All[which(All$animal < 1000),c("sample","ID")])
samp.num <- (gsub("([0-9]+).*$", "\\1", unlist(ew.samp.fam.id[,1])))

# Check correct order of assignment
identical(as.vector(samp.num[match(sample.id[1:288],samp.num)]),sample.id[1:288])
ew.samp.fam.id <- ew.samp.fam.id$ID[match(sample.id[1:288],samp.num)]

lgep.samp.fam.id <- unique(All[which(All$animal > 1000),c("sample","ID")])
samp.num <- (gsub("([0-9]+).*$", "\\1", unlist(lgep.samp.fam.id[,1])))

identical(as.vector(samp.num[match(sample.id[-c(1:288)],samp.num)]),sample.id[-c(1:288)])
lgep.samp.fam.id <- lgep.samp.fam.id$ID[match(sample.id[-c(1:288)],samp.num)]

fam.id <- c(ew.samp.fam.id,lgep.samp.fam.id)
all.crosses <- grep("x", (fam.id),value=TRUE)
par1 <- fam.id
par2 <- rep("NA",length(par1))
lapply(strsplit(all.crosses,split = "x"),function(x) x[1])
par2[which(par1 %in% all.crosses)] <- unlist(lapply(strsplit(all.crosses,split = "x"),function(x) x[2]))
par1[which(par1 %in% all.crosses)] <- unlist(lapply(strsplit(all.crosses,split = "x"),function(x) x[1]))

## Create final data set #### 

exp.dat <- data.frame("animal_id"=animal.id,"fam_id"=fam.id,"p1"=par1,"p2"=par2, "batch"=batch.name,"lane"=lane.names,
                      "index_seq"=tech.rep.index[,1], "fl_mean"=fl.lengths,"index_let"=index.let,
                      "sample_id"=sample.id, "salmon_input_path"=tech.rep.index[,2])

# This data set is ready to be used for calling counts in Salmon ####
# Next add phenotypes and location information ####
# Write out families:
unique.fams <- unique(as.character(exp.dat$fam_id))
unique.fams <- unique(unlist(strsplit(x = unique.fams,split = "x")))
unique.fams <- data.frame("ID"=unique.fams)
# Accessed: tracker.cnr.ncsu.edu/tipdb/index.php
# >> Pedigree Information >> Pedigrees & Relatives
write_csv(unique.fams,col_names = T,path ="/media/disk6/ARF/RNASEQ/ufams.csv")
# Pedigree csv file was uploaded to tiproot and Ancestor info was downloaded:
pedigree <- read_excel("/media/disk6/ARF/RNASEQ/pedigree_2018_04_07_21_08_39.xlsx")

unique.crosses <- unique(as.character(exp.dat$fam_id))
unique.crosses <- grep("x", unique.crosses,value=TRUE)
unique.crosses <- strsplit(x = unique.crosses,split = "x")
unique.crosses <- data.frame(do.call(rbind,unique.crosses))
colnames(unique.crosses) <- c("parent1", "parent2")
unique.crosses <- unique.crosses[-c(3,4,6),]
write_delim(unique.crosses,col_names = T,append = F,delim = "\t",path ="/media/disk6/ARF/RNASEQ/ucrosses.csv")
# N22037	N22010 had to be manually reversed to be found in genetic_vals
# UC values do not exsit; Same for:
# N22029 B48:: B48 doesn't exist
# N10010 B21 :: B21 doesn't exist
# 
 breeding_values <- read_excel("/media/disk6/ARF/RNASEQ/breeding_values_2018_04_07_21_17_15.xlsx")
 genetic_values <- read_excel("/media/disk6/ARF/RNASEQ/genetic_values_2018_04_07_21_57_58.xlsx")
 
 phenotypes <- rbind(breeding_values[,-1],genetic_values[,-c(1,2)])
par1 <- c(breeding_values$Parent1,genetic_values$Parent1)
par2 <- c(rep(NA,nrow(breeding_values)),genetic_values$Parent2)

phenotypes <- cbind("par1"=par1,"par2"=par2,phenotypes)
phenotypes[58,1] <- "N22010"
phenotypes[58,2] <- "N22037"

missing.phenos <- which(exp.dat$p2 %in% c("UC","B21","B48"))

exp.dat$p_info = "A"
exp.dat$p_info[which(exp.dat$p2 %in% c("UC"))] <- "B"
exp.dat$p_info[which(exp.dat$p2 %in% c("B21","B48"))] <- "C"
exp.dat$p_info[which(is.na(exp.dat$Volume))] <- "D"
exp.dat$p2[missing.phenos] <- "NA"

match.phenos <- phenotypes[match(paste0(exp.dat$p1,exp.dat$p2),paste0(phenotypes$par1,phenotypes$par2)),]

exp.dat <- cbind(exp.dat,match.phenos[,-c(1,2)])

out.name <-  paste0("Sample_",exp.dat$sample_id,"_",exp.dat$lane)
 true.lane <- gsub(pattern = "ew_",replacement = "",x = exp.dat$lane[1:288])
exp.dat$salmon_output_path = paste0("media/disk6/ARF/RNASEQ/counts/86kSalmon/EW/", true.lane,"/",out.name)
exp.dat$salmon_output_path[-c(1:288)] = paste0("media/disk6/ARF/RNASEQ/counts/86kSalmon/LGEP/", true.lane[-c(1:288)],"/",out.name[-c(1:288)])
out.name <- paste0("Sample_",exp.dat$sample_id[each.file],"_",exp.dat$lane[each.file])
  true.lane <- gsub(pattern = "ew_",replacement = "",x = exp.dat$lane[each.file])
exp.dat.p <- exp.dat[-which(is.na(exp.dat$Volume)),]
