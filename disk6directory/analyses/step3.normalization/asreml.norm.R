load(file="/media/disk6/ARF/RNASEQ/counts/86kSalmon/Step2_Load_Counts.RData")
u.fams <- unique(c(as.character(load.counts$phenos$p1),
                   as.character(load.counts$phenos$p2)))
unique.fams <- unique(as.character(load.counts$phenos$fam_id))
unique.fams <- unique(unlist(strsplit(x = unique.fams,split = "x")))
unique.fams <- data.frame("ID"=unique.fams)
write_csv(unique.fams,col_names = T,path ="/media/disk6/ARF/RNASEQ/Breeding-Value-Prediction/disk6directory/resources/phenos/ufams.v2.csv")
##Accessed tip root
##pedigree infomration --> pedigree and relatives
##  upload csv written up and generate ancestors
##  download file and put in /resources/pedigree/
##  Open file and delete the comment column
##  Load pedigree
pedigree <- readxl::read_xlsx("/media/disk6/ARF/RNASEQ/Breeding-Value-Prediction/disk6directory/resources/pedigree/pedigree_2018_04_19_20_36_15.xlsx")
load.counts$phenos$fam_id
as.character(unique.fams[,1])[which(!as.character(unique.fams[,1]) %in% as.character(pedigree$Id))]
