#RUN SALMON ####
load("/media/disk6/ARF/RNASEQ/Breeding-Value-Prediction/disk6directory/resources/expt.dat.720.RData")
setwd("/media/disk6/ARF/RNASEQ")

for(each.file in 1:288){
  file <- paste0("/media/disk6/ARF/RNASEQ/trimmedfiltreads/86k/EW/",exp.dat$salmon_input_path[each.file])
  out.name <- paste0("Sample_",exp.dat$sample_id[each.file],"_",exp.dat$lane[each.file])
  true.lane <- gsub(pattern = "ew_",replacement = "",x = exp.dat$lane[each.file])
  cmd <- paste0("./software/Salmon/bin/salmon quant --seqBias -r ",file, " -i ./references/Pita.86k -o ./counts/86kSalmon/EW/", true.lane,"/",out.name," -p 40 -l U --numBootstraps 10 --fldMean ", exp.dat$fl_mean[each.file]," --fldSD 21")
  system(cmd)
}  


for(each.file in 289:720){
  file <- paste0("/media/disk6/ARF/RNASEQ/trimmedfiltreads/86k/LGEP/",exp.dat$salmon_input_path[each.file])
  out.name <- paste0("Sample_",exp.dat$sample_id[each.file],"_",exp.dat$lane[each.file])
  cmd <- paste0("./software/Salmon/bin/salmon quant --seqBias -r ",file, " -i ./references/Pita.86k -o ./counts/86kSalmon/LGEP/", exp.dat$lane[each.file],"/",out.name," -p 40 -l U --numBootstraps 10 --fldMean ", exp.dat$fl_mean[each.file]," --fldSD 21")
  system(cmd)
}  
