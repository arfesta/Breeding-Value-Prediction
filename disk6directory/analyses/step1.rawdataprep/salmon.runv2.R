#RUN SALMON ####
load("/media/disk6/ARF/RNASEQ/Breeding-Value-Prediction/disk6directory/resources/expt.dat.720.RData")
setwd("/media/disk6/ARF/RNASEQ")
ew.animals <- unique(as.character(expt.dat.720$animal_id[which(expt.dat.720$batch == "EW")]))
for(each.file in 2:96){
 select.rows <-  which(as.character(expt.dat.720$animal_id) %in% ew.animals[each.file])

  file.cat.names <- (unlist( lapply(1:length(select.rows),function(x) { paste0("/media/disk6/ARF/RNASEQ/trimmedfiltreads/86k/EW/",expt.dat.720$salmon_input_path[select.rows[x]])}) ))
  
  file.cat.names <- paste(file.cat.names[1],file.cat.names[2],file.cat.names[3],sep = " ")
  
  out.name <- paste0("Sample_",ew.animals[each.file])
 FL.value <- mean(expt.dat.720$fl_mean[select.rows])
  cmd <- paste0("./software/Salmon/bin/salmon quant --seqBias -r ",file.cat.names, " -i ./references/Pita.86k -o ./counts/86kSalmon/bio_EW/",out.name," -p 50 -l U --numBootstraps 10 --fldMean ", FL.value," --fldSD 21")
  system(cmd)
}  

each.file=28
lgep.animals <- unique(as.character(expt.dat.720$animal_id[which(expt.dat.720$batch == "LGEP")]))
for(each.file in 1:144){
 select.rows <-  which(as.character(expt.dat.720$animal_id) %in% lgep.animals[each.file])

  file.cat.names <- (unlist( lapply(1:length(select.rows),function(x) { paste0("/media/disk6/ARF/RNASEQ/trimmedfiltreads/86k/LGEP/",expt.dat.720$salmon_input_path[select.rows[x]])}) ))
  
  file.cat.names <- paste(file.cat.names[1],file.cat.names[2],file.cat.names[3],sep = " ")
  
  out.name <- paste0("Sample_",lgep.animals[each.file])
 FL.value <- mean(expt.dat.720$fl_mean[select.rows])
  cmd <- paste0("./software/Salmon/bin/salmon quant --seqBias -r ",file.cat.names, " -i ./references/Pita.86k -o ./counts/86kSalmon/bio_LGEP/",out.name," -p 55 -l U --numBootstraps 10 --fldMean ", FL.value," --fldSD 21")
  system(cmd)
}  

make.int <- as.character(c(54,19,33,35,36,41,47,114,116))
for(each.file in 2:length(make.int)){
 select.rows <-  which(as.character(expt.dat.720$animal_id) %in% make.int[each.file])

  file.cat.names <- (unlist( lapply(1:length(select.rows),function(x) { paste0("/media/disk6/ARF/RNASEQ/trimmedfiltreads/86k/LGEP/",expt.dat.720$salmon_input_path[select.rows[x]])}) ))
  
  file.cat.names <- paste(file.cat.names[1],file.cat.names[2],file.cat.names[3],sep = " ")
  
  out.name <- paste0("Sample_",make.int[each.file])
 FL.value <- as.integer(mean(expt.dat.720$fl_mean[select.rows]))
  cmd <- paste0("./software/Salmon/bin/salmon quant --seqBias -r ",file.cat.names, " -i ./references/Pita.86k -o ./counts/86kSalmon/bio_LGEP/",out.name," -p 55 -l U --numBootstraps 10 --fldMean ", FL.value," --fldSD 21")
  system(cmd)
}  
