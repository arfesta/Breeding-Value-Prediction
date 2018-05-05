library(edgeR)
edgeRUsersGuide()

load("/media/disk6/ARF/RNASEQ/counts/86kSalmon/Step2_TECH_Load_Counts.RData")
names(load.counts_tech)
new.rnames <- paste0(load.counts_tech$phenos$animal_id,".",load.counts_tech$phenos$lane,".",load.counts_tech$phenos$index_let)
rownames(load.counts_tech$phenos) <- new.rnames

cts <- load.counts_tech$txi_object$counts
normMat <- load.counts_tech$txi_object$length
normMat <- normMat/exp(rowMeans(log(normMat)))
library(edgeR)
o <- log(calcNormFactors(cts/normMat)) + log(colSums(cts/normMat))
y <- DGEList(cts)
y$offset <- t(t(log(normMat)) + o)
phenos <- load.counts_tech$phenos
design <- model.matrix(~phenos$batch + phenos$lane + phenos$index_seq + phenos$animal_id)
y <- estimateDisp(y,design)
