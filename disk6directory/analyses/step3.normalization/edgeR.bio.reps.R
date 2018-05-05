library(edgeR)
edgeRUsersGuide()

load("/media/disk6/ARF/RNASEQ/counts/86kSalmon/Step2_Load_Counts.RData")
names(load.counts)

cts <- load.counts$txi_object$counts
normMat <- load.counts$txi_object$length
normMat <- normMat/exp(rowMeans(log(normMat)))
library(edgeR)
o <- log(calcNormFactors(cts/normMat)) + log(colSums(cts/normMat))
y <- DGEList(cts)
y$offset <- t(t(log(normMat)) + o)
phenos <- load.counts$phenos
design <- model.matrix(~phenos$batch + phenos$index_seq + phenos$fam_id)
y <- estimateDisp(y,design)
