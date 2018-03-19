#Compare asreml and sommer normalization methods:

####
load("/mnt/ARF/step2.normalization/asreml/asreml.counts.RData")
load("/mnt/ARF/step2.normalization/asreml/192.expt.data.RData")
load("/mnt/ARF/step2.normalization/sommer/sommer.counts.RData")
####

###
rownames(asreml.counts)
rownames(sommer.counts)
expt.dat.192$animal.id
###

# What is the correlation between all of the counts?
library(parallel)
salmon.vs.asreml.logcor <- unlist(mclapply(1:ncol(asreml.counts),function(each.count){
  cor(asreml.counts[,each.count],sommer.counts[,each.count])
},mc.cores=25))
hist(salmon.vs.asreml.logcor)
plot(salmon.vs.asreml.logcor)

salmon.vs.asreml.expcor <- unlist(mclapply(1:ncol(asreml.counts),function(each.count){
  cor(2^asreml.counts[,each.count],2^sommer.counts[,each.count])
},mc.cores=25))
hist(salmon.vs.asreml.expcor)
plot(salmon.vs.asreml.expcor)

plot(salmon.vs.asreml.expcor,salmon.vs.asreml.logcor)
t.test(salmon.vs.asreml.expcor,salmon.vs.asreml.logcor)

# What is the distribution of counts?
library(fitdistrplus)
descdist(2^asreml.counts[,2], discrete = TRUE,method = "sample")
fit.weibull <- fitdist(2^asreml.counts[,1], "weibull")
plot(fit.weibull)
library(logspline)

