load("~/Desktop/prediction.paper/3.4/ht.snp.and.de.or.txtrr.comp.RDA.RData")
load("~/Desktop/prediction.paper/3.4/volume.snp.and.de.or.txtrr.comp.RDA.RData")

se.98 <- seq(.01,.98,.01)
se.99 <- seq(.01,.99,.01)

pdf(file = "~/Desktop/snp.vs.txt.pdf",width = 4,height = 4)
plot(se.99,cor.list,xlab = "", ylab="", cex.axis=.7, cex.lab=1, cex=.5)
title(xlab = "Weight of snps vs. DE transcripts", ylab="R^2", mgp=c(2,1,0),cex.lab=.8)
plot(se.98,cor.list2,xlab = "", ylab="", cex.axis=.7, cex.lab=1, cex=.5)
title(xlab = "Weight of snps vs. LOO txpts", ylab="R^2", mgp=c(2,1,0),cex.lab=.8)
plot(se.99,cor.list3,xlab = "", ylab="", cex.axis=.7, cex.lab=1, cex=.5)
title(xlab = "Weight of snps vs. DE txpts", ylab="R^2", mgp=c(2,1,0),cex.lab=.8)
plot(se.98,cor.list4,xlab = "", ylab="", cex.axis=.7, cex.lab=1, cex=.5)
title(xlab = "Weight of snps vs. LOO txpts", ylab="R^2", mgp=c(2,1,0),cex.lab=.8)
dev.off()

plot(pred,true,xlab = "", ylab="", cex.axis=.8, cex.lab=1.1, cex=.5,
     ylim = c(90,145), xlim=c(95,140))
title(xlab = "prediction", ylab="true value", mgp=c(2,1,0),cex.lab=1)
lml <- summary(lm(pred~true))
legend(x='topleft', legend=paste("R^2 =",round(lml$r.squared,2),sep=""),bty = "n",cex = .75)
dev.off()