?glmnet
## Standardize variables: (need to use n instead of (n-1) as denominator)
mysd <- function(y) sqrt(sum((y-mean(y))^2)/length(y))
sx <- scale(log2(tmp.dat.fam[-1,]), scale = apply(log2(tmp.dat.fam[-1,]), 2, mysd))
sx <- as.matrix(sx, ncol = ncol(tmp.dat.fam), nrow = 56)
y <- fam.phenos$vol[2:57]
sy <- as.vector(scale(y, scale = mysd(y)))

n=56
## Calculate lambda path (first get lambda_max):
lambda_max <- max(abs(colSums(sx*y)))/(n)
epsilon <- .01
K <- 100
lambdapath <- round(exp(seq(log(lambda_max), log(lambda_max*epsilon), 
                            length.out = K)), digits = 10)

plot(lambdapath[1:92],gl.mod$lambda)


plot(gl.mod$dev.ratio)
