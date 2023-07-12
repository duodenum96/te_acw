# Calculate TE from data from matlab
rm(list = ls())
library(future)
library(RTransferEntropy)
library(R.matlab)
plan(multicore)

setwd("/BICNAS2/ycatal/te_acw/modeling/.gitignore/results")
task <- "defaultnetw"
tesavename <- paste(task, "_te_plain.mat", sep = "")

qfile <- paste("greedy_", task, ".mat", sep = "")
qdata <- readMat(qfile)
datamatrix <- qdata$datamatrix
q_xx <- qdata$q.xx
q_xy <- qdata$q.xy
q_yx <- qdata$q.yx
q_yy <- qdata$q.yy

nsim <- dim(datamatrix)[3]
nroi <- dim(datamatrix)[1]

te  <- array(data = 0, dim = c(nroi, nroi, nsim))
for (i in 1:nsim) {
    for (p in 1:nroi) {
        for (q in (p+1):nroi) {
            if (q <= nroi) {
                x <- datamatrix[p, , i]
                y <- datamatrix[q, , i]
                lx_1  <- q_xy[p, q, i] # From x to y
                ly_1  <- q_yy[p, q, i] # From y to y
                lx_2  <- q_yx[p, q, i] # From y to x
                ly_2  <- q_xx[p, q, i] # From x to x

                te[p, q, i] <- calc_te(x, y, lx = lx_1, ly = ly_1)

                # Do y to x
                te[q, p, i] <- calc_te(y, x, lx = lx_2, ly = ly_2)
            }
        }
    }
    print(paste("Simulation ", i, " finished!!", sep = ""))
}
writeMat(tesavename, te = te)
print("Everything finished, yayyyyyyyy!!!!!!!!!!!!!!!!!!!!!!!!!!!!")