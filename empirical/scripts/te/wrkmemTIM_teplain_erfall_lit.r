#!/home/ycatal/bin/ Rscript
# Calculate transfer entropy with given data from matlab
rm(list = ls())
library(future)
library(RTransferEntropy)
library(R.matlab)
plan(multicore)

setwd("/BICNAS2/ycatal/sourcereconstruction/data")

task <- "Wrkmem"
taskext <- "TIM_noprestim"
subj <- readMat("Wrkmem_subjs.mat")
subjs <- subj$subjs
nsubj <- dim(subjs)[1]

for (i in 1:nsubj) {
    isubj  <- subjs[i]
    # import data
    tryCatch( {
    datafolder  <- paste("/BICNAS2/ycatal/sourceparcellated/", isubj, "/", task, "/icamne/source", sep = "")
    setwd(datafolder) 
    }, 
    error = function(e){
        datafolder  <- paste("/BICNAS2/ycatal/sourceparcellated/", isubj, "/", task, "/source", sep = "")
        setwd(datafolder)
    } )
    qfile  <- paste(isubj, "_MEG", "-", task, "_", taskext, "_qs_erfall_lit.mat", sep = "")
    qdata  <- readMat(qfile)
    datamatrix  <- qdata$datamatrix
    q_xx  <- qdata$q.xx
    q_xy  <- qdata$q.xy
    q_yx  <- qdata$q.yx
    q_yy  <- qdata$q.yy
    nroi  <- dim(datamatrix)[1]

    # initialize te matrix
    te  <- array(data = 0, dim = c(nroi, nroi))
    for (p in 1:nroi) {
        for (q in (p + 1):nroi) {
            if (q <= nroi){
            x  <- datamatrix[p, ]
            y  <- datamatrix[q, ]
            lx_1  <- q_xy[p, q] # From x to y
            ly_1  <- q_yy[p, q] # From y to y
            lx_2  <- q_yx[p, q] # From y to x
            ly_2  <- q_xx[p, q] # From x to x

            # Do x to y
            te[p, q] <- calc_te(x, y, lx = lx_1, ly = ly_1)
            
            te[q, p] <- calc_te(y, x, lx = lx_2, ly = ly_2)
            }
        }
    }
    tesavename  <- paste(isubj, "_MEG", "-", task, "_", taskext, "_teplain_erfall_lit.mat", sep = "")
    writeMat(tesavename, te = te)
    print(paste("Finished: Subject ", i, sep = ""))
}

print("Everything finished, yayyyyyyyy!!!!!!!!!!!!!!!!!!!!!!!!!!!!")