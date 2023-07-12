library(boot)
array2tibble  <- function(ndarray, columnnames, pagenames, pagecol) {
    # Turn a 3D array into a tibble. Array has to be organized such that
    # first dim is observations (don't have a name), second and third 
    # dimensions have a significance. In the end, the tibble will have 
    # pagenames in the page column pagecol. 
    # Input:
    #   ndarray: a 3d array
    #   columnnames: a vector of strings that match the second dim of ndarray
    #   pagenames: a vector of strings that match the third dim
    #   pagecol: the name of the column where pagenames are stored (string)
    # Output:
    #   nicetibble: a nice tibble

    npages  <- dim(ndarray)[3]
    nrows  <- dim(ndarray)[1]
    # Initialize tibble
    matx <- ndarray[, , 1]
    rownames(matx) <- rep(pagenames[1], nrows)
    colnames(matx)  <-  columnnames
    nicetibble  <- as_tibble(matx, rownames = pagecol)

    for (i in 2:npages) {
        matx <- ndarray[, , i]
        rownames(matx) <- rep(pagenames[i], nrows)
        colnames(matx)  <-  columnnames
        itibble  <- as_tibble(matx, rownames = pagecol)
        nicetibble <- full_join(nicetibble, itibble)
    }
    return(nicetibble)
}
temap2tibble <- function(ndarray, rois, pagenames, pagecol, value) {
    # Turn a 4D TE structure to a nice tibble
    # TE structure has to be rois x rois x manipulation x simulations
    # Input: 
    #   ndarray: 4D array (rois x rois x manipulation x simulations)
    #   rois: name of rois in 1st and 2nd dims (string vector)
    #   pagenames: names of the third dim (string vector)
    #   pagecol: the name of the column where pagenames are stored (string)
    #   value: name of the value stored (string) (for example, "TE")
    # Output:
    #   A not-so-nice tibble. For the purposes of the paper, I set Incoming and Outgoing for column names of TE.
    #   Felt too lazy to finesse it. R makes me tired. 
    nrois <- dim(ndarray)[1]
    tex <- apply(ndarray, c(1,2,3), mean) # Average across simulations
    nmanip  <- dim(tex)[3]
    # Initialize tibble
    matx <- tex[, , 1]
    rownames(matx)  <- rois
    colnames(matx)  <- rois
    nicedf  <- data.frame(Outgoing=rownames(matx)[row(matx)], Incoming=colnames(matx)[col(matx)], baddesign=c(matx))
    names(nicedf)[names(nicedf) == "baddesign"] <- value
    nicedf <- add_column(nicedf, baddesign=pagenames[1])
    names(nicedf)[names(nicedf) == "baddesign"]  <- pagecol

    nicetibble  <- as_tibble(nicedf)

    for (i in 2:nmanip) {
        matx <- tex[, , i]
        rownames(matx)  <- rois
        colnames(matx)  <- rois
        nicedf  <- data.frame(Outgoing=rownames(matx)[row(matx)], Incoming=colnames(matx)[col(matx)], baddesign=c(matx))
        names(nicedf)[names(nicedf) == "baddesign"] <- value
        nicedf <- add_column(nicedf, baddesign=pagenames[i])
        names(nicedf)[names(nicedf) == "baddesign"]  <- pagecol
        itibble  <- as_tibble(nicedf)
        nicetibble <- full_join(nicetibble, itibble)
    }
    return(nicetibble)
}
stat2str <- function(stattibble) {
    # Add a string column to a stattibble that contains the statistic that 
    # you probably want to report (for correlation only) like r = 0.46, p = 0.02
    ps <- stattibble$p.value * dim(stattibble)[1]
    pstring <- c()
    for (i in 1:length(ps)) {
        if (ps[i] < 0.001) {
            # pstring[i] <- sprintf(ps[i], fmt = "%#.1e")
            pstring[i] <- "p < 0.001"
        } else {
            pstring[i] <- paste0("p = ", sprintf(ps[i], fmt = "%#.3f"))
        }
    }
    eststring <- sprintf(stattibble$estimate, fmt = "%#.3f")
    statstring <- paste0(r"( $\rho = $)", eststring, ", ", pstring)
    stattibble2 <- add_column(stattibble, statstr = statstring)
    return(stattibble2)
}
stat2strwilcox <- function(stattibble) {
    # Add a string column to a stattibble that contains the statistic that 
    # you probably want to report (for correlation only) like r = 0.46, p = 0.02
    ps <- stattibble$p.value
    pstring <- c()
    for (i in 1:length(ps)) {
        if (ps[i] < 0.001) {
            # pstring[i] <- sprintf(ps[i], fmt = "%#.1e")
            pstring[i] <- "p < 0.001"
        } else {
            pstring[i] <- sprintf(ps[i], fmt = "%#.3f")
        }
    }
    eststring <- sprintf(stattibble$estimate, fmt = "%#.3f")
    statstring <- paste0("W = ", eststring, ", p = ", pstring)
    stattibble2 <- add_column(stattibble, statstr = statstring)
    return(stattibble2)
}
se <- function(x) {
        b <- boot(data = x, statistic = function(y, i) median(y[i]), R = 10000)
        ci <- boot.ci(b, type = "perc")
        sevals <- c(ci$percent[4], ci$percen[5])
        if (is.null(sevals)) {
            damntibble <- tibble(lowprc = median(x), highprc = median(x))
        } else {
        damntibble <- tibble(lowprc = sevals[1], highprc = sevals[2])
        }
    return(damntibble)
}
