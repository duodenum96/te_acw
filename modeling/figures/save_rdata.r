###### Save MATLAB style data to R tibbles for ggplot
rm(list = ls())
library(R.matlab)
library(broom)
library(ggthemes)
library(ggpubr)
library(corrplot)
library(patchwork)
library(effectsize)
library(multcomp)
library(coin)
library(tidyverse)
library(latex2exp)
# source("C:\\Users\\user\\Desktop\\brain_stuff\\ggplot_theme_Publication-master\\ggplot_theme_Publication-2.R")
setwd("/BICNAS2/ycatal/te_acw/modeling/src")
source("rhelpers.r")

setwd("/BICNAS2/ycatal/te_acw/modeling/.gitignore/teplain_results")
datasavelocation <- "/BICNAS2/ycatal/te_acw/modeling/figures/figuredata"

rois  <- c("V1", "V2", "V4", "DP", "MT", "8m", "5", "8l", "TEO", "2", "F1", "STPc", "7A", 
    "46d", "10", "9/46v", "9/46d", "F5", "TEpd", "PBr", "7m", "7B", "F2", 
    "STPi", "ProM", "F7", "8B", "STPr", "24c")
visual  <-  c(1, 2, 3, 5, 4, 9, 6, 19, 13, 12)
irois  <- rois[visual]

# Do batch stuff
netws <- c("w_EE", "w_IE", "w_EI", "w_II", "beta_E", "beta_I", "mu_EE", "mu_IE", "tau_E", "tau_I", "eta")
manips <- list(c(0.5, 0.6, 0.7, 0.8, 0.9, 1), 
               c(1, 1.1, 1.2, 1.3, 1.4, 1.5),
               c(1, 1.1, 1.2),
               c(0.8, 0.9, 1),
               c(0.5, 0.6, 0.7, 0.8, 0.9, 1),
               c(1, 1.2, 1.4, 1.6, 1.8, 2),
               c(0.9, 0.95, 1, 1.01),
               c(0.99, 1, 1.05, 1.1, 1.2),
               c(0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2),
               c(0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2),
               c(0.8, 0.82, 0.84, 0.86, 0.88, 0.9, 0.92, 0.94, 0.96, 0.98, 1))
defvals <- c(24.3, 12.2, 19.7, 12.5, 0.06, 0.3510, 33.7, 25.3, 20, 10, 0.68)

nnetw <- length(netws)
for (zaaxd in 1:nnetw) {
    # General function for LaTeX
    appender <- function(string) {
        TeX(paste0(string))
    }
    # Import data
    netw  <- netws[zaaxd]
    netwsave <- netw # Add variable for figure saving, netw will be changed to LaTeX
    # File name to be saved
    savename <- paste0(datasavelocation, "\\", netwsave, "_rformat.RData")

    tes  <- readMat(paste(netw, "netw_task_te_plain.mat", sep = ""))
    te  <- tes$te
    acws_rest  <- readMat(paste("acws_", netw, "netw_rest.mat", sep = ""))
    acw_rest  <- acws_rest$acw0

    # Make the name LaTeX style
    if (netw != "eta") {
        netw <- paste0(r"( $\ )", netw)
        netw_X <- strsplit(netw, "_")
        netw_X[[1]][2] <- paste0("{", netw_X[[1]][2], "}")
        netw <- paste0(netw_X[[1]][1], "_", netw_X[[1]][2], "$")
    } else {
        netw <- paste0(r"( $\ )", netw, "$")
    }

    manip  <- manips[[zaaxd]]
    defval <- defvals[zaaxd]

    dim3names  <- paste0(netw, " = ", as.character(manip * defval))
    defi <- manip == 1
    dim3names[defi] <- paste0("## ", netw, " = ", as.character(manip[defi] * defval), " ##")
    dim3names <- factor(dim3names, level = dim3names)

    nsim  <- dim(te)[3]
    nroi  <- dim(te)[1]
    nmanip  <- dim(te)[4]

    # ACWs : First dim: ROIs, second dim: simulations, 3rd dim: manipulations
    # TEs: roi x roi x sim x manipulation

    # %% For TEs, calculate node degree
    nodedegree  <- array(data = NA, dim = c(nroi, nsim, nmanip))
    for (m in 1:nmanip) {
        for (i in 1:nsim) {
            for (j in 1:nroi) {
                nodedegree[j, i, m]  <- sum(te[j, , i, m]) + sum(te[ , j, i, m])
            }
        }
    }
    
    # %% Turn arrays into a tibble
    # Now, TEs are roi x sim x manipulation
    # First, convert both ACWs and TE calculations to sim x roi x man
    nodedegreep  <- aperm(nodedegree, c(2, 1, 3)) # Node degree permuted
    acw_restp  <- aperm(acw_rest, c(2, 1, 3))
    # Node Degree tibble
    ndt  <- nodedegreep %>% 
        array2tibble(columnnames = irois, pagenames = dim3names, pagecol = "Manipulation") %>%
        pivot_longer(cols = !Manipulation, names_to = "ROI", values_to = "Degree")
    
    ndt <- ndt %>% 
        add_column(ROIf = as.numeric(factor(ndt$ROI, levels = irois))) %>%  # Add factorized column
        group_by(Manipulation)


    acwrt <- acw_restp %>% 
        array2tibble(columnnames = irois, pagenames = dim3names, pagecol = "Manipulation")  %>% 
        pivot_longer(cols = !Manipulation, names_to = "ROI", values_to = "ACW_0")

    acwrt <- acwrt %>% add_column(ROIf = as.numeric(factor(acwrt$ROI, levels = irois))) %>%  # Add factorized column
        group_by(Manipulation)

    # Turn 4d TE array to 3D by averaging
    tex <- aperm(te, c(1,2,4,3))
    # Then turn it into a tibble, this is for the heatmap
    tet <- temap2tibble(tex, irois, dim3names, "Manipulation", "TE")

    # Get summary statistics: median and standard error
    tesummarytibble <- ndt %>% 
        group_by(Manipulation, ROI, ROIf) %>%
        summarise(median = median(Degree), se = se(Degree))
    
    acwsummarytibble <- acwrt %>% 
        group_by(Manipulation, ROI, ROIf) %>% 
        summarise(ACW_0_median = median(ACW_0), ACW_0_se = se(ACW_0))

    # Tibble that contains both
    acwandte <- left_join(acwsummarytibble, tesummarytibble) # Note that TE is encoded as median and se, not TE_median

    # Nest for correlation (figure 1)
    nestedacwrt <- acwrt %>% select(Manipulation, ACW_0, ROIf) %>%  nest(data = c("ACW_0", "ROIf"))
    acwcorrres <- nestedacwrt %>% mutate(
        test = map(data, ~cor.test(.x$ROIf, .x$ACW_0, method = "spearman")),
        tidied = map(test, tidy)) %>% 
        unnest(cols = tidied) %>% 
        select(-data, -test) %>% 
        stat2str()
    
    # Get the statistics str
    statstr1 <- c()
    statstr2 <- c()
    for (i in 1:length(acwcorrres$statstr)) {
        statstr1[i] <- strsplit(acwcorrres$statstr[i], ", ")[[1]][1]
        statstr2[i] <- strsplit(acwcorrres$statstr[i], ", ")[[1]][2]
    }
    acwcorrres <- add_column(acwcorrres, statstr1 = statstr1)
    acwcorrres <- add_column(acwcorrres, statstr2 = statstr2)

    # For figure 1, use acwrt and acwcorrres
    # For heatmap: use tet

    # Do the stats for TE (rank of ROI vs TE)
    idata <- group_by(ndt, Manipulation)
    nestedidata <- idata %>% select(Manipulation, Degree, ROIf) %>% nest(data = c("Degree", "ROIf"))
    tecorrres <- nestedidata %>% mutate(
        test = map(data, ~cor.test(.x$ROIf, .x$Degree, method = "spearman")),
        tidied = map(test, tidy)) %>% 
        unnest(cols = tidied) %>% 
        select(-data, -test) %>% 
        stat2str()
    
    # Again, extract the strings for various correlation results
    statstr1 <- c()
    statstr2 <- c()
    for (i in 1:length(tecorrres$statstr)) {
        statstr1[i] <- strsplit(tecorrres$statstr[i], ", ")[[1]][1]
        statstr2[i] <- strsplit(tecorrres$statstr[i], ", ")[[1]][2]
    }
    tecorrres <- add_column(tecorrres, statstr1 = statstr1)
    tecorrres <- add_column(tecorrres, statstr2 = statstr2)

    # For violinplot of TE, use idata, tecorres

    ######## ACW and TE correlation #############
    nestedacwte <- acwandte %>% 
        group_by(Manipulation) %>% 
        select(Manipulation, ACW_0_median, median) %>% 
        nest(data = c("ACW_0_median", "median"))

    acwtecorrres <- nestedacwte %>% mutate(
        test = map(data, ~cor.test(.x$ACW_0_median, .x$median, method = "spearman")),
        tidied = map(test, tidy)) %>% 
        unnest(cols = tidied) %>% 
        select(-data, -test) %>% stat2str()

    statstr1 <- c()
    statstr2 <- c()
    for (i in 1:length(acwtecorrres$statstr)) {
        statstr1[i] <- strsplit(acwtecorrres$statstr[i], ", ")[[1]][1]
        statstr2[i] <- strsplit(acwtecorrres$statstr[i], ", ")[[1]][2]
    }
    acwtecorrres <- add_column(acwtecorrres, statstr1 = statstr1)
    acwtecorrres <- add_column(acwtecorrres, statstr2 = statstr2)

    # For acw-te correlation: use acwandte and acwtecorres

    ####### Hier vs No Hier
    hierindices <- (acwcorrres$p.value * dim(acwcorrres)[1]) < 0.05 & acwcorrres$estimate > 0
    if (length(unique(hierindices)) > 1) {
        acwcorrres <- add_column(acwcorrres, hierindices)
        hiertibble <- right_join(acwcorrres, acwrt)
        fulldata <- add_column(hiertibble, Degree = idata$Degree)
        fulldata <- add_column(fulldata, hierindicesnum = factor(as.numeric(fulldata$hierindices)))

        fulldata <- add_column(fulldata, ROIf2 = factor(fulldata$ROIf))

        model <- anova(lm(Degree ~ ROIf2 + hierindicesnum + ROIf2 : hierindicesnum, data = fulldata)) # 2-way ANOVA
        eta2 <- sprintf(eta_squared(model)$Eta2_partial[2], fmt = "%#.3f")
        Fval <- sprintf(model$"F value"[2], fmt = "%#.3f")
        df1 <- model$Df[2]
        df_res <- model$Df[4]
        p_anova <- model$"Pr(>F)"[2]
        if (p_anova < 0.001) {p_anova <- "p < 0.001"} else {p_anova <- paste0("p = ", sprintf(p_anova, fmt = "%#.3f"))}
        anovareport <- paste0("F(", df1, ", ", df_res, ") = ", Fval, ", ", 
        p_anova, ", ", 
        r"( $\eta^2$ = )", eta2)

        ws <- c()
        ps <- c()
        zs <- c()
        pstring <- c()
        # Get Wilcoxon results
        for (i in 1:length(unique(fulldata$ROIf))) {
            itibble <- filter(fulldata, ROIf == i)
            res <- wilcox.test(filter(itibble, hierindicesnum == 0)$Degree, filter(itibble, hierindicesnum == 1)$Degree, 
            conf.int = TRUE)
            
            ws[i] <- res$statistic
            ps[i] <- res$p.value * length(unique(fulldata$ROIf))
            zs[i] <- res$estimate

            if (ps[i] < 0.001) {
                pstring[i] <- "p < 0.001"
            } else {
                pstring[i] <- paste0("p = ", sprintf(ps[i], fmt = "%#.3f"))
            }
        }

        compstr <- paste0("W = ", ws, ", ", pstring, ", dil = ", sprintf(zs, fmt = "%#.3f"))

        comptibble <- tibble(label = paste(compstr, collapse = ""))
        
        # Hier vs No Hier
        fulldata_medians <- fulldata %>% group_by(ROIf, hierindicesnum) %>% summarise(Medians = median(Degree), se = se(Degree))
        fulldata_medians <- add_column(fulldata_medians, ROI = factor(irois[fulldata_medians$ROIf], levels = irois))
    } else {
        fulldata_medians <- c()
        anovareport <- c()
        compstr <- c()
    }
    ### For hier plot, use fulldata_medians, anovareport, compstr


    ################## SAVE THE RESULTS #####################
    # Recap on what to be used:
    # For figure 1, use acwrt and acwcorrres
    # For heatmap: use tet
    # For acw-te correlation: use acwandte and acwtecorres
    # For hier plot, use fulldata_medians, anovareport, compstr
    # In general, use netwsave, irois and hierindices for checking if there is nohier
    
    save(acwrt, acwcorrres, tet, acwandte, acwtecorrres, fulldata_medians, anovareport, compstr, hierindices, 
    file = savename)
}