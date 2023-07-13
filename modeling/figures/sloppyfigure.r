# %% Do a pipeline to plot all relevant results from simulation manipulations
# Don't forget to add default value of network parameter

rm(list = ls())
library(tidyverse)
library(R.matlab)
library(broom)
library(ggthemes)
library(ggpubr)
library(corrplot)
library(patchwork)
library(rstatix)
library(latex2exp)
setwd("C:\\Users\\user\\Desktop\\brain_stuff\\te_acw_simulations_v2\\src")
source("rhelpers.r")
source("C:\\Users\\user\\Desktop\\brain_stuff\\ggplot_theme_Publication-master\\ggplot_theme_Publication-2.R")

setwd("C:\\Users\\user\\Desktop\\brain_stuff\\te_acw_simulations_v2\\.gitignore\\teplain_results")
textsize  <- 30

rois  <- c("V1", "V2", "V4", "DP", "MT", "8m", "5", "8l", "TEO", "2", "F1", "STPc", "7A", 
    "46d", "10", "9/46v", "9/46d", "F5", "TEpd", "PBr", "7m", "7B", "F2", 
    "STPi", "ProM", "F7", "8B", "STPr", "24c")
visual  <-  c(1, 2, 3, 5, 4, 9, 6, 19, 13, 12)
irois  <- rois[visual]

# Do batch stuff
netws <- c("w_EE", "w_IE", "w_EI", "w_II", "beta_E", "beta_I", "mu_EE", "mu_IE", "tau_E", "tau_I", "eta")
manips <- list(c(0.5, 0.6, 0.7, 0.8, 0.9, 1),           # 1
               c(1, 1.1, 1.2, 1.3, 1.4, 1.5),           # 2
               c(1, 1.1, 1.2),                          # 3  
               c(0.8, 0.9, 1),                          # 4  
               c(0.5, 0.6, 0.7, 0.8, 0.9, 1),           # 5
               c(1, 1.2, 1.4, 1.6, 1.8, 2),             # 6
               c(0.9, 0.95, 1, 1.01),                   # 7
               c(0.99, 1, 1.05, 1.1, 1.2),              # 8
               c(0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2),   # 9 
               c(0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2),   # 10
               c(0.8, 0.82, 0.84, 0.86, 0.88, 0.9, 0.92, 0.94, 0.96, 0.98, 1)) # 11
defvals <- c(24.3, 12.2, 19.7, 12.5, 0.06, 0.3510, 33.7, 25.3, 20, 10, 0.68)

manipdummy <- c(rep("w_EE", 6), rep("w_IE", 5), rep("w_EI", 2), rep("w_II", 2), rep("beta_E", 5),
                rep("beta_I", 5), rep("mu_EE", 3), rep("mu_IE", 4), rep("tau_E", 7), rep("tau_I", 7), rep("eta", 10))

nnetw <- length(netws)


acwcorresults <- list()
teresults <- list()
teslope <- list()
# Start batch figure generation
## These three lines change in batch figure generation
for (zaaxd in 1:nnetw) {
    netw  <- netws[zaaxd]

    manip  <- manips[[zaaxd]]
    defval <- defvals[zaaxd]

    dim3names  <- paste(netw, " = ", as.character(manip * defval, sep = ""))
    dim3names <- factor(dim3names, level = dim3names)

    tes  <- readMat(paste(netw, "netw_task_te_plain.mat", sep = ""))
    te  <- tes$te
    acws_rest  <- readMat(paste("acws_", netw, "netw_rest.mat", sep = ""))
    acw_rest  <- acws_rest$acw0

    nsim  <- dim(te)[3]
    nroi  <- dim(te)[1]
    nmanip  <- dim(te)[4]

    # ACWs : First dim: ROIs, second dim: simulations, 3rd dim: manipulations
    # TEs: roi x roi x sim x manipulation

    # %% For TEs, calculate node degree, incoming TE, outgoing TE per ROI
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
    nodedegreep  <- aperm(nodedegree, c(2, 1, 3))
    acw_restp  <- aperm(acw_rest, c(2, 1, 3))
    

    ndt  <- nodedegreep %>% array2tibble(columnnames = irois, pagenames = dim3names, pagecol = "Manipulation") %>%
        pivot_longer(cols = !Manipulation, names_to = "ROI", values_to = "Degree")

    acwrt <- acw_restp %>% array2tibble(columnnames = irois, pagenames = dim3names, pagecol = "Manipulation")  %>% 
        pivot_longer(cols = !Manipulation, names_to = "ROI", values_to = "ACW_0")

    # acwallt  <- add_column(acwrt, ACW_0_Task = acwtt$ACW_0)
    # Turn 4d TE array to 3D by averaging
    tex <- aperm(te, c(1,2,4,3))
    # Then turn it into a tibble
    tet <- temap2tibble(tex, irois, dim3names, "Manipulation", "TE")

    summarytibble <- ndt %>% group_by(Manipulation, ROI) %>%
        summarise(degree_median = median(Degree))

    acwsummary <- acwrt %>% group_by(Manipulation, ROI) %>% 
        summarise(acw_median = median(ACW_0))

    tesummary <- summarytibble

    acwandte <- full_join(acwsummary, tesummary)

    # Drop c = 1
    if ((zaaxd == 1)) {
        acwtibble <- acwrt
        tetibble <- ndt
    } else {
        acwtibble <- filter(acwrt, Manipulation != (dim3names[manip == 1]))
        tetibble <- filter(ndt, Manipulation != (dim3names[manip == 1]))
    }

    remmanip <- unique(acwtibble$Manipulation) # Remaining manipulations after dropping default

    acwcorresults[[zaaxd]] <- vector(mode = "numeric", length = length(remmanip))
    teresults[[zaaxd]] <- vector(mode = "numeric", length = length(remmanip))
    teslope[[zaaxd]] <- vector(mode = "numeric", length = length(remmanip))

    # Calculate correlation coefficients
    for (i in 1:(nmanip - 1)) {
        imanip <- filter(acwtibble, Manipulation == remmanip[i])
        imanip <- add_column(imanip, ROIf = as.numeric(factor(imanip$ROI, levels = irois)))
        linfit <- lm(imanip$ACW_0 ~ imanip$ROIf)[[1]][[2]] # Extract the slope
        acwcorresults[[zaaxd]][i] <- linfit

        imanipte <- filter(tetibble, Manipulation == remmanip[i])
        imanipte <- add_column(imanipte, ROIf = as.numeric(factor(imanip$ROI, levels = irois)))
        teresults[[zaaxd]][i] <- sum(imanipte$Degree)
        teslope[[zaaxd]][i] <- lm(imanipte$Degree ~ imanipte$ROIf)[[1]][2]
    }
}
finalte <- unlist(teresults)
finalacw <- unlist(acwcorresults)
finalteslope <- unlist(teslope)

manipdummy <- paste0(r"( $\ )", manipdummy)
netw_X <- strsplit(manipdummy, "_")
manipdummy_tex <- c()
for (i in 1:length(netw_X)) {
    if (length(netw_X[[i]]) == 2) {
        netw_X[[i]][2] <- paste0("{", netw_X[[i]][2], "}")
        manipdummy_tex[i] <- paste0(netw_X[[i]][1], "_", netw_X[[i]][2], "$")
    } else {
        manipdummy_tex[i] <- paste0(netw_X[[i]][1], "$")
    }
}
appender <- function(string) {
    TeX(paste0(string))
}

teacwtbl <- tibble(te=finalte, acw=finalacw, teslope=finalteslope, manip = manipdummy_tex)
ggplot(teacwtbl, mapping = aes(x = acw, y = te)) + 
    geom_point(size = 7, aes(x = acw, y = te, color = manip)) + 
    stat_smooth(method = "lm", color = "black") + 
    stat_cor(method = "spearman", 
    label.sep = "\n", size = 12, cor.coef.name = "rho",
    p.accuracy = 0.001, label.x.npc = 0.7) +
    labs(x = "Slope of linear fit to ROIs vs ACW", y = "Total TE across ROIs", color="") + 
    scale_color_discrete(labels = TeX(manipdummy_tex)) + 
    theme_Publication() + theme(aspect.ratio = 1, 
    text = element_text(size = textsize))
savename <- "figures/simulation_sloppy.png"
ggsave(savename)

ggplot(teacwtbl, aes(x = te)) + 
    geom_histogram(bins = 15, color = "black", fill = "darkolivegreen4") + 
    labs(y = "Count", x = "Total TE across ROIs") + 
    theme_Publication() + theme(text = element_text(size = textsize), aspect.ratio = 0.5)
savename <- "figures/simulation_sloppy_tehist.png"
ggsave(savename)

ggplot(teacwtbl, aes(x = acw)) + 
    geom_histogram(bins = 15, color = "black", fill = "deeppink4") + 
    labs(x = "Slope of linear fit to ROIs vs ACW", y = "Count") + 
    theme_Publication() + theme(text = element_text(size = textsize), aspect.ratio = 0.5)
savename <- "figures/simulation_sloppy_acwhist.png"
ggsave(savename)
