# %% Do a pipeline to plot all relevant results from simulation manipulations
# Don't forget to add default value of network parameter

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
source("C:\\Users\\user\\Desktop\\brain_stuff\\ggplot_theme_Publication-master\\ggplot_theme_Publication-2.R")
setwd("C:\\Users\\user\\Desktop\\brain_stuff\\te_acw_simulations_v2\\src")
source("rhelpers.r")

setwd("C:\\Users\\user\\Desktop\\brain_stuff\\te_acw_simulations_v2\\.gitignore\\teplain_results")

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

# Start batch figure generation
zaaxd <- 1
    ## These three lines change in batch figure generation
    netw  <- netws[zaaxd]
    netwsave <- netw

    tes  <- readMat(paste(netw, "netw_task_te_plain.mat", sep = ""))
    te  <- tes$te

    qs <- readMat(paste0("greedy_", netw, "netw_task.mat"))
    q_xx <- qs$q.xx # roi x roi x sim x manip
    q_xy <- qs$q.xy
    q_yy <- qs$q.yy
    q_yx <- qs$q.yx

    betw <- q_xy + aperm(q_yx, c(2, 1, 3, 4))
    with <- q_xx + aperm(q_yy, c(2, 1, 3, 4))

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

    # %% For TEs, calculate node degree, incoming TE, outgoing TE per ROI
    nodedegree  <- array(data = NA, dim = c(nroi, nsim, nmanip))
    qbetw  <- array(data = NA, dim = c(nroi, nsim, nmanip))
    qwith  <- array(data = NA, dim = c(nroi, nsim, nmanip))
    for (m in 1:nmanip) {
        for (i in 1:nsim) {
            for (j in 1:nroi) {
                nodedegree[j, i, m]  <- sum(te[j, , i, m]) + sum(te[ , j, i, m])
                qbetw[j, i, m]  <- sum(betw[j, , i, m]) + sum(betw[ , j, i, m])
                qwith[j, i, m] <- sum(with[j, , i, m]) + sum(with[, j, i, m])
            }
        }
    }

    # %% Turn arrays into a tibble
    # Now, TEs are roi x sim x manipulation
    # First, convert both ACWs and TE calculations to sim x roi x man
    nodedegreep  <- aperm(nodedegree, c(2, 1, 3))
    qbetwp <- aperm(qbetw, c(2, 1, 3))
    qwithp <- aperm(qwith, c(2, 1, 3))

    ndt  <- nodedegreep %>% array2tibble(columnnames = irois, pagenames = dim3names, pagecol = "Manipulation") %>%
        pivot_longer(cols = !Manipulation, names_to = "ROI", values_to = "Degree")
    bigtibble <- ndt

    qbt <- qbetwp %>% array2tibble(columnnames = irois, pagenames = dim3names, pagecol = "Manipulation") %>% 
        pivot_longer(cols = !Manipulation, names_to = "ROI", values_to = "q_betw")
    qwt <- qwithp %>% array2tibble(columnnames = irois, pagenames = dim3names, pagecol = "Manipulation") %>% 
        pivot_longer(cols = !Manipulation, names_to = "ROI", values_to = "q_with")

    summarytibble <- bigtibble %>% group_by(Manipulation, ROI) %>%
        summarise(mean = median(Degree), se = se(Degree))
    qws <- qwt %>% group_by(Manipulation, ROI) %>%
        summarise(qw_mean = median(q_with), qw_se = se(q_with))
    qbs <- qbt %>% group_by(Manipulation, ROI) %>%
        summarise(qb_mean = median(q_betw), qb_se = se(q_betw))

    funslist <- list(median, se)
    names(funslist) <- c("mean", "se")

    tesummary <- summarytibble #  %>% pivot_wider(names_from = testyle, values_from = c(mean, se))

    finaltibble <- left_join(tesummary, qws)
    finaltibble <- left_join(finaltibble, qbs)

    finaltibble <- filter(finaltibble, Manipulation == "##  $\\ w_{EE}$ = 24.3 ##")

    # %% First, do the within q

    corresults <- finaltibble %>% 
                    select(Manipulation, mean, qw_mean) %>% 
                    nest(data = c("mean", "qw_mean")) %>% 
                    mutate(
                    test = map(data, ~cor.test(.x$mean, .x$qw_mean, method = "spearman")),
                    tidied = map(test, tidy)
                    ) %>% 
                    unnest(cols = tidied) %>% 
                    select(-data, -test) %>% 
                    stat2str()

    
    statstr1 <- c()
    statstr2 <- c()
    for (i in 1:length(corresults$statstr)) {
        statstr1[i] <- strsplit(corresults$statstr[i], ", ")[[1]][1]
        statstr2[i] <- strsplit(corresults$statstr[i], ", ")[[1]][2]
    }
    corresults <- add_column(corresults, statstr1 = statstr1)
    corresults <- add_column(corresults, statstr2 = statstr2)

    finaltibble_selected <- finaltibble
    corresults_selected <- corresults # an artifact

    ggplot(finaltibble_selected) + 
        geom_point(mapping = aes(x = mean, y = qw_mean, color = factor(ROI, level = irois)), size = 5) + 
        geom_smooth(mapping = aes(x = mean, y = qw_mean), method = "lm", color = "black") + 
        geom_errorbar(mapping = aes(x = mean, y = qw_mean, 
        ymin = qw_se$lowprc, ymax = qw_se$highprc)) + 
        geom_errorbarh(mapping = aes(x = mean, y = qw_mean, 
        xmin = se$lowprc, xmax = se$highprc)) + 
        geom_text(data = corresults_selected, parse = TRUE, # rho
        aes(x = -Inf, y = Inf, label = TeX(statstr1, output = "character"), color = "black"), 
        vjust = 1, hjust = -0.1, size = 8, color = "black") + 
        geom_text(data = corresults_selected, parse = TRUE, # p
        aes(x = -Inf, y = Inf, label = TeX(statstr2, output = "character"), color = "black"), 
        vjust = 2, hjust = -0.2, size = 8, color = "black") + 
        labs(x = "q_xx + q_yy", y = "Total TE", color = "") + 
        theme_Publication() + theme(aspect.ratio=1,
        text = element_text(size = 20), axis.text.x = element_text(angle = 90),
        strip.text = element_text(size=20))
    savename <- paste("figures/","defaultnetwork_qxx_te.png", sep = "")
    ggsave(savename)

    # Now, do the same for xy / yx

    corresults <- finaltibble %>% 
                select(Manipulation, mean, qb_mean) %>% 
                nest(data = c("mean", "qb_mean")) %>% 
                mutate(
                test = map(data, ~cor.test(.x$mean, .x$qb_mean, method = "spearman")),
                tidied = map(test, tidy)
                ) %>% 
                unnest(cols = tidied) %>% 
                select(-data, -test) %>% 
                stat2str()

    
    statstr1 <- c()
    statstr2 <- c()
    for (i in 1:length(corresults$statstr)) {
        statstr1[i] <- strsplit(corresults$statstr[i], ", ")[[1]][1]
        statstr2[i] <- strsplit(corresults$statstr[i], ", ")[[1]][2]
    }
    corresults <- add_column(corresults, statstr1 = statstr1)
    corresults <- add_column(corresults, statstr2 = statstr2)

    finaltibble_selected <- finaltibble
    corresults_selected <- corresults # an artifact

    ggplot(finaltibble_selected) + 
        geom_point(mapping = aes(x = mean, y = qb_mean, color = factor(ROI, level = irois)), size = 5) + 
        geom_smooth(mapping = aes(x = mean, y = qb_mean), method = "lm", color = "black") + 
        geom_errorbar(mapping = aes(x = mean, y = qb_mean, 
        ymin = qb_se$lowprc, ymax = qb_se$highprc)) + 
        geom_errorbarh(mapping = aes(x = mean, y = qb_mean, 
        xmin = se$lowprc, xmax = se$highprc)) + 
        geom_text(data = corresults_selected, parse = TRUE, # rho
        aes(x = -Inf, y = Inf, label = TeX(statstr1, output = "character"), color = "black"), 
        vjust = 1, hjust = -0.1, size = 8, color = "black") + 
        geom_text(data = corresults_selected, parse = TRUE, # p
        aes(x = -Inf, y = Inf, label = TeX(statstr2, output = "character"), color = "black"), 
        vjust = 2, hjust = -0.2, size = 8, color = "black") + 
        labs(x = "q_xy + q_yx", y = "Total TE", color = "") + 
        theme_Publication() + theme(aspect.ratio=1,
        text = element_text(size = 20), axis.text.x = element_text(angle = 90),
        strip.text = element_text(size=20))
    savename <- paste("figures/", "defaultnetwork_qxy_te.png", sep = "")
    ggsave(savename)
