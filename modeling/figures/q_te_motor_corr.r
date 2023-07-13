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

setwd("C:\\Users\\user\\Desktop\\brain_stuff\\tepaper_figuresv2\\data")
textsize <- 20


tasksubjs <- readMat("Motort_subjs.mat")$subjs


roilabels <-  c("L_V1", "R_V1", "L_V2", "R_V2", "L_6d", "R_6d", "L_6mp", "R_6mp", "L_6v", "R_6v", "L_4", "R_4")
irois <- roilabels

    te <- readMat("te_motor.mat")$te # roi x roi x subject
    qs <- readMat("qs_Motort.mat")
    q_xx <- qs$q.xx # roi x roi x sim x manip
    q_xy <- qs$q.xy
    q_yy <- qs$q.yy
    q_yx <- qs$q.yx

    betw <- q_xy + aperm(q_yx, c(2, 1, 3))
    with <- q_xx + aperm(q_yy, c(2, 1, 3))


    nsim  <- dim(te)[3] # nsubj
    nroi  <- dim(te)[1]

    # ACWs : First dim: ROIs, second dim: simulations, 3rd dim: manipulations
    # TEs: roi x roi x sim x manipulation

    # %% For TEs, calculate node degree, incoming TE, outgoing TE per ROI
    nodedegree  <- array(data = NA, dim = c(nroi, nsim))
    qbetw  <- array(data = NA, dim = c(nroi, nsim))
    qwith  <- array(data = NA, dim = c(nroi, nsim))
        for (i in 1:nsim) {
            for (j in 1:nroi) {
                nodedegree[j, i]  <- sum(te[j, , i]) + sum(te[ , j, i])
                qbetw[j, i]  <- sum(betw[j, , i]) + sum(betw[ , j, i])
                qwith[j, i] <- sum(with[j, , i]) + sum(with[, j, i])
            }
        }

    # %% Turn arrays into a tibble
    rownames(nodedegree) <- roilabels
    colnames(nodedegree) <- tasksubjs
    rownames(qbetw) <- roilabels
    colnames(qbetw) <- tasksubjs
    rownames(qwith) <- roilabels
    colnames(qwith) <- tasksubjs
    
    ndt <- nodedegree %>% data.frame() %>% rownames_to_column(var = "ROI") %>% 
        as_tibble() %>% pivot_longer(cols = !ROI, names_to = "Subject", values_to = "Degree")
    bigtibble <- ndt
    qbt <- qbetw %>% data.frame() %>% rownames_to_column(var = "ROI") %>% 
        as_tibble() %>% pivot_longer(cols = !ROI, names_to = "Subject", values_to = "q_betw")
    qwt <- qwith %>% data.frame() %>% rownames_to_column(var = "ROI") %>% 
        as_tibble() %>% pivot_longer(cols = !ROI, names_to = "Subject", values_to = "q_with")

    summarytibble <- bigtibble %>% group_by(ROI) %>%
        summarise(mean = median(Degree), se = se(Degree))
    qws <- qwt %>% group_by(ROI) %>%
        summarise(qw_mean = median(q_with), qw_se = se(q_with))
    qbs <- qbt %>% group_by(ROI) %>%
        summarise(qb_mean = median(q_betw), qb_se = se(q_betw))

    tesummary <- summarytibble #  %>% pivot_wider(names_from = testyle, values_from = c(mean, se))

    finaltibble <- left_join(tesummary, qws)
    finaltibble <- left_join(finaltibble, qbs)

    # %% First, do the within q

    corresults <- finaltibble %>% 
                    select(mean, qw_mean) %>% 
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
    savename <- paste("figures/", "empirical_qxx_te_motor.png", sep = "")
    ggsave(savename)

    # Now, do the same for xy / yx

    corresults <- finaltibble %>% 
                select(mean, qb_mean) %>% 
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
    savename <- paste("figures/", "empirical_qxy_te_motor.png", sep = "")
    ggsave(savename)
