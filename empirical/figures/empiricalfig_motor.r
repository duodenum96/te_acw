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
library(matrixStats)
setwd("C:\\Users\\user\\Desktop\\brain_stuff\\te_acw_simulations_v2\\src")
source("rhelpers.r")
source("C:\\Users\\user\\Desktop\\brain_stuff\\ggplot_theme_Publication-master\\ggplot_theme_Publication-2.R")

setwd("C:\\Users\\user\\Desktop\\brain_stuff\\tepaper_figuresv2\\data")
textsize <- 30
tes <- readMat("te_motor.mat")$te # roi x roi x subject

tasksubjs <- readMat("Motort_subjs.mat")$subjs
restsubjs <- readMat("Restin_subjs.mat")$subjs

restacws <- readMat("acw_restmatched_MotortTFLA_round2.mat")$acw0
restacws2 <- apply(restacws, c(1, 2), mean) # Mean across scans

roilabels <-  c("L_V1", "R_V1", "L_V2", "R_V2", "L_6d", "R_6d", "L_6mp", "R_6mp", "L_6v", "R_6v", "L_4", "R_4")
irois <- roilabels

# Get ROI factor
roif <- roilabels[order(rowMedians(restacws2))]

# Take medians acws
restacwm <- rowMedians(restacws2)

nroi <- dim(tes)[1]
nsubj <- dim(tes)[3]

# Make a nice TE tibble

matx <- apply(tes, c(1, 2), mean)
rownames(matx)  <- roilabels
colnames(matx)  <- roilabels
nicedf  <- data.frame(Outgoing=rownames(matx)[row(matx)], Incoming=colnames(matx)[col(matx)], baddesign=c(matx))
names(nicedf)[names(nicedf) == "baddesign"] <- "TE"

tetibble  <- as_tibble(nicedf)

# Calculate node degree etc.
nodedegree  <- array(data = NA, dim = c(nroi, nsubj))
for (i in 1:nsubj) {
    for (j in 1:nroi) {
        nodedegree[j, i]  <- sum(tes[j, , i]) + sum(tes[ , j, i])
    }
}

rownames(nodedegree) <- roilabels
colnames(nodedegree) <- tasksubjs

ndt <- nodedegree %>% data.frame() %>% rownames_to_column(var = "ROI") %>% 
    as_tibble() %>% pivot_longer(cols = !ROI, names_to = "Subject", values_to = "nodedegree")

# Make ACW tibble
rownames(restacws2) <- roilabels
colnames(restacws2) <- restsubjs

acwrt <- restacws2 %>% data.frame() %>% rownames_to_column(var = "ROI") %>% 
    as_tibble() %>% pivot_longer(cols = !ROI, names_to = "Subject", values_to = "rest_ACW")

bigt <- full_join(acwrt, ndt)
# summary tibble
funslist <- list(median, se)
names(funslist) <- c("median", "se")

summarytibble <- bigt %>% group_by(ROI) %>% drop_na() %>% 
    summarise_at(.vars = c("rest_ACW", "nodedegree"), .funs = funslist)

setwd("C:\\Users\\user\\Desktop\\brain_stuff\\te_acw_simulations_v2\\.gitignore\\teplain_results")
# %% Start plotting
# Rest ACW
acwrt <- add_column(acwrt, ROIf = as.numeric(factor(acwrt$ROI, levels = roif)))

nestedacwrt <- acwrt %>% select(ROIf, rest_ACW) %>% nest(data = c("rest_ACW", "ROIf"))
acwcorrres <- nestedacwrt %>% mutate(
    test = map(data, ~cor.test(.x$ROIf, .x$rest_ACW, method = "spearman")),
    tidied = map(test, tidy)) %>% 
    unnest(cols = tidied) %>% 
    select(-data, -test) %>% 
    stat2str()

ggplot(data = bigt, 
    mapping = aes(x = factor(ROI, level = roif), y = rest_ACW, color = ROI)) + 
    geom_violin() + geom_point() + stat_summary(fun = "median", geom = "crossbar", width = 0.3, color = "black") + 
    geom_text(data = acwcorrres, parse = TRUE,
    aes(x = -Inf, y = Inf, label = TeX(statstr, output = "character"), color = "black"), 
        vjust = 1, hjust = -0.1, size = 12, color = "black") + 
    labs(x = "", y = "ACW-0") + 
    theme_Publication() + 
    theme(legend.position = "none", axis.text.x = element_text(angle = 90), aspect.ratio=1,
    text = element_text(size = textsize))
savename <- paste("figures/", "empirical_acw_rest_motor.png", sep = "")
ggsave(savename)

# Heatmap of TE

ggplot(data=tetibble, aes(x=factor(Incoming, level=roif), y = factor(Outgoing, level=rev(roif)))) + 
    geom_tile(aes(fill=TE)) + 
    scale_fill_gradient(low = "white", high = "red") +
    labs(x = "To", y = "From") + 
    coord_equal() + theme_Publication() + 
    theme(legend.text = element_text(angle = 45), axis.text.x = element_text(angle = 90),
    text = element_text(size = textsize), legend.key.width = unit(1.5, "cm"))
savename <- paste("figures/", "empirical_te_heatmap_motor.png", sep = "")
ggsave(savename)

# Violin plot of TE stats
idata <- select(bigt, ROI, nodedegree)
idata <- add_column(idata, ROIf = as.numeric(factor(idata$ROI, levels = roif)))
nestedidata <- idata %>% select(nodedegree, ROIf) %>% nest(data = c("nodedegree", "ROIf"))
tecorrres <- nestedidata %>% mutate(
    test = map(data, ~cor.test(.x$ROIf, .x$nodedegree, method = "spearman")),
    tidied = map(test, tidy)) %>% 
    unnest(cols = tidied) %>% 
    select(-data, -test) %>% stat2str()

ggplot(bigt, mapping = aes(x = factor(ROI, level = roif), y = nodedegree, color = ROI)) + 
    geom_violin() + 
    geom_point() + 
    stat_summary(fun = "median", geom = "crossbar", width = 0.3, color = "black") + 
    geom_text(data = tecorrres, parse = TRUE,
    aes(x = -Inf, y = Inf, label = TeX(statstr, output = "character"), color = "black"), 
        vjust = 1, hjust = -0.1, size = 12, color = "black") + 
    labs(x = "", y = "Total TE", color = "") + 
    theme_Publication() + 
    theme(legend.position = "none", axis.text.x = element_text(angle = 90), aspect.ratio=1,
    text = element_text(size = textsize))
savename <- paste("figures/", "empirical_testats_motor.png", sep = "")
ggsave(savename)

# ACW - TE correlation
ggplot(summarytibble) + 
    geom_point(mapping = aes(x = rest_ACW_median, y = nodedegree_median, color = ROI), size = 10) + 
    stat_cor(mapping = aes(x = rest_ACW_median, y = nodedegree_median), method = "spearman", 
    size = 10, cor.coef.name = "rho", label.x.npc = 0.7, label.sep = "\n") + 
    geom_smooth(aes(x = rest_ACW_median, y = nodedegree_median), method = "lm", color = "black") + 
    geom_errorbar(aes(x = rest_ACW_median, y = nodedegree_median, 
    ymin = nodedegree_se$lowprc, ymax = nodedegree_se$highprc)) + 
    geom_errorbarh(aes(x = rest_ACW_median, y = nodedegree_median,
    xmin = rest_ACW_se$lowprc, xmax = rest_ACW_se$highprc)) + 
    labs(x = "ACW-0", y = "Total TE") + 
    theme_Publication() + theme(aspect.ratio=1,
    text = element_text(size = textsize))
savename <- paste("figures/", "empirical_acw_te_motor.png", sep = "")
ggsave(savename)
