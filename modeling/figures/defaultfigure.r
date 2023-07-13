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
datasavelocation <- "C:\\Users\\user\\Desktop\\brain_stuff\\te_acw_simulations_v2\\figures\\figuredata"

axtextsize <- 30

zaaxd <- 1
    # General function for LaTeX
    appender <- function(string) {
        TeX(paste0(string))
    }
    ## These three lines change in batch figure generation
    netw  <- netws[zaaxd]
    netwsave <- "defaultnetw"
    # File name to be saved
    savename <- paste0(datasavelocation, "\\", netwsave, "_rformat.RData")
    load(savename)

    ###### ACW hierarchy
    acwrt <- filter(acwrt, Manipulation == "##  $\\ w_{EE}$ = 24.3 ##")
    acwcorrres <- filter(acwcorrres, Manipulation == "##  $\\ w_{EE}$ = 24.3 ##")
    ggplot(data = right_join(acwrt, acwcorrres), 
        mapping = aes(x = factor(ROI, level = irois), y = ACW_0, color = ROI)) + 
        geom_violin() + geom_point() + 
        geom_text(data = acwcorrres, parse = TRUE, # rho
        aes(x = -Inf, y = Inf, label = TeX(statstr1, output = "character"), color = "black"), 
        vjust = 1, hjust = -0.1, size = 12, color = "black") + 
        geom_text(data = acwcorrres, parse = TRUE, # p
        aes(x = -Inf, y = Inf, label = TeX(statstr2, output = "character"), color = "black"), 
        vjust = 2, hjust = -0.2, size = 12, color = "black") + 
        stat_summary(fun = "median", geom = "crossbar", width = 0.3, color = "black") + 
        labs(x = "", y = "ACW-0") + 
        theme_Publication() + 
        theme(legend.position = "none", axis.text.x = element_text(angle = 90), aspect.ratio=1,
        text = element_text(size = axtextsize))
    savename <- paste("figures/", netwsave, "_acw_rest.png", sep = "")
    ggsave(savename)

    tet <- filter(tet, Manipulation == "##  $\\ w_{EE}$ = 24.3 ##")
    #### TE heatmap
    ggplot(data = tet, aes(x=factor(Incoming, level=irois), y = factor(Outgoing, level=rev(irois)))) + 
        geom_tile(aes(fill=TE)) + 
        scale_fill_gradient(low = "white", high = "red") +
        labs(x = "To", y = "From") + 
        coord_equal() + theme_Publication() + 
        theme(legend.text = element_text(angle = 45), axis.text.x = element_text(angle = 90),
        text = element_text(size = axtextsize),
        strip.text = element_text(size = axtextsize), legend.key.width = unit(1.5, "cm"))
    savename <- paste("figures/", netwsave, "te_heatmap.png", sep = "")
    ggsave(savename)

    idata <- filter(idata, Manipulation == "##  $\\ w_{EE}$ = 24.3 ##")
    tecorrres <- filter(tecorrres, Manipulation == "##  $\\ w_{EE}$ = 24.3 ##")
    ##### TE hierarchy
    ggplot(data = idata, mapping = aes(x = factor(ROI, level = irois), y = Degree, color = ROI)) + 
        geom_violin() + geom_point() + 
        stat_summary(fun = "median", geom = "crossbar", width = 0.3, color = "black") + 
        geom_text(data = tecorrres, parse = TRUE, # rho
        aes(x = -Inf, y = Inf, label = TeX(statstr1, output = "character"), color = "black"), 
        vjust = 1, hjust = -0.1, size = 12, color = "black") + 
        geom_text(data = tecorrres, parse = TRUE, # p
        aes(x = -Inf, y = Inf, label = TeX(statstr2, output = "character"), color = "black"), 
        vjust = 2, hjust = -0.2, size = 12, color = "black") + 
        labs(x = "", y = "Total TE") + 
        theme_Publication() + 
        theme(legend.position = "none", axis.text.x = element_text(angle = 90), aspect.ratio = 1,
        text = element_text(size = axtextsize))
    savename <- paste("figures/", netwsave, "testats.png", sep = "")
    ggsave(savename)

    acwandte <- filter(acwandte, Manipulation == "##  $\\ w_{EE}$ = 24.3 ##")
    acwtecorrres <- filter(acwtecorrres, Manipulation == "##  $\\ w_{EE}$ = 24.3 ##")
    ## ACW - TE correlation ##################
    ggplot(acwandte) + 
        geom_point(mapping = aes(x = ACW_0_median, y = median, color = factor(ROI, level = irois)), size = 8) + 
        geom_smooth(mapping = aes(x = ACW_0_median, y = median), method = "lm", color = "black") + 
        geom_errorbar(mapping = aes(x = ACW_0_median, y = median, 
        ymin = se$lowprc, ymax = se$highprc)) + 
        geom_errorbarh(mapping = aes(x = ACW_0_median, y = median, 
        xmin = ACW_0_se$lowprc, xmax = ACW_0_se$highprc)) + 
        geom_text(data = acwtecorrres, parse = TRUE, # rho
        aes(x = -Inf, y = Inf, label = TeX(statstr1, output = "character"), color = "black"), 
        vjust = 1, hjust = -0.1, size = 12, color = "black") + 
        geom_text(data = acwtecorrres, parse = TRUE, # p
        aes(x = -Inf, y = Inf, label = TeX(statstr2, output = "character"), color = "black"), 
        vjust = 2, hjust = -0.2, size = 12, color = "black") + 
        labs(x = "ACW-0", y = "Total TE", color = "") + 
        theme_Publication() + theme(aspect.ratio=1,
        text = element_text(size = axtextsize), axis.text.x = element_text(angle = 90),
        strip.text = element_text(size = axtextsize))
    savename <- paste("figures/", netwsave, "acw_te.png", sep = "")
    ggsave(savename)
