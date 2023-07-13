# %% Supplementary
# %% Do a pipeline to plot all relevant results from simulation manipulations

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
               c(0.8, 0.9, 1, 1.1, 1.2))
defvals <- c(24.3, 12.2, 19.7, 12.5, 0.06, 0.3510, 33.7, 25.3, 20, 10, 0.68)

nnetw <- length(netws)

# Start batch figure generation
zaaxd <- 7
    ## These three lines change in batch figure generation
    netw  <- netws[zaaxd]
    
    tes  <- readMat(paste(netw, "netw_task_te_plain.mat", sep = ""))
    te  <- tes$te
    acws_rest  <- readMat(paste("acws_", netw, "netw_rest.mat", sep = ""))
    acw_rest  <- acws_rest$acw0
    acws_task  <- readMat(paste("acws_", netw, "netw_task.mat", sep = ""))
    acw_task  <- acws_task$acw0

    netw <- paste0(r"( $\ )", netw)
    netw_X <- strsplit(netw, "_")
    netw_X[[1]][2] <- paste0("{", netw_X[[1]][2], "}")
    netw <- paste0(netw_X[[1]][1], "_", netw_X[[1]][2], "$")
    
    manip  <- manips[[zaaxd]]
    defval <- defvals[zaaxd]

    dim3names  <- paste0(netw, " = ", as.character(manip * defval))
    defi <- manip == 1
    dim3names[defi] <- paste0("## ", netw, " = ", as.character(manip[defi] * defval), " ##")
    dim3names <- factor(dim3names, level = dim3names)


    # get rid of nans at eta
    if (netw == "eta") {
        te <- te[, , , 1:3]
        acw_rest <- acw_rest[, , 1:3]
        acw_task <- acw_task[, , 1:3]
        dim3names <- dim3names[1:3]
        manip <- manip[1:3]
    }

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
    #indegree  <- colSums(te)
    #outdegree  <- rowSums(te)
    outdegree <- apply(te, c(1, 3, 4), sum) # Sum along rows
    indegree <- apply(te, c(2, 3, 4), sum) # Sum along columns

    # %% Turn arrays into a tibble
    # Now, TEs are roi x sim x manipulation
    # First, convert both ACWs and TE calculations to sim x roi x man
    nodedegreep  <- aperm(nodedegree, c(2, 1, 3))
    indegreep  <- aperm(indegree, c(2, 1, 3))
    outdegreep  <- aperm(outdegree, c(2, 1, 3))
    acw_restp  <- aperm(acw_rest, c(2, 1, 3))
    acw_taskp  <- aperm(acw_task, c(2, 1, 3))

    ndt  <- nodedegreep %>% array2tibble(columnnames = irois, pagenames = dim3names, pagecol = "Manipulation") %>%
        pivot_longer(cols = !Manipulation, names_to = "ROI", values_to = "Degree")

    indt <- indegreep  %>% array2tibble(columnnames = irois, pagenames = dim3names, pagecol = "Manipulation")  %>% 
        pivot_longer(cols = !Manipulation, names_to = "ROI", values_to = "Degree")

    outdt <- outdegreep  %>% array2tibble(columnnames = irois, pagenames = dim3names, pagecol = "Manipulation")  %>%
        pivot_longer(cols = !Manipulation, names_to = "ROI", values_to = "Degree")

    acwrt <- acw_restp %>% array2tibble(columnnames = irois, pagenames = dim3names, pagecol = "Manipulation")  %>% 
        pivot_longer(cols = !Manipulation, names_to = "ROI", values_to = "ACW_0") 
        # add_column(RT = rep("Rest", prod(dim(acw_restp))))

    acwtt <- acw_taskp %>% array2tibble(columnnames = irois, pagenames = dim3names, pagecol = "Manipulation")  %>% 
        pivot_longer(cols = !Manipulation, names_to = "ROI", values_to = "ACW_0")
        # add_column(RT = rep("Task", prod(dim(acw_taskp))))

    acwallt  <- add_column(acwrt, ACW_0_Task = acwtt$ACW_0)
    # Turn 4d TE array to 3D by averaging
    tex <- aperm(te, c(1,2,4,3))
    # Then turn it into a tibble
    tet <- temap2tibble(tex, irois, dim3names, "Manipulation", "TE")

    # Merge tibbles ndt, indt, outdt
    indt <- add_column(indt, testyle = "in_degree")
    outdt <- add_column(outdt, testyle = "out_degree")
    ndt <- add_column(ndt, testyle = "node_degree")


    bigtibble <- rbind(ndt, indt)
    bigtibble <- rbind(bigtibble, outdt)

    summarytibble <- bigtibble %>% group_by(Manipulation, ROI, testyle) %>%
        summarise(median = median(Degree), se = sd(Degree) / length(Degree))

    # se <- function(x) sd(x) / length(x)
    funslist <- list(median, se)
    names(funslist) <- c("median", "se")

    acwsummary <- acwallt %>% group_by(Manipulation, ROI) %>% 
        summarise_at(.vars = c("ACW_0", "ACW_0_Task"), .funs = funslist)

    tesummary <- summarytibble %>% pivot_wider(names_from = testyle, values_from = c(median, se))

    acwandte <- left_join(acwsummary, tesummary)

    # %% Use ggplot2 for stuff
    # Rest ACW
    # %% Different analyses
    #############################################
    acwrt <- add_column(acwrt, ROIf = as.numeric(factor(acwrt$ROI, levels = irois)))
    acwrt <- group_by(acwrt, Manipulation)

    nestedacwrt <- acwrt %>% select(Manipulation, ACW_0, ROIf) %>%  nest(data = c("ACW_0", "ROIf"))
    anan <- nestedacwrt %>% mutate(
        test = map(data, ~cor.test(.x$ROIf, .x$ACW_0, method = "spearman")),
        tidied = map(test, tidy)) %>% 
        unnest(cols = tidied) %>% 
        select(-data, -test) %>% stat2str()

    idata <- filter(bigtibble, testyle == "node_degree")
    idata <- add_column(idata, ROIf = as.numeric(factor(idata$ROI, levels = irois)))
    idata <- group_by(idata, Manipulation)
    nestedidata <- idata %>% select(Manipulation, Degree, ROIf) %>% nest(data = c("Degree", "ROIf"))
    anan2 <- nestedidata %>% mutate(
        test = map(data, ~cor.test(.x$ROIf, .x$Degree, method = "spearman")),
        tidied = map(test, tidy)) %>% 
        unnest(cols = tidied) %>% 
        select(-data, -test) %>% stat2str()
    
        ananzaaxd <- right_join(anan, acwrt)
        fulldata <- add_column(ananzaaxd, Degree = idata$Degree)

        fulldata <- add_column(fulldata, ROIf2 = factor(fulldata$ROIf))
        moi <- filter(fulldata, (Manipulation == "##  $\\ mu_{EE}$ = 33.7 ##") |  (Manipulation == " $\\ mu_{EE}$ = 32.015"))

        model <- anova(lm(Degree ~ ROIf2 + Manipulation + ROIf2 : Manipulation, data = moi))
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
        for (i in 1:length(unique(moi$ROIf))) {
            itibble <- filter(moi, ROIf == i)
            res <- wilcox.test(filter(itibble, Manipulation == "##  $\\ mu_{EE}$ = 33.7 ##")$Degree, 
                filter(itibble, Manipulation == " $\\ mu_{EE}$ = 32.015")$Degree, 
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
        
        
        fulldata_medians <- moi %>% group_by(ROIf, Manipulation) %>% summarise(Medians = median(Degree), se = se(Degree))
        fulldata_medians <- add_column(fulldata_medians, ROI = factor(irois[fulldata_medians$ROIf], levels = irois))


        p <- ggplot(fulldata_medians, 
            aes(x = Manipulation, y = Medians, colour = ROI, group = ROI)) +
            geom_point(size = 8) + geom_line() + 
            geom_errorbar(mapping = aes(x = Manipulation, y = Medians, 
            ymin = se$lowprc, ymax = se$highprc), width = 0.1) 

            for (i in 1:10) {
                p <- p + annotate(geom = "text", 
                label = paste0(irois[i], ": ", compstr[i]), x = Inf, y = Inf, hjust = 1, vjust = (i)*1.1, size = 7)
            }
            p <- p + annotate(geom = "text", label = TeX(anovareport, output = "character"), 
            x = -Inf, y = -Inf, hjust = 0, vjust = -1, size = 8, parse = TRUE)
            p + labs(x = "", y = "Median of Total TE", color = "") + 
            scale_x_discrete(labels=parse(text = c(TeX("$\\ mu_{EE}$ = 32.015", output = "character"), 
            TeX("$\\ mu_{EE}$ = 33.7"), output = "character"))) + 
            theme_Publication() + theme(aspect.ratio = 1,
            text = element_text(size = 25))
        
        savename <- "figures/mu_EE_supp_mu_EE.png"
        ggsave(savename)
