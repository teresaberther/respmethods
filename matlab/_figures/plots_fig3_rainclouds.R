###' Plotting script: raincloud plot for surrogate comparison
#'
#' This script reads in the phase-locking values of surrogate
#' signals with the original signal for different surrogate 
#' generation methods and creates a nice raincloud plot 
#' (Manuscript Fig. 3d)
#' 
#' Input:
#'   - phase-locking values (calculation see "plots_figure3.mat")
#'
#' Required R Packages:
#'   - dplyr, tidyr, ggplot2, ggdist
#'   
#' Copyright (C) 2026 Teresa Berther, University of Münster, Germany

# setup
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggdist)

setwd("~/respmethods/matlab/_figures")

# some nice colors 
cols = c(rgb(217, 216, 134, maxColorValue = 255),
         rgb(159, 193, 120, maxColorValue = 255),
         rgb(135, 198, 194, maxColorValue = 255),
         rgb(103, 150, 189, maxColorValue = 255))
cols_long = rep(cols, 15)

# load PLV data
plvs <- read.csv("_data4plotting/allplv.csv", sep=";", stringsAsFactors=TRUE)

# reshape
plvs_long <- as.data.frame(plvs) %>%
  mutate(sub = 1:n()) %>%
  pivot_longer(cols = -sub, names_to = "group", values_to = "V1")
plvs_long$group <- factor(plvs_long$group, levels = colnames(plvs))

# setup for nicer spacing
spacing <- 1.3                                                          # increase for more space between clouds
plvs_long$group_num <- as.numeric(factor(plvs_long$group)) * spacing

spacing_violin  <- 0.40   
spacing_boxplot <- 0.25
label_offset <- (spacing_boxplot + spacing_violin) / 2  

# make nice labels 
surr_labels <- c('iaaft', 'circular shifting', 'segment shuffling', 'random shuffling')

# create plot 
p <- ggplot(plvs_long, aes(x = group_num, y = V1, fill = group)) +
  
  stat_slab(
    side     = "right",
    scale    = 0.5,
    width    = 1.3,
    alpha    = 0.7,
    color    = NA,
    position = position_nudge(x = spacing_violin)
  ) +
  
  geom_boxplot(
    width         = 0.2,
    alpha         = 0.9,
    outlier.shape = NA,
    linewidth     = 0.8,
    color         = "black",
    position      = position_nudge(x = spacing_boxplot)
  ) +
  
  geom_point(
    stroke   = 0,
    shape    = 19,
    alpha    = 0.65,
    size     = 5,
    color    = cols_long,
    position = position_jitter(width = 0.04, height = 0)
  ) +
  
  scale_x_continuous(
    breaks = sort(unique(plvs_long$group_num)) + label_offset,
    labels = surr_labels
  ) +
  
  scale_fill_manual(values = cols, labels = surr_labels) +
  theme_classic(base_size = 18) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x  = element_text(size = 18),
    axis.ticks.x = element_blank(),
    axis.text.y  = element_text(size = 20),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.3),
    axis.line    = element_blank()
  ) +
  ylim(0, 0.07)

# show plot
p
