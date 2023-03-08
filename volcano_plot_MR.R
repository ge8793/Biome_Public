rm(list = ls())

# environment ====
library(stringr)
library(remotes)
library(data.table)
library(calibrate)
library(ggrepel)
library(ggthemes)
library(devtools)
library(plyr) 
library(dplyr)
library(ggplot2)
library(png)
library(tidyr)
library(fuzzyjoin)
library(cowplot)

# colour palette ====
library(wesanderson)
colours <- names(wes_palettes)
discrete_palette <- wes_palette(colours[8], type = "discrete")
discrete_palette2 <- wes_palette(colours[3], type = "discrete")

# axis limits ====
a <- readxl::read_xlsx("data/Supplement_S4_forMatt.xlsx", skip = 1)
a <- a[c(1:414),]
table(a$`MR Method`)
b <- subset(a, `MR Method` == "Inverse variance weighted")
c <- subset(a, `MR Method` == "Wald ratio")
a <- rbind(b,c)

# x 
x_min <- min(a$Beta, na.rm = T)
x_max <- max(a$Beta, na.rm = T)

# y 
y_min <- min(-log10(a$P), na.rm = T)
y_max <- max(-log10(a$P), na.rm = T)

# overall ====
plot_data <- subset(a, Outcome == "Overall Breast Cancer")
plot_data$sig[plot_data$P < 8.47E-4] <- TRUE
plot_data$sig[plot_data$P > 8.47E-4] <- FALSE
plot_data$colour <- discrete_palette[3]
label <- subset(plot_data, sig == T)
label$label <- label$Exposure
label$colour <- discrete_palette[1]
plot_data <- subset(plot_data, sig == F)
plot_data$label <- NA
plot_data <- rbind(plot_data, label)
labels_of_interest <- label$Exposure

pdf("figures/volcano/MR_overall.pdf",
    height = 6, width = 6)

ggplot(plot_data, aes(x = Beta, y = -log10(P))) +
  geom_point(colour=plot_data$colour) +
  
  geom_text_repel(data = subset(plot_data, label %in% labels_of_interest), 
                  aes(label = label),
                  min.segment.length = 0, box.padding = 1,
                  xlim = 1, point.padding = NA, hjust = 0, # this will force the labels in the centre and line them up
                  force = 1,
                  colour = "black", 
                  max.iter = 10000, max.overlaps = 10,
                  seed = 821) +
  
  labs(x = "Association between metabolites and overall breast cancer (beta)", 
       y = "P value (-log10)") +
  
  theme_cowplot() + 
  theme(axis.title.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        legend.position = "none") 

dev.off()

# TNBC ====
plot_data <- subset(a, Outcome == "Meta-analysed TNBC")
plot_data$sig[plot_data$P < 8.06E-4] <- TRUE
plot_data$sig[plot_data$P > 8.06E-4] <- FALSE
plot_data$colour <- discrete_palette[3]
label <- subset(plot_data, sig == T)
label$label <- label$Exposure
label$colour <- discrete_palette[1]
plot_data <- subset(plot_data, sig == F)
plot_data$label <- NA
plot_data <- rbind(plot_data, label)
labels_of_interest <- label$Exposure

pdf("figures/volcano/MR_TNBC.pdf",
    height = 6, width = 6)

ggplot(plot_data, aes(x = Beta, y = -log10(P))) +
  geom_point(colour=plot_data$colour) +
  
  geom_text_repel(data = subset(plot_data, label %in% labels_of_interest), 
                  aes(label = label),
                  min.segment.length = 0, box.padding = 1,
                  xlim = 1, point.padding = NA, hjust = 0, # this will force the labels in the centre and line them up
                  force = 1,
                  colour = "black", 
                  max.iter = 10000, max.overlaps = 10,
                  seed = 821) +
  
  labs(x = "Association between metabolites and triple negative breast cancer (beta)", 
       y = "P value (-log10)") +
  
  theme_cowplot() + 
  theme(axis.title.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        legend.position = "none") 

dev.off()

# HER2negative ====
plot_data <- subset(a, Outcome == "HER2 negative Breast Cancer")
plot_data$sig[plot_data$P < 7.24E-4] <- TRUE
plot_data$sig[plot_data$P > 7.24E-4] <- FALSE
plot_data$colour <- discrete_palette[3]
label <- subset(plot_data, sig == T)
label$label <- label$Exposure
label$colour <- discrete_palette[1]
plot_data <- subset(plot_data, sig == F)
plot_data$label <- NA
plot_data <- rbind(plot_data, label)
labels_of_interest <- label$Exposure

pdf("figures/volcano/MR_HER2negative.pdf",
    height = 6, width = 6)

ggplot(plot_data, aes(x = Beta, y = -log10(P))) +
  geom_point(colour=plot_data$colour) +
  
  geom_text_repel(data = subset(plot_data, label %in% labels_of_interest), 
                  aes(label = label),
                  min.segment.length = 0, box.padding = 1,
                  xlim = 1, point.padding = NA, hjust = 0, # this will force the labels in the centre and line them up
                  force = 1,
                  colour = "black", 
                  max.iter = 10000, max.overlaps = 10,
                  seed = 821) +
  
  labs(x = "Association between metabolites and HER2 negative breast cancer (beta)", 
       y = "P value (-log10)") +
  
  theme_cowplot() + 
  theme(axis.title.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        legend.position = "none") 

dev.off()
# HER2positive ====
plot_data <- subset(a, Outcome == "HER2 positive Breast Cancer")
plot_data$sig[plot_data$P < 7.24E-4] <- TRUE
plot_data$sig[plot_data$P > 7.24E-4] <- FALSE
plot_data$colour <- discrete_palette[3]
label <- subset(plot_data, sig == T)
label$label <- label$Exposure
label$colour <- discrete_palette[1]
plot_data <- subset(plot_data, sig == F)
plot_data$label <- NA
plot_data <- rbind(plot_data, label)
labels_of_interest <- label$Exposure

pdf("figures/volcano/MR_HER2positive.pdf",
    height = 6, width = 6)

ggplot(plot_data, aes(x = Beta, y = -log10(P))) +
  geom_point(colour=plot_data$colour) +
  
  geom_text_repel(data = subset(plot_data, label %in% labels_of_interest), 
                  aes(label = label),
                  min.segment.length = 0, box.padding = 1,
                  xlim = 1, point.padding = NA, hjust = 0, # this will force the labels in the centre and line them up
                  force = 1,
                  colour = "black", 
                  max.iter = 10000, max.overlaps = 10,
                  seed = 821) +
  
  labs(x = "Association between metabolites and HER2 positive breast cancer (beta)", 
       y = "P value (-log10)") +
  
  theme_cowplot() + 
  theme(axis.title.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        legend.position = "none") 

dev.off()

# luminal A ====
plot_data <- subset(a, Outcome == "Luminal A Breast Cancer")
plot_data$sig[plot_data$P < 7.24E-4] <- TRUE
plot_data$sig[plot_data$P > 7.24E-4] <- FALSE
plot_data$colour <- discrete_palette[3]
label <- subset(plot_data, sig == T)
label$label <- label$Exposure
label$colour <- discrete_palette[1]
plot_data <- subset(plot_data, sig == F)
plot_data$label <- NA
plot_data <- rbind(plot_data, label)
labels_of_interest <- label$Exposure

pdf("figures/volcano/MR_luminalA.pdf",
    height = 6, width = 6)

ggplot(plot_data, aes(x = Beta, y = -log10(P))) +
  geom_point(colour=plot_data$colour) +
  
  geom_text_repel(data = subset(plot_data, label %in% labels_of_interest), 
                  aes(label = label),
                  min.segment.length = 0, box.padding = 1,
                  xlim = 1, point.padding = NA, hjust = 0, # this will force the labels in the centre and line them up
                  force = 1,
                  colour = "black", 
                  max.iter = 10000, max.overlaps = 10,
                  seed = 821) +
  
  labs(x = "Association between metabolites and luminal A breast cancer (beta)", 
       y = "P value (-log10)") +
  
  theme_cowplot() + 
  theme(axis.title.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        legend.position = "none") 

dev.off()

# luminal B ====
plot_data <- subset(a, Outcome == "Luminal B Breast Cancer")
plot_data$sig[plot_data$P < 7.24E-4] <- TRUE
plot_data$sig[plot_data$P > 7.24E-4] <- FALSE
plot_data$colour <- discrete_palette[3]
label <- subset(plot_data, sig == T)
label$label <- label$Exposure
label$colour <- discrete_palette[1]
plot_data <- subset(plot_data, sig == F)
plot_data$label <- NA
plot_data <- rbind(plot_data, label)
labels_of_interest <- label$Exposure

pdf("figures/volcano/MR_luminalB.pdf",
    height = 6, width = 6)

ggplot(plot_data, aes(x = Beta, y = -log10(P))) +
  geom_point(colour=plot_data$colour) +
  
  geom_text_repel(data = subset(plot_data, label %in% labels_of_interest), 
                  aes(label = label),
                  min.segment.length = 0, box.padding = 1,
                  xlim = 1, point.padding = NA, hjust = 0, # this will force the labels in the centre and line them up
                  force = 1,
                  colour = "black", 
                  max.iter = 10000, max.overlaps = 10,
                  seed = 821) +
  
  labs(x = "Association between metabolites and luminal B breast cancer (beta)", 
       y = "P value (-log10)") +
  
  theme_cowplot() + 
  theme(axis.title.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        legend.position = "none") 

dev.off()
