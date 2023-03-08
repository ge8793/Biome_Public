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
a <- readxl::read_xlsx("data/Supplement_S5_forMatt.xlsx", skip = 1)
a <- a[c(1:3147),]

# x 
x_min <- min(a$Beta, na.rm = T)
x_max <- max(a$Beta, na.rm = T)

# y 
y_min <- min(-log10(a$P), na.rm = T)
y_max <- max(-log10(a$P), na.rm = T)

# G rum ====
plot_data <- subset(a, Exposure == "G. Ruminococcus (presence vs. absence)")
plot_data$sig[plot_data$P < 3.08E-5] <- TRUE # old = 1.96E-6; new = 3.08E-5
plot_data$sig[plot_data$P > 3.08E-5] <- FALSE # old = 1.96E-6; new = 3.08E-5
plot_data$colour <- discrete_palette[3]
label <- subset(plot_data, sig == T)
label$label <- label$Outcome
label$colour <- discrete_palette[1]
plot_data <- subset(plot_data, sig == F)
plot_data$label <- NA
plot_data <- rbind(plot_data, label)
labels_of_interest <- label$Outcome

pdf("figures/volcano/obs_GRuminococcus.pdf",
    height = 10, width = 10)

ggplot(plot_data, aes(x = Beta, y = -log10(P))) +
  geom_point(colour=plot_data$colour) +
  
  geom_text_repel(data = subset(plot_data, label %in% labels_of_interest), 
                  aes(label = label),
                  min.segment.length = 0, box.padding = 1,
                  point.padding = NA, hjust = 0, # this will force the labels in the centre and line them up
                  force = 100,
                  colour = "black", 
                  max.iter = 10000, max.overlaps = 1000,
                  seed = 821) +
  
  labs(x = "Association between G. Ruminococcus and metabolites (beta)", 
       y = "P value (-log10)") +
  
  theme_cowplot() + 
  theme(axis.title.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        legend.position = "none") 

dev.off()

# G. Parabacteroides (abundance) ====
plot_data <- subset(a, Exposure == "G. Parabacteroides (abundance)")
plot_data$sig[plot_data$P < 3.08E-5] <- TRUE # old = 1.96E-6; new = 3.08E-5
plot_data$sig[plot_data$P > 3.08E-5] <- FALSE # old = 1.96E-6; new = 3.08E-5
plot_data$colour <- discrete_palette[3]
label <- subset(plot_data, sig == T)
label$label <- label$Outcome
label$colour <- discrete_palette[1]
plot_data <- subset(plot_data, sig == F)
plot_data$label <- NA
plot_data <- rbind(plot_data, label)
labels_of_interest <- label$Outcome

pdf("figures/volcano/obs_GParabacteroides.pdf",
    height = 10, width = 10)

ggplot(plot_data, aes(x = Beta, y = -log10(P))) +
  geom_point(colour=plot_data$colour) +
  
  geom_text_repel(data = subset(plot_data, label %in% labels_of_interest), 
                  aes(label = label),
                  min.segment.length = 0, box.padding = 1,
                  point.padding = NA, hjust = 0, # this will force the labels in the centre and line them up
                  force = 1,
                  colour = "black", 
                  max.iter = 10000, max.overlaps = 10,
                  seed = 821) +
  
  labs(x = "Association between G. Parabacteroides and metabolites (beta)", 
       y = "P value (-log10)") +
  
  theme_cowplot() + 
  theme(axis.title.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        legend.position = "none") 

dev.off()

# G. unclassified, and O.  Bacteroidales ====
plot_data <- subset(a, Exposure == "G. unclassified, O.  Bacteroidales (presence vs. absence)")
plot_data$sig[plot_data$P < 3.08E-5] <- TRUE # old = 1.96E-6; new = 3.08E-5
plot_data$sig[plot_data$P > 3.08E-5] <- FALSE # old = 1.96E-6; new = 3.08E-5
plot_data$colour <- discrete_palette[3]
label <- subset(plot_data, sig == T)
label$label <- label$Outcome
label$colour <- discrete_palette[1]
plot_data <- subset(plot_data, sig == F)
plot_data$label <- NA
plot_data <- rbind(plot_data, label)
labels_of_interest <- label$Outcome

pdf("figures/volcano/obs_GUnclassified-OBacteroidales.pdf",
    height = 10, width = 10)

ggplot(plot_data, aes(x = Beta, y = -log10(P))) +
  geom_point(colour=plot_data$colour) +
  
  geom_text_repel(data = subset(plot_data, label %in% labels_of_interest), 
                  aes(label = label),
                  min.segment.length = 0, box.padding = 1,
                  point.padding = NA, hjust = 0, # this will force the labels in the centre and line them up
                  force = 10,
                  colour = "black", 
                  max.iter = 10000, max.overlaps = 1000,
                  seed = 821) +
  
  labs(x = "Association between G. unclassified, and O.  Bacteroidales and metabolites (beta)", 
       y = "P value (-log10)") +
  
  theme_cowplot() + 
  theme(axis.title.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        legend.position = "none") 

dev.off()
