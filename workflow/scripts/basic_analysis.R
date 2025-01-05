rm(list = ls(all.names = TRUE))


########## Load libraries ##########
library("tidyverse")
library("edgeR")
library("dplyr")
library("readr")
library("ggfortify")
library("ggplot2")
library("limma")
library("reshape2")
########## Load libraries ##########

########## Import data as raw_counts and libs' specification as spec ##########
raw_counts <- readr::read_tsv(file = "src/raw_counts.tsv", col_names = TRUE)
spec <- readr::read_tsv(file = "src/spec.tsv", col_names = TRUE)
########## Import data as raw_counts and phenodata as lib_spec ##########


########## Prep lib_spec ##########
spec <- spec %>% 
  tibble::column_to_rownames("Sample")
########## Prep lib_spec ##########


########## Prep df ##########
df <- raw_counts
anyDuplicated(df)
df <- data.frame(df, row.names = 1)
all(rownames(spec) == colnames(df))
########## Prep df ##########

########## Decide factors ##########
spec$Timepoint <- factor(spec$Timepoint, levels = c("t0", "t1"))
spec$Animal <- factor(spec$Animal, levels = c("01", "02", "03"))
spec$Class <- factor(spec$Class, levels = c("control", "treat"))
########## Decide factors ##########


########## Basic analysis ##########
# All outputs will be saved on this folder
dirOut <- paste0("workflow/results/basics/")
if(!file.exists(dirOut)) { dir.create(path = dirOut, recursive = T) }

###### 1. CPM ######
cpm <- cpm(df, log = TRUE, prior.count = 2, normalized.lib.sizes = TRUE)
min_cpm <- min(cpm)
max_cpm <- max(cpm)

png(filename = paste0(dirOut, "/log2cpm_boxplot.png"))
boxplot(cpm, xlab="", ylab="Log2 counts per million",las=2)
abline(h=median(cpm),col="blue")
title("Boxplots of log2CPMs")
dev.off()


cpm_normalized <- normalizeQuantiles(cpm)
cpm_melted <- melt(cpm_normalized)


p <- ggplot(cpm_melted, aes(x = Var2, y = value)) +
  geom_boxplot() +
  labs(title = "Boxplot of Normalized CPM Data", x = "Samples", y = "Normalized Expression Level") +
  theme_minimal() + 
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5))
p

ggsave(filename = paste0(dirOut, "/", "log2cpm_normalized.png"),
       width = 5, height = 5)
ggsave(filename = paste0(dirOut, "/", "log2cpm_normalized.pdf"),
       width = 5, height = 5)

###### 2. PCA ######
pca <- prcomp(t(cpm))

# Without clustering
pcaPlot <- autoplot(pca, 
                    data = spec, 
                    colour = "Timepoint", 
                    shape = "Animal", 
                    size = 3) +
  ggtitle("PCA") + 
  theme(plot.title = element_text(hjust = 0.5))

# With clustering (frame)
pcaPlot2 <- autoplot(pca, 
                     data = spec, 
                     colour = "Timepoint", 
                     shape = "Animal", 
                     size = 4, 
                     frame = TRUE) + 
  ggtitle("PCA - clustering") + 
  theme(plot.title = element_text(hjust = 0.5)) 


ggsave(filename = paste0(dirOut, "/", "PCA_no_cluster.png"), plot = pcaPlot)
ggsave(filename = paste0(dirOut, "/", "PCA_no_cluster.pdf"), plot = pcaPlot)

ggsave(filename = paste0(dirOut, "/", "PCA_cluster.png"), plot = pcaPlot2)
ggsave(filename = paste0(dirOut, "/", "PCA_cluster.png"), plot = pcaPlot2)
########## Basic analysis ##########