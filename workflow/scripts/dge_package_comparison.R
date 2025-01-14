rm(list = ls(all.names = TRUE))

########## Load libraries ##########
library("tidyverse")
library("edgeR")
library("DESeq2")
library("dplyr")
library("readr")
library("ggplot2")
library("RColorBrewer")
library("ashr")
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


########## Directory Output  ##########
dirOut <- paste0("workflow/results/DGE-comparison/")
if(!file.exists(dirOut)) { dir.create(path = dirOut, recursive = T) }
########## Directory Output  ##########


########## DGE - edgeR ##########
########## 1. Factors Preparing ##########
tp <- factor(spec$Timepoint, levels = c("t0", "t1"))
animal <- factor(spec$Animal, levels = c("01", "02", "03"))
design_edger <- model.matrix(~ 0 + tp + animal) 
########## 1. Factors Preparing ##########


###### 2. Create DGE List ######
y <- DGEList(counts = df, group = tp)
y <- normLibSizes(object = y, method = "TMM")
y <- estimateDisp(y, design_edger)
###### 2. Matrix Preparing ######


###### 3. Contrasts ######
colnames(design_edger) <- base::gsub(pattern = "tp", 
                                     replacement = "", 
                                     x = colnames(design_edger))

fit <- glmQLFit(y = y, design = design_edger)

contrasts <- makeContrasts(t01_vs_t0 = t1 - t0,
                           levels = colnames(design_edger))
###### 3. Contrasts ######


###### 4. DEG ######
qlf <- glmQLFTest(glmfit = fit, contrast = contrasts)
tTags <- as.data.frame(topTags(qlf, n = nrow(y), adjust.method = "BH"))

tTags$DirectionPvalue <- ifelse(test = tTags$logFC > 1 & 
                                  tTags$PValue < 0.05, yes = "Up", no = 
                                  ifelse(test = tTags$PValue < -1 & 
                                           tTags$FDR < 0.05, 
                                         yes = "Down", no = "None"))

tTags$Direction <- ifelse(test = tTags$logFC > 1 & 
                            tTags$FDR < 0.05, yes = "Up", no = 
                            ifelse(test = tTags$logFC < -1 & 
                                     tTags$FDR < 0.05, 
                                   yes = "Down", no = "None"))

###### 4. DEG ######
########## DGE - edgeR ##########


########## DGE - DESeq2 ##########
########## 1. Create a DESeqDataSet ##########
dds <- DESeqDataSetFromMatrix(countData = df, 
                              colData = spec, 
                              design = ~ Animal + Timepoint) 
########## 1. Prep ##########


########## 2. DGE  ##########
dds <- DESeq(dds)
resultsNames(dds)
resLFC <- as.data.frame(lfcShrink(dds, coef = "Timepoint_t1_vs_t0", type = "ashr"))

resLFC$DirectionPvalue <- ifelse(test = resLFC$log2FoldChange > 1 &
                                   resLFC$pvalue < 0.05 &
                                   resLFC$baseMean > 0,
                                 yes = "Up",
                                 no = ifelse(test = resLFC$log2FoldChange < -1 &
                                               resLFC$pvalue < 0.05 &
                                               resLFC$baseMean > 0,
                                             yes = "Down",
                                             no = "None"))

resLFC$Direction <- ifelse(test = resLFC$log2FoldChange > 1 &
                             resLFC$padj < 0.05 &
                             resLFC$baseMean > 0,
                           yes = "Up",
                           no = ifelse(test = resLFC$log2FoldChange < -1 &
                                         resLFC$padj < 0.05 &
                                         resLFC$baseMean > 0,
                                       yes = "Down",
                                       no = "None"))

resLFC$Direction <- replace(resLFC$Direction, is.na(resLFC$Direction), "None")
########## 2. DGE  ##########
########## DGE - DESeq2 ##########


########## Comparison ##########
########## Prepare table  ##########
dge_edgeR <- tTags %>% dplyr::select(logFC, PValue, FDR, Direction, DirectionPvalue)
dge_edgeR <- rownames_to_column(dge_edgeR)
names(dge_edgeR)[1] <- "Gene"

dge_deseq2 <- resLFC %>% dplyr::select(log2FoldChange, pvalue, padj, Direction, DirectionPvalue)
dge_deseq2 <- rownames_to_column(dge_deseq2)
names(dge_deseq2)[1] <- "Gene"

dge_inner <- inner_join(x = dge_edgeR, y = dge_deseq2, by = "Gene")


dge_edgeR_down <- sum(dge_edgeR$DirectionPvalue == "Down")
dge_edgeR_up <- sum(dge_edgeR$DirectionPvalue == "Up")
dge_edgeR_none <- sum(dge_edgeR$DirectionPvalue == "None")

dge_deseq2_down <- sum(dge_deseq2$DirectionPvalue == "Down")
dge_deseq2_up <- sum(dge_deseq2$DirectionPvalue == "Up")
dge_deseq2_none <- sum(dge_deseq2$DirectionPvalue == "None")


edger <- c(dge_edgeR_down, dge_edgeR_up, dge_edgeR_none)
deseq <- c(dge_deseq2_down, dge_deseq2_up, dge_deseq2_none)
overall <- rbind(edger, deseq)

colnames(overall) <- c("down", "up", "none")
########## Prepare table  ##########


########## Plot  ##########
p <- ggplot() + 
  geom_point(data = dge_inner, aes(x = logFC, y = log2FoldChange), size = 1) +
  ggtitle("t1 - t0") + xlab("log2FC - edgeR") + ylab("log2FC - DESeq2") +
  xlim(c(-10,10)) + 
  ylim(c(-7,7))

p

ggsave(filename = paste0(dirOut, "/", "dge_comparison.png"),
       width = 5, height = 5)
ggsave(filename = paste0(dirOut, "/", "dge_comparison.pdf"),
       width = 5, height = 5)
########## Plot  ##########
########## Comparison ##########