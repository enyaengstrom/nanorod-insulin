# Installing Packages -----------------------------------------------------

library(DESeq2)
library(readxl)
library(MASS)
library(data.table)
library(tximport)
library(ggplot2)
library(ggfortify)
library(reshape2)
library(ComplexHeatmap)
library(clusterProfiler)
library(circlize)
library(RColorBrewer)
library(org.Mm.eg.db) # Mouse
library(ggnewscale)

# DE Analysis -----------------------------------------

# Load quantification files and collapse to gene-level
setwd("/Volumes/projects/Engström/Insulin Project/2022-02-10 Sequencing of NR-1 and NR-7")

# NOTE! Change path, working directory and number in "names(fls)" depending on which files we want.
fls <- list.files(path = "Salmon Quantification Files", pattern = "quant.sf", full.names = T, recursive = T)
names(fls) <- sapply(strsplit(fls, "/"), "[[", 2)

# Load reference 
t2g <- fread("https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.metadata.MGI.gz", header = F)

# Import quantification and link to gene ID
txi <- tximport(fls, type = "salmon", tx2gene = t2g[,1:2], ignoreAfterBar = TRUE)

# Construct file for metadata 
construct.meta <- function(files){
  meta <- data.table(sample=names(files))
  meta[, c("experiment", "condition", "time", "replicate") := tstrsplit(sample, "_")] 
  meta[, condition := factor(gsub("-", ".", condition, fixed = T))] 
  control = paste(meta$condition[grep("Ctrl", meta)])
  meta[, condition := relevel(condition, ref=control)]
  
  return(meta)
}
meta <- construct.meta(fls)

## DIFFERENTIAL EXPRESSION
dds.input <- DESeqDataSetFromTximport(txi, meta, design = ~condition) 

# LRT 
dds.LRT <- DESeq(dds.input, test = "LRT", reduced = ~1) 
res.LRT <- as.data.table(as.data.frame(results(dds.LRT)), keep.rownames = T)
res.LRT$diffExpressed <- "NO"
res.LRT$diffExpressed[res.LRT$log2FoldChange > 0.58 & res.LRT$padj < 1e-3] <- "UP"
res.LRT$diffExpressed[res.LRT$log2FoldChange < -0.58 & res.LRT$padj < 1e-3] <- "DOWN"
DE.LRT <- res.LRT[res.LRT$diffExpressed != "NO"]
colnames(res.LRT)[1] <- "Symbol"
colnames(DE.LRT)[1] <- "Symbol"

# Wald Test
dds.wald <- DESeq(dds.input)
res.wald <- as.data.table(as.data.frame(results(dds.wald)), keep.rownames = T)
res.wald$diffExpressed <- "NO"
res.wald$diffExpressed[res.wald$log2FoldChange > 0.58 & res.wald$padj < 1e-3] <- "UP"
res.wald$diffExpressed[res.wald$log2FoldChange < -0.58 & res.wald$padj < 1e-3] <- "DOWN"
DE.wald = res.wald[res.wald$diffExpressed != "NO"]
colnames(res.wald)[1] <- "Symbol"
colnames(DE.wald)[1] <- "Symbol"

# Run Wald testing DESeq analysis for all samples individually 
res.condition <- function(dds, sample, control){
  sample <- paste(sample)
  control <- paste(control)
  
  res <- as.data.table(as.data.frame(results(dds, contrast = c("condition", sample, control))), keep.rownames = T)
  res$diffExpressed <- "NO"
  res$diffExpressed[res$log2FoldChange > 0.58 & res$padj < 0.001] <- "UP"
  res$diffExpressed[res$log2FoldChange < - 0.58 & res$padj < 0.001] <- "DOWN"
  DE <- res[res$diffExpressed != "NO"]
  colnames(DE)[1] <- "Symbol"
  resList <- list(res, DE)
  return(resList)
}

# SECOND ROUND OF SEQUENCING - NR SAMPLES
# NR-7, 10nM
all.NR7.10 <- res.condition(dds.wald, meta$condition[15], meta$condition[1])
res.NR7.10 <- as.data.table(all.NR7.10[1])
DE.NR7.10 <- as.data.table(all.NR7.10[2])


# NR-7, 100nM
all.NR7.100 <- res.condition(dds.wald, meta$condition[13], meta$condition[1])
res.NR7.100  <- as.data.table(all.NR7.100[1])
DE.NR7.100 <- as.data.table(all.NR7.100[2])

# NR-1, 10nm
all.NR1.10 <- res.condition(dds.wald, meta$condition[10], meta$condition[1])
res.NR1.10 <- as.data.table(all.NR1.10[1])
DE.NR1.10 <- as.data.table(all.NR1.10[2])

# Free insulin, 10nM
all.ins.10 <- res.condition(dds.wald, meta$condition[7], meta$condition[1])
res.ins.10 <- as.data.table(all.ins.10[1])
DE.ins.10 <- as.data.table(all.ins.10[2])

# Free insulin, 100nM
all.ins.100 <- res.condition(dds.wald, meta$condition[4], meta$condition[1])
res.ins.100 <- as.data.table(all.ins.100[1])
DE.ins.100 <- as.data.table(all.ins.100[2])

DE.list <- list("All Genes" = DE.LRT, 
                   "Genes Insulin 10nM" = DE.ins.10,
                   "Genes Insulin 100nM" = DE.ins.100, 
                   "Genes NR-7 10nM" = DE.NR7.10, 
                   "Genes NR-7 100nM" = DE.NR7.100, 
                   "Genes NR-1 10nM" = DE.NR1.10)


# GO-Term GSEA --------------------------------------------
# Generate GSEA based on ranked list of log-fold change in genes
# DE-expression based on Wald-testing + contrast for each condition

ranked.list <- function(res){
  ranked.list <- res$log2FoldChange
  names(ranked.list) <- res$rn
  ranked.list <- na.omit(ranked.list)
  sorted.ranked.list <- sort(ranked.list, decreasing = TRUE)
  return(sorted.ranked.list)
}

ranked.ins.100 <- ranked.list(res.ins.100)
ranked.NR1.10 <- ranked.list(res.NR1.10)
ranked.NR7.10 <- ranked.list(res.NR7.10)
ranked.NR7.100 <- ranked.list(res.NR7.100)

# Replace name and gene list depending on sample analyzed 
gse.ins.100 <- gseGO(geneList=ranked.ins.100, 
             ont ="BP", 
             keyType = "SYMBOL", 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.1, 
             verbose = TRUE, 
             OrgDb = org.Mm.eg.db, 
             pAdjustMethod = "fdr",
             nPermSimple = 10000)

gsea.ins.100 <- as.data.table(gse.ins.100@result)
gsea.NR7.10 <- as.data.table(gse.NR7.10@result)
gsea.NR7.100 <- as.data.table(gse.NR7.100@result)
gsea.all <- list("NR7 10nM" = gsea.NR7.10, "NR7 100nM" = gsea.NR7.100, "Ins 100nM" = gsea.ins.100)

# Combine GSEA results 
combine.gsea <- function(gsea.list,number.of.pathways){
  # Create combined+reduced list
  gsea <- rbindlist(l = gsea.list, idcol = TRUE)
  names(gsea)[1] <- "Condition"
  
  # Only keep pathways occurring in X conditions
  pathways.freq <- as.data.frame(table(gsea$Description))
  common.pathways <- pathways.freq[pathways.freq$Freq <= number.of.pathways,]
  common.gsea <- gsea[gsea$Description %in% common.pathways$Var1]
  
  # Include -log(FDR) p-values
  common.gsea$logFDR <- -log10(common.gsea$p.adjust)
  
  return(common.gsea)
}

gsea <- combine.gsea(gsea.all,3)
gsea <- gsea[order(-NES),]
gsea.plot <- gsea[1:70,] # change depending on which to plot

# Make Dotplot (DP)
pdf(file = "GSEA Heatmap All.1 GO-Terms NES.pdf")
HM <- ggplot(gsea.plot, aes(x= factor(Condition, levels = c("Insulin 100nM", "NR7 10nM", "NR7 100nM")), 
                       y=reorder(Description,NES,mean), 
                       group= factor(Condition, levels = c("Insulin 100nM", "NR7 10nM", "NR7 100nM")))) + 
  geom_tile(aes(fill = gsea.plot$NES)) +
  scale_fill_distiller(palette = "RdBu", limits=c(min(gsea$NES),max(gsea$NES))) +
  ylab("Pathway") + xlab("Condition") + labs(fill ="NES") +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  theme_classic()
HM
dev.off()


# KEGG GSEA ------------------------------------------------

ranked.kegg <- function(res, ranked.list){
  ids = bitr(names(ranked.list), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
  dedup.ids = ids[!duplicated(ids[c("SYMBOL")]),]
  colnames(dedup.ids) = c("rn", "EntrezID")
  df = merge(res, dedup.ids, by = "rn")
  
  
  kegg.gene.list = df$log2FoldChange    # Create a vector of the gene universe
  names(kegg.gene.list) = df$EntrezID   # Name vector with ENTREZ ids
  kegg.gene.list = na.omit(kegg.gene.list)  # omit any NA values
  kegg.gene.list = sort(kegg.gene.list, decreasing = TRUE)  # sort the list in decreasing order (required for clusterProfiler)
  
  return(kegg.gene.list)
}

kegg.ins.100 <- ranked.kegg(res.ins.100, ranked.ins.100)
kegg.NR7.10 <- ranked.kegg(res.NR7.10, ranked.NR7.10)
kegg.NR7.100 <- ranked.kegg(res.NR7.100, ranked.NR7.100)
kegg.NR1.10 <- ranked.kegg(res.NR1.10, ranked.NR1.10)

# Change for each sample
gse.kegg.NR7.10 <- gseKEGG(geneList = kegg.NR7.10,
                    organism = "mmu",
                    minGSSize    = 3,
                    maxGSSize    = 800,
                    pvalueCutoff = 0.05,
                    pAdjustMethod = "fdr",
                    keyType       = "ncbi-geneid")

gsea.kegg.ins.100 <- gse.kegg.ins.100@result
gsea.kegg.NR7.10 <- gse.kegg.NR7.10@result
gsea.kegg.NR7.100 <- gse.kegg.NR7.100@result
gsea.kegg.NR1.10 <- gse.kegg.NR1.10@result
gsea.kegg.all <- list("NR7 10nM" = gsea.kegg.NR7.10,"NR7 100nM" = gsea.kegg.NR7.100, "NR1 10nM" = gsea.kegg.NR1.10)

gsea.kegg <- rbindlist(l = gsea.kegg.all, idcol = TRUE)
names(gsea.kegg)[1] <- "Condition"

pdf(file = "GSEA KEGG Dotplot NES.pdf")
KDP <- ggplot(gsea.kegg , aes(x= factor(Condition, levels = c("NR7 10nM", "NR7 100nM", "NR1 10nM")), 
                             y=reorder(Description,NES,mean), 
                             size=-log10(p.adjust), 
                             color=NES, 
                             group= factor(Condition, levels = c("NR7 10nM", "NR7 100nM", "NR1 10nM")))) + 
  geom_point(alpha = 0.9) +
  scale_colour_distiller(palette = "RdBu", limits=c(min(gsea.kegg$NES),max(gsea.kegg$NES))) +
  ylab("Pathway") + xlab("Condition") + labs(size = '-log FDR') +
  theme_classic()
KDP
dev.off()


# Plots -------------------------------------------------------------------

# Show dispersion and fit
plotDispEsts(dds.LRT)
dev.off()

# MA Plot
plotMA(results(dds.LRT))
dev.off()

# PCA plot
normalized.counts = counts(dds.LRT, normalized = T)
pca.input <- normalized.counts[rowSums(normalized.counts[])>0,]
pca <- prcomp(t(pca.input), scale. = TRUE)

pdf("PCA.pdf")
autoplot(pca, data = meta, colour = "condition") + scale_color_brewer(palette="Paired", direction = -1) + cowplot::theme_cowplot()
dev.off()

# Heatmap
ncol <- 6 # number of conditions for colors 
colOrder <- c(1,2,3,6,7,8,4,5,9,10,11,14,15,16,12,13) # Change to plot columns in a specific order

# Imported from http://www.informatics.jax.org/vocab/gene_ontology/GO:0050796
insulin.genes <- read_excel("/Users/enya.engstrom/OneDrive - Karolinska Institutet/quant/GO_term_summary_20220127_060749.xlsx")

mat <- normalized.counts[res.LRT[padj < 1e-3, Symbol],]
rha <- rowAnnotation(gene = anno_mark(at = which(DE.LRT$Symbol %in% insulin.genes$Symbol), labels = DE.LRT$Symbol[which(DE.LRT$Symbol %in% insulin.genes$Symbol)]))
cha <- HeatmapAnnotation(Condition = sapply(strsplit(colnames(mat), "_"), "[[", 2), 
                         col = list(Condition = setNames(brewer.pal(ncol, "Paired"), 
                                        meta[rev(order(sample)), unique(gsub(".","-",condition, fixed = T))])))

pdf(file = "Heatmap.pdf")
Heatmap(
  matrix = t(apply(mat, 1, scale)),
  col = colorRamp2(seq(-2,2, length.out = 11), rev(brewer.pal(11, "RdBu"))),
  right_annotation = rha,
  top_annotation = cha,
  name = "Z-score", 
  border = "black",
  column_order = colOrder,
  # column_split = rep(c("Insulin", "NR"), each = 8),
  show_row_names = F,
  show_row_dend = F,
  use_raster = F
)
dev.off()

# Upset plot
up.set.input <- list("Insulin 100nM" = DE.ins.100$Symbol,
                  "NR-7 10nM" = DE.NR7.10$Symbol,
                  "NR-7 100nM" = DE.NR7.100$Symbol,
                  "NR-1 10nM" = DE.NR1.10$Symbol)
a <- make_comb_mat(up.set.input, mode = "intersect")
UpSet(a, comb_order = order(comb_size(a), decreasing = TRUE))
