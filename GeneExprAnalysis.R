library(DESeq2)
library(ggplot2)
library(clusterProfiler)
library(biomaRt)
library(ReactomePA)
library(DOSE)
library(KEGGREST)
library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(pheatmap)
library(genefilter)
library(RColorBrewer)
library(GO.db)
library(topGO)
library(dplyr)
library(gage)
library(ggsci)

setwd("Desktop/scRNA/new_workflow")
countdata <- read.table("example/final_counts.txt", header=TRUE, skip=1, row.names=1)
colnames(countdata) <- gsub(".bam", "", colnames(countdata), fixed=T)
colnames(countdata) <- gsub("..",   "", colnames(countdata), fixed=T)
countdata <- countdata[,c(-1:-5)]
countdata <- countdata[,c(3,4,1,2)]
countdata <- countdata[,c(1,2,4,3)]

metadata <- read.delim("example/metadata.txt", row.names = 1)
metadata$sampleid <- row.names(metadata)

ddsMat <- DESeqDataSetFromMatrix(countData=countdata, colData=metadata, design = ~Group)

ddsAnalysis <- DESeq(ddsMat)

results_table <- results(ddsAnalysis, pAdjustMethod = "fdr", alpha=0.05)

results_table$description <- mapIds(x=org.Mm.eg.db, keys=row.names(results_table),
                                    column="GENENAME", keytype="SYMBOL", multiVals="first")
results_table$symbol <- mapIds(x=org.Mm.eg.db, keys=row.names(results_table), 
                               column="ENTREZID", keytype="SYMBOL", multiVals = "first")
results_table$ensemble <- mapIds(x=org.Mm.eg.db, keys=row.names(results_table),
                                 column="ENSEMBL", keytype="SYMBOL", multiVals="first")

sig_results <- subset(results_table,padj<0.05)

ddsAnalysis_rlog <- rlog(ddsAnalysis)

plotPCA(ddsAnalysis_rlog, intgroup="Group", ntop=100) +
  scale_y_continuous(limits=c(-5,5))

topsig_genes <- assay(ddsAnalysis_rlog[row.names(sig_results)])[1:100,]

annotation_col <- data.frame(Group=factor(colData(ddsAnalysis_rlog)$Group),
                             Replicate=factor(colData(ddsAnalysis_rlog)$Replicate),
                             row.names = colData(ddsAnalysis_rlog)$sampleid)
pheatmap(topsig_genes,annotation_col = annotation_col)









# Gather Log-fold change and FDR-corrected pvalues from DESeq2 results
## - Change pvalues to -log10 (1.3 = 0.05)
data <- data.frame(gene = row.names(results_table),
                   padj = -log10(results_table$padj), 
                   log2FoldChange = results_table$log2FoldChange)

# Remove any rows that have NA as an entry
data <- na.omit(results_table)

# Color the points which are up or down
## If fold-change > 0 and pvalue > 1.3 (Increased significant)
## If fold-change < 0 and pvalue > 1.3 (Decreased significant)
results_table <- mutate(results_table, color = case_when(data$log2FoldChange > 0 & data$padj > 1.3 ~ "Increased",
                                       data$log2FoldChange < 0 & data$padj > 1.3 ~ "Decreased",
                                       data$padj < 1.3 ~ "nonsignificant"))

# Make a basic ggplot2 object with x-y values
vol <- ggplot(data, aes(x = log2FoldChange, y = padj, color = color))

# Add ggplot2 layers
vol +   
  ggtitle(label = "Volcano Plot", subtitle = "Colored by fold-change direction") +
  geom_point(size = 2.5, alpha = 0.8, na.rm = T) +
  scale_color_manual(name = "Directionality",
                     values = c(Increased = "#008B00", Decreased = "#CD4F39", nonsignificant = "darkgray")) +
  theme_bw(base_size = 14) + # change overall theme
  theme(legend.position = "right") + # change the legend
  xlab(expression(log[2]("LoGlu" / "HiGlu"))) + # Change X-Axis label
  ylab(expression(-log[10]("adjusted p-value"))) + # Change Y-Axis label
  geom_hline(yintercept = 1.3, colour = "darkgrey") + # Add p-adj value cutoff line
  scale_y_continuous(trans = "log1p") # Scale yaxis due to large p-values
