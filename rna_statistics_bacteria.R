# === CLIENT CONFIGURATION =====================================================
SPECIES_DB <- "org.EcK12.eg.db"      
PRIMARY_ID_TYPE <- "SYMBOL"         
# ==============================================================================
suppressPackageStartupMessages({
    library(DESeq2); library(tidyverse); library(EnhancedVolcano); library(pheatmap)
    library(RColorBrewer); library(ggplot2); library(SPECIES_DB, character.only=TRUE)
    library(clusterProfiler); library(AnnotationDbi)
})

dir.create("R_analysis/Plots", showWarnings = FALSE, recursive = TRUE)
dir.create("R_analysis/GO_Analysis", showWarnings = FALSE, recursive = TRUE)

# FATAL CRASH FIX: Prevents hyphenated sample names from breaking DESeq2
counts <- read.delim("R_analysis/counts.tsv", row.names=1, check.names=FALSE)
metadata <- read.delim("R_analysis/metadata.tsv", row.names=1, check.names=FALSE)

# FATAL CRASH FIX: Prevents model matrix rank failure 
if (length(unique(metadata$type)) < 2) {
    stop("FATAL ERROR: 'metadata.tsv' contains only 1 condition ('normal'). You must edit 'R_analysis/metadata.tsv' to define at least two groups (e.g., 'control' vs 'treated') before running DESeq2.")
}

dds = DESeqDataSetFromMatrix(counts, metadata, ~type)
rld <- vst(dds)
dds$type <- relevel(factor(dds$type), unique(metadata$type)[1])
dds <- DESeq(dds)
res <- results(dds)

pdf("R_analysis/Plots/01_PCA_Plot.pdf")
print(plotPCA(rld, intgroup = c("type")))
dev.off()

pdf("R_analysis/Plots/02_Volcano_Plot.pdf")
print(EnhancedVolcano(res, lab=rownames(res), x='log2FoldChange', y='padj', pCutoff=0.05, FCcutoff=1.5, title='Bacterial Differential Expression'))
dev.off()

resOrdered <- res[order(res$padj),]
write.csv(as.data.frame(resOrdered), file="R_analysis/Final_Differential_Expression.csv")

sig_genes <- rownames(subset(resOrdered, padj < 0.05 & abs(log2FoldChange) > 1))

if(length(sig_genes) > 0) {
    ego <- enrichGO(gene=sig_genes, OrgDb=get(SPECIES_DB), keyType=PRIMARY_ID_TYPE, ont="BP", pAdjustMethod="BH", pvalueCutoff=0.05)
    if(!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
        write.csv(as.data.frame(ego), "R_analysis/GO_Analysis/GO_Enrichment_Results.csv")
        pdf("R_analysis/Plots/03_GO_Dotplot.pdf")
        print(dotplot(ego, showCategory=15, title="Top Biological Processes"))
        dev.off()
    }
}
print("Bacterial Analysis Complete.")
