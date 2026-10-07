# === CLIENT CONFIGURATION =====================================================
SPECIES_DB <- "org.Hs.eg.db"
MSIGDBR_SPECIES <- "Homo sapiens"
PRIMARY_ID_TYPE <- "ENSEMBL"
# ==============================================================================
suppressPackageStartupMessages({
    library(DESeq2); library(tidyverse); library(EnhancedVolcano); library(pheatmap)
    library(RColorBrewer); library(ggplot2); library(SPECIES_DB, character.only=TRUE)
    library(clusterProfiler); library(msigdbr); library(AnnotationDbi); library(enrichplot)
})

dir.create("R_analysis/Plots", showWarnings = FALSE, recursive = TRUE)

counts <- read.delim("R_analysis/counts.tsv", row.names=1, check.names=FALSE)
metadata <- read.delim("R_analysis/metadata.tsv", row.names=1, check.names=FALSE)

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
print(EnhancedVolcano(res, lab=rownames(res), x='log2FoldChange', y='padj', pCutoff=10e-5, FCcutoff=2.0, title='Human Differential Expression'))
dev.off()

res$symbol <- tryCatch({ mapIds(get(SPECIES_DB), keys=row.names(res), column="SYMBOL", keytype=PRIMARY_ID_TYPE, multiVals="first") }, error=function(e){rep(NA,nrow(res))})
resOrdered <- res[order(res$padj),]
write.csv(as.data.frame(resOrdered), file="R_analysis/Final_Differential_Expression.csv")

dir.create("R_analysis/GSEA", showWarnings = FALSE)
df_gsea <- as.data.frame(resOrdered) %>% rownames_to_column("GeneID") %>% filter(!is.na(symbol)) %>% distinct(symbol, .keep_all=TRUE)
lfc <- sort(setNames(df_gsea$log2FoldChange, df_gsea$symbol), decreasing=TRUE)
pathways <- msigdbr(species=MSIGDBR_SPECIES, category="H")
gsea_res <- GSEA(geneList=lfc, pvalueCutoff=0.05, pAdjustMethod="BH", TERM2GENE=dplyr::select(pathways, gs_name, gene_symbol))

# FATAL CRASH FIX: Checking S4 object safely
if (!is.null(gsea_res) && nrow(as.data.frame(gsea_res)) > 0) {
    gseaResTidy <- as.data.frame(gsea_res) %>% arrange(desc(NES))
    readr::write_tsv(gseaResTidy, "R_analysis/GSEA/GSEA_results.tsv")

    pdf("R_analysis/Plots/03_GSEA_Summary.pdf")
    print(ggplot(gseaResTidy, aes(reorder(ID, NES), NES)) + geom_col(aes(fill=p.adjust<0.01)) + coord_flip())
    dev.off()

    pdf("R_analysis/Plots/04_GSEA_Top_Pathway.pdf")
    print(enrichplot::gseaplot(gsea_res, geneSetID=gseaResTidy$ID[1], title=gseaResTidy$ID[1]))
    dev.off()
}
print("Human Analysis Complete.")
