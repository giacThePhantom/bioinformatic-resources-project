crc<-load("Colorectal_Cancer.RData")
library(biomaRt)
library(edgeR)
library(dplyr)
library(ggrepel)


# annotate and retrieve protein coding genes
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
r_c<-raw_counts_df
pc_genes<-getBM(attributes = c("ensembl_gene_id"), filters = c("transcript_biotype"),values = list("protein_coding"), mart = ensembl)
pc_raw_counts_df<-raw_counts_df[rownames(raw_counts_df) %in% pc_genes[,],]
gene_length_symbol<-r_anno_df[r_anno_df$'gene_id' %in% pc_genes[,],]

##### edge R

# Setup
gene_sc <- c_anno_df
raw_counts <- pc_raw_counts_df

# Punto b
raw_count <- raw_counts[rowSums(raw_counts > 20) >= 5,]


dgeFull <- DGEList(raw_counts, group=gene_sc$condition)
dgeFull <- calcNormFactors(dgeFull, method = "TMM")
normCounts <- cpm(dgeFull)
dgeFull <- estimateCommonDisp(dgeFull)
dgeFull <- estimateTagwiseDisp(dgeFull)

dgeTest <- exactTest(dgeFull)

# Filtrattio per PValue e logCPM
deg <- filter(dgeTest$table, logCPM > 1 & PValue < 0.01)

# Up and down regulated
upreg <- filter(deg, logFC > 1.5)
downreg <- filter(deg, logFC < -1.5)


# Punto c
deg$diffexpressed <- "NO"
deg$diffexpressed[deg$logFC > 1.5] <- "UP"
deg$diffexpressed[deg$logFC < -1.5] <- "DOWN"
ggplot(data=deg, aes(x=logFC, y=-log10(PValue), col=diffexpressed, label="")) +
    geom_point() + 
    theme_minimal() +
    geom_text_repel() +
    scale_color_manual(values=c("blue", "black", "red")) +
    geom_vline(xintercept=c(-1.5, 1.5), col="red") +
    geom_hline(yintercept=-log10(0.01), col="red")


#Punto d
library(pheatmap)
pheatmap(upreg$logCPM)
