crc<-load("Colorectal_Cancer.RData")
library(biomaRt)
library(edgeR)
library(dplyr)
library(ggrepel)
library(fgsea)
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)
library(ggnewscale)
library(DOSE)
library(pathview)
library(tidyverse)



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
degs <- filter(dgeTest$table, logCPM > 1 & PValue < 0.01)
degs$id = rownames(degs)

### Serve per punto 4
ensembl2 <- useEnsembl(biomart = "ensembl",dataset = "hsapiens_gene_ensembl")
convert <- getBM(attributes=c("ensembl_gene_id","entrezgene_id","external_gene_name"),
                 filters=c("ensembl_gene_id"), 
                 values=degs$id,
                 mart = ensembl)

degs <- merge(degs,convert,by.x="id",by.y="ensembl_gene_id")
degs <- degs[which(!is.na(degs$entrezgene_id)),]
degs <- degs[-which(duplicated(degs$entrezgene_id)),]


# Up and down regulated
upreg <- filter(degs, logFC > 1.5)
downreg <- filter(degs, logFC < -1.5)


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


upego_BP <- enrichGO(gene = upreg$external_gene_name,
                   OrgDb = org.Hs.eg.db,
                   keyType = 'SYMBOL',
                   ont = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.05)

downego_BP <- enrichGO(gene = upreg$external_gene_name,
                   OrgDb = org.Hs.eg.db,
                   keyType = 'SYMBOL',
                   ont = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.05)


upego_MF <- enrichGO(gene = upreg$external_gene_name,
                     OrgDb = org.Hs.eg.db,
                     keyType = 'SYMBOL',
                     ont = "MF",
                     pAdjustMethod = "BH",
                     pvalueCutoff = 0.05,
                     qvalueCutoff = 0.05)

downego_MF <- enrichGO(gene = downreg$external_gene_name,
                       OrgDb = org.Hs.eg.db,
                       keyType = 'SYMBOL',
                       ont = "MF",
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.05)

upkegg <- enrichKEGG(gene = upreg$entrezgene_id,
                    organism = 'human',
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.05)

downkegg <- enrichKEGG(gene = downreg$entrezgene_id,
                    organism = 'human',
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.05)



