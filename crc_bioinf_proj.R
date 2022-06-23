crc<-load("data/Colorectal_Cancer.RData")
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
pc_genes<-getBM(attributes = c("ensembl_gene_id","entrezgene_id","external_gene_name"), filters = c("transcript_biotype"),values = list("protein_coding"), mart = ensembl)

pc_raw_counts_df<-raw_counts_df[rownames(raw_counts_df) %in% pc_genes$ensembl_gene_id,]

#Per il punto b
raw_counts <- pc_raw_counts_df[rowSums(pc_raw_counts_df > 20) >= 5,]

gene_length_symbol<-r_anno_df[r_anno_df$gene_id %in% rownames(raw_counts),]

# create a DGRList object
edge_c <- DGEList(counts=raw_counts,group=c_anno_df$condition,samples=c_anno_df,genes=gene_length_symbol)

# normalization with the edgeR package (TMM method)
edge_n <- calcNormFactors(edge_c,method="TMM")

# display normalization factors (also accounting for library size)
norm_factors <- mean(edge_n$samples$lib.size*edge_n$samples$norm.factors)/(edge_n$samples$lib.size*edge_n$samples$norm.factors)
names(norm_factors) <- edge_n$samples$sample

# create a cpm-rpkm table (normalized expression values)
cpm_table <- as.data.frame(round(cpm(edge_n),2))
head(cpm_table)

# define the experimental design matrix
design <- model.matrix(~0+group, data=edge_n$samples)
colnames(design) <- levels(edge_n$samples$group)
rownames(design) <- edge_n$samples$sample

# calculate dispersion and fit with edgeR (necessary for differential expression analysis)
edge_d <- estimateDisp(edge_n,design)
edge_f <- glmQLFit(edge_d,design)

# definition of the contrast (conditions to be compared)
contro <- makeContrasts("Case-Control", levels=design)

# fit the model with generalized linear models
edge_t <- glmQLFTest(edge_f,contrast=contro)

degs <- edge_t$table

degs['gene_id'] <- rownames(degs)
 degs <- merge(degs,pc_genes,by.x='gene_id',by.y="ensembl_gene_id")
 degs <- degs[which(!is.na(degs$entrezgene_id)),]
 degs <- degs[-which(duplicated(degs$entrezgene_id)),]


# Up and down regulated
upreg <- filter(degs, logCPM > 1 & PValue < 0.01 & logFC > 1.5 )
downreg <- filter(degs, logCPM > 1 & PValue < 0.01 & logFC < -1.5 )

# Volcano plot
degs$diffexpressed <- "NO"
degs$diffexpressed[degs$logFC > 1.5] <- "UP"
degs$diffexpressed[degs$logFC < -1.5] <- "DOWN"
ggplot(data=degs, aes(x=logFC, y=-log10(PValue), col=diffexpressed, label="")) +
    geom_point() +
    theme_minimal() +
    scale_color_manual(values=c("blue", "black", "red")) +
    geom_vline(xintercept=c(-1.5, 1.5), col="red") +
    geom_hline(yintercept=-log10(0.01), col="red")


#Heatmap
heatmap(as.matrix(cpm_table[head(rbind(upreg, downreg), 10)$gene_id,]))


#Gene set enrichment
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

head(upego_BP, 10)
head(downego_BP, 10)

head(upego_MF, 10)
head(downego_MF, 10)

head(upkegg, 10)
head(downkegg, 10)

#Plot the pathway with DEG genes, need to choose the pathway
log2FC <- degs$logFC
names(log2FC) <- degs$entrezgene_id
pathview(gene.data = log2FC, pathway.id="hsa05210")

