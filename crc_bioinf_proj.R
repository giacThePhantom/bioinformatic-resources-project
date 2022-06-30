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
ensembl <- useMart(biomart="ensembl",dataset="hsapiens_gene_ensembl")

pc_genes <- getBM(attributes=c("ensembl_gene_id","gene_biotype"),
          filters=c("ensembl_gene_id"), 
          values=rownames(raw_counts_df),
          mart = ensembl)

pc_genes <- pc_genes[pc_genes$gene_biotype == "protein_coding",]

raw_counts_df <- raw_counts_df[which(rownames(raw_counts_df)%in%pc_genes$ensembl_gene_id),]
r_anno_df <- r_anno_df[which(r_anno_df$gene_id%in%pc_genes$ensembl_gene_id),]

#3) Now first we filter genes with raw reads above 20 in at least 1 normal and 1 tumor

# count threshold
count_thr <- 20

# number of replicates with more counts than the count threshold
repl_thr <- 5

#we convert conditions in factors for the following function
c_anno_df$condition <- as.factor(c_anno_df$condition)

#we create a filtering vector
filter_vec <- apply(raw_counts_df,1,function(y) min(by(y, c_anno_df$condition, function(x) sum(x>count_thr))))

#now we select only genes with filtering characteristics and obtain a new matrix
filter_counts_df <- raw_counts_df[filter_vec>=repl_thr,]

# apply the filter on gene annotation
filtered_anno_df <- r_anno_df[rownames(filter_counts_df),]

# create a DGRList object
edge_c <- DGEList(counts=filter_counts_df,group=c_anno_df$condition,samples=c_anno_df,genes=filter_anno_df) 

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
deg$diffexpressed <- "NO"
deg$diffexpressed[deg$logFC > 1.5] <- "UP"
deg$diffexpressed[deg$logFC < -1.5] <- "DOWN"
ggplot(data=deg, aes(x=logFC, y=-log10(PValue), col=diffexpressed, label="")) +
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






