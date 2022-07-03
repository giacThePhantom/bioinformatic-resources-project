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
library(Biostrings)
library(enrichTF)
library(MotifDb)
library(seqLogo)
library(PWMEnrich)
library(PWMEnrich.Hsapiens.background)
library(igraph)

# Task 1. Load the RData file. The following three data-frames are available:
#       a) raw_counts_df = contains the raw RNA-seq counts
#       b) c_anno_df = contains sample name and condition (Case and Control)
#       c) r_ anno_df = contains the ENSEMBL genes ids, the length of the genes and
#       the genes symbols

crc<-load("data/Colorectal_Cancer.RData")


# Task 2. Update raw_count_df and r_anno_df extracting only protein coding genes.
#       a) Use biomaRt package to retrieve the needed information
#       b) Next tasks should use the new data-frames you have created


# connect to biomart
ensembl <- useMart(biomart="ensembl",dataset="hsapiens_gene_ensembl")

# get protein coding gens among the genes present in raw_counts_df
pc_genes <- getBM(attributes=c("ensembl_gene_id","gene_biotype","entrezgene_id","external_gene_name", "description"),
                  filters=c("ensembl_gene_id", "biotype"),
                  values=list(rownames(raw_counts_df), c("protein_coding")),
                  mart = ensembl)
# collapse records with same ensembl_gene_id into one
pc_genes <- pc_genes %>% distinct(ensembl_gene_id, .keep_all = TRUE)

# filter for protein coding genes
raw_counts_df <- raw_counts_df[which(rownames(raw_counts_df)%in%pc_genes$ensembl_gene_id),]
r_anno_df <- r_anno_df[which(r_anno_df$gene_id%in%pc_genes$ensembl_gene_id),]

# Task 3. Perform differential expression analysis using edgeR package and select up- and
#   down-regulated genes using a p-value cutoff of 0.01, a log fold change ratio >1.5 for
#   up-regulated genes and < (-1.5) for down-regulated genes and a log CPM >1. Relax
#   the thresholds if no or few results are available.
#     a) Use the workflow we developed during the course
#     b) Filter raw counts data retaining only genes with a raw count >20 in at least
#       5 Cases or 5 Control samples
#     c) Create a volcano plot of your results
#     d) Create an annotated heatmap focusing only on up- and downregulated
#       genes

# Count threshold
count_thr <- 20

# Number of replicates with more counts than the count threshold
repl_thr <- 5

# Convert conditions in factors for the following function
c_anno_df$condition <- as.factor(c_anno_df$condition)

# Create a filtering vector
filter_vec <- apply(raw_counts_df,1,function(y) min(by(y, c_anno_df$condition, function(x) sum(x>count_thr))))

# Select only genes with filtering characteristics and obtain a new matrix
filter_counts_df <- raw_counts_df[filter_vec>=repl_thr,]

# Apply the filter on gene annotation
filter_anno_df <- r_anno_df[match(rownames(filter_counts_df),r_anno_df$gene_id),]

# Create a DGEList object
edge_c <- DGEList(counts=filter_counts_df,group=c_anno_df$condition,samples=c_anno_df,genes=filter_anno_df)

# Normalization with the edgeR package (TMM method)
edge_n <- calcNormFactors(edge_c,method="TMM")

# Display normalization factors (also accounting for library size)
norm_factors <- mean(edge_n$samples$lib.size*edge_n$samples$norm.factors)/(edge_n$samples$lib.size*edge_n$samples$norm.factors)
names(norm_factors) <- edge_n$samples$sample

symbols <- edge_n$genes$symbol
symbols <- lapply(symbols, paste, collapse = ',')

# Create a cpm-rpkm table (normalized expression values)
cpm_table <- as.data.frame(round(cpm(edge_n),2))

# Define the experimental design matrix
design <- model.matrix(~0+group, data=edge_n$samples)
colnames(design) <- levels(edge_n$samples$group)
rownames(design) <- edge_n$samples$sample

# Calculate dispersion and fit with edgeR (necessary for differential expression analysis)
edge_d <- estimateDisp(edge_n,design)
edge_f <- glmQLFit(edge_d,design)

# Definition of the contrast (conditions to be compared)
contro <- makeContrasts("Case-Control", levels=design)

# Fit the model with generalized linear models
edge_t <- glmQLFTest(edge_f,contrast=contro)

# Retrieve data from the edge class
degs <- edge_t$table

# Merge edgeR-derived data with ENSEMBL-derived gene information, such that there is an entrezgene_id and it is unique (needed for KEGG)
degs['gene_id'] <- rownames(degs)
degs <- merge(degs,pc_genes,by.x='gene_id',by.y="ensembl_gene_id")
degs <- degs[which(!is.na(degs$entrezgene_id)),]
degs <- degs[-which(duplicated(degs$entrezgene_id)),]


# Split degs into up and down regulated genes
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
#ggsave("volcano_plot.pdf",path="plots")

#Heatmap representing the top 10 differentially expressed genes 5 from the up set and 5 from the down
upreg <- upreg[order(-upreg$logFC),]
downreg <- downreg[order(downreg$logFC),]
cols <- c(rep("chartreuse4",50),rep("green",50))
pal <- c("blue","white","red") 
pal <- colorRampPalette(pal)(100)
heatmap_matrix <- as.matrix(cpm_table[rbind(head(upreg, 5), head(downreg, 5))$gene_id,]) 
heatmap(heatmap_matrix, ColSideColors = cols, cexCol = 0.05,margins = c(4,4),col=pal,cexRow = 1)

        
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






