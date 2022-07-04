# Bioinformatics Resources Project A.Y. 2021/2022
# Group: Ilaria Cherchi, Giacomo Fantoni, Elisa Pettinà, Alessandro Polignano

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
filter_anno_df <- r_anno_df[rownames(filter_counts_df),]

# Create a DGEList object
edge_c <- DGEList(counts=filter_counts_df,group=c_anno_df$condition,samples=c_anno_df,genes=filter_anno_df)

# Normalization with the edgeR package (TMM method)
edge_n <- calcNormFactors(edge_c,method="TMM")

# Display normalization factors (also accounting for library size)
norm_factors <- mean(edge_n$samples$lib.size*edge_n$samples$norm.factors)/(edge_n$samples$lib.size*edge_n$samples$norm.factors)
names(norm_factors) <- edge_n$samples$sample

# Create a cpm-rpkm table (normalized expression values)
cpm_table <- as.data.frame(round(cpm(edge_n),2))
head(cpm_table)

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
ggsave("volcano_plot.pdf",path="plots")

#Heatmap representing the top 10 differentially expressed genes 5 from the up set and 5 from the down
upreg <- upreg[order(-upreg$logFC),]
downreg <- downreg[order(downreg$logFC),]

pdf(file="plots/heatmap.pdf")
cols <- c(rep("chartreuse4",50),rep("green",50))
pal <- c("blue","white","red") 
pal <- colorRampPalette(pal)(100)

ordered_anno <- transform(c_anno_df, n=nchar(as.character(condition)))
ordered_anno <- ordered_anno[with(ordered_anno, order(n, condition)), ]
ordered_anno <- subset(ordered_anno, select = -c(n))

heatmap_matrix <- as.matrix(cpm_table[rbind(head(upreg, 10), head(downreg, 10))$gene_id,])[, ordered_anno$sample]
heatmap_matrix <- as.data.frame(heatmap_matrix)
heatmap_matrix$gene_id <- rownames(heatmap_matrix)
hm_with_genename <- merge(degs, heatmap_matrix)
hm_with_genename <- subset(hm_with_genename, select = -c(logFC, logCPM, F, PValue, gene_biotype, entrezgene_id, description, diffexpressed, gene_id))
row.names(hm_with_genename) <- hm_with_genename$external_gene_name 
hm_with_genename <- subset(hm_with_genename, select = -c(external_gene_name))
heatmap(as.matrix(hm_with_genename), ColSideColors = cols, cexCol = 1,margins = c(4,4),col = pal, cexRow = 1, Colv = NA, Rowv = TRUE)
dev.off()

# Task 4. Perform gene set enrichment analysis using clusterProfiler R package.
#       a) Perform both GO (BP and MF) and KEGG analysis
#       b) Report the top 10 enriched GO terms and the top 10 enriched KEGG
#         pathways resulting from both up- and down-regulated gene lists

#Gene set enrichment
#From GO BP
upego_BP <- enrichGO(gene = upreg$external_gene_name,
                     OrgDb = org.Hs.eg.db,
                     keyType = 'SYMBOL',
                     ont = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff = 0.05,
                     qvalueCutoff = 0.05)

downego_BP <- enrichGO(gene = downreg$external_gene_name,
                       OrgDb = org.Hs.eg.db,
                       keyType = 'SYMBOL',
                       ont = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.05)

#From GO MF
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

# From KEGG
upkegg <- enrichKEGG(gene = upreg$entrezgene_id,
                     organism = 'human',
                     pvalueCutoff = 0.05,
                     qvalueCutoff = 0.05)

downkegg <- enrichKEGG(gene = downreg$entrezgene_id,
                       organism = 'human',
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.05)

# Top 10 enriched terms
head(upego_BP, 10)
head(downego_BP, 10)

head(upego_MF, 10)
head(downego_MF, 10)

head(upkegg, 10)
head(downkegg, 10)

get_enrich_plots <- function(df){
    #show(barplot(df, showCategory = 10))
    show(dotplot(df, showCategory = 10))
    #heatplot(df, showCategory = 4)
}
pdf("plots/uprego_BP.pdf")
get_enrich_plots(upego_BP)
dev.off()
pdf("plots/uprego_MF.pdf")
get_enrich_plots(upego_MF)
dev.off()
pdf("plots/upkegg.pdf")
get_enrich_plots(upkegg)
dev.off()
pdf("plots/downrego_BP.pdf")
get_enrich_plots(downego_BP)
dev.off()
pdf("plots/downego_MF.pdf")
get_enrich_plots(downego_MF)
dev.off()
pdf("plots/downkegg.pdf")
get_enrich_plots(downkegg)
dev.off()

# Task 5. Use the pathview R package to visualize one pathway you find enriched using the
#     upregulated gene list. 

#Plot the pathway with DEG genes, hsa05210 chosen as the "Colorectal cancer - Homo sapiens (human)" from KEGG
log2FC <- degs$logFC
names(log2FC) <- degs$entrezgene_id
setwd("plots/")
pathview(gene.data = log2FC, pathway.id="hsa05210",species = "human")
setwd("../")

# Task 6. Identify which transcription factors (TFs) have enriched scores in the promoters of all
#     up-regulated (or down-regulated if you prefer) genes.
#       a) use a window of 500 nucleotides upstream each gene

#Get sequences
seqs <- getSequence(id = upreg$gene_id,
                    type = "ensembl_gene_id",
                    seqType = "gene_flank",
                    upstream = 500,
                    mart = ensembl)
# 666 DNA sequences
DNA_set<- lapply(seqs$gene_flank,function(x) DNAString(x))
data(PWMLogn.hg19.MotifDb.Hsap)
enrich <- motifEnrichment(DNA_set,PWMLogn.hg19.MotifDb.Hsap,score = "affinity")
report <- groupReport(enrich,by.top.motifs = T)

#2287 genes
sig_report <- report[report$p.value < 0.05]
#643 with a p-value lower than 0.05

# Task 7. Select one among the top enriched TFs, compute the empirical distributions of scores
#     for all PWMs that you find in MotifDB for the selected TF and determine for all of
#     them the distribution (log2) threshold cutoff at 99.75%.

# we chose JUND
mdb.human.jund <- subset(MotifDb, organism=='Hsapiens' & geneSymbol=="JUND")
PWM <- toPWM(as.list(mdb.human.jund))
names(PWM) <- sapply(names(PWM),function(x) strsplit(x,"-")[[1]][3])

# raw.scores = T list of raw score values before cutoff
scores<- motifScores(DNA_set,PWM, raw.scores = T, verbose=T)
ecdf <- motifEcdf(PWM,organism = "hg19",quick=TRUE)


# Task 8. Identify which up-regulated (or down-regulated depending on the choice you made
#     at point 7) genes have a region in their promoter (defined as previously) with binding
#     scores above the computed thresholds for any of the previously selected PWMs.
#       a) Use pattern matching as done during the course

thresholds <- lapply(ecdf,function(x) quantile(x,0.9975))
scores <-  motifScores(DNA_set,PWM,raw.score=FALSE,cutoff=unlist(thresholds))
scores_sign <- which(apply(scores,1,sum)>0)
enriched_jund <- upreg[scores_sign,7]


# Task 9. Use STRING database to find PPI interactions among differentially expressed genes
#     and export the network in TSV format.
upreg_string <- upreg[order(upreg$PValue),]
downreg_string <- downreg[order(downreg$PValue),]
diff_string<-rbind(upreg_string,downreg_string)
# All differentially expressed genes
write(diff_string$gene_id,"data/degs_IDs.txt")
# Choose the 100 upregulated and 100 downregulated genes with lowest p-value
write(upreg_string$gene_id[1:100], "data/upreg_IDs.txt")
write(downreg_string$gene_id[1:100], "data/downreg_IDs.txt")


# Task 10. Import the network in R and using igraph package and identify and plot the largest
#     connected component 
# load STRING .tsv files
links_degs <-  read.delim("data/string_interactions_degs.tsv")
links_up <- read.delim("data/string_interactions_up.tsv")
links_down <- read.delim("data/string_interactions_down.tsv")

draw_largest_comp <- function(links){
    nodes<-union(links[,1],links[,2])
    net <- graph_from_data_frame(d=links,vertices=nodes,directed=FALSE)
    comp <- components(net, mode = "strong")
    # find largest component
    biggest_comp <- which.max(comp$csize)
    # print largest component size ## 1676
    print(comp$csize[biggest_comp])
    # isolate largest component
    first_c<-induced_subgraph(net,V(net)[comp$membership == biggest_comp])
    plot(first_c, 
         # edge proportional to combined score
         edge.width=E(first_c)$combined_score*3,
         vertex.color="orange",
         vertex.size=10,
         vertex.frame.color="darkgray",
         vertex.label.color="black", 
         vertex.label.cex=0.7,
         edge.curved=0.1)
}

pdf("plots/string_degs.pdf")
draw_largest_comp(links_degs)
dev.off()
pdf("plots/string_up.pdf")
draw_largest_comp(links_up)
dev.off()
pdf("plots/string_down.pdf")
draw_largest_comp(links_down)
dev.off()