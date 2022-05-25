crc<-load("./Downloads/RData files-20220510/Colorectal_Cancer.RData")
library(biomaRt)
library(edgeR)

# annotate and retrieve protein coding genes
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
r_c<-raw_counts_df
pc_genes<-getBM(attributes = c("ensembl_gene_id"), filters = c("transcript_biotype"),values = list("protein_coding"), mart = ensembl)
pc_raw_counts_df<-raw_counts_df[rownames(raw_counts_df) %in% pc_genes[,],]
pc_r_anno_df<-r_anno_df[r_anno_df$'gene_id' %in% pc_genes[,],]

#edgeR


