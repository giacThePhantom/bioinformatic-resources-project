library(MotifDb)
library(seqLogo)
library(PWMEnrich)
BiocManager::install("PWMEnrich.Hsapiens.background")
library(PWMEnrich.Hsapiens.background)
#1)first we load data
crc<-load("data/Colorectal_Cancer.RData")

#crc<-load("Colorectal_Cancer.RData")

#2)Using biomart, I can extract only protein coding genes
library('biomaRt')
ensembl <- useMart(biomart="ensembl",dataset="hsapiens_gene_ensembl")

pr_cod_df <- 
    getBM(attributes=c("ensembl_gene_id","gene_biotype"),
          filters=c("ensembl_gene_id"), 
          values=rownames(raw_counts_df),
          mart = ensembl)
pr_cod_df <- pr_cod_df[pr_cod_df$gene_biotype == "protein_coding",]
# now I update both raw_counts_df and r_anno_df
raw_counts_df <- raw_counts_df[which(rownames(raw_counts_df)%in%pr_cod_df$ensembl_gene_id),]
