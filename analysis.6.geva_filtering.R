###GET READY######

library(RcppCNPy)
library(tidyverse)
library(data.table)
library(qvalue)
library(wesanderson)

setwd(dir = "~/Desktop/Working/Sea_Age/")

#read in datasets
MNL_norm_sig_gene_overlap <- fread("MNL_norm_sig_gene_overlap") %>% select(V5, V6, V7, V4, V8, V9, V10, V11) %>%  arrange(V9) %>% distinct(V4, .keep_all = TRUE)
colnames(MNL_norm_sig_gene_overlap) <- c("CHROM", "POS", "BPtoo", "gene",  "top_Xpcent_wins_by_frac_gt_2", "frac_sites_xpnsl_gt_2", "top_Xpcent_of_wins_by_frac_lt_2", "frac_sites_xpnsl_lt_2")

FNL_norm_sig_gene_overlap <- fread("FNL_norm_sig_gene_overlap") %>% select(V5, V6, V7, V4, V8, V9, V10, V11) %>%  arrange(V9) %>% distinct(V4, .keep_all = TRUE)
colnames(FNL_norm_sig_gene_overlap) <- c("CHROM", "POS", "BPtoo", "gene",  "top_Xpcent_wins_by_frac_gt_2", "frac_sites_xpnsl_gt_2", "top_Xpcent_of_wins_by_frac_lt_2", "frac_sites_xpnsl_lt_2")

MQN_norm_sig_gene_overlap <- fread("MQN_norm_sig_gene_overlap") %>% select(V5, V6, V7, V4, V8, V9, V10, V11) %>%  arrange(V9) %>% distinct(V4, .keep_all = TRUE)
colnames(MNL_norm_sig_gene_overlap) <- c("CHROM", "POS", "BPtoo", "gene",  "top_Xpcent_wins_by_frac_gt_2", "frac_sites_xpnsl_gt_2", "top_Xpcent_of_wins_by_frac_lt_2", "frac_sites_xpnsl_lt_2")

FQN_norm_sig_gene_overlap <- fread("FQN_norm_sig_gene_overlap") %>% select(V5, V6, V7, V4, V8, V9, V10, V11) %>%  arrange(V9) %>% distinct(V4, .keep_all = TRUE)
colnames(MNL_norm_sig_gene_overlap) <- c("CHROM", "POS", "BPtoo", "gene",  "top_Xpcent_wins_by_frac_gt_2", "frac_sites_xpnsl_gt_2", "top_Xpcent_of_wins_by_frac_lt_2", "frac_sites_xpnsl_lt_2")

PC_FST_GWA_PCcor_MSW_1SW_gene_overlap <- fread("All_2PCcorr/PC_FST_GWA_PCcor_MSW1SW_gene_overlap_stats.tsv")

#get XPnsl outliers
OL_gene_all <- unique(c(MNL_norm_sig_gene_overlap$gene,  FNL_norm_sig_gene_overlap$gene, MQN_norm_sig_gene_overlap$gene,FQN_norm_sig_gene_overlap$gene))

#get coordinates
PC_FST_GWA_PCcor_MSW_1SW_gene_overlap_noOL <- PC_FST_GWA_PCcor_MSW_1SW_gene_overlap %>% filter(!gene %in% OL_gene_all)
PC_FST_GWA_PCcor_MSW_1SW_gene_overlap_OL <- PC_FST_GWA_PCcor_MSW_1SW_gene_overlap %>% filter(gene %in% OL_gene_all)
XPnsl_sig_genes <- unique(PC_FST_GWA_PCcor_MSW_1SW_gene_overlap_OL$gene)

chroms <- unique(PC_FST_GWA_PCcor_MSW_1SW_gene_overlap_noOL$CHROM)

#get random positions for each chromosome
geva_chrom_filt <- function(chrom) {
  chrom_noOL_rando_1K <- PC_FST_GWA_PCcor_MSW_1SW_gene_overlap_noOL %>% filter(CHROM %in% chrom)
  chrom_noOL_rando_1K <-  chrom_noOL_rando_1K %>%  select(POS) %>% 
    sample_n(1000) %>% arrange(POS)
  write.table(chrom_noOL_rando_1K, paste0("rand1k", "_", chrom, ".txt"), 
              quote = F, row.names = F,
              col.names = F, sep = "\t")
  }

map(chroms, geva_chrom_filt)

#get positions for each gene in XPnsl outliers
geva_gene_filt <- function(genes) {
  XPnsl_gene <- PC_FST_GWA_PCcor_MSW_1SW_gene_overlap_OL %>% filter(gene %in% genes)
  XPnsl_gene_pos <-  XPnsl_gene %>%  select(POS, CHROM) %>% 
  arrange(POS)
  CHROM <- unique(XPnsl_gene_pos$CHROM)
  POS <- XPnsl_gene_pos$POS
  write.table(POS, paste0("positions", "_", genes, "_", CHROM, ".txt"), 
              quote = F, row.names = F,
              col.names = F, sep = "\t")
}

map(XPnsl_sig_genes, geva_gene_filt)





#get positions for max scoring gene in PCA 
PC_maxsigPC1_gene_overlap<- PC_FST_GWA_PCcor_MSW_1SW_gene_overlap %>%  filter(qvals_PC1 < 0.05) %>%  group_by(gene) %>%
  filter(qvals_PC1 == min(qvals_PC1))
PC_maxsigPC1_gene_overlap <- data.frame(PC_maxsigPC1_gene_overlap)

PC1_gene_overlap <- PC_maxsigPC1_gene_overlap$gene


PCA_geva_gene_filt <- function(genes) {
  PCA_gene <- PC_maxsigPC1_gene_overlap %>% filter(gene %in% genes)
  PCA_gene_pos <- PCA_gene  %>%  select(POS, CHROM) %>% 
    arrange(POS)
  CHROM <- unique(PCA_gene_pos$CHROM)
  POS <- PCA_gene_pos$POS
  write.table(POS, paste0("PCAOL_positions", "_", genes, "_", CHROM, ".txt"), 
              quote = F, row.names = F,
              col.names = F, sep = "\t")
}

map(PC1_gene_overlap , PCA_geva_gene_filt )



