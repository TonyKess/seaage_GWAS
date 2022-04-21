###GET READY######

library(RcppCNPy)
library(tidyverse)
library(data.table)
library(qvalue)
library(wesanderson)

setwd(dir = "~/Desktop/Working/Sea_Age/")
#Read in metadata #####
#write IDs for running ANGSD genotyping and GL
Salmon_Metadata_Bam_Unique_HD  <- fread("Salmon_Metadata_Bam_Unique_HD.txt")
#Recode
Salmon_Metadata_Bam_Unique_HD$AgeClass[Salmon_Metadata_Bam_Unique_HD$AgeClass=="1SW"] <- 0
Salmon_Metadata_Bam_Unique_HD$AgeClass[Salmon_Metadata_Bam_Unique_HD$AgeClass=="MSW"] <- 1

#Set pop groupss
QCNB <- c("MSW", "MUN", "RTG", "TRI")
NL <- c("ENG", "MBB", "SAN", "RVF")

#Separate by region and sex for 
Male_1SW_MSW_QN <- Salmon_Metadata_Bam_Unique_HD %>% filter(PCR_Sex =="M", RiverCode %in% QCNB)
Male_1SW_MSW_NL <-  Salmon_Metadata_Bam_Unique_HD %>% filter(PCR_Sex =="M", RiverCode %in% NL)
Female_1SW_MSW_QN <- Salmon_Metadata_Bam_Unique_HD %>% filter(PCR_Sex =="F", RiverCode %in% QCNB)
Female_1SW_MSW_NL <- Salmon_Metadata_Bam_Unique_HD %>% filter(PCR_Sex =="F", RiverCode %in% NL)

#FQN
#list files
FQN_norm_files <- list.files(path = "./", pattern = ".*Female_1SWMSW_QN.xpnsl.*windows")

# get chrom IDs
chr_IDs <- str_replace(FQN_norm_files,  "_CIGENE.*", "") %>%  str_replace("Salmo_", "")

#read in files as list
All_FQN_norm <- map(.x = FQN_norm_files, .f = fread)

#change chrom names
names(All_FQN_norm) <- chr_IDs

FQN_norm <- map_df(All_FQN_norm , ~as.data.frame(.x), .id="CHROM")
colnames(FQN_norm)[2:10] <- c("win_start", "win_end", "num_sites", "frac_sites_xpnsl_gt_2",
                              "frac_sites_xpnsl_lt_2", "top_Xpcent_wins_by_frac_gt_2", 
                              "top_Xpcent_of_wins_by_frac_lt_2", "max_xpnsl", "min_xpnsl")
FQN_norm_sig_MSW <-  FQN_norm %>% filter(top_Xpcent_wins_by_frac_gt_2 %in% 1)
FQN_norm_sig_1SW <-  FQN_norm %>% filter(top_Xpcent_of_wins_by_frac_lt_2 %in% 1)
FQN_norm_sig <- bind_rows(FQN_norm_sig_MSW, FQN_norm_sig_1SW)
FQN_norm_sig_bed <- FQN_norm_sig  %>%  select(CHROM, win_start, win_end, top_Xpcent_wins_by_frac_gt_2, frac_sites_xpnsl_gt_2, top_Xpcent_of_wins_by_frac_lt_2, frac_sites_xpnsl_lt_2)
fwrite(FQN_norm_sig_bed, "FQN_norm_sig.bed", col.names = F, row.names = F, sep = "\t", quote = F)
system("bedtools intersect -b FQN_norm_sig.bed -a SSA_genes.bed -wb > FQN_norm_sig_gene_overlap")
FQN_norm_sig_gene_overlap <- fread("FQN_norm_sig_gene_overlap") %>% select(V5, V6, V7, V4, V8, V9, V10, V11) %>%  arrange(V9) %>% distinct(V4, .keep_all = TRUE)
colnames(FQN_norm_sig_gene_overlap) <- c("CHROM", "POS", "BPtoo", "gene",  "top_Xpcent_wins_by_frac_gt_2", "frac_sites_xpnsl_gt_2", "top_Xpcent_of_wins_by_frac_lt_2", "frac_sites_xpnsl_lt_2")


#MQN
#list files
MQN_norm_files <- list.files(path = "./", pattern = ".*Male_1SWMSW_QN.xpnsl.*windows")

#read in files as list
All_MQN_norm <- map(.x = MQN_norm_files, .f = fread)

#change chrom names
names(All_MQN_norm) <- chr_IDs

MQN_norm <- map_df(All_MQN_norm , ~as.data.frame(.x), .id="CHROM")
colnames(MQN_norm)[2:10] <- c("win_start", "win_end", "num_sites", "frac_sites_xpnsl_gt_2",
                              "frac_sites_xpnsl_lt_2", "top_Xpcent_wins_by_frac_gt_2", 
                              "top_Xpcent_of_wins_by_frac_lt_2", "max_xpnsl", "min_xpnsl")

MQN_norm_sig_MSW <-  MQN_norm %>% filter(top_Xpcent_wins_by_frac_gt_2 %in% 1)
MQN_norm_sig_1SW <-  MQN_norm %>% filter(top_Xpcent_of_wins_by_frac_lt_2 %in% 1)
MQN_norm_sig <- bind_rows(MQN_norm_sig_MSW, MQN_norm_sig_1SW)

MQN_norm_sig_bed <- MQN_norm_sig  %>%  select(CHROM, win_start, win_end, top_Xpcent_wins_by_frac_gt_2, frac_sites_xpnsl_gt_2, top_Xpcent_of_wins_by_frac_lt_2, frac_sites_xpnsl_lt_2)
fwrite(MQN_norm_sig_bed, "MQN_norm_sig.bed", col.names = F, row.names = F, sep = "\t", quote = F)
system("bedtools intersect -b MQN_norm_sig.bed -a SSA_genes.bed -wb > MQN_norm_sig_gene_overlap")
MQN_norm_sig_gene_overlap <- fread("MQN_norm_sig_gene_overlap") %>% select(V5, V6, V7, V4, V8, V9, V10, V11) %>%  arrange(V9) %>% distinct(V4, .keep_all = TRUE)
colnames(MQN_norm_sig_gene_overlap) <- c("CHROM", "POS", "BPtoo", "gene",  "top_Xpcent_wins_by_frac_gt_2", "frac_sites_xpnsl_gt_2", "top_Xpcent_of_wins_by_frac_lt_2", "frac_sites_xpnsl_lt_2")



#FNL
#list files
FNL_norm_files <- list.files(path = "./", pattern = ".*Female_1SWMSW_NL.xpnsl.*windows")

#read in files as list
All_FNL_norm <- map(.x = FNL_norm_files, .f = fread)

#change chrom names
names(All_FNL_norm) <- chr_IDs

FNL_norm <- map_df(All_FNL_norm , ~as.data.frame(.x), .id="CHROM")
colnames(FNL_norm)[2:10] <- c("win_start", "win_end", "num_sites", "frac_sites_xpnsl_gt_2",
                              "frac_sites_xpnsl_lt_2", "top_Xpcent_wins_by_frac_gt_2", 
                              "top_Xpcent_of_wins_by_frac_lt_2", "max_xpnsl", "min_xpnsl")
FNL_norm_sig_MSW <-  FNL_norm %>% filter(top_Xpcent_wins_by_frac_gt_2 %in% 1)
FNL_norm_sig_1SW <-  FNL_norm %>% filter(top_Xpcent_of_wins_by_frac_lt_2 %in% 1)
FNL_norm_sig <- bind_rows(FNL_norm_sig_MSW, FNL_norm_sig_1SW)
FNL_norm_sig_bed <- FNL_norm_sig  %>%  select(CHROM, win_start, win_end, top_Xpcent_wins_by_frac_gt_2, frac_sites_xpnsl_gt_2, top_Xpcent_of_wins_by_frac_lt_2, frac_sites_xpnsl_lt_2)
fwrite(FNL_norm_sig_bed, "FNL_norm_sig.bed", col.names = F, row.names = F, sep = "\t", quote = F)
system("bedtools intersect -b FNL_norm_sig.bed -a SSA_genes.bed -wb > FNL_norm_sig_gene_overlap")
FNL_norm_sig_gene_overlap <- fread("FNL_norm_sig_gene_overlap") %>% select(V5, V6, V7, V4, V8, V9, V10, V11) %>%  arrange(V9) %>% distinct(V4, .keep_all = TRUE)
colnames(FNL_norm_sig_gene_overlap) <- c("CHROM", "POS", "BPtoo", "gene",  "top_Xpcent_wins_by_frac_gt_2", "frac_sites_xpnsl_gt_2", "top_Xpcent_of_wins_by_frac_lt_2", "frac_sites_xpnsl_lt_2")


#MNL
#list files
MNL_norm_files <- list.files(path = "./", pattern = ".*Male_1SWMSW_NSL.xpnsl.*windows")

#read in files as list
All_MNL_norm <- map(.x = MNL_norm_files, .f = fread)

#change chrom names
names(All_MNL_norm) <- chr_IDs

MNL_norm <- map_df(All_MNL_norm , ~as.data.frame(.x), .id="CHROM")
colnames(MNL_norm)[2:10] <- c("win_start", "win_end", "num_sites", "frac_sites_xpnsl_gt_2",
                              "frac_sites_xpnsl_lt_2", "top_Xpcent_wins_by_frac_gt_2", 
                              "top_Xpcent_of_wins_by_frac_lt_2", "max_xpnsl", "min_xpnsl")

MNL_norm_sig_MSW <-  MNL_norm %>% filter(top_Xpcent_wins_by_frac_gt_2 %in% 1)
MNL_norm_sig_1SW <-  MNL_norm %>% filter(top_Xpcent_of_wins_by_frac_lt_2 %in% 1)
MNL_norm_sig <- bind_rows(MNL_norm_sig_MSW, MNL_norm_sig_1SW)
MNL_norm_sig_bed <- MNL_norm_sig  %>%  select(CHROM, win_start, win_end, top_Xpcent_wins_by_frac_gt_2, frac_sites_xpnsl_gt_2, top_Xpcent_of_wins_by_frac_lt_2, frac_sites_xpnsl_lt_2)
fwrite(MNL_norm_sig_bed, "MNL_norm_sig.bed", col.names = F, row.names = F, sep = "\t", quote = F)
system("bedtools intersect -b MNL_norm_sig.bed -a SSA_genes.bed -wb > MNL_norm_sig_gene_overlap")
MNL_norm_sig_gene_overlap <- fread("MNL_norm_sig_gene_overlap") %>% select(V5, V6, V7, V4, V8, V9, V10, V11) %>%  arrange(V9) %>% distinct(V4, .keep_all = TRUE)
colnames(MNL_norm_sig_gene_overlap) <- c("CHROM", "POS", "BPtoo", "gene",  "top_Xpcent_wins_by_frac_gt_2", "frac_sites_xpnsl_gt_2", "top_Xpcent_of_wins_by_frac_lt_2", "frac_sites_xpnsl_lt_2")

OL_gene_all <- unique(c(MNL_norm_sig_gene_overlap$gene,  FNL_norm_sig_gene_overlap$gene, MQN_norm_sig_gene_overlap$gene,FQN_norm_sig_gene_overlap$gene))

#get GWAS PCA FST stats
PC_FST_GWA_PCcor_MSW_1SW_gene_overlap <- fread("All_2PCcorr/PC_FST_GWA_PCcor_MSW1SW_gene_overlap_stats.tsv")

FQN_PC_FST_GWA_PCcorr_gene_overlap <- fread("F_QN_GWAS/FQN_PC_FST_GWA_PCcor_gene_overlap_stats.tsv")

FNL_PC_FST_GWA_PCcorr_gene_overlap <- fread("F_NL_GWAS/FNL_PC_FST_GWA_PCcor_gene_overlap_stats.tsv")

MQN_PC_FST_GWA_PCcorr_gene_overlap <- fread("M_QN_GWAS/MQN_PC_FST_GWA_PCcor_gene_overlap_stats.tsv")

MNL_PC_FST_GWA_PCcorr_gene_overlap <- fread("M_NL_GWAS/MNL_PC_FST_GWA_PCcor_gene_overlap_stats.tsv")

Overlap_genes <- unique(PC_FST_GWA_PCcor_MSW_1SW_gene_overlap$gene)

#This function resamples the number of outlier genes from each group randomly from the gene list
#then plots the distribution of expected overlap, and places the 99% boundary in red and the observed proportion in blue
parfunc <- function(Group1, Group2){
  G1L <- length(Group1$gene)
  G2L <- length(Group2$gene)
  
set.seed(121) ## for reproducibility
nsim <- 10000
res <- numeric(nsim) ## set aside space for results
for (i in 1:nsim) {
  ## standard approach: scramble response value
  x <- sample(Overlap_genes, G1L)
  y <- sample(Overlap_genes, G2L)
  ## compute & store difference in means; store the value
  res[i] <- length(y[y %in% x])/length(unique(c(x,y)))}

quantile(res) * 100
res <- data.frame(res)
ggplot() + geom_density(data = res, aes(x = res)) +
  geom_vline(xintercept = quantile(res$res, 0.999), colour = "red") +
  geom_vline(xintercept =  length(Group1$gene[Group1$gene %in% Group2$gene])/length(unique(c(Group1$gene, Group2$gene))),colour = "blue") + theme_classic()
}


#in regions
parfunc(Group1 = MNL_norm_sig_gene_overlap , Group2 = FNL_norm_sig_gene_overlap)
length(MNL_norm_sig_gene_overlap$gene[MNL_norm_sig_gene_overlap$gene %in% FNL_norm_sig_gene_overlap$gene])/
  length(unique(c(MNL_norm_sig_gene_overlap$gene,FNL_norm_sig_gene_overlap$gene)))

parfunc(Group1 = MQN_norm_sig_gene_overlap , Group2 = FQN_norm_sig_gene_overlap)
length(MQN_norm_sig_gene_overlap$gene[MQN_norm_sig_gene_overlap$gene %in% FQN_norm_sig_gene_overlap$gene])/
  length(unique(c(MQN_norm_sig_gene_overlap$gene,FQN_norm_sig_gene_overlap$gene)))

#in sex
parfunc(Group1 = FQN_norm_sig_gene_overlap , Group2 = FNL_norm_sig_gene_overlap)
length(FQN_norm_sig_gene_overlap$gene[FQN_norm_sig_gene_overlap$gene %in% FNL_norm_sig_gene_overlap$gene])/
  length(unique(c(FQN_norm_sig_gene_overlap$gene,FNL_norm_sig_gene_overlap$gene)))

parfunc(Group1 = MNL_norm_sig_gene_overlap , Group2 = MQN_norm_sig_gene_overlap)
length(MNL_norm_sig_gene_overlap$gene[MNL_norm_sig_gene_overlap$gene %in% MQN_norm_sig_gene_overlap$gene])/
  length(unique(c(MNL_norm_sig_gene_overlap$gene,MQN_norm_sig_gene_overlap$gene)))


#plot for chr 25 and 9, illustrious homes of vgll3 and six6.

ggplot() + geom_point(data = FQN_norm %>% filter(CHROM %in% "ssa25",top_Xpcent_wins_by_frac_gt_2 %in% 1), aes(y = max_xpnsl, x = win_start, colour = "FQN")) + 
  geom_point(data = FNL_norm %>% filter(CHROM %in% "ssa25", top_Xpcent_wins_by_frac_gt_2 %in% 1), aes(y = max_xpnsl, x = win_start, colour = "FNL")) +
  geom_point(data = MNL_norm %>% filter(CHROM %in% "ssa25", top_Xpcent_wins_by_frac_gt_2 %in% 1), aes(y = max_xpnsl, x = win_start, colour = "MNL")) +
  geom_point(data = MQN_norm %>% filter(CHROM %in% "ssa25", top_Xpcent_wins_by_frac_gt_2 %in% 1), aes(y = max_xpnsl, x = win_start, colour = "MQN")) +
  geom_smooth(data = FQN_norm %>% filter(CHROM %in% "ssa25",  !top_Xpcent_wins_by_frac_gt_2 %in% 1), aes(y = max_xpnsl, x = win_start, colour = "FQN"), method = "loess", span = 0.03, se = F) + 
  geom_smooth(data = FNL_norm %>% filter(CHROM %in% "ssa25",  !top_Xpcent_wins_by_frac_gt_2 %in% 1), aes(y = max_xpnsl, x = win_start, colour = "FNL"), method = "loess", span = 0.03, se = F) +
  geom_smooth(data = MNL_norm %>% filter(CHROM %in% "ssa25",  !top_Xpcent_wins_by_frac_gt_2 %in% 1), aes(y = max_xpnsl, x = win_start, colour = "MNL"), method = "loess", span = 0.03, se = F) +
  geom_smooth(data = MQN_norm %>% filter(CHROM %in% "ssa25",  !top_Xpcent_wins_by_frac_gt_2 %in% 1), aes(y = max_xpnsl, x = win_start, colour = "MQN"), method = "loess", span = 0.03, se = F) +
  geom_vline(xintercept =  28654947) +geom_vline(xintercept =  28659019) + theme_classic() 


ggplot() + geom_point(data = FQN_norm %>% filter(CHROM %in% "ssa09",top_Xpcent_wins_by_frac_gt_2 %in% 1), aes(y = max_xpnsl, x = win_start, colour = "FQN")) + 
  geom_point(data = FNL_norm %>% filter(CHROM %in% "ssa09", top_Xpcent_wins_by_frac_gt_2 %in% 1), aes(y = max_xpnsl, x = win_start, colour = "FNL")) +
  geom_point(data = MNL_norm %>% filter(CHROM %in% "ssa09", top_Xpcent_wins_by_frac_gt_2 %in% 1), aes(y = max_xpnsl, x = win_start, colour = "MNL")) +
  geom_point(data = MQN_norm %>% filter(CHROM %in% "ssa09", top_Xpcent_wins_by_frac_gt_2 %in% 1), aes(y = max_xpnsl, x = win_start, colour = "MQN")) +
  geom_smooth(data = FQN_norm %>% filter(CHROM %in% "ssa09",  !top_Xpcent_wins_by_frac_gt_2 %in% 1), aes(y = max_xpnsl, x = win_start, colour = "FQN"), method = "loess", span = 0.03, se = F) + 
  geom_smooth(data = FNL_norm %>% filter(CHROM %in% "ssa09",  !top_Xpcent_wins_by_frac_gt_2 %in% 1), aes(y = max_xpnsl, x = win_start, colour = "FNL"), method = "loess", span = 0.03, se = F) +
  geom_smooth(data = MNL_norm %>% filter(CHROM %in% "ssa09",  !top_Xpcent_wins_by_frac_gt_2 %in% 1), aes(y = max_xpnsl, x = win_start, colour = "MNL"), method = "loess", span = 0.03, se = F) +
  geom_smooth(data = MQN_norm %>% filter(CHROM %in% "ssa09",  !top_Xpcent_wins_by_frac_gt_2 %in% 1), aes(y = max_xpnsl, x = win_start, colour = "MQN"), method = "loess", span = 0.03, se = F) +
  geom_vline(xintercept =  24902777) +
  geom_vline(xintercept =  24905552) + theme_classic()

##plot all just to see what's up
ggplot() + geom_point(data = FQN_norm %>% filter(top_Xpcent_wins_by_frac_gt_2 %in% 1), aes(y = max_xpnsl, x = win_start, colour = "FQN")) + 
  geom_point(data = FNL_norm %>% filter(top_Xpcent_wins_by_frac_gt_2 %in% 1), aes(y = max_xpnsl, x = win_start, colour = "FNL")) +
  geom_point(data = MNL_norm %>% filter(top_Xpcent_wins_by_frac_gt_2 %in% 1), aes(y = max_xpnsl, x = win_start, colour = "MNL")) +
  geom_point(data = MQN_norm %>% filter(top_Xpcent_wins_by_frac_gt_2 %in% 1), aes(y = max_xpnsl, x = win_start, colour = "MQN")) +
  geom_smooth(data = FQN_norm %>% filter(!top_Xpcent_wins_by_frac_gt_2 %in% 1), aes(y = max_xpnsl, x = win_start, colour = "FQN"), method = "loess", span = 0.03, se = F) + 
  geom_smooth(data = FNL_norm %>% filter(!top_Xpcent_wins_by_frac_gt_2 %in% 1), aes(y = max_xpnsl, x = win_start, colour = "FNL"), method = "loess", span = 0.03, se = F) +
  geom_smooth(data = MNL_norm %>% filter(!top_Xpcent_wins_by_frac_gt_2 %in% 1), aes(y = max_xpnsl, x = win_start, colour = "MNL"), method = "loess", span = 0.03, se = F) +
  geom_smooth(data = MQN_norm %>% filter(!top_Xpcent_wins_by_frac_gt_2 %in% 1), aes(y = max_xpnsl, x = win_start, colour = "MQN"), method = "loess", span = 0.03, se = F) +
  theme_classic() + facet_wrap(~CHROM, scales = "free_x")


