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
#lsit files
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
#lsit files
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
#lsit files
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
#lsit files
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

#six6 vgll3 etc
SCW_loci <- c("LOC106603085", "phs", "calhm3", "LOC106607716", "LOC106610932", "LOC106611671", "LOC106563476", "scamp1", "ptpa", "LOC106570862", "LOC106598444", "LOC106600468", "LOC106602561", "LOC106602634", "LOC106602724", "acsl1", "LOC106602846", "LOC106602938", "LOC106602949", "atg4a", "LOC106603569", "LOC106603552", "LOC106603548", "LOC106603550", "LOC106603415", "LOC106603634", "LOC106603677", "LOC106603676", "opcml", "LOC106605010", "LOC106605123", "LOC106605325", "LOC106606975", "LOC106606972", "LOC106606978", "LOC106609940", "LOC106610586", "LOC106610702", "LOC106611044", "six6", "LOC106610951", "ypel5", "LOC106611117", "LOC106611159", "LOC106611209", "LOC106611243", "ocrl", "LOC106611537", "LOC106611736", "LOC106611993", "LOC106612445", "pid1", "LOC106613053", "LOC106613255", "scaf4", "LOC106613903", "LOC100194696", "LOC106560912", "LOC106560912", "mppd2", "LOC106562713", "LOC106565279", "LOC106565354", "LOC106566502", "LOC106566519", "LOC106566728", "LOC106566890", "LOC106567043", "LOC106567439", "LOC106567659", "fa76b", "metk1", "LOC106568331", "LOC106568538", "LOC106568853", "LOC106569180", "tep1", "LOC106569569", "LOC106569685", "eomes", "LOC106570041", "LOC106570240", "LOC106570256", "LOC106570261", "ndufs4", "rad51b", "LOC106571738", "LOC106571870", "LOC106571940", "LOC106572154", "LOC106572115", "LOC106573506", "LOC106577101", "LOC106577101", "LOC106581384", "LOC106581613", "fgf14", "LOC106581859", "LOC106582292", "LOC106582431", "LOC106582845", "LOC106584946", "LOC106585937", "LOC106586064", "LOC106586268", "LOC106586269", "LOC106586317", "LOC106586310", "LOC106586456", "LOC106586520", "LOC106586514", "LOC106586625", "LOC106586653", "LOC106586753", "LOC106586749", "LOC106586824", "LOC106586862", "nrxn1", "LOC106590053", "LOC106590047")
MNL_norm_sig_gene_overlap %>% filter(gene %in% c("magi2", "six6", "vgll3", "picalma","LOC106563683",  "akap11", "ndufs4", "rora","cntn4", "LOC106581831", "LOC106581727", "LOC106581830", "tead3"))
FNL_norm_sig_gene_overlap %>% filter(gene %in% c("magi2", "six6", "vgll3", "picalma","LOC106563683",  "akap11", "ndufs4", "rora","cntn4", "LOC106581831", "LOC106581727", "LOC106581830", "tead3"))
MQN_norm_sig_gene_overlap %>% filter(gene %in% c("magi2", "six6", "vgll3", "picalma","LOC106563683",  "akap11", "ndufs4", "rora","cntn4", "LOC106581831", "LOC106581727", "LOC106581830","tead3"))
FQN_norm_sig_gene_overlap %>% filter(gene %in% c("magi2", "six6", "vgll3", "picalma","LOC106563683",  "akap11", "ndufs4", "rora","cntn4", "LOC106581831", "LOC106581727", "LOC106581830","tead3"))


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


#get top loci for each group, by GWAS score - get top 1000 first to filter beagle via command line, can downsample in R from there

All_PCcor_top100_LRT_genes <- PC_FST_GWA_PCcor_MSW_1SW_gene_overlap %>%
  group_by(gene) %>%
  filter(LRT == max(LRT))  
All_PCcor_top100_LRT_genes <- data.frame(All_PCcor_top100_LRT_genes) %>%  slice_max(LRT, n = 100)
sig_ALL <- PC_FST_GWA_PCcor_MSW_1SW_gene_overlap %>% filter(GWAS_PCcor_q < 0.05)
Sig_All_gene <- inner_join(sig_ALL, All_PCcor_top100_LRT_genes)
fwrite(Sig_All_gene, "Sig_All_gene", col.names = T, row.names = F, sep = "\t", quote = F)


FQN_PC_FST_GWA_PCcorr_topLRT_genes <- FQN_PC_FST_GWA_PCcorr_gene_overlap %>% 
  group_by(gene) %>%
  filter(LRT == max(LRT)) 
FQN_sig <- FQN_PC_FST_GWA_PCcorr_gene_overlap %>% filter(GWAS_PCcor_q < 0.05) %>%  distinct(gene)

FQN_top100_LRT <- data.frame(FQN_PC_FST_GWA_PCcorr_topLRT_genes) %>%  slice_max(LRT, n = 100)


MQN_PC_FST_GWA_PCcorr_topLRT_genes <- MQN_PC_FST_GWA_PCcorr_gene_overlap %>% 
  group_by(gene) %>%
  filter(LRT == max(LRT)) 
MQN_top100_LRT <-  data.frame(MQN_PC_FST_GWA_PCcorr_topLRT_genes) %>%  slice_max(LRT, n = 100)


FNL_PC_FST_GWA_PCcorr_topLRT_genes <- FNL_PC_FST_GWA_PCcorr_gene_overlap %>% 
  group_by(gene) %>%
  filter(LRT == max(LRT)) 
FNL_top100_LRT <-  data.frame(FNL_PC_FST_GWA_PCcorr_topLRT_genes) %>%  slice_max(LRT, n = 100)


MNL_PC_FST_GWA_PCcorr_topLRT_genes <- MNL_PC_FST_GWA_PCcorr_gene_overlap %>% 
  group_by(gene) %>%
  filter(LRT == max(LRT)) 
MNL_top100_LRT <-  data.frame(MNL_PC_FST_GWA_PCcorr_topLRT_genes) %>%  slice_max(LRT, n = 100)



#check gene overlap
MNL_LRT100_gene <- MNL_PC_FST_GWA_PCcorr_gene_overlap  %>%  slice_max(LRT, n = 100)
FNL_LRT100_gene <- FNL_PC_FST_GWA_PCcorr_gene_overlap  %>%  slice_max(LRT, n = 100)
MQN_LRT100_gene <- MQN_PC_FST_GWA_PCcorr_gene_overlap  %>%  slice_max(LRT, n = 100)
FQN_LRT100_gene <- FQN_PC_FST_GWA_PCcorr_gene_overlap  %>%  slice_max(LRT, n = 100)

FNL_LRT100_gene$gene %in% FQN_LRT100_gene$gene


MNL_LRT100_gene[MNL_LRT100_gene$gene %in%  MNL_norm_sig_gene_overlap$gene, ] %>%  distinct(gene, .keep_all = T)
FNL_LRT100_gene[FNL_LRT100_gene$gene %in%  FNL_norm_sig_gene_overlap$gene, ] %>%  distinct(gene, .keep_all = T)
MQN_LRT100_gene[MQN_LRT100_gene$gene %in%  MQN_norm_sig_gene_overlap$gene, ] %>%  distinct(gene, .keep_all = T)
FQN_LRT100_gene[FQN_LRT100_gene$gene %in%  FQN_norm_sig_gene_overlap$gene, ] %>%  distinct(gene, .keep_all = T)

MNL_LRT100_gene[MNL_LRT100_gene$gene %in%  MNL_norm_sig_gene_overlap$gene, ] %>%  distinct(gene, .keep_all = T) %>% filter(gene %in% c("magi2", "six6", "vgll3", "picalma", "akap11", "ndufs4", "rora","cntn4", "LOC106581831", "LOC106581727", "LOC106581830", "tead3"))
FNL_LRT100_gene[FNL_LRT100_gene$gene %in%  FNL_norm_sig_gene_overlap$gene, ] %>%  distinct(gene, .keep_all = T) %>% filter(gene %in% c("magi2", "six6", "vgll3", "picalma", "akap11", "ndufs4", "rora","cntn4", "LOC106581831", "LOC106581727", "LOC106581830", "tead3"))
MQN_LRT100_gene[MQN_LRT100_gene$gene %in%  MQN_norm_sig_gene_overlap$gene, ] %>%  distinct(gene, .keep_all = T) %>% filter(gene %in% c("magi2", "six6", "vgll3", "picalma", "akap11", "ndufs4", "rora","cntn4", "LOC106581831", "LOC106581727", "LOC106581830","tead3"))
FQN_LRT100_gene[FQN_LRT100_gene$gene %in%  FQN_norm_sig_gene_overlap$gene, ] %>%  distinct(gene, .keep_all = T) %>% filter(gene %in% c("magi2", "six6", "vgll3", "picalma", "akap11", "ndufs4", "rora","cntn4", "LOC106581831", "LOC106581727", "LOC106581830","tead3"))



 
MNL_LRT_XPnsl_OL_bed <-  MNL_LRT100_gene[MNL_LRT100_gene$gene %in%  MNL_norm_sig_gene_overlap$gene, ] %>% select(CHROM, POS, BPtoo, gene )  
FNL_LRT_XPnsl_OL_bed <-  FNL_LRT100_gene[FNL_LRT100_gene$gene %in%  FNL_norm_sig_gene_overlap$gene, ] %>% select(CHROM, POS, BPtoo, gene )  
MQN_LRT_XPnsl_OL_bed <-  MQN_LRT100_gene[MQN_LRT100_gene$gene %in%  MQN_norm_sig_gene_overlap$gene, ] %>% select(CHROM, POS, BPtoo, gene )  
FQN_LRT_XPnsl_OL_bed <-  FQN_LRT100_gene[FQN_LRT100_gene$gene %in%  FQN_norm_sig_gene_overlap$gene, ] %>% select(CHROM, POS, BPtoo, gene ) 


fwrite(MNL_LRT_XPnsl_OL_bed, "MNL_LRT_XPnsl_OL.bed", col.names = F, row.names = F, sep = "\t", quote = F)
fwrite(FNL_LRT_XPnsl_OL_bed, "FNL_LRT_XPnsl_OL.bed", col.names = F, row.names = F, sep = "\t", quote = F)
fwrite(MQN_LRT_XPnsl_OL_bed, "MQN_LRT_XPnsl_OL.bed", col.names = F, row.names = F, sep = "\t", quote = F)
fwrite(FQN_LRT_XPnsl_OL_bed, "FQN_LRT_XPnsl_OL.bed", col.names = F, row.names = F, sep = "\t", quote = F)

system("sed -i.bun 's/ssa/chr/' *LRT_XPnsl_OL.bed")
system("mv *LRT_XPnsl_OL.bed ~/Desktop/Working/Salmon_GO_terms")


MNL_LRT_XPnsl_OL <-  inner_join(MNL_LRT100_gene, MNL_norm_sig_gene_overlap, by = "gene")
FNL_LRT_XPnsl_OL_bed <-  FNL_LRT100_gene[FNL_LRT100_gene$gene %in%  FNL_norm_sig_gene_overlap$gene, ] %>% select(CHROM, POS, BPtoo, gene )  
MQN_LRT_XPnsl_OL_bed <-  MQN_LRT100_gene[MQN_LRT100_gene$gene %in%  MQN_norm_sig_gene_overlap$gene, ] %>% select(CHROM, POS, BPtoo, gene )  
FQN_LRT_XPnsl_OL_bed <-  FQN_LRT100_gene[FQN_LRT100_gene$gene %in%  FQN_norm_sig_gene_overlap$gene, ] %>% select(CHROM, POS, BPtoo, gene ) 

  geom_vline(xintercept =  28659019) + theme_classic()

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


ggplot() + geom_point(data = FQN_norm %>% filter(top_Xpcent_wins_by_frac_gt_2 %in% 1), aes(y = max_xpnsl, x = win_start, colour = "FQN")) + 
  geom_point(data = FNL_norm %>% filter(top_Xpcent_wins_by_frac_gt_2 %in% 1), aes(y = max_xpnsl, x = win_start, colour = "FNL")) +
  geom_point(data = MNL_norm %>% filter(top_Xpcent_wins_by_frac_gt_2 %in% 1), aes(y = max_xpnsl, x = win_start, colour = "MNL")) +
  geom_point(data = MQN_norm %>% filter(top_Xpcent_wins_by_frac_gt_2 %in% 1), aes(y = max_xpnsl, x = win_start, colour = "MQN")) +
  geom_smooth(data = FQN_norm %>% filter(!top_Xpcent_wins_by_frac_gt_2 %in% 1), aes(y = max_xpnsl, x = win_start, colour = "FQN"), method = "loess", span = 0.03, se = F) + 
  geom_smooth(data = FNL_norm %>% filter(!top_Xpcent_wins_by_frac_gt_2 %in% 1), aes(y = max_xpnsl, x = win_start, colour = "FNL"), method = "loess", span = 0.03, se = F) +
  geom_smooth(data = MNL_norm %>% filter(!top_Xpcent_wins_by_frac_gt_2 %in% 1), aes(y = max_xpnsl, x = win_start, colour = "MNL"), method = "loess", span = 0.03, se = F) +
  geom_smooth(data = MQN_norm %>% filter(!top_Xpcent_wins_by_frac_gt_2 %in% 1), aes(y = max_xpnsl, x = win_start, colour = "MQN"), method = "loess", span = 0.03, se = F) +
  theme_classic() + facet_wrap(~CHROM, scales = "free_x")


MNL_LRT100_gene <- MNL_PC_FST_GWA_PCcorr_gene_overlap  %>%  slice_max(LRT, n = 100)
FNL_LRT100_gene <- FNL_PC_FST_GWA_PCcorr_gene_overlap  %>%  slice_max(LRT, n = 100)
MQN_LRT100_gene <- MQN_PC_FST_GWA_PCcorr_gene_overlap  %>%  slice_max(LRT, n = 100)
FQN_LRT100_gene <- FQN_PC_FST_GWA_PCcorr_gene_overlap  %>%  slice_max(LRT, n = 100)


 write.table(unique(MNL_LRT100_gene$gene), "MNL_LRT100_gene.txt", col.names = F, row.names = F, sep = "\t", quote = F) 
 write.table(unique(FNL_LRT100_gene$gene), "FNL_LRT100_gene.txt", col.names = F, row.names = F, sep = "\t", quote = F) 
 write.table(unique(MQN_LRT100_gene$gene), "MQN_LRT100_gene.txt", col.names = F, row.names = F, sep = "\t", quote = F) 
 write.table(unique(FQN_LRT100_gene$gene), "FQN_LRT100_gene.txt", col.names = F, row.names = F, sep = "\t", quote = F) 

quantile(FNL_PC_FST_GWA_PCcorr_gene_overlap$FST, 0.99)


#revigo plot

pal <- wes_palette("Zissou1", 100, type = "continuous")
GO_plot <- function(REVIGO) {
  GO <- REVIGO
  GO$PlotX <- as.numeric(GO$PlotX) 
  GO$PlotY <- as.numeric(GO$PlotY) 
  
  PLOT <- ggplot(GO) + geom_point(aes(x = PlotX, y = PlotY, colour = Uniqueness, size = LogSize)) +
    labs (y = "semantic space x", x = "semantic space y") + 
    geom_text(aes(PlotX, PlotY, label = Name), colour = I(alpha("black", 0.85)), size = 3 ) + theme_classic() 
  PLOT +  scale_colour_gradientn(colours = pal) }


All_GO <- fread("Revigo_tiny.csv")

GO_plot(All_GO)

FQN_GO$Name[FQN_GO$Name %in% MQN_GO$Name]

