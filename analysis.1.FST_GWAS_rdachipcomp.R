library(RcppCNPy)
library(tidyverse)
library(data.table)
library(qvalue)

setwd(dir = "~/Desktop/Working/Sea_Age/")
#Read in metadata
#write IDs for running ANGSD genotyping and GL
Salmon_Metadata_Bam_Unique_HD  <- fread("Salmon_Metadata_Bam_Unique_HD.txt")
#Recode
Salmon_Metadata_Bam_Unique_HD$AgeClass[Salmon_Metadata_Bam_Unique_HD$AgeClass=="1SW"] <- 0
Salmon_Metadata_Bam_Unique_HD$AgeClass[Salmon_Metadata_Bam_Unique_HD$AgeClass=="MSW"] <- 1

QCNB <- c("MSW", "MUN", "RTG", "TRI")
NL <- c("ENG", "MBB", "SAN", "RVF")

#Separate by region and sex for 
Male_1SW_MSW_QN <- Salmon_Metadata_Bam_Unique_HD %>% filter(PCR_Sex =="M", RiverCode %in% QCNB)
Male_1SW_MSW_QN %>% group_by(SiteName) %>%  count(AgeClass)

Male_1SW_MSW_NL <-  Salmon_Metadata_Bam_Unique_HD %>% filter(PCR_Sex =="M", RiverCode %in% NL)
Male_1SW_MSW_NL %>% group_by(SiteName) %>%  count(AgeClass)


Female_1SW_MSW_QN <- Salmon_Metadata_Bam_Unique_HD %>% filter(PCR_Sex =="F", RiverCode %in% QCNB)
Female_1SW_MSW_QN %>% group_by(SiteName) %>%  count(AgeClass)

Female_1SW_MSW_NL <- Salmon_Metadata_Bam_Unique_HD %>% filter(PCR_Sex =="F", RiverCode %in% NL)
Female_1SW_MSW_NL %>% group_by(SiteName) %>%  count(AgeClass)

#pop structure first

covmat <-  read.table("Salmo_CIGENE_80geno_HD.cov", quote="\"",
                      comment.char="", stringsAsFactors=FALSE)
PCA <- eigen(covmat)

#var explained per axis
eigenvalues <- PCA$values
(eigenvalues[1]/sum(eigenvalues))*100
(eigenvalues[2]/sum(eigenvalues))*100
(eigenvalues[3]/sum(eigenvalues))*100

PCAscores <- data.frame(cbind(Salmon_Metadata_Bam_Unique_HD,
                              PCA$vectors[,1:10]), stringsAsFactors = F)

varpart(as.numeric(PCAscores$AgeClass),  ~PCAscores$V1, ~PCAscores$V2, ~PCAscores$V3)
AC_rda <- rda(as.numeric(PCAscores$AgeClass)  ~PCAscores$V1 + PCAscores$V2 + PCAscores$V3)
anova.cca(AC_rda, parallel= 16,by = "terms")
RsquareAdj(AC_rda)

ggplot() + geom_point(data = PCAscores, aes(x = V1, y = V2, colour
                                            =AgeClass))
ggplot() + geom_point(data = PCAscores, aes(x = V1, y = V2, colour
                                            =RiverCode))+ theme_classic()

ggplot() + geom_point(data =  PCAscores, aes(x = V3, y = V4, colour
                                             =AgeClass))
ggplot() + geom_point(data =  PCAscores, aes(x = V4, y = V5, colour
                                             =RiverCode)) + theme_classic()

#Get PCA axes for covariates in GWAS
write.table(PCAscores[7:8], "All_HD_PCA_12.txt", col.names = F, row.names = F, sep = "\t", quote = F)

###Make an admixture plot too because wooohoo
All_admix <- fread("Salmo_CIGENE_80geno_HD.admix.6.Q")
Salmo_admix_meta <- data.frame(cbind(Salmon_Metadata_Bam_Unique_HD, All_admix))
colnames(Salmo_admix_meta)
Admix_table <- Salmo_admix_meta
rownames(Admix_table ) <- Salmo_admix_meta$FGL_ID

plot_data <-  Admix_table  %>% 
  gather('pop', 'prob', V1:V6) %>% 
  group_by(RiverCode)

ggplot(plot_data, aes(FGL_ID, prob, fill = pop)) +
  geom_col() +
  facet_grid(~RiverCode, scales = 'free', space = 'free')

#Read in PCA results per region now
#F_QN

F_QN_covmat <-  read.table("Salmo_CIGENE_80geno_F_QN.cov", quote="\"",
                           comment.char="", stringsAsFactors=FALSE)
F_QN_PCA <- eigen(F_QN_covmat)
F_QN_PCAscores <- data.frame(cbind(Female_1SW_MSW_QN,
                                   F_QN_PCA$vectors[,1:10]), stringsAsFactors = F)

F_QN_eigenvalues <- F_QN_PCA$values
(F_QN_eigenvalues[1]/sum(F_QN_eigenvalues))*100
(F_QN_eigenvalues[2]/sum(F_QN_eigenvalues))*100

ggplot() + geom_point(data = F_QN_PCAscores, aes(x = V1, y = V2, colour
                                                 =AgeClass))
ggplot() + geom_point(data = F_QN_PCAscores, aes(x = V1, y = V2, colour
                                                 =RiverCode))+ theme_classic()

write.table(F_QN_PCAscores[7:8], "F_QN_PCA_12.txt", col.names = F, row.names = F, sep = "\t", quote = F)

#F_NL_covmat 
F_NL_covmat <-  read.table("Salmo_CIGENE_80geno_F_NL.cov", quote="\"",
                           comment.char="", stringsAsFactors=FALSE)
F_NL_PCA <- eigen(F_NL_covmat)
F_NL_PCAscores <- data.frame(cbind(Female_1SW_MSW_NL,
                                   F_NL_PCA$vectors[,1:10]), stringsAsFactors = F)

F_NL_eigenvalues <- F_NL_PCA$values
(F_NL_eigenvalues[1]/sum(F_NL_eigenvalues))*100
(F_NL_eigenvalues[2]/sum(F_NL_eigenvalues))*100

ggplot() + geom_point(data = F_NL_PCAscores, aes(x = V1, y = V2, colour
                                                 =AgeClass))
ggplot() + geom_point(data = F_NL_PCAscores, aes(x = V1, y = V2, colour
                                                 =RiverCode))+ theme_classic()

write.table(F_NL_PCAscores[7:8], "F_NL_PCA_12.txt", col.names = F, row.names = F, sep = "\t", quote = F)

#M QN
M_QN_covmat <-  read.table("Salmo_CIGENE_80geno_M_QN.cov", quote="\"",
                           comment.char="", stringsAsFactors=FALSE)
M_QN_PCA <- eigen(M_QN_covmat)
M_QN_PCAscores <- data.frame(cbind(Male_1SW_MSW_QN,
                                   M_QN_PCA$vectors[,1:10]), stringsAsFactors = F)

M_QN_eigenvalues <- M_QN_PCA$values
(M_QN_eigenvalues[1]/sum(M_QN_eigenvalues))*100
(M_QN_eigenvalues[2]/sum(M_QN_eigenvalues))*100

ggplot() + geom_point(data = M_QN_PCAscores, aes(x = V1, y = V2, colour
                                                 =AgeClass))
ggplot() + geom_point(data = M_QN_PCAscores, aes(x = V1, y = V2, colour
                                                 =RiverCode))+ theme_classic()

write.table(M_QN_PCAscores[7:8], "M_QN_PCA_12.txt", col.names = F, row.names = F, sep = "\t", quote = F)

#M NL
M_NL_covmat <-  read.table("Salmo_CIGENE_80geno_M_NL.cov", quote="\"",
                           comment.char="", stringsAsFactors=FALSE)
M_NL_PCA <- eigen(M_NL_covmat)
M_NL_PCAscores <- data.frame(cbind(Male_1SW_MSW_NL,
                                   M_NL_PCA$vectors[,1:10]), stringsAsFactors = F)

M_NL_eigenvalues <- M_NL_PCA$values
(M_NL_eigenvalues[1]/sum(M_NL_eigenvalues))*100
(M_NL_eigenvalues[2]/sum(M_NL_eigenvalues))*100

ggplot() + geom_point(data = M_NL_PCAscores, aes(x = V1, y = V2, colour
                                                 =AgeClass))
ggplot() + geom_point(data = M_NL_PCAscores, aes(x = V1, y = V2, colour
                                                 =RiverCode))+ theme_classic()

write.table(M_NL_PCAscores[7:8], "M_NL_PCA_12.txt", col.names = F, row.names = F, sep = "\t", quote = F)

#PCANGSD Selscan
#import pcangsd fastPCA loadings per site
All_PCLoadings <- data.frame(npyLoad("Salmo_CIGENE_80geno_HD.selection.npy")) 

#import pcangsd sites output
All_sites <- fread("Salmo_CIGENE_80geno_HD.sites", data.table = F,
                   stringsAsFactors = F, header = F)

Loc_pos <- fread("Beagle_loc_pos")
setnames(Loc_pos, new = c("CHROM", "POS"))

All_Loc_pos <- bind_cols(Loc_pos, All_sites)
All_Loc_pos <- All_Loc_pos %>% filter(V1 %in% "1")

All_sites_pos <- bind_cols(All_Loc_pos, All_PCLoadings)

colnames(All_sites_pos)[3:8] <- c("KEEP", "PC1", "PC2", "PC3", "PC4", "PC5")
All_sites_pos$pvals_PC1 <- 1- pchisq(q =All_sites_pos$PC1, df = 1)
qvals <- qvalue(All_sites_pos$pvals_PC1)
All_sites_pos$qvals_PC1 <- qvals$qvalues
All_sites_sigsites_PC1 <- All_sites_pos %>%  filter(qvals_PC1 < 0.01)

All_sites_six6 <-  All_sites_pos %>% filter(CHROM %in% "ssa09", POS > 24902777, POS <24905552)
All_sites_vgll3 <-  All_sites_pos %>% filter(CHROM %in% "ssa25", POS > 28654947, POS <28659019)

All_sites_six6 <- All_sites_six6 %>%  mutate(gene = "six6")
All_sites_vgll3 <- All_sites_vgll3 %>%  mutate(gene = "vgll3")

six6_vgll3_WG <- bind_rows(All_sites_six6, All_sites_vgll3) %>%  
  mutate(SNP = paste0(CHROM, "_", POS)) %>% 
  select(SNP, gene)
  

ggplot() + geom_point(data = All_sites_sigsites_PC1, aes(x = POS, y =-log10(qvals_PC1))) + facet_wrap(~CHROM, scales = "free_x") + 
  theme_classic()  + geom_point(data = All_sites_six6 %>%  filter(qvals_PC1 < 0.01), aes(x = POS, y =-log10(qvals_PC1), colour = "six6"))

All_sites_sigsites_PC1$SNP <- paste0(All_sites_sigsites_PC1$CHROM, "_", All_sites_sigsites_PC1$POS)

ALL_PCload <- ggman(gwas = All_sites_sigsites_PC1, snp = "SNP", bp = "POS", chrom = "CHROM", pvalue = "qvals_PC1", 
      xlabel = "Chromosome", ylabel = "-log10(qvalue)", pointSize = 1, 
      logTransform = T) + theme_classic()

ggmanHighlightGroup(ALL_PCload, highlightDfm = six6_vgll3_WG, snp = "SNP", group = "gene")

All_sites_sigsites_PC1_bed <- All_sites_sigsites_PC1 %>%  select(CHROM, POS) %>% 
  mutate(POSJR = POS + 1)


#Check SV overlap
All_sites_sigsites_PC1_bed %>%  filter(CHROM %in% "ssa01", POS > 44000000,  POS < 53000000)
All_sites_sigsites_PC1_bed  %>%  filter(CHROM %in% "ssa23",POS < 8000000)



fwrite(All_sites_sigsites_PC1_bed, "All_sites_sigsites_PC1.bed", sep = "\t", quote = F, row.names = F, col.names = F)
system("bedtools intersect -a SSA_GENES.bed -b All_sites_sigsites_PC1.bed > All_sites_sigsites_PC1_genes.txt ")

distinct_PC1_WG <- fread("../All_sites_sigsites_PC1_genes.txt") %>%  select(V4) %>%  distinct
fwrite(distinct_PC1_WG, "distinct_PC1_WG", col.names = F, row.names = F, sep = "\t", quote = F)

#F_QN_PCA
#import pcangsd fastPCA loadings per site
F_QN_PCLoadings <- data.frame(npyLoad("Salmo_CIGENE_80geno_F_QN.selection.npy")) 

#import pcangsd sites output

F_QN_sites <- fread("Salmo_CIGENE_80geno_F_QN.sites", data.table = F,
                    stringsAsFactors = F, header = F)
F_QN_pos <- bind_cols(Loc_pos,F_QN_sites)
F_QN_pos <- F_QN_pos %>% filter(V1 %in% "1")

F_QN_pos_loadings <- bind_cols(F_QN_pos, F_QN_PCLoadings)

colnames(F_QN_pos_loadings)[3:4] <- c("KEEP", "PC1")
F_QN_pos_loadings$pvals_PC1 <- 1- pchisq(q =F_QN_pos_loadings$PC1, df = 1)
qvals <- qvalue(F_QN_pos_loadings$pvals_PC1)
F_QN_pos_loadings$qvals_PC1 <- qvals$qvalues
F_QN_sigsites_PC1 <- F_QN_pos_loadings %>%  filter(qvals_PC1 < 0.01)

ggplot() + geom_point(data = F_QN_sigsites_PC1, aes(x = POS, y =-log10(qvals_PC1))) + facet_wrap(~CHROM, scales = "free_x") + 
  theme_classic()


#F_NL

F_NL_PCLoadings <- data.frame(npyLoad("Salmo_CIGENE_80geno_F_NL.selection.npy")) 

#import pcangsd sites output

F_NL_sites <- fread("Salmo_CIGENE_80geno_F_NL.sites", data.table = F,
                    stringsAsFactors = F, header = F)
F_NL_pos <- bind_cols(Loc_pos,F_NL_sites)
F_NL_pos <- F_NL_pos %>% filter(V1 %in% "1")

F_NL_pos_loadings <- bind_cols(F_NL_pos, F_NL_PCLoadings)

colnames(F_NL_pos_loadings)[3:4] <- c("KEEP", "PC1")
F_NL_pos_loadings$pvals_PC1 <- 1- pchisq(q =F_NL_pos_loadings$PC1, df = 1)
qvals <- qvalue(F_NL_pos_loadings$pvals_PC1)
F_NL_pos_loadings$qvals_PC1 <- qvals$qvalues
F_NL_sigsites_PC1 <- F_NL_pos_loadings %>%  filter(qvals_PC1 < 0.01)


ggplot() + geom_point(data = F_NL_sigsites_PC1, aes(x = POS, y =-log10(qvals_PC1))) + facet_wrap(~CHROM, scales = "free_x") + 
  theme_classic()


#M_QN

M_QN_PCLoadings <- data.frame(npyLoad("Salmo_CIGENE_80geno_M_QN.selection.npy")) 

#import pcangsd sites output

M_QN_sites <- fread("Salmo_CIGENE_80geno_M_QN.sites", data.table = F,
                    stringsAsFactors = F, header = F)
M_QN_pos <- bind_cols(Loc_pos,M_QN_sites)
M_QN_pos <- M_QN_pos %>% filter(V1 %in% "1")

M_QN_pos_loadings <- bind_cols(M_QN_pos, M_QN_PCLoadings)

colnames(M_QN_pos_loadings)[3:4] <- c("KEEP", "PC1")
M_QN_pos_loadings$pvals_PC1 <- 1- pchisq(q =M_QN_pos_loadings$PC1, df = 1)
qvals <- qvalue(M_QN_pos_loadings$pvals_PC1)
M_QN_pos_loadings$qvals_PC1 <- qvals$qvalues
M_QN_sigsites_PC1 <- M_QN_pos_loadings %>%  filter(qvals_PC1 < 0.01)

ggplot() + geom_point(data = M_QN_sigsites_PC1, aes(x = POS, y =-log10(qvals_PC1))) + facet_wrap(~CHROM, scales = "free_x") + 
  theme_classic()


#M_NL

M_NL_PCLoadings <- data.frame(npyLoad("Salmo_CIGENE_80geno_M_NL.selection.npy")) 

#import pcangsd sites output

M_NL_sites <- fread("Salmo_CIGENE_80geno_M_NL.sites", data.table = F,
                    stringsAsFactors = F, header = F)
M_NL_pos <- bind_cols(Loc_pos,M_NL_sites)
M_NL_pos <- M_NL_pos %>% filter(V1 %in% "1")

M_NL_pos_loadings <- bind_cols(M_NL_pos, M_NL_PCLoadings)

colnames(M_NL_pos_loadings)[3:4] <- c("KEEP", "PC1")
M_NL_pos_loadings$pvals_PC1 <- 1- pchisq(q =M_NL_pos_loadings$PC1, df = 1)
qvals <- qvalue(M_NL_pos_loadings$pvals_PC1)
M_NL_pos_loadings$qvals_PC1 <- qvals$qvalues
M_NL_sigsites_PC1 <- M_NL_pos_loadings %>%  filter(qvals_PC1 < 0.01)


ggplot() + geom_point(data = M_NL_sigsites_PC1, aes(x = POS, y =-log10(qvals_PC1))) + facet_wrap(~CHROM, scales = "free_x") + 
  theme_classic()


#FST comparisons 

MSW_1SW_FST <- fread("MSW_1SW_FST")
colnames(MSW_1SW_FST)[3] <- "FST"
MSW_1SW_FST$FST <- as.numeric(MSW_1SW_FST$FST)
MSW_1SW_FST <- MSW_1SW_FST %>% drop_na()
MSW_1SW_FST$POS <- as.integer(MSW_1SW_FST$POS)

PC_FST_MSW_1SW<- inner_join(MSW_1SW_FST, All_sites_pos)

#F_QN
F_QN_FST <- fread("F_QN_FST")
colnames(F_QN_FST)[3] <- "FST"
F_QN_FST$FST <- as.numeric(F_QN_FST$FST)
F_QN_FST <- F_QN_FST %>% drop_na()
F_QN_FST$POS <- as.integer(F_QN_FST$POS)
PC_FST_F_QN<- inner_join(F_QN_FST, F_QN_pos_loadings)

#F_NL
F_NL_FST <- fread("F_NL_FST")
colnames(F_NL_FST)[3] <- "FST"
F_NL_FST$FST <- as.numeric(F_NL_FST$FST)
F_NL_FST <- F_NL_FST %>% drop_na()
F_NL_FST$POS <- as.integer(F_NL_FST$POS)
PC_FST_F_NL<- inner_join(F_NL_FST, F_NL_pos_loadings)

#M_QN
M_QN_FST <- fread("M_QN_FST")
colnames(M_QN_FST)[3] <- "FST"
M_QN_FST$FST <- as.numeric(M_QN_FST$FST)
M_QN_FST <- M_QN_FST %>% drop_na()
M_QN_FST$POS <- as.integer(M_QN_FST$POS)
PC_FST_M_QN<- inner_join(M_QN_FST, M_QN_pos_loadings)

#M_NL
M_NL_FST <- fread("M_NL_FST")
colnames(M_NL_FST)[3] <- "FST"
M_NL_FST$FST <- as.numeric(M_NL_FST$FST)
M_NL_FST <- M_NL_FST %>% drop_na()
M_NL_FST$POS <- as.integer(M_NL_FST$POS)
PC_FST_M_NL<- inner_join(M_NL_FST, M_NL_pos_loadings)

#GWAS 
setwd("All_nocor/")

All_nocorr<- 
  list.files(pattern = "All_hd.*lrt.*gz") %>% 
  map_df(~fread(.))

colnames(All_nocorr)[1:2] <- c("CHROM", "POS")

All_nocorr <- All_nocorr[!(All_nocorr$LRT %in% "-999"),]
All_nocorr$GWAS_nocor_pvals <- 1- pchisq(q = All_nocorr$LRT, df = 1)
PC_FST_GWA_MSW_1SW<- inner_join(PC_FST_MSW_1SW, All_nocorr)

allnocorrrq <- qvalue(PC_FST_GWA_MSW_1SW$GWAS_nocor_pvals)
PC_FST_GWA_MSW_1SW$GWAS_uncor_q <- allnocorrrq$qvalues

fwrite(PC_FST_GWA_MSW_1SW, "PC_FST_GWA_MSW_1SW.tsv", col.names = T, row.names = F, sep = "\t", quote = F)


#comparisons of FST, frequency, etc.
ggplot() + geom_point(data =PC_FST_GWA_MSW_1SW %>%  sample_n(50000), aes(y = FST, x = -log10(GWAS_uncor_q))) + theme_classic()
ggplot() + geom_point(data =PC_FST_GWA_MSW_1SW %>%  sample_n(50000), aes(y = FST, x = LRT)) + theme_classic()
ggplot() + geom_point(data =PC_FST_GWA_MSW_1SW %>%  sample_n(50000), aes(y = FST, x = beta)) + theme_classic()
ggplot() + geom_point(data =PC_FST_GWA_MSW_1SW %>% sample_n(50000), aes(x = Frequency, y = beta)) + theme_classic()
cor.test(PC_FST_GWA_MSW_1SW$FST, PC_FST_GWA_MSW_1SW$LRT)
#PC corrected
setwd("../All_2PCcorr/")
All_PCcorr<- 
  list.files(pattern = "All_hd.*lrt.*gz") %>% 
  map_df(~fread(.))

colnames(All_PCcorr)[1:2] <- c("CHROM", "POS")

All_PCcorr <- All_PCcorr[!(All_PCcorr$LRT %in% "-999"),]
All_PCcorr$GWAS_PCcor_pvals <- 1- pchisq(q = All_PCcorr$LRT, df = 1)
PC_FST_GWA_PCcor_MSW_1SW<- inner_join(PC_FST_MSW_1SW, All_PCcorr)

allcorrrq <- qvalue(PC_FST_GWA_PCcor_MSW_1SW$GWAS_PCcor_pvals)

PC_FST_GWA_PCcor_MSW_1SW$GWAS_PCcor_q <- allcorrrq$qvalues

fwrite(PC_FST_GWA_PCcor_MSW_1SW, "PC_FST_GWA_PCcor_MSW_1SW.tsv", col.names = T, row.names = F, sep = "\t", quote = F)

PC_FST_GWA_PCcor_MSW_1SW$POSEND <- PC_FST_GWA_PCcor_MSW_1SW$POS + 1 
write.table(cbind(PC_FST_GWA_PCcor_MSW_1SW$CHROM, as.integer(PC_FST_GWA_PCcor_MSW_1SW$POS), as.integer(PC_FST_GWA_PCcor_MSW_1SW$POSEND)), "PC_FST_GWA_PCcor_MSW_1SW.bed", sep = "\t",
            col.names = F, row.names = F, quote = F)

system("bedtools intersect -b PC_FST_GWA_PCcor_MSW_1SW.bed -a ../SSA_genes.bed -wb > PC_FST_GWA_PCcor_MSW_1SW_gene_overlap")
PC_FST_GWA_PCcor_MSW_1SW_gene_overlap <- fread("PC_FST_GWA_PCcor_MSW_1SW_gene_overlap") %>% select(V5, V6, V7, V4)
colnames(PC_FST_GWA_PCcor_MSW_1SW_gene_overlap) <- c("CHROM", "POS", "BPtoo", "gene")

PC_FST_GWA_PCcor_MSW1SW_gene_overlap_stats <- inner_join(PC_FST_GWA_PCcor_MSW_1SW_gene_overlap, PC_FST_GWA_PCcor_MSW_1SW)
fwrite(PC_FST_GWA_PCcor_MSW1SW_gene_overlap_stats, "PC_FST_GWA_PCcor_MSW1SW_gene_overlap_stats.tsv", col.names = T, row.names = F,
       quote = F, sep = "\t")

ggplot() + geom_point(data =PC_FST_GWA_PCcor_MSW_1SW %>%  sample_n(50000), aes(y = FST, x = -log10(GWAS_PCcor_q))) + theme_classic()
ggplot() + geom_point(data =PC_FST_GWA_PCcor_MSW_1SW %>%  sample_n(50000), aes(y = FST, x = LRT)) + theme_classic()
ggplot() + geom_point(data =PC_FST_GWA_PCcor_MSW_1SW %>%  sample_n(50000), aes(y = FST, x = beta)) + theme_classic()
ggplot() + geom_point(data =PC_FST_GWA_PCcor_MSW_1SW %>% filter(Frequency <= 0.5) %>% sample_n(50000), aes(x = Frequency, y = beta)) + theme_classic()

cor.test(PC_FST_GWA_PCcor_MSW1SW_gene_overlap_stats$FST, PC_FST_GWA_PCcor_MSW1SW_gene_overlap_stats$LRT)

All_PCcorr_sig <- PC_FST_GWA_PCcor_MSW_1SW %>% filter(GWAS_PCcor_q < 0.05)
All_PCcorr_sig_gene <- inner_join(All_PCcorr_sig, PC_FST_GWA_PCcor_MSW_1SW_gene_overlap_stats)
fwrite(All_PCcorr_sig, "All_PCcorr_sig_gene", col.names = T, row.names = F, sep = "\t", quote = F)

#F_QN PC corrected
setwd("../F_QN_GWAS//")
FQN_PCcorr<- 
  list.files(pattern = "*lrt.*gz") %>% 
  map_df(~fread(.))

colnames(FQN_PCcorr)[1:2] <- c("CHROM", "POS")

FQN_PCcorr<- FQN_PCcorr[!(FQN_PCcorr$LRT %in% "-999"),]
FQN_PCcorr$GWAS_PCcor_pvals <- 1- pchisq(q = FQN_PCcorr$LRT, df = 1)
FQN_PC_FST_GWA_PCcor <- inner_join(PC_FST_F_QN, FQN_PCcorr)

FQN_PCcorrq <- qvalue(FQN_PC_FST_GWA_PCcor$GWAS_PCcor_pvals, pi0 = 1)

FQN_PC_FST_GWA_PCcor$GWAS_PCcor_q <- FQN_PCcorrq$qvalues

fwrite(FQN_PC_FST_GWA_PCcor, "FQN_PC_FST_GWA_PCcorr.tsv", col.names = T, row.names = F, sep = "\t", quote = F)

FQN_PC_FST_GWA_PCcor$POSEND <- FQN_PC_FST_GWA_PCcor$POS + 1 
write.table(cbind(FQN_PC_FST_GWA_PCcor$CHROM, as.integer(FQN_PC_FST_GWA_PCcor$POS), as.integer(FQN_PC_FST_GWA_PCcor$POSEND)), "FQN_PC_FST_GWA_PCcor.bed", sep = "\t",
            col.names = F, row.names = F, quote = F)

system("bedtools intersect -b FQN_PC_FST_GWA_PCcor.bed -a ../SSA_genes.bed -wb > FQN_PC_FST_GWA_PCcor_gene_overlap")
FQN_PC_FST_GWA_PCcor_gene_overlap <- fread("FQN_PC_FST_GWA_PCcor_gene_overlap") %>% select(V5, V6, V7, V4)
colnames(FQN_PC_FST_GWA_PCcor_gene_overlap) <- c("CHROM", "POS", "BPtoo", "gene")

FQN_PC_FST_GWA_PCcor_gene_overlap_stats <- inner_join(FQN_PC_FST_GWA_PCcor_gene_overlap, FQN_PC_FST_GWA_PCcor)
fwrite(FQN_PC_FST_GWA_PCcor_gene_overlap_stats, "FQN_PC_FST_GWA_PCcor_gene_overlap_stats.tsv", col.names = T, row.names = F,
       quote = F, sep = "\t")


ggplot() + geom_point(data =FQN_PC_FST_GWA_PCcor %>%  sample_n(50000), aes(y = FST, x = -log10(GWAS_PCcor_q))) + theme_classic()
ggplot() + geom_point(data =FQN_PC_FST_GWA_PCcor %>%  sample_n(50000), aes(y = FST, x = LRT)) + theme_classic()
ggplot() + geom_point(data =FQN_PC_FST_GWA_PCcor %>%  sample_n(50000), aes(y = LRT, x = beta)) + theme_classic()
ggplot() + geom_point(data =FQN_PC_FST_GWA_PCcor %>% filter(Frequency <= 0.5) %>% sample_n(50000), aes(x = Frequency, y = beta)) + theme_classic()

cor.test(FQN_PC_FST_GWA_PCcor_gene_overlap_stats$LRT, FQN_PC_FST_GWA_PCcor_gene_overlap_stats$FST)

FQN_PC_FST_GWA_PCcor_gene_overlap_stats %>% slice_max(LRT, n =100) %>%  select(gene)
FQN_PC_FST_GWA_PCcor_gene_overlap_stats %>% slice_max(FST, n =100) %>%  select(gene)


FQN_sig <- FQN_PC_FST_GWA_PCcorr %>%  filter(GWAS_PCcor_q < 0.05)
FQN_sig_gene <- inner_join(FQN_sig, FQN_PC_FST_GWA_PCcorr_gene_overlap) %>%  arrange(GWAS_PCcor_q)
fwrite(FQN_sig_gene, "FQN_LRT_sig_gene")
FQN_sig_gene %>%  distinct(gene)
#F_NL PC corrected
setwd("../F_NL_GWAS//")
FNL_PCcorr<- 
  list.files(pattern = "*lrt.*gz") %>% 
  map_df(~fread(.))

colnames(FNL_PCcorr)[1:2] <- c("CHROM", "POS")

FNL_PCcorr<- FNL_PCcorr[!(FNL_PCcorr$LRT %in% "-999"),]
FNL_PCcorr$GWAS_PCcor_pvals <- 1- pchisq(q = FNL_PCcorr$LRT, df = 1)
FNL_PC_FST_GWA_PCcor <- inner_join(PC_FST_F_NL, FNL_PCcorr)

FNL_PCcorrq <- qvalue(FNL_PC_FST_GWA_PCcor$GWAS_PCcor_pvals, pi0 = 1)

FNL_PC_FST_GWA_PCcor$GWAS_PCcor_q <- FNL_PCcorrq$qvalues

fwrite(FNL_PC_FST_GWA_PCcor, "FNL_PC_FST_GWA_PCcorr.tsv", col.names = T, row.names = F, sep = "\t", quote = F)

FNL_PC_FST_GWA_PCcor$POSEND <- FNL_PC_FST_GWA_PCcor$POS + 1 
write.table(cbind(FNL_PC_FST_GWA_PCcor$CHROM, as.integer(FNL_PC_FST_GWA_PCcor$POS), as.integer(FNL_PC_FST_GWA_PCcor$POSEND)), "FNL_PC_FST_GWA_PCcor.bed", sep = "\t",
            col.names = F, row.names = F, quote = F)

system("bedtools intersect -b FNL_PC_FST_GWA_PCcor.bed -a ../SSA_genes.bed -wb > FNL_PC_FST_GWA_PCcor_gene_overlap")
FNL_PC_FST_GWA_PCcor_gene_overlap <- fread("FNL_PC_FST_GWA_PCcor_gene_overlap") %>% select(V5, V6, V7, V4)
colnames(FNL_PC_FST_GWA_PCcor_gene_overlap) <- c("CHROM", "POS", "BPtoo", "gene")

FNL_PC_FST_GWA_PCcor_gene_overlap_stats <- inner_join(FNL_PC_FST_GWA_PCcor_gene_overlap, FNL_PC_FST_GWA_PCcor)
fwrite(FNL_PC_FST_GWA_PCcor_gene_overlap_stats, "FNL_PC_FST_GWA_PCcor_gene_overlap_stats.tsv", col.names = T, row.names = F,
       quote = F, sep = "\t")


ggplot() + geom_point(data =FNL_PC_FST_GWA_PCcor %>%  sample_n(50000), aes(y = FST, x = -log10(GWAS_PCcor_q))) + theme_classic()
ggplot() + geom_point(data =FNL_PC_FST_GWA_PCcor %>%  sample_n(50000), aes(y = FST, x = LRT)) + theme_classic()
ggplot() + geom_point(data =FNL_PC_FST_GWA_PCcor %>%  sample_n(50000), aes(y = LRT, x = beta)) + theme_classic()
ggplot() + geom_point(data =FNL_PC_FST_GWA_PCcor %>% filter(Frequency <= 0.5) %>% sample_n(50000), aes(x = Frequency, y = beta)) + theme_classic()


cor.test(FNL_PC_FST_GWA_PCcor_gene_overlap_stats$FST, FNL_PC_FST_GWA_PCcor_gene_overlap_stats$LRT)

#M_QN PC corrected
setwd("../M_QN_GWAS//")
MQN_PCcorr<- 
  list.files(pattern = "*lrt.*gz") %>% 
  map_df(~fread(.))

colnames(MQN_PCcorr)[1:2] <- c("CHROM", "POS")

MQN_PCcorr<- MQN_PCcorr[!(MQN_PCcorr$LRT %in% "-999"),]
MQN_PCcorr$GWAS_PCcor_pvals <- 1- pchisq(q = MQN_PCcorr$LRT, df = 1)
MQN_PC_FST_GWA_PCcor <- inner_join(PC_FST_M_QN, MQN_PCcorr)

MQN_PCcorrq <- qvalue(MQN_PC_FST_GWA_PCcor$GWAS_PCcor_pvals, pi0 = 1)

MQN_PC_FST_GWA_PCcor$GWAS_PCcor_q <- MQN_PCcorrq$qvalues

fwrite(MQN_PC_FST_GWA_PCcor, "MQN_PC_FST_GWA_PCcorr.tsv", col.names = T, row.names = F, sep = "\t", quote = F)

MQN_PC_FST_GWA_PCcor$POSEND <- MQN_PC_FST_GWA_PCcor$POS + 1 
write.table(cbind(MQN_PC_FST_GWA_PCcor$CHROM, as.integer(MQN_PC_FST_GWA_PCcor$POS), as.integer(MQN_PC_FST_GWA_PCcor$POSEND)), "MQN_PC_FST_GWA_PCcor.bed", sep = "\t",
            col.names = F, row.names = F, quote = F)

system("bedtools intersect -b MQN_PC_FST_GWA_PCcor.bed -a ../SSA_genes.bed -wb > MQN_PC_FST_GWA_PCcor_gene_overlap")
MQN_PC_FST_GWA_PCcor_gene_overlap <- fread("MQN_PC_FST_GWA_PCcor_gene_overlap") %>% select(V5, V6, V7, V4)
colnames(MQN_PC_FST_GWA_PCcor_gene_overlap) <- c("CHROM", "POS", "BPtoo", "gene")

MQN_PC_FST_GWA_PCcor_gene_overlap_stats <- inner_join(MQN_PC_FST_GWA_PCcor_gene_overlap, MQN_PC_FST_GWA_PCcor)
fwrite(MQN_PC_FST_GWA_PCcor_gene_overlap_stats, "MQN_PC_FST_GWA_PCcor_gene_overlap_stats.tsv", col.names = T, row.names = F,
       quote = F, sep = "\t")

cor.test(MQN_PC_FST_GWA_PCcor_gene_overlap_stats$FST, MQN_PC_FST_GWA_PCcor_gene_overlap_stats$LRT)
ggplot() + geom_point(data =MQN_PC_FST_GWA_PCcor %>%  sample_n(50000), aes(y = FST, x = -log10(GWAS_PCcor_q))) + theme_classic()
ggplot() + geom_point(data =MQN_PC_FST_GWA_PCcor %>%  sample_n(50000), aes(y = FST, x = LRT)) + theme_classic()
ggplot() + geom_point(data =MQN_PC_FST_GWA_PCcor %>%  sample_n(50000), aes(y = LRT, x = beta)) + theme_classic()
ggplot() + geom_point(data =MQN_PC_FST_GWA_PCcor %>% filter(Frequency <= 0.5) %>% sample_n(50000), aes(x = Frequency, y = beta)) + theme_classic()


#M_NL PC corrected
setwd("../M_NL_GWAS//")
MNL_PCcorr<- 
  list.files(pattern = "*lrt.*gz") %>% 
  map_df(~fread(.))

colnames(MNL_PCcorr)[1:2] <- c("CHROM", "POS")

MNL_PCcorr<- MNL_PCcorr[!(MNL_PCcorr$LRT %in% "-999"),]
MNL_PCcorr$GWAS_PCcor_pvals <- 1- pchisq(q = MNL_PCcorr$LRT, df = 1)
MNL_PC_FST_GWA_PCcor <- inner_join(PC_FST_M_NL, MNL_PCcorr)

MNL_PCcorrq <- qvalue(MNL_PC_FST_GWA_PCcor$GWAS_PCcor_pvals, pi0 = 1)

MNL_PC_FST_GWA_PCcor$GWAS_PCcor_q <- MNL_PCcorrq$qvalues

fwrite(MNL_PC_FST_GWA_PCcor, "MNL_PC_FST_GWA_PCcorr.tsv", col.names = T, row.names = F, sep = "\t", quote = F)

MNL_PC_FST_GWA_PCcor$POSEND <- MNL_PC_FST_GWA_PCcor$POS + 1 
write.table(cbind(MNL_PC_FST_GWA_PCcor$CHROM, as.integer(MNL_PC_FST_GWA_PCcor$POS), as.integer(MNL_PC_FST_GWA_PCcor$POSEND)), "MNL_PC_FST_GWA_PCcor.bed", sep = "\t",
            col.names = F, row.names = F, quote = F)

system("bedtools intersect -b MNL_PC_FST_GWA_PCcor.bed -a ../SSA_genes.bed -wb > MNL_PC_FST_GWA_PCcor_gene_overlap")
MNL_PC_FST_GWA_PCcor_gene_overlap <- fread("MNL_PC_FST_GWA_PCcor_gene_overlap") %>% select(V5, V6, V7, V4)
colnames(MNL_PC_FST_GWA_PCcor_gene_overlap) <- c("CHROM", "POS", "BPtoo", "gene")

MNL_PC_FST_GWA_PCcor_gene_overlap_stats <- inner_join(MNL_PC_FST_GWA_PCcor_gene_overlap, MNL_PC_FST_GWA_PCcor)
fwrite(MNL_PC_FST_GWA_PCcor_gene_overlap_stats, "MNL_PC_FST_GWA_PCcor_gene_overlap_stats.tsv", col.names = T, row.names = F,
       quote = F, sep = "\t")


ggplot() + geom_point(data =MNL_PC_FST_GWA_PCcor %>%  sample_n(50000), aes(y = FST, x = -log10(GWAS_PCcor_q))) + theme_classic()
ggplot() + geom_point(data =MNL_PC_FST_GWA_PCcor %>%  sample_n(50000), aes(y = FST, x = LRT)) + theme_classic()
ggplot() + geom_point(data =MNL_PC_FST_GWA_PCcor %>%  sample_n(50000), aes(y = LRT, x = beta)) + theme_classic()
ggplot() + geom_point(data =MNL_PC_FST_GWA_PCcor %>% filter(Frequency <= 0.5) %>% sample_n(50000), aes(x = Frequency, y = beta)) + theme_classic()

cor.test(MNL_PC_FST_GWA_PCcor_gene_overlap_stats$LRT, MNL_PC_FST_GWA_PCcor_gene_overlap_stats$FST )

#Plotting and overlap with six6 or vggll3 
PC_FST_GWA_PCcor_MSW1SW_gene_overlap_stats <- PC_FST_GWA_PCcor_MSW1SW_gene_overlap_stats %>%  mutate(SNP =  paste0(CHROM, "_", POS))
MNL_PC_FST_GWA_PCcor_gene_overlap_stats <- MNL_PC_FST_GWA_PCcor_gene_overlap_stats %>%  mutate(SNP =  paste0(CHROM, "_", POS))
FNL_PC_FST_GWA_PCcor_gene_overlap_stats <- FNL_PC_FST_GWA_PCcor_gene_overlap_stats %>%  mutate(SNP =  paste0(CHROM, "_", POS))
MQN_PC_FST_GWA_PCcor_gene_overlap_stats <- MQN_PC_FST_GWA_PCcor_gene_overlap_stats %>%  mutate(SNP =  paste0(CHROM, "_", POS))
FQN_PC_FST_GWA_PCcor_gene_overlap_stats <- FQN_PC_FST_GWA_PCcor_gene_overlap_stats %>%  mutate(SNP =  paste0(CHROM, "_", POS))

All_PCcor_top2000 <- PC_FST_GWA_PCcor_MSW1SW_gene_overlap_stats %>%
  top_n(2000, LRT)

FQN_PC_FST_GWA_PCcor_top500 <- FQN_PC_FST_GWA_PCcor_gene_overlap_stats %>% 
  top_n(500, LRT)

FNL_PC_FST_GWA_PCcor_top500 <-FNL_PC_FST_GWA_PCcor_gene_overlap_stats  %>% 
  top_n(500, LRT)

MQN_PC_FST_GWA_PCcor_top500 <- MQN_PC_FST_GWA_PCcor_gene_overlap_stats  %>% 
  top_n(500, LRT)

MNL_PC_FST_GWA_PCcor_top500 <- FQN_PC_FST_GWA_PCcor_gene_overlap_stats  %>% 
  top_n(500, LRT)


ALL_PCload <- ggman(gwas = All_sites_sigsites_PC1, snp = "SNP", bp = "POS", chrom = "CHROM", pvalue = "qvals_PC1", 
                    xlabel = "Chromosome", ylabel = "-log10(qvalue)", pointsize = 5, 
                    logTransform = T) + theme_classic()

ggmanHighlightGroup(ALL_PCload, highlightDfm = six6_vgll3_WG, snp = "SNP", group = "gene")

All_PC_GWAS_man <- ggman(gwas = All_PCcor_top500, snp = "SNP", chrom = "CHROM", bp = "POS", pvalue = "GWAS_PCcor_q",
                         xlabel = "Chromosome", ylabel = "-log10(qvalue)", pointSize = 1, 
                         logTransform = T, sigLine = -log10(0.05)) + theme_classic()
ggmanHighlightGroup(All_PC_GWAS_man, highlightDfm = six6_vgll3_WG, snp = "SNP", group = "gene")


All_PCcor_top500 %>%  filter(GWAS_PCcor_q < 0.05)
FQN_PC_FST_GWA_PCcor_top500 %>%  filter(GWAS_PCcor_q < 0.05)


FQN_PC_FST_GWA_PCcor_top500 %>%  top_n(FST, n = 10)
FQN_PC_GWAS_man <- ggman(gwas = FQN_PC_FST_GWA_PCcor_top500, snp = "SNP", chrom = "CHROM", bp = "POS", pvalue = "GWAS_PCcor_q",
                         xlabel = "Chromosome", ylabel = "-log10(qvalue)", pointSize = 2, 
                         logTransform = T, sigLine = -log10(0.05)) + theme_classic()
ggmanHighlightGroup(FQN_PC_GWAS_man, highlightDfm = six6_vgll3_WG, snp = "SNP", group = "gene")


MNL_LRT100 <- MNL_PC_FST_GWA_PCcor %>%  slice_max(LRT, n = 100)
FNL_LRT100 <- FNL_PC_FST_GWA_PCcor  %>%  slice_max(LRT, n = 100)
MQN_LRT100 <- MQN_PC_FST_GWA_PCcor  %>%  slice_max(LRT, n = 100)
FQN_LRT100 <- FQN_PC_FST_GWA_PCcor  %>%  slice_max(LRT, n = 100)

mean_check <- function(allDF,OLdf){
  mean_100 = mean(OLdf$PC1)
  mean_all = mean(allDF$PC1)
  return(data.frame(mean_100, mean_all))}
mean_check(OLdf = FQN_LRT100, allDF = FQN_PC_FST_GWA_PCcor)
wilcox.test(abs(FQN_LRT100$beta), abs(FQN_PC_FST_GWA_PCcor$beta))

#for topGO
MNL_LRT100_bed <- MNL_LRT100 %>% select(CHROM, POS) %>%  mutate(POS2 = POS + 1)
FNL_LRT100_bed <- FNL_LRT100 %>% select(CHROM, POS) %>%  mutate(POS2 = POS + 1)
MQN_LRT100_bed <- MQN_LRT100 %>% select(CHROM, POS) %>%  mutate(POS2 = POS + 1)
FQN_LRT100_bed <- FQN_LRT100 %>% select(CHROM, POS) %>%  mutate(POS2 = POS + 1)
fwrite(MNL_LRT100_bed, "MNL_LRT100.bed", col.names = F, sep = "\t", row.names = F, quote = F)
fwrite(FNL_LRT100_bed, "FNL_LRT100.bed", col.names = F, sep = "\t", row.names = F, quote = F)
fwrite(MQN_LRT100_bed, "MQN_LRT100.bed", col.names = F, sep = "\t", row.names = F, quote = F)
fwrite(FQN_LRT100_bed, "FQN_LRT100.bed", col.names = F, sep = "\t", row.names = F, quote = F)

ssalar_genes <- fread("ssalar_gene_result.txt")
ssalar_bed <- ssalar_genes

