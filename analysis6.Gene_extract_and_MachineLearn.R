###GET READY######

library(RcppCNPy)
library(tidyverse)
library(data.table)
library(qvalue)
library(randomForest)
library(pROC)
library(VSURF)

setwd(dir = "~/Desktop/Working/Sea_Age/")
#Read in metadata #####
Salmon_Metadata_Bam_Unique_HD  <- fread("Salmon_Metadata_Bam_Unique_HD.txt")
#Recode
Salmon_Metadata_Bam_Unique_HD$AgeClass[Salmon_Metadata_Bam_Unique_HD$AgeClass=="1SW"] <- 0
Salmon_Metadata_Bam_Unique_HD$AgeClass[Salmon_Metadata_Bam_Unique_HD$AgeClass=="MSW"] <- 1

#Set pop groups
QCNB <- c("MSW", "MUN", "RTG", "TRI")
NL <- c("ENG", "MBB", "SAN", "RVF")

#Separate by region and sex
Male_1SW_MSW_QN <- Salmon_Metadata_Bam_Unique_HD %>% filter(PCR_Sex =="M", RiverCode %in% QCNB)
Male_1SW_MSW_NL <-  Salmon_Metadata_Bam_Unique_HD %>% filter(PCR_Sex =="M", RiverCode %in% NL)
Female_1SW_MSW_QN <- Salmon_Metadata_Bam_Unique_HD %>% filter(PCR_Sex =="F", RiverCode %in% QCNB)
Female_1SW_MSW_NL <- Salmon_Metadata_Bam_Unique_HD %>% filter(PCR_Sex =="F", RiverCode %in% NL)

#get GWAS PCA FST stats
PC_FST_GWA_PCcor_MSW_1SW <- fread("All_2PCcorr/PC_FST_GWA_PCcor_MSW_1SW.tsv")
PC_FST_GWA_PCcor_MSW_1SW_gene_overlap <- fread("All_2PCcorr/PC_FST_GWA_PCcor_MSW1SW_gene_overlap_stats.tsv")

PC_FST_GWA_PCcor_MSW_1SW %>% filter(GWAS_PCcor_q < 0.05) %>%  arrange(desc(LRT)) 

PC_FST_GWA_PCcor_MSW_1SW$beagle_coord <- paste0(PC_FST_GWA_PCcor_MSW_1SW$CHROM, "_", PC_FST_GWA_PCcor_MSW_1SW$POS)
PC_FST_GWA_PCcor_MSW_1SW_gene_overlap$beagle_coord <- paste0(PC_FST_GWA_PCcor_MSW_1SW_gene_overlap$CHROM, "_", PC_FST_GWA_PCcor_MSW_1SW_gene_overlap$POS)

FQN_PC_FST_GWA_PCcor<- fread("F_QN_GWAS/FQN_PC_FST_GWA_PCcorr.tsv")
FQN_PC_FST_GWA_PCcorr_gene_overlap <- fread("F_QN_GWAS/FQN_PC_FST_GWA_PCcor_gene_overlap_stats.tsv")

FQN_PC_FST_GWA_PCcor %>% filter(GWAS_PCcor_q < 0.05) %>%  arrange(desc(-log10(GWAS_PCcor_q)))
FQN_PC_FST_GWA_PCcorr_gene_overlap %>% filter(GWAS_PCcor_q < 0.05) %>%  arrange(desc(-log10(GWAS_PCcor_q)))

FQN_PC_FST_GWA_PCcorr_gene_overlap$beagle_coord <- paste0(FQN_PC_FST_GWA_PCcorr_gene_overlap$CHROM, "_", FQN_PC_FST_GWA_PCcorr_gene_overlap$POS)
FQN_PC_FST_GWA_PCcor$beagle_coord <- paste0(FQN_PC_FST_GWA_PCcor$CHROM, "_", FQN_PC_FST_GWA_PCcor$POS)

FNL_PC_FST_GWA_PCcor<- fread("F_NL_GWAS/FNL_PC_FST_GWA_PCcorr.tsv")
FNL_PC_FST_GWA_PCcorr_gene_overlap <- fread("F_NL_GWAS/FNL_PC_FST_GWA_PCcor_gene_overlap_stats.tsv")

FNL_PC_FST_GWA_PCcor %>% filter(GWAS_PCcor_q < 0.05) %>%  arrange(desc(-log10(GWAS_PCcor_q)))
FNL_PC_FST_GWA_PCcorr_gene_overlap %>% filter(GWAS_PCcor_q < 0.05) %>%  arrange(desc(-log10(GWAS_PCcor_q)))

FNL_PC_FST_GWA_PCcorr_gene_overlap$beagle_coord <- paste0(FNL_PC_FST_GWA_PCcorr_gene_overlap$CHROM, "_", FNL_PC_FST_GWA_PCcorr_gene_overlap$POS)
FNL_PC_FST_GWA_PCcor$beagle_coord <- paste0(FNL_PC_FST_GWA_PCcor$CHROM, "_", FNL_PC_FST_GWA_PCcor$POS)

MQN_PC_FST_GWA_PCcor<- fread("M_QN_GWAS/MQN_PC_FST_GWA_PCcorr.tsv")
MQN_PC_FST_GWA_PCcorr_gene_overlap <- fread("M_QN_GWAS/MQN_PC_FST_GWA_PCcor_gene_overlap_stats.tsv")

MQN_PC_FST_GWA_PCcorr_gene_overlap$beagle_coord <- paste0(MQN_PC_FST_GWA_PCcorr_gene_overlap$CHROM, "_", MQN_PC_FST_GWA_PCcorr_gene_overlap$POS)
MQN_PC_FST_GWA_PCcor$beagle_coord <- paste0(MQN_PC_FST_GWA_PCcor$CHROM, "_", MQN_PC_FST_GWA_PCcor$POS)

MQN_PC_FST_GWA_PCcor %>% filter(GWAS_PCcor_q < 0.05) %>%  arrange(desc(-log10(GWAS_PCcor_q)))
MQN_PC_FST_GWA_PCcorr_gene_overlap %>% filter(GWAS_PCcor_q < 0.05) %>%  arrange(desc(-log10(GWAS_PCcor_q)))

MNL_PC_FST_GWA_PCcor<- fread("M_NL_GWAS/MNL_PC_FST_GWA_PCcorr.tsv")
MNL_PC_FST_GWA_PCcorr_gene_overlap <- fread("M_NL_GWAS/MNL_PC_FST_GWA_PCcor_gene_overlap_stats.tsv")

MNL_PC_FST_GWA_PCcorr_gene_overlap$beagle_coord <- paste0(MNL_PC_FST_GWA_PCcorr_gene_overlap$CHROM, "_", MNL_PC_FST_GWA_PCcorr_gene_overlap$POS)
MNL_PC_FST_GWA_PCcor$beagle_coord <- paste0(MNL_PC_FST_GWA_PCcor$CHROM, "_", MNL_PC_FST_GWA_PCcor$POS)

MNL_PC_FST_GWA_PCcor %>% filter(GWAS_PCcor_q < 0.05) %>%  arrange(desc(-log10(GWAS_PCcor_q)))
MNL_PC_FST_GWA_PCcorr_gene_overlap %>% filter(GWAS_PCcor_q < 0.05) %>%  arrange(desc(-log10(GWAS_PCcor_q)))

#get top loci for each group, by GWAS score - get top 500 first to filter beagle via command line, can downsample in R from there


All_PCcor_top500 <- PC_FST_GWA_PCcor_MSW_1SW %>%
  top_n(500, LRT)

FQN_PC_FST_GWA_PCcor_top500 <- FQN_PC_FST_GWA_PCcor %>% 
  top_n(500, LRT)

FNL_PC_FST_GWA_PCcor_top500 <- FNL_PC_FST_GWA_PCcor %>% 
  top_n(500, LRT)

MQN_PC_FST_GWA_PCcor_top500 <- MQN_PC_FST_GWA_PCcor %>% 
  top_n(500, LRT)

MNL_PC_FST_GWA_PCcor_top500 <- MNL_PC_FST_GWA_PCcor %>% 
  top_n(500, LRT)

#output for beagle file filtering

All_top500 <- unique(c(All_PCcor_top500$beagle_coord, 
                       FQN_PC_FST_GWA_PCcor_top500$beagle_coord, 
                       FNL_PC_FST_GWA_PCcor_top500$beagle_coord,
                       MQN_PC_FST_GWA_PCcor_top500$beagle_coord,
                       MNL_PC_FST_GWA_PCcor_top500$beagle_coord))


Rando_1000 <- PC_FST_GWA_PCcor_MSW_1SW %>% 
  filter(!beagle_coord %in% All_top500) %>% 
  sample_n(1000) %>% 
  select(beagle_coord)

GWAS_seaage_All_top500 <- PC_FST_GWA_PCcor_MSW_1SW %>%  filter(beagle_coord %in% All_top500)

Beagle_dose_header <- c("beagle_coord", "allele1", "allele2", Salmon_Metadata_Bam_Unique_HD$FGL_ID)

#beagle filtering - filter all individuals and loci first, and then DS here 

write.table(All_top500, "All_top500", col.names = F, row.names = F, sep = "\t", quote = F)
write.table(Rando_1000$beagle_coord, "Rando_1000",  col.names = F, row.names = F, sep = "\t", quote = F)

All_top1000_dose_allgroups <-fread("All_top1000.beagle.dose")
colnames(All_top1000_dose_allgroups ) <- Beagle_dose_header

#all
All_PCcor_top_500_beagcoord <- All_PCcor_top500$beagle_coord
AllGWAS_top500_dose <-  All_top1000_dose_allgroups %>%  filter(beagle_coord %in% All_PCcor_top500$beagle_coord) 
AllGWAS_top500_dose  <- AllGWAS_top500_dose  %>%   select(-c(1,2,3))

AllGWAS_top500_dose_t <- data.frame(t(AllGWAS_top500_dose)) %>% mutate(FGL_ID = colnames(AllGWAS_top500_dose))  %>% 
  relocate(FGL_ID, .before = X1)

colnames(AllGWAS_top500_dose_t)[2:501] <-All_PCcor_top_500_beagcoord
AllGWAS_top500_dose_phenos <- inner_join(Salmon_Metadata_Bam_Unique_HD, AllGWAS_top500_dose_t)


#split data
table(AllGWAS_top500_dose_phenos$AgeClass)/2
ALL_sample_size <- c(140,140)

All_PCcor_top_500_beagcoord <- All_PCcor_top500$beagle_coord
All_PCcor_top_100_beagcoord<-  All_PCcor_top500  %>%  top_n(100, LRT) %>% select(beagle_coord)
All_PCcor_top_200_beagcoord<-  All_PCcor_top500  %>%  top_n(200, LRT) %>% select(beagle_coord)
All_PCcor_top_300_beagcoord<-  All_PCcor_top500  %>%  top_n(300, LRT) %>% select(beagle_coord)
All_PCcor_top_400_beagcoord<-  All_PCcor_top500  %>%  top_n(400, LRT) %>% select(beagle_coord)

All_PCcor_top_100_beagcoord <- All_PCcor_top_100_beagcoord$beagle_coord
All_PCcor_top_200_beagcoord <- All_PCcor_top_200_beagcoord$beagle_coord
All_PCcor_top_300_beagcoord <- All_PCcor_top_300_beagcoord$beagle_coord
All_PCcor_top_400_beagcoord <- All_PCcor_top_400_beagcoord$beagle_coord


ALL_rf_500 <- randomForest(x = AllGWAS_top500_dose_phenos[,..All_PCcor_top_500_beagcoord], y = as.factor(AllGWAS_top500_dose_phenos$AgeClass), importance=TRUE ,proximity=TRUE, mtry=400, ntree=25000, strata=as.factor(AllGWAS_top500_dose_phenos$AgeClass), sampsize=ALL_sample_size)
ALL_rf_400 <- randomForest(x = AllGWAS_top500_dose_phenos[,..All_PCcor_top_400_beagcoord], y = as.factor(AllGWAS_top500_dose_phenos$AgeClass), importance=TRUE ,proximity=TRUE, mtry=320, ntree=25000, strata=as.factor(AllGWAS_top500_dose_phenos$AgeClass), sampsize=ALL_sample_size)
ALL_rf_300 <- randomForest(x = AllGWAS_top500_dose_phenos[,..All_PCcor_top_300_beagcoord ], y = as.factor(AllGWAS_top500_dose_phenos$AgeClass), importance=TRUE ,proximity=TRUE, mtry=240, ntree=25000, strata=as.factor(AllGWAS_top500_dose_phenos$AgeClass), sampsize=ALL_sample_size)
ALL_rf_200 <- randomForest(x = AllGWAS_top500_dose_phenos[,..All_PCcor_top_200_beagcoord ], y = as.factor(AllGWAS_top500_dose_phenos$AgeClass), importance=TRUE ,proximity=TRUE, mtry=160, ntree=25000, strata=as.factor(AllGWAS_top500_dose_phenos$AgeClass), sampsize=ALL_sample_size)
ALL_rf_100 <- randomForest(x = AllGWAS_top500_dose_phenos[,..All_PCcor_top_100_beagcoord ], y = as.factor(AllGWAS_top500_dose_phenos$AgeClass), importance=TRUE ,proximity=TRUE, mtry=80, ntree=25000, strata=as.factor(AllGWAS_top500_dose_phenos$AgeClass), sampsize=ALL_sample_size)

All_err_rate <- data.frame(cbind(rbind(ALL_rf_500$err.rate[25000], ALL_rf_400$err.rate[25000], ALL_rf_300$err.rate[25000], ALL_rf_200$err.rate[25000], ALL_rf_100$err.rate[25000]), c(500, 400, 300, 200, 100), 'AllGWAS_top_SNPs'))
colnames(All_err_rate) <- c("OOB", "SNPs", "panel")

#Get ROC
ALL_rf_500_ROC <- roc(AllGWAS_top500_dose_phenos$AgeClass, ALL_rf_500$votes[,2])
plot(ALL_rf_500_ROC)
auc(ALL_rf_500_ROC)

#break down into groups by region and sex 
#FQN
#filter by GWAS outlier
F_QN_dose <-  All_top500_dose %>% filter(beagle_coord %in% FQN_PC_FST_GWA_PCcor_top500$beagle_coord) 

F_QN_top500_beagcoord <- F_QN_dose$beagle_coord

F_QN_dose <- F_QN_dose %>%   select(-c(1,2,3))

#transpose
F_QN_dose_t  <- data.frame(t(F_QN_dose)) %>% mutate(FGL_ID = colnames(F_QN_dose))  %>% 
  relocate(FGL_ID, .before = X1)

colnames(F_QN_dose_t)[2:501] <- F_QN_top500_beagcoord
F_QN_dose_phenos <- inner_join(Female_1SW_MSW_QN, F_QN_dose_t)


table(F_QN_dose_phenos$AgeClass)/2
FQN_sample_size <- c(29,29)

F_QN_top100_beagcoord<-  FQN_PC_FST_GWA_PCcor_top500  %>%  top_n(100, LRT) %>% select(beagle_coord)
F_QN_top200_beagcoord<-  FQN_PC_FST_GWA_PCcor_top500  %>%  top_n(200, LRT) %>% select(beagle_coord)
F_QN_top300_beagcoord<-  FQN_PC_FST_GWA_PCcor_top500  %>%  top_n(300, LRT) %>% select(beagle_coord)
F_QN_top400_beagcoord<-  FQN_PC_FST_GWA_PCcor_top500  %>%  top_n(400, LRT) %>% select(beagle_coord)

F_QN_top400_beagcoord <- F_QN_top400_beagcoord$beagle_coord
F_QN_top300_beagcoord <- F_QN_top300_beagcoord$beagle_coord
F_QN_top200_beagcoord <- F_QN_top200_beagcoord$beagle_coord
F_QN_top100_beagcoord <- F_QN_top100_beagcoord$beagle_coord


FQN_rf_500 <- randomForest(x = F_QN_dose_phenos[,..F_QN_top500_beagcoord ], y = as.factor(F_QN_dose_phenos$AgeClass), importance=TRUE ,proximity=TRUE, mtry=400, ntree=25000, strata=as.factor(F_QN_dose_phenos$AgeClass), sampsize=FQN_sample_size)
FQN_rf_400 <- randomForest(x = F_QN_dose_phenos[,..F_QN_top400_beagcoord], y = as.factor(F_QN_dose_phenos$AgeClass), importance=TRUE ,proximity=TRUE, mtry=320, ntree=25000, strata=as.factor(F_QN_dose_phenos$AgeClass), sampsize=FQN_sample_size)
FQN_rf_300 <- randomForest(x = F_QN_dose_phenos[,..F_QN_top300_beagcoord ], y = as.factor(F_QN_dose_phenos$AgeClass), importance=TRUE ,proximity=TRUE, mtry=240, ntree=25000, strata=as.factor(F_QN_dose_phenos$AgeClass), sampsize=FQN_sample_size)
FQN_rf_200 <- randomForest(x = F_QN_dose_phenos[,..F_QN_top200_beagcoord ], y = as.factor(F_QN_dose_phenos$AgeClass), importance=TRUE ,proximity=TRUE, mtry=160, ntree=25000, strata=as.factor(F_QN_dose_phenos$AgeClass), sampsize=FQN_sample_size)
FQN_rf_100 <- randomForest(x = F_QN_dose_phenos[,..F_QN_top100_beagcoord ], y = as.factor(F_QN_dose_phenos$AgeClass), importance=TRUE ,proximity=TRUE, mtry=80, ntree=25000, strata=as.factor(F_QN_dose_phenos$AgeClass), sampsize=FQN_sample_size)

F_QN_err_rate <- data.frame(cbind(rbind(FQN_rf_500$err.rate[25000], FQN_rf_400$err.rate[25000], FQN_rf_300$err.rate[25000], FQN_rf_200$err.rate[25000], FQN_rf_100$err.rate[25000]), c(500, 400, 300, 200, 100), 'FQN_top_SNPs'))
colnames(F_QN_err_rate ) <- c("OOB", "SNPs", "panel")

#F NL

F_NL_dose <-  All_top500_dose %>% filter(beagle_coord %in% FNL_PC_FST_GWA_PCcor_top500$beagle_coord) 

F_NL_top500_beagcoord <- F_NL_dose$beagle_coord

F_NL_dose <- F_NL_dose %>%   select(-c(1,2,3))

#transpose
F_NL_dose_t  <- data.frame(t(F_NL_dose)) %>% mutate(FGL_ID = colnames(F_NL_dose))  %>% 
  relocate(FGL_ID, .before = X1)

colnames(F_NL_dose_t)[2:501] <- F_NL_top500_beagcoord
F_NL_dose_phenos <- inner_join(Female_1SW_MSW_NL, F_NL_dose_t)


table(F_NL_dose_phenos$AgeClass)/2
FNL_sample_size <- c(41,41)

F_NL_top100_beagcoord <-  FNL_PC_FST_GWA_PCcor_top500  %>%  top_n(100, LRT) %>% select(beagle_coord)
F_NL_top200_beagcoord <-  FNL_PC_FST_GWA_PCcor_top500  %>%  top_n(200, LRT) %>% select(beagle_coord)
F_NL_top300_beagcoord <-  FNL_PC_FST_GWA_PCcor_top500  %>%  top_n(300, LRT) %>% select(beagle_coord)
F_NL_top400_beagcoord <-  FNL_PC_FST_GWA_PCcor_top500  %>%  top_n(400, LRT) %>% select(beagle_coord)

F_NL_top400_beagcoord <- F_NL_top400_beagcoord$beagle_coord
F_NL_top300_beagcoord <- F_NL_top300_beagcoord$beagle_coord
F_NL_top200_beagcoord <- F_NL_top200_beagcoord$beagle_coord
F_NL_top100_beagcoord <- F_NL_top100_beagcoord$beagle_coord


FNL_rf_500 <- randomForest(x = F_NL_dose_phenos[,..F_NL_top500_beagcoord ], y = as.factor(F_NL_dose_phenos$AgeClass), importance=TRUE ,proximity=TRUE, mtry=400, ntree=25000, strata=as.factor(F_NL_dose_phenos$AgeClass), sampsize=FNL_sample_size)
FNL_rf_400 <- randomForest(x = F_NL_dose_phenos[,..F_NL_top400_beagcoord], y = as.factor(F_NL_dose_phenos$AgeClass), importance=TRUE ,proximity=TRUE, mtry=320, ntree=25000, strata=as.factor(F_NL_dose_phenos$AgeClass), sampsize=FNL_sample_size)
FNL_rf_300 <- randomForest(x = F_NL_dose_phenos[,..F_NL_top300_beagcoord ], y = as.factor(F_NL_dose_phenos$AgeClass), importance=TRUE ,proximity=TRUE, mtry=240, ntree=25000, strata=as.factor(F_NL_dose_phenos$AgeClass), sampsize=FNL_sample_size)
FNL_rf_200 <- randomForest(x = F_NL_dose_phenos[,..F_NL_top200_beagcoord ], y = as.factor(F_NL_dose_phenos$AgeClass), importance=TRUE ,proximity=TRUE, mtry=160, ntree=25000, strata=as.factor(F_NL_dose_phenos$AgeClass), sampsize=FNL_sample_size)
FNL_rf_100 <- randomForest(x = F_NL_dose_phenos[,..F_NL_top100_beagcoord ], y = as.factor(F_NL_dose_phenos$AgeClass), importance=TRUE ,proximity=TRUE, mtry=80, ntree=25000, strata=as.factor(F_NL_dose_phenos$AgeClass), sampsize=FNL_sample_size)


F_NL_err_rate <- data.frame(cbind(rbind(FNL_rf_500$err.rate[25000], FNL_rf_400$err.rate[25000], FNL_rf_300$err.rate[25000], FNL_rf_200$err.rate[25000], FNL_rf_100$err.rate[25000]), c(500, 400, 300, 200, 100), 'FNL_top_SNPs'))
colnames(F_NL_err_rate ) <- c("OOB", "SNPs", "panel")
hmm[250000]
#break down into groups by region and sex 
#filter by GWAS outlier
#MQN

M_QN_dose <-  All_top500_dose %>% filter(beagle_coord %in% MQN_PC_FST_GWA_PCcor_top500$beagle_coord) 

M_QN_top500_beagcoord <- M_QN_dose$beagle_coord

M_QN_dose <- M_QN_dose %>%   select(-c(1,2,3))

#transpose
M_QN_dose_t  <- data.frame(t(M_QN_dose)) %>% mutate(FGL_ID = colnames(M_QN_dose))  %>% 
  relocate(FGL_ID, .before = X1)

colnames(M_QN_dose_t)[2:501] <- M_QN_top500_beagcoord
M_QN_dose_phenos <- inner_join(Male_1SW_MSW_QN, M_QN_dose_t)


table(M_QN_dose_phenos$AgeClass)/2
MQN_sample_size <- c(34,34)

M_QN_top100_beagcoord<-  MQN_PC_FST_GWA_PCcor_top500  %>%  top_n(100, LRT) %>% select(beagle_coord)
M_QN_top200_beagcoord<-  MQN_PC_FST_GWA_PCcor_top500  %>%  top_n(200, LRT) %>% select(beagle_coord)
M_QN_top300_beagcoord<-  MQN_PC_FST_GWA_PCcor_top500  %>%  top_n(300, LRT) %>% select(beagle_coord)
M_QN_top400_beagcoord<-  MQN_PC_FST_GWA_PCcor_top500  %>%  top_n(400, LRT) %>% select(beagle_coord)

M_QN_top400_beagcoord <- M_QN_top400_beagcoord$beagle_coord
M_QN_top300_beagcoord <- M_QN_top300_beagcoord$beagle_coord
M_QN_top200_beagcoord <- M_QN_top200_beagcoord$beagle_coord
M_QN_top100_beagcoord <- M_QN_top100_beagcoord$beagle_coord


MQN_rf_500 <- randomForest(x = M_QN_dose_phenos[,..M_QN_top500_beagcoord ], y = as.factor(M_QN_dose_phenos$AgeClass), importance=TRUE ,proximity=TRUE, mtry=400, ntree=25000, strata=as.factor(M_QN_dose_phenos$AgeClass), sampsize=MQN_sample_size)
MQN_rf_400 <- randomForest(x = M_QN_dose_phenos[,..M_QN_top400_beagcoord], y = as.factor(M_QN_dose_phenos$AgeClass), importance=TRUE ,proximity=TRUE, mtry=320, ntree=25000, strata=as.factor(M_QN_dose_phenos$AgeClass), sampsize=MQN_sample_size)
MQN_rf_300 <- randomForest(x = M_QN_dose_phenos[,..M_QN_top300_beagcoord ], y = as.factor(M_QN_dose_phenos$AgeClass), importance=TRUE ,proximity=TRUE, mtry=240, ntree=25000, strata=as.factor(M_QN_dose_phenos$AgeClass), sampsize=MQN_sample_size)
MQN_rf_200 <- randomForest(x = M_QN_dose_phenos[,..M_QN_top200_beagcoord ], y = as.factor(M_QN_dose_phenos$AgeClass), importance=TRUE ,proximity=TRUE, mtry=160, ntree=25000, strata=as.factor(M_QN_dose_phenos$AgeClass), sampsize=MQN_sample_size)
MQN_rf_100 <- randomForest(x = M_QN_dose_phenos[,..M_QN_top100_beagcoord ], y = as.factor(M_QN_dose_phenos$AgeClass), importance=TRUE ,proximity=TRUE, mtry=80, ntree=25000, strata=as.factor(M_QN_dose_phenos$AgeClass), sampsize=MQN_sample_size)

M_QN_err_rate <- data.frame(cbind(rbind(MQN_rf_500$err.rate[25000], MQN_rf_400$err.rate[25000], MQN_rf_300$err.rate[25000], MQN_rf_200$err.rate[25000], MQN_rf_100$err.rate[25000]), c(500, 400, 300, 200, 100), 'MQN_top_SNPs'))
colnames(M_QN_err_rate ) <- c("OOB", "SNPs", "panel")

###
M_NL_dose <-  All_top500_dose %>% filter(beagle_coord %in% MNL_PC_FST_GWA_PCcor_top500$beagle_coord) 

M_NL_top500_beagcoord <- M_NL_dose$beagle_coord

M_NL_dose <- M_NL_dose %>%   select(-c(1,2,3))

#transpose
M_NL_dose_t  <- data.frame(t(M_NL_dose)) %>% mutate(FGL_ID = colnames(M_NL_dose))  %>% 
  relocate(FGL_ID, .before = X1)

colnames(M_NL_dose_t)[2:501] <- M_NL_top500_beagcoord
M_NL_dose_phenos <- inner_join(Male_1SW_MSW_NL, M_NL_dose_t)


table(M_NL_dose_phenos$AgeClass)/2
MNL_sample_size <- c(18,18)

M_NL_top100_beagcoord<-  MNL_PC_FST_GWA_PCcor_top500  %>%  top_n(100, LRT) %>% select(beagle_coord)
M_NL_top200_beagcoord<-  MNL_PC_FST_GWA_PCcor_top500  %>%  top_n(200, LRT) %>% select(beagle_coord)
M_NL_top300_beagcoord<-  MNL_PC_FST_GWA_PCcor_top500  %>%  top_n(300, LRT) %>% select(beagle_coord)
M_NL_top400_beagcoord<-  MNL_PC_FST_GWA_PCcor_top500  %>%  top_n(400, LRT) %>% select(beagle_coord)

M_NL_top400_beagcoord <- M_NL_top400_beagcoord$beagle_coord
M_NL_top300_beagcoord <- M_NL_top300_beagcoord$beagle_coord
M_NL_top200_beagcoord <- M_NL_top200_beagcoord$beagle_coord
M_NL_top100_beagcoord <- M_NL_top100_beagcoord$beagle_coord

set.seed(1347238923)
MNL_rf_500 <- randomForest(x = M_NL_dose_phenos[,..M_NL_top500_beagcoord ], y = as.factor(M_NL_dose_phenos$AgeClass), importance=TRUE ,proximity=TRUE, mtry=400, ntree=25000, strata=as.factor(M_NL_dose_phenos$AgeClass), sampsize=MNL_sample_size)
MNL_rf_400 <- randomForest(x = M_NL_dose_phenos[,..M_NL_top400_beagcoord], y = as.factor(M_NL_dose_phenos$AgeClass), importance=TRUE ,proximity=TRUE, mtry=320, ntree=25000, strata=as.factor(M_QN_dose_phenos$AgeClass), sampsize=MNL_sample_size)
MNL_rf_300 <- randomForest(x = M_NL_dose_phenos[,..M_NL_top300_beagcoord ], y = as.factor(M_NL_dose_phenos$AgeClass), importance=TRUE ,proximity=TRUE, mtry=240, ntree=25000, strata=as.factor(M_QN_dose_phenos$AgeClass), sampsize=MNL_sample_size)
MNL_rf_200 <- randomForest(x = M_NL_dose_phenos[,..M_NL_top200_beagcoord ], y = as.factor(M_NL_dose_phenos$AgeClass), importance=TRUE ,proximity=TRUE, mtry=160, ntree=25000, strata=as.factor(M_QN_dose_phenos$AgeClass), sampsize=MNL_sample_size)
MNL_rf_100 <- randomForest(x = M_NL_dose_phenos[,..M_NL_top100_beagcoord ], y = as.factor(M_NL_dose_phenos$AgeClass), importance=TRUE ,proximity=TRUE, mtry=80, ntree=25000, strata=as.factor(M_QN_dose_phenos$AgeClass), sampsize=MNL_sample_size)

M_NL_err_rate <- data.frame(cbind(rbind(MNL_rf_500$err.rate[25000], MNL_rf_400$err.rate[25000], MNL_rf_300$err.rate[25000], MNL_rf_200$err.rate[25000], MNL_rf_100$err.rate[25000]), c(500, 400, 300, 200, 100), 'MNL_top_SNPs'))
colnames(M_NL_err_rate ) <- c("OOB", "SNPs", "panel")



###RANDOM######

#now do random loci

All_1000rando_dose <-fread("All_random1000.beagle.dose")
colnames(All_1000rando_dose) <- Beagle_dose_header
All_1000rando_dose_beagle_coord <- All_1000rando_dose$beagle_coord
  
All_1000rando_dose <- All_1000rando_dose %>%   select(-c(1,2,3))


#transpose
All_1000rando_dose_dose_t  <- data.frame(t(All_1000rando_dose)) %>% mutate(FGL_ID = colnames(All_1000rando_dose))  %>% 
  relocate(FGL_ID, .before = X1)
colnames(All_1000rando_dose_dose_t)[2:1026] <- All_1000rando_dose_beagle_coord

All_dose_genes_phenos <- inner_join(Salmon_Metadata_Bam_Unique_HD, All_1000rando_dose_dose_t)



All_random_500_beagle_coord <- sample(x = All_1000rando_dose_beagle_coord, 500)
All_random_400_beagle_coord <- sample(x = All_1000rando_dose_beagle_coord, 400)
All_random_300_beagle_coord <- sample(x = All_1000rando_dose_beagle_coord, 300)
All_random_200_beagle_coord <- sample(x = All_1000rando_dose_beagle_coord, 200)
All_random_100_beagle_coord <- sample(x = All_1000rando_dose_beagle_coord, 100)


ALL_random_rf_1000 <- randomForest(x = All_dose_genes_phenos[,..All_1000rando_dose_beagle_coord], y = as.factor(All_dose_genes_phenos$AgeClass), importance=TRUE ,proximity=TRUE, mtry=800, ntree=25000, strata=as.factor(All_dose_genes_phenos$AgeClass), sampsize=ALL_sample_size)
ALL_random_rf_500 <- randomForest(x = All_dose_genes_phenos[,..All_random_500_beagle_coord], y = as.factor(All_dose_genes_phenos$AgeClass), importance=TRUE ,proximity=TRUE, mtry=400, ntree=25000, strata=as.factor(All_dose_genes_phenos$AgeClass), sampsize=ALL_sample_size)
ALL_random_rf_400 <- randomForest(x = All_dose_genes_phenos[,..All_random_400_beagle_coord], y = as.factor(All_dose_genes_phenos$AgeClass), importance=TRUE ,proximity=TRUE, mtry=320, ntree=25000, strata=as.factor(All_dose_genes_phenos$AgeClass), sampsize=ALL_sample_size)
ALL_random_rf_300 <- randomForest(x = All_dose_genes_phenos[,..All_random_300_beagle_coord], y = as.factor(All_dose_genes_phenos$AgeClass), importance=TRUE ,proximity=TRUE, mtry=240, ntree=25000, strata=as.factor(All_dose_genes_phenos$AgeClass), sampsize=ALL_sample_size)
ALL_random_rf_200 <- randomForest(x = All_dose_genes_phenos[,..All_random_200_beagle_coord], y = as.factor(All_dose_genes_phenos$AgeClass), importance=TRUE ,proximity=TRUE, mtry=160, ntree=25000, strata=as.factor(All_dose_genes_phenos$AgeClass), sampsize=ALL_sample_size)
ALL_random_rf_100 <- randomForest(x = All_dose_genes_phenos[,..All_random_100_beagle_coord], y = as.factor(All_dose_genes_phenos$AgeClass), importance=TRUE ,proximity=TRUE, mtry=80, ntree=25000, strata=as.factor(All_dose_genes_phenos$AgeClass), sampsize=ALL_sample_size)

All_random_err_rate <- data.frame(cbind(rbind(ALL_random_rf_500$err.rate[25000], ALL_random_rf_400$err.rate[25000], ALL_random_rf_300$err.rate[25000], ALL_random_rf_200$err.rate[25000], ALL_random_rf_100$err.rate[25000]), c(500, 400, 300, 200, 100), 'All_sample_random_SNPs'))
colnames(All_random_err_rate) <- c("OOB", "SNPs", "panel")

I_AM_ERROR <- rbind(All_random_err_rate, All_err_rate, F_QN_err_rate, M_QN_err_rate, F_NL_err_rate, M_NL_err_rate)
I_AM_ERROR$OOB <- as.numeric(I_AM_ERROR$OOB) 
I_AM_ERROR$SNPs <- as.numeric(I_AM_ERROR$SNPs) 

#all this to get this figure!
ggplot() + geom_line(data = I_AM_ERROR, aes(x = SNPs, y = OOB, colour = panel)) +
geom_point(data = I_AM_ERROR, aes(x = SNPs, y = OOB, colour = panel)) +
  theme_classic()

