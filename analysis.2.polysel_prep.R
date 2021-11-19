###GET READY######

library(RcppCNPy)
library(tidyverse)
library(data.table)
library(qvalue)


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


#get GWAS PCA FST stats
PC_FST_GWA_PCcor_MSW_1SW_gene_overlap <- fread("All_2PCcorr/PC_FST_GWA_PCcor_MSW1SW_gene_overlap_stats.tsv") %>%  
  mutate(Symbol = gene ) %>% select(-gene)

FQN_PC_FST_GWA_PCcorr_gene_overlap <- fread("F_QN_GWAS/FQN_PC_FST_GWA_PCcor_gene_overlap_stats.tsv") %>% 
  mutate(Symbol = gene ) %>% select(-gene)

FNL_PC_FST_GWA_PCcorr_gene_overlap <- fread("F_NL_GWAS/FNL_PC_FST_GWA_PCcor_gene_overlap_stats.tsv") %>% 
  mutate(Symbol = gene )%>% select(-gene)
  
MQN_PC_FST_GWA_PCcorr_gene_overlap <- fread("M_QN_GWAS/MQN_PC_FST_GWA_PCcor_gene_overlap_stats.tsv")  %>%  
  mutate(Symbol = gene ) %>% select(-gene)

MNL_PC_FST_GWA_PCcorr_gene_overlap <- fread("M_NL_GWAS/MNL_PC_FST_GWA_PCcor_gene_overlap_stats.tsv")  %>%  
  mutate(Symbol = gene ) %>% select(-gene)


#get top loci for each group, by GWAS score - get top 1000 first to filter beagle via command line, can downsample in R from there
All_PCcor_topLRT_genes <- PC_FST_GWA_PCcor_MSW_1SW_gene_overlap %>%
  group_by(Symbol) %>%
  filter(LRT == max(LRT)) 
All_PCcor_topLRT_genes <- as.data.frame(unique(setDT(All_PCcor_topLRT_genes), by = "Symbol"))

All_PCcor_topPC1_genes <- PC_FST_GWA_PCcor_MSW_1SW_gene_overlap %>%
  group_by(Symbol) %>%
  filter(PC1 == max(PC1))
All_PCcor_topPC1_genes <- as.data.frame(unique(setDT(All_PCcor_topPC1_genes), by = "Symbol"))

FQN_PC_FST_GWA_PCcorr_topLRT_genes <- FQN_PC_FST_GWA_PCcorr_gene_overlap %>% 
  group_by(Symbol) %>%
  filter(LRT == max(LRT)) 
FQN_PC_FST_GWA_PCcorr_topLRT_genes <- as.data.frame(unique(setDT(FQN_PC_FST_GWA_PCcorr_topLRT_genes), by = "Symbol"))

FQN_PC_FST_GWA_PCcorr_topPC1_genes <- FQN_PC_FST_GWA_PCcorr_gene_overlap %>% 
  group_by(Symbol) %>%
  filter(PC1 == max(PC1)) 
FQN_PC_FST_GWA_PCcorr_topPC1_genes<- as.data.frame(unique(setDT(FQN_PC_FST_GWA_PCcorr_topPC1_genes), by = "Symbol"))

FNL_PC_FST_GWA_PCcorr_topLRT_genes <- FNL_PC_FST_GWA_PCcorr_gene_overlap %>% 
  group_by(Symbol) %>%
  filter(LRT == max(LRT)) 
FNL_PC_FST_GWA_PCcorr_topLRT_genes <- as.data.frame(unique(setDT(FNL_PC_FST_GWA_PCcorr_topLRT_genes), by = "Symbol"))

FNL_PC_FST_GWA_PCcorr_topPC1_genes <- FNL_PC_FST_GWA_PCcorr_gene_overlap %>% 
  group_by(Symbol) %>%
  filter(PC1 == max(PC1)) 
FNL_PC_FST_GWA_PCcorr_topPC1_genes <- as.data.frame(unique(setDT(FNL_PC_FST_GWA_PCcorr_topPC1_genes), by = "Symbol"))

MQN_PC_FST_GWA_PCcorr_topLRT_genes <- MQN_PC_FST_GWA_PCcorr_gene_overlap %>% 
  group_by(Symbol) %>%
  filter(LRT == max(LRT)) 
MQN_PC_FST_GWA_PCcorr_topLRT_genes <- as.data.frame(unique(setDT(MQN_PC_FST_GWA_PCcorr_topLRT_genes), by = "Symbol"))

MQN_PC_FST_GWA_PCcorr_topPC1_genes <- MQN_PC_FST_GWA_PCcorr_gene_overlap %>% 
  group_by(Symbol) %>%
  filter(PC1 == max(PC1)) 
MQN_PC_FST_GWA_PCcorr_topPC1_genes <- as.data.frame(unique(setDT(MQN_PC_FST_GWA_PCcorr_topPC1_genes), by = "Symbol"))

MNL_PC_FST_GWA_PCcorr_topLRT_genes <- MNL_PC_FST_GWA_PCcorr_gene_overlap %>% 
  group_by(Symbol) %>%
  filter(LRT == max(LRT)) 
MNL_PC_FST_GWA_PCcorr_topLRT_genes <- as.data.frame(unique(setDT(MNL_PC_FST_GWA_PCcorr_topLRT_genes), by = "Symbol"))

MNL_PC_FST_GWA_PCcorr_topPC1_genes <- MNL_PC_FST_GWA_PCcorr_gene_overlap %>% 
  group_by(Symbol) %>%
  filter(PC1 == max(PC1)) 
MNL_PC_FST_GWA_PCcorr_topPC1_genes <- as.data.frame(unique(setDT(MNL_PC_FST_GWA_PCcorr_topPC1_genes), by = "Symbol"))

#get gene info
ssalar_genes <- fread("ssalar_gene_result.txt")
Obj_Info <- inner_join(MQN_PC_FST_GWA_PCcorr_topPC1_genes, ssalar_genes)

Obj_Info$SNPcount <- 1
colnames(Obj_Info)
Obj_Info$GeneLength <- Obj_Info$end_position_on_the_genomic_accession - Obj_Info$start_position_on_the_genomic_accession
#Make an object info file with charr positions and FST, but used humann gene IDs
Obj_Info <- data.frame(cbind(Obj_Info$GeneID, Obj_Info$LRT, Obj_Info$Symbol, Obj_Info$GeneLength, Obj_Info$CHROM, Obj_Info$start_position_on_the_genomic_accession, Obj_Info$end_position_on_the_genomic_accession, Obj_Info$orientation), stringsAsFactors = F)
colnames(Obj_Info) <- c("objID",	"objStat",	"objName",	"GeneLength",	"chr",	"startpos",	"endpos",	"strand")
Obj_Info$strand <- NA

#SetInfo - TBD

SetInfo <- read.csv("ssalar_kegg.csv")

colnames(SetInfo) <- c("setID", "setName", "setSource")


#SetObj
#turns into SetObj after filtering by GeneIDs in Salmon dataset
SetObj <- read.delim("biosystems_gene", header=FALSE)
colnames(SetObj) <- c("setID",	"objID")
SetObj <- SetObj[1:2]

SetObj <- SetObj[SetObj$setID %in% SetInfo$setID,]
colnames(SetObj)
colnames(Obj_Info)


#ensure consistency of genes, pathways across comparisons
Obj_Info <- Obj_Info[Obj_Info$objID  %in% SetObj$objID,]
SetObj <- SetObj[SetObj$objID  %in% Obj_Info$objID,]
SetInfo <- SetInfo[SetInfo$setID %in% SetObj$setID,]

#Just go through and ensure class consistency for all datasets compared to PolyLinkR example
Obj_Info$objID <- as.integer(Obj_Info$objID)
Obj_Info$objStat <- as.numeric(Obj_Info$objStat)
Obj_Info$GeneLength <- as.integer(Obj_Info$GeneLength)
Obj_Info$startpos <- as.integer(Obj_Info$startpos)
Obj_Info$endpos <- as.integer(Obj_Info$endpos)
setDT(Obj_Info)

SetInfo$setID <- as.integer(SetInfo$setID)
setDT(SetInfo)

SetObj$setID 

class(SetObj$objID)
setDT(SetObj)
class(SetObj)

system("mkdir ~/Desktop/Software/polysel/data/salmon_MQN_PC1/ ")

write.table(Obj_Info, "~/Desktop/Software/polysel/data/salmon_MQN_PC1/ObjInfo.txt", col.names = T, row.names = F, sep = "\t", quote = F)
write.table(SetInfo, "~/Desktop/Software/polysel/data/salmon_MQN_PC1/SetInfo.txt", col.names = T, row.names = F, sep = "\t", quote = F)
write.table(SetObj, "~/Desktop/Software/polysel/data/salmon_MQN_PC1/SetObj.txt", col.names = T, row.names = F, sep = "\t", quote = F)



