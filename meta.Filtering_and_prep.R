library(RcppCNPy)
library(tidyverse)
library(data.table)
setwd(dir = "~/Desktop/Working/Sea_Age/")
#write IDs for running ANGSD genotyping and GL
Salmon_Metadata <- read.delim("Salmon_AOM_meta_simple.txt",
                              stringsAsFactors = F)
Salmon_Metadata$FGL_ID
Salmo_bamo <- read.table("Salmo_bamo.txt",
                         quote="\"", comment.char="", stringsAsFactors=FALSE)
colnames(Salmo_bamo) <- "BAM"
Salmo_bamo <- Salmo_bamo %>% mutate(FGL_ID = str_replace(BAM, "NS.*i5.", ""),
                                    FGL_ID = str_replace(FGL_ID, ".realigned.bam", ""))

Salmon_Metadata_Bam <- inner_join(Salmo_bamo, Salmon_Metadata)

Salmon_Metadata_Bam_Unique <- Salmon_Metadata_Bam[!(duplicated(Salmon_Metadata_Bam$FGL_ID)),]

write.table(Salmon_Metadata_Bam_Unique$BAM, "Salmo_bams_all_nodup.txt",
            col.names = F, row.names = F, sep = "\t", quote = F)

#get IDs after depth analysis in ANGSD/vcftools for phasing and GWAS - get vcf order using vcf-query -l on genotyped bcf > vcf file
#drop individuals with low SNP call rates (inferred from mega low depths)
VCF_ID_order<- fread("VCF_ID_order.txt", header = F)
colnames(VCF_ID_order) <- "BAM"

Salmon_Metadata_Bam_Unique_HD <- inner_join(VCF_ID_order, Salmon_Metadata_Bam_Unique)

#Recode
Salmon_Metadata_Bam_Unique_HD$AgeClass[Salmon_Metadata_Bam_Unique_HD$AgeClass=="1SW"] <- 0
Salmon_Metadata_Bam_Unique_HD$AgeClass[Salmon_Metadata_Bam_Unique_HD$AgeClass=="MSW"] <- 1

#Set pop groups
QCNB <- c("MSW", "MUN", "RTG", "TRI")
NL <- c("ENG", "MBB", "SAN", "RVF")

#Write out all pheno data
write.table(Salmon_Metadata_Bam_Unique_HD$AgeClass,
            "All_Ageclass_HD.txt", sep = "\t",
            col.names = F, row.names = F, quote = F)

#Separate by region and sex for VCF filtering, FST in VCFtools, and GWAS
Male_1SW_MSW_QN <- Salmon_Metadata_Bam_Unique_HD %>% filter(PCR_Sex =="M", RiverCode %in% QCNB)
Male_1SW_MSW_NL <-  Salmon_Metadata_Bam_Unique_HD %>% filter(PCR_Sex =="M", RiverCode %in% NL)
Female_1SW_MSW_QN <- Salmon_Metadata_Bam_Unique_HD %>% filter(PCR_Sex =="F", RiverCode %in% QCNB)
Female_1SW_MSW_NL <- Salmon_Metadata_Bam_Unique_HD %>% filter(PCR_Sex =="F", RiverCode %in% NL)


#write out for filtering VCF
write.table(Male_1SW_MSW_QN$BAM,
            "Male_1SW_MSW_QN_bam.txt", sep = "\t",
            col.names = F, row.names = F, quote = F)

write.table(Male_1SW_MSW_NL$BAM,
            "Male_1SW_MSW_NL_bam.txt", sep = "\t",
            col.names = F, row.names = F, quote = F)

write.table(Female_1SW_MSW_QN$BAM, "Female_1SW_MSW_QN_bam.txt", sep =
              "\t", col.names = F, row.names = F, quote = F)

write.table(Female_1SW_MSW_NL$BAM, "Female_1SW_MSW_NL_bam.txt", sep =
              "\t", col.names = F, row.names = F, quote = F)


#check correct orders
Male_1SW_MSW_QN_VCF_order <- fread("Male_1SW_MSW_QN_VCF_order.txt", header = F) %>% select(V1) %>% rename(BAM = V1)
Male_1SW_MSW_QN_VCF_order$BAM==Male_1SW_MSW_QN$BAM
Male_1SW_MSW_NL_VCF_order <- fread("Male_1SW_MSW_NL_VCF_order.txt", header = F) %>% select(V1) %>% rename(BAM = V1)
Male_1SW_MSW_NL_VCF_order$BAM==Male_1SW_MSW_NL$BAM
Female_1SW_MSW_QN_VCF_order <- fread("Female_1SW_MSW_QN_VCF_order.txt", header = F) %>% select(V1) %>% rename(BAM = V1)
Female_1SW_MSW_QN_VCF_order$BAM==Female_1SW_MSW_QN$BAM
Female_1SW_MSW_QN_VCF_order <- fread("Female_1SW_MSW_QN_VCF_order.txt", header = F) %>% select(V1) %>% rename(BAM = V1)
Female_1SW_MSW_QN_VCF_order$BAM==Female_1SW_MSW_QN$BAM


#write out pheno data
write.table(Male_1SW_MSW_QN$AgeClass,
            "Male_1SW_MSW_QN_AgeClass.txt", sep =
              "\t", col.names = F, row.names = F, quote = F)


write.table(Male_1SW_MSW_NL$AgeClass,
            "Male_1SW_MSW_NL_AgeClass.txt", sep =
              "\t", col.names = F, row.names = F, quote = F)

write.table(Female_1SW_MSW_QN$AgeClass,
            "Female_1SW_MSW_QN_AgeClass.txt", sep = "\t", col.names = F, row.names
            = F, quote = F)

write.table(Female_1SW_MSW_NL$AgeClass,
            "Female_1SW_MSW_NL_AgeClass.txt", sep = "\t", col.names = F, row.names
            = F, quote = F)

#Write samples for FST comparisons
Male_1SW_QN <- Male_1SW_MSW_QN %>%  filter(AgeClass %in% 0)
write.table(Male_1SW_QN$BAM, "Male_1SW_QN_bam.txt", sep =
              "\t", col.names = F, row.names = F, quote = F)


Male_MSW_QN <- Male_1SW_MSW_QN %>%  filter(AgeClass %in% 1)
write.table(Male_MSW_QN$BAM, "Male_MSW_QN_bam.txt", sep =
              "\t", col.names = F, row.names = F, quote = F)


Male_1SW_NL <- Male_1SW_MSW_NL %>%  filter(AgeClass %in% 0)
write.table(Male_1SW_NL$BAM, "Male_1SW_NL_bam.txt", sep =
              "\t", col.names = F, row.names = F, quote = F)


Male_MSW_NL <- Male_1SW_MSW_NL %>%  filter(AgeClass %in% 1)
write.table(Male_MSW_NL$BAM, "Male_MSW_NL_bam.txt", sep =
              "\t", col.names = F, row.names = F, quote = F)


Female_1SW_QN <- Female_1SW_MSW_QN %>%  filter(AgeClass %in% 0)
write.table(Female_1SW_QN$BAM, "Female_1SW_QN_bam.txt", sep =
              "\t", col.names = F, row.names = F, quote = F)


Female_MSW_QN <- Female_1SW_MSW_QN %>%  filter(AgeClass %in% 1)
write.table(Female_MSW_QN$BAM, "Female_MSW_QN_bam.txt", sep =
              "\t", col.names = F, row.names = F, quote = F)



Female_1SW_NL <- Female_1SW_MSW_NL %>%  filter(AgeClass %in% 0)
write.table(Female_1SW_NL$BAM, "Female_1SW_NL_bam.txt", sep =
              "\t", col.names = F, row.names = F, quote = F)


Female_MSW_NL <- Female_1SW_MSW_NL %>%  filter(AgeClass %in% 1)
write.table(Female_MSW_NL$BAM, "Female_MSW_NL_bam.txt", sep =
              "\t", col.names = F, row.names = F, quote = F)


SingleSW <- Salmon_Metadata_Bam_Unique_HD %>%  filter(AgeClass %in% 0)
write.table(SingleSW$BAM, "1SW_bam.txt", sep =
              "\t", col.names = F, row.names = F, quote = F)

MSW <- Salmon_Metadata_Bam_Unique_HD %>%  filter(AgeClass %in% 1)
write.table(MSW$BAM, "MSW_bam.txt", sep =
              "\t", col.names = F, row.names = F, quote = F)

#set up files for filtering PCANGSD data
Rep_each <- rep(Salmon_Metadata_Bam_Unique$BAM, each=3)

Beag_Colnames <-  data.frame(t(fread("Beag_Colnames", header = F)))

Inds <- Beag_Colnames[4:1767,]
Ind_BAM <- data.frame(cbind(Inds, Rep_each))
colnames(Ind_BAM) <- c("Ind", "BAM")

Cut_HD <- which(Ind_BAM$BAM %in% Salmon_Metadata_Bam_Unique_HD$BAM) + 3
First3 <-  paste0(rep(1:3))

write.table(c(First3, Cut_HD), "CutHD.csv", col.names = F, row.names = F, quote = F, sep = ",")

Cut_Male_1SW_MSW_QN <- which(Ind_BAM$BAM %in% Male_1SW_MSW_QN$BAM) + 3
write.table(c(First3, Cut_Male_1SW_MSW_QN), "Cut_Male_1SW_MSW_QN.csv", col.names = F, row.names = F, quote = F, sep = ",")

Cut_Male_1SW_MSW_NL <- which(Ind_BAM$BAM %in% Male_1SW_MSW_NL$BAM) + 3
write.table(c(First3, Cut_Male_1SW_MSW_NL), "Cut_Male_1SW_MSW_NL.csv", col.names = F, row.names = F, quote = F, sep = ",")

Cut_Female_1SW_MSW_QN <- which(Ind_BAM$BAM %in% Female_1SW_MSW_QN$BAM ) + 3
write.table(c(First3, Cut_Female_1SW_MSW_QN), "Cut_Female_1SW_MSW_QN.csv", col.names = F, row.names = F, quote = F, sep = ",")

Cut_Female_1SW_MSW_NL<- which(Ind_BAM$BAM %in% Female_1SW_MSW_NL$BAM ) + 3
write.table(c(First3, Cut_Female_1SW_MSW_NL), "Cut_Female_1SW_MSW_NL.csv", col.names = F, row.names = F, quote = F, sep = ",")


#for demographic reconstructions on oh let's say 20 individuals per region 
QNB20 <- Salmon_Metadata_Bam_Unique_HD %>% filter(RiverCode %in% QCNB) %>%  sample_n(20) %>%  select(BAM)



