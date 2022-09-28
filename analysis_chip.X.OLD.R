library(tidyverse)
library(data.table)
library(vegan)
library(qvalue)
setwd("~/Desktop/Working/Sea_Age/Chip_AOM/")

#Read in metadata
Samples_meta_simple <- read.delim("Samples_meta_simple")
Combined_CAN_Inds <- read.table("Salmo_220K_May_Recluster_June_Baseline_2021_nodup_corrected.fam", quote="\"", comment.char="")

Combined_CAN_Inds_wmeta  <- Combined_CAN_Inds[!(Combined_CAN_Inds$V1 %in% Samples_meta_simple$Code),]

Pheno_all <- fread("Sea_Age_220K_Hutchings_data.txt", data.table = F)
Pheno_all <- inner_join(Pheno_all, Samples_meta_simple)
Pheno_grilse_prop <- Pheno_all[!(is.na(Pheno_all$Grilse_prop)),]

Combined_CAN_Inds_wmeta  <- Combined_CAN_Inds[!(Combined_CAN_Inds$V1 %in% Samples_meta_simple$Code),]

#Map
library(marmap)

bathydata <- getNOAA.bathy(-71,-43,40,57, res=3,keep=T)
plot(bathydata)
map=autoplot(bathydata, geom=c("r", "c"), colour="grey", size=0.1) +
  scale_fill_etopo(guide = FALSE) +
  geom_contour(aes(z=z), breaks=c(-100, -200, -500, -1000, -2000, -4000), colour="grey90", size=0.01) +
  xlab("Degrees Longitude") +
  ylab("Degrees Latitude") +
  theme(legend.position = "none")
library(tidyverse)
Samples_meta_simple <- arrange(Samples_meta_simple, Lat, Long)

Samples_meta_simple$SiteName <- factor(Samples_meta_simple$SiteName , levels = Samples_meta_simple$SiteName)
Combined_CAN_Inds_wmeta  <- Combined_CAN_Inds[!(Combined_CAN_Inds$V1 %in% Samples_meta_simple$Code),]

colfunc <- colorRampPalette(c("red", "blue"))
col_pal <- colfunc(26)
map + geom_point(data = Pheno_grilse_prop, aes(x = Long, y =Lat, col = Grilse_prop), size = 3, inherit.aes = F) + theme(legend.position="right") +
  scale_colour_gradient(low="blue", high="red") #+ geom_text(data = Pheno_grilse_prop, aes(x = Long, y = Lat, label =FID))

#bring in genomic and phenotype data
write.table(Pheno_grilse_prop$FID, "Pops_with_grilse_data.txt", col.names = F, row.names = F, sep = "\t", quote = F)
system("~/Desktop/Software/plink_mac_20200219/plink --bfile Salmo_220K_May_Recluster_June_Baseline_2021_nodup_corrected --keep-fam Pops_with_grilse_data.txt --chr-set 30 no-xy  --recode A --maf 0.05 --out Ssa_NAsiteswith_grilse_prop ")
system("~/Desktop/Software/plink_mac_20200219/plink --bfile Salmo_220K_May_Recluster_June_Baseline_2021_nodup_corrected --keep-fam Pops_with_grilse_data.txt --chr-set 30 no-xy --make-bed --maf 0.05 --out Ssa_NAsiteswith_grilse_prop ")

Geno_all <- fread("Ssa_NAsiteswith_grilse_prop.raw", data.table = F)
Pheno_grilse_prop$FID <- Pheno_grilse_prop$Code
Geno_Pheno_grilse_prop <- inner_join(Pheno_grilse_prop, Geno_all)
Grilse_prop_geno_matrix <- Geno_Pheno_grilse_prop[,19:97703]
Grilse_prop_pheno_matrix <- Geno_Pheno_grilse_prop[,1:18]
Grilse_geno_IMP <- apply(Grilse_prop_geno_matrix, 2, function(x) replace(x, is.na(x), as.numeric(names(which.max(table(x))))))


Grilse_prop_table1 <- Grilse_prop_pheno_matrix %>%  select(SiteName, Region, Lat, Long, Grilse_prop) %>% 
  count(SiteName, Grilse_prop, Region, Lat, Long) %>% 
  arrange(Region) 
fwrite(Grilse_prop_table1, "Table1_Chip_GrilseProportion_Sites.txt", col.names = T, row.names = F, sep = "\t", quote = F)

#First we should do a PCA with pcadapt to make sure we are estimating pop structure in LFMM properly
library(pcadapt)
#Entire dataset
Salmo_grilse <- read.pcadapt("Ssa_NAsiteswith_grilse_prop.bed", type = "bed")
PCs <- pcadapt(Salmo_grilse, K =5, min.maf = 0.05)
FAM <- read.delim("Ssa_NAsiteswith_grilse_prop.fam", sep = "", header = F)
colnames(FAM)[1:2] <- c("FID", "IID")
Grilse_meta <- inner_join(FAM[1:2], Grilse_prop_pheno_matrix)
PC_scores_pop <- as.data.frame(cbind(Grilse_meta, PCs$scores))

plot(PCs, option = "manhattan" ) + theme_classic()
plot(PCs, option = "screeplot" ) + theme_classic()
#K5 
PC_scores_pop <- as.data.frame(cbind(Grilse_meta, PCs$scores))

MAP <- fread("Ssa_NAsiteswith_grilse_prop.bim") %>%  select(V1, V4, V2)
colnames(MAP) <- c("CHROM", "BP", "SNP")
#plot by multiple themes
ggplot() + geom_point(data = PC_scores_pop, aes(x = `1`, y = `2`, colour = Grilse_prop)) + theme_classic() + scale_colour_gradient(low="blue", high="red")
ggplot() + geom_point(data = PC_scores_pop, aes(x = `1`, y = `2`, colour = Region)) + theme_classic()
ggplot() + geom_point(data = PC_scores_pop, aes(x = `1`, y = `2`, colour = SiteName)) + theme_classic()

#Get var explained
(PCs$singular.values^2) * 100

#Get outliers
PCA_qvals <- qvalue(PCs$pvalues)

PCA_qvals_map <- data.frame(cbind(MAP, PCA_qvals$qvalues))
colnames(PCA_qvals_map)[4] <- "qval"
PCA_qvals_map_OL <- PCA_qvals_map %>%  filter(qval < 0.05)
PCA_qvals_map_OL %>%  filter(CHROM %in% "1")

six6_SNPs <- PCA_qvals_map %>% filter(CHROM %in% 9, BP > 24902777 , BP <	24905552)

Vgll3_SNPs_100Kwindo<- PCA_qvals_map %>% filter(CHROM %in% 25,BP > (28654947 - 50000)	, BP < (28659019 + 50000))
six6_SNPs_100Kwindo <- PCA_qvals_map %>% filter(CHROM %in% 9, BP > (24902777 - 50000), BP <	(24905552 + 50000))


PCA_qvals_map_DS20K  <- PCA_qvals_map %>%  filter(!SNP %in% c(Vgll3_SNPs_100Kwindo$SNP, six6_SNPs_100Kwindo$SNP, PCA_qvals_map_OL$SNP )) %>% 
  sample_n(20000)

PCA_qvals_map_DS20K_withOLs <-  bind_rows(PCA_qvals_map_DS20K, Vgll3_SNPs_100Kwindo,six6_SNPs_100Kwindo , PCA_qvals_map_OL) %>%  
  distinct()

Vgll3_SNPs_list <- Vgll3_SNPs_100Kwindo %>%  mutate(gene = "vgll3") %>% 
  select(SNP, gene)

six6_SNPs_list <- six6_SNPs_100Kwindo %>%  mutate(gene = "six6") %>% 
  select(SNP, gene)

SNPs_highlight <- bind_rows(Vgll3_SNPs_list, six6_SNPs_list)
#plot
library(ggman)

grilseplot_PCA <- ggman(PCA_qvals_map_DS20K_withOLs, 
                        chrom = "CHROM", pvalue = "qval", snp = "SNP", bp="BP", 
                        pointSize = 1, title = "Proportion grilse per river",
                        xlabel = "Chromosome", ylabel = "-log10(qvalue)", 
                        logTransform = T,  sigLine = 1.3) + theme_classic()
ggmanHighlightGroup(grilseplot_PCA, highlightDfm = SNPs_highlight, snp = "SNP", group = "gene")


#get gene overlap
ssa_chromconvert <- fread("ssa_chromconvert.txt")
PCA_qvals_map$PC1 <- PCs$loadings[,1]
PCA_qvals_map<- inner_join(PCA_qvals_map, ssa_chromconvert)
PCA_PC1_map_bed <- PCA_qvals_map %>%  select(CHROM_nam, BP, PC1) %>% 
  mutate(BP_jr = BP + 1) %>%  select(CHROM_nam, BP, BP_jr, PC1)

fwrite(PCA_PC1_map_bed, "PCA_PC1_map.bed",
      col.names = F, row.names = F, sep = "\t", quote = F )

system("bedtools intersect -a ../SSA_genes.bed -b PCA_PC1_map.bed -wb > PCA_PC1_chipgrilse_geneoverlap")

PCA_qvals_map_OL
PCA_qvals_map_OL<- inner_join(PCA_qvals_map_OL, ssa_chromconvert)
PCA_qvals_map_OL_bed <- PCA_qvals_map_OL %>%  select(CHROM_nam, BP, qval) %>% 
  mutate(BP_jr = BP + 1) %>%  select(CHROM_nam, BP, BP_jr, qval)

fwrite(PCA_qvals_map_OL_bed , "PCA_qvals_map_OL.bed ",
       col.names = F, row.names = F, sep = "\t", quote = F )

system("bedtools intersect -a ../SSA_genes.bed -b PCA_qvals_map_OL.bed -wb > PCA_qvals_map_OL_geneoverlap")

#Check SV overlap
PCA_qvals_map_OL %>%  filter(CHROM %in% "1", BP > 44000000,  BP < 53000000)
PCA_qvals_map_OL %>%  filter(CHROM %in% "23", BP < 8000000)


PCA_qvals_map_OL_geneoverlap <- fread("PCA_qvals_map_OL_geneoverlap") %>%  select(V4) %>%  distinct()
fwrite(PCA_qvals_map_OL_geneoverlap, "PCA_qvals_map_OL_geneoverlap_genes_unique", sep = "\t", col.names = F, row.names = F, quote = F)
PCA_qvals_map_OL_geneoverlap_distinct <- fread("~/Desktop/Working/Sea_Age/Chip_AOM/PCA_qvals_map_OL_geneoverlap") %>%  arrange(V4) %>% 
distinct(V4)

fwrite(PCA_qvals_map_OL_geneoverlap_distinct, "Supplementary_Table_1_Axiom_pcadapt_genes_unique.txt", sep = "\t", col.names = F, row.names = F, )
#LFMM
library(LEA)
library(parallel)
###SNMF #ALL  ###
system("~/Desktop/Software/plink_mac_20200219/plink --bfile ~/Desktop/Working/Salmon_Introgression_2021/CIGENE_2021_Reculster/Salmo_220K_May_Recluster_June_Baseline_2021_nodup_corrected --keep-fam Pops_with_grilse_data.txt --chr-set 30 no-xy  --recode 12 --maf 0.05 --out Ssa_NAsiteswith_grilse_prop_12 ")
ped2lfmm(input.file = "Ssa_NAsiteswith_grilse_prop_12.ped")
pc = pca("Ssa_NAsiteswith_grilse_prop_12.lfmm") 
tc = tracy.widom(pc)
plot(tc$percentage)
project2 = NULL
project2 = snmf(input.file = "Ssa_NAsiteswith_grilse_prop_12.lfmm", K = 1:15, entropy = TRUE, project = "new", CPU = 16 )
plot(project2, col = "blue", pch = 19, cex = 1.2)
best = which.min(cross.entropy(project2, K = 6))


Qval <- data.frame(cbind(Grilse_meta, Q(object = project2, K = 15, run = 1)), stringsAsFactors = F)
Qval_ord <- Qval[order(Qval$Region),]
tbl = Qval_ord
rownames(tbl) <- Qval_ord$ID

plot_data <-  tbl %>% 
  gather('pop', 'prob', V1:V15) %>% 
  group_by(IID)

ggplot(plot_data, aes(IID, prob, fill = pop)) +
  geom_col() +
  facet_grid(~Region, scales = 'free', space = 'free')

#RDA - grilse proportion
Grilse_rda <- rda(Grilse_geno_IMP ~ Grilse_prop_pheno_matrix$Grilse_prop ,scale=T)

#variamce partitioning check with PCs
Grilse_prop_pheno_matrix_PCs <- inner_join(Grilse_prop_pheno_matrix,PC_scores_pop)
colnames(Grilse_prop_pheno_matrix_PCs)[19:23] <- paste0("PC", rep(1:5))
RsquareAdj(Grilse_rda)
anova.cca(Grilse_rda, parallel= 16)

VEEP <- varpart(Grilse_prop_pheno_matrix$Grilse_prop, ~ Grilse_prop_pheno_matrix_PCs$PC1, ~ Grilse_prop_pheno_matrix_PCs$PC2, ~ Grilse_prop_pheno_matrix_PCs$PC5)
plot(VEEP)
VP_sigtest <-rda(Grilse_prop_pheno_matrix$Grilse_prop ~ Grilse_prop_pheno_matrix_PCs$PC1 + Grilse_prop_pheno_matrix_PCs$PC2 + Grilse_prop_pheno_matrix_PCs$PC5)
anova.cca(VP_sigtest, parallel= 16,by = "terms")
RsquareAdj(VP_sigtest)

ecos <- as.factor(Grilse_prop_pheno_matrix$Grilse_prop)
bg <- colfunc(20)
plot(Grilse_rda, type="n", scaling=3, choices = c(1,2))
points(Grilse_rda, display="species", pch=20, cex=0.7, col="gray32", scaling=3, choices = c(1,2))           # the SNPs
points(Grilse_rda, display="sites", pch=21, cex=1.3, col="gray32", scaling=3, bg=bg[ecos], choices = c(1,2)) # the fish
text(Grilse_rda, scaling=3, display="bp", col="red", cex=1)                           # the predictors
legend("bottomright", legend=levels(ecos), bty="n", col="gray32", pch=21, cex=1, pt.bg=bg)


#Get SNP scores
RDA_SNPSCORES <- data.frame(Grilse_rda$CCA$v, stringsAsFactors = F)
RDA_SNPSCORES$SNP <- MAP$SNP

RDA_SNPSCORES_MAPPED <- inner_join(MAP, RDA_SNPSCORES)
RDA_SNPSCORES_MAPPED$RDA1_abs <-  abs(as.numeric(as.character(RDA_SNPSCORES$RDA1)))

#RDA outliers
RDA_grilse_OL_99 <- RDA_SNPSCORES_MAPPED$SNP[which(RDA_SNPSCORES_MAPPED$RDA1_abs > quantile(RDA_SNPSCORES_MAPPED$RDA1_abs, 0.99))]

RDA_SNPSCORES_MAPPED_OL <- RDA_SNPSCORES_MAPPED %>%  filter(SNP %in%RDA_grilse_OL_99)
#gene overlap
RDA_SNPSCORES_MAPPED_OL  <- inner_join(RDA_SNPSCORES_MAPPED_OL, ssa_chromconvert)
RDA_SNPSCORES_MAPPED_OL_bed   <-  RDA_SNPSCORES_MAPPED_OL   %>%  select(CHROM_nam, BP, RDA1_abs) %>% 
  mutate(BP_jr = BP + 1) %>%  select(CHROM_nam, BP, BP_jr, RDA1_abs)

fwrite(RDA_SNPSCORES_MAPPED_OL_bed, "RDA_SNPSCORES_MAPPED_OL.bed",
       col.names = F, row.names = F, sep = "\t", quote = F )

system("bedtools intersect -a ../SSA_genes.bed -b RDA_SNPSCORES_MAPPED_OL.bed -wb > RDA_SNPSCORES_MAPPED_OL_chip_geneoverlap")
RDA_SNPSCORES_MAPPED_OL_chip_geneoverlap <- fread("RDA_SNPSCORES_MAPPED_OL_chip_geneoverlap", header = F) %>% select(V1, V2, V3, V4, V8)
colnames(RDA_SNPSCORES_MAPPED_OL_chip_geneoverlap) <- c("CHROM", "POS", "BPtoo", "Symbol", "RDA1")
arrange(RDA_SNPSCORES_MAPPED_OL_chip_geneoverlap, desc(RDA1))

OL_grileK2_chipgrilse_geneoverlap %>% select(Symbol) %>% distinct()




##RIDGE REGRESSION LFMM - do K 1, 2, 5
library(lfmm)
library(LEA)
## ## Fit an LFMM, i.e, compute B, U, V estimates:
mod.lfmm <- lfmm_ridge(Y = Grilse_geno_IMP, 
                       X = Grilse_prop_pheno_matrix$Grilse_prop, 
                       K = 1)
## performs association testing using the fitted model:
pv <- lfmm_test(Y = Grilse_geno_IMP, 
                X = Grilse_prop_pheno_matrix$Grilse_prop, 
                lfmm = mod.lfmm, 
                calibrate = "gif")

pvalues <- pv$calibrated.pvalue 
qvalues <- qvalue(pv$calibrated.pvalue)
grilse_LFMM_K1 <- data.frame(cbind(rownames(pvalues), pvalues, qvalues$qvalues))
colnames(grilse_LFMM_K1) <- c("SNP", "pval", "qval")
grilse_LFMM_K1$qval <- as.numeric(as.character((grilse_LFMM_K1$qval)))
grilse_LFMM_K1$pval <- as.numeric(as.character((grilse_LFMM_K1$pval)))
grilse_LFMM_K1  <- grilse_LFMM_K1 %>%  mutate(SNP = str_replace(SNP, "_.*", ""))
grilse_LFMM_K1 <- inner_join(grilse_LFMM_K1, MAP)

grilse_LFMM_OL_K1 <- grilse_LFMM_K1 %>% filter(qval < 0.05)
Vgll3_SNPs_100Kwindo<- grilse_LFMM_K1  %>% filter(CHROM %in% 25,BP > (28654947 - 50000)	, BP < (28659019 + 50000))
six6_SNPs_100Kwindo <- grilse_LFMM_K1  %>% filter(CHROM %in% 9, BP > (24902777 - 50000), BP <	(24905552 + 50000))

qqplot(rexp(length(pvalues), rate = log(10)),
       -log10(pvalues), xlab = "Expected quantile",
       pch = 19, cex = .4)
abline(0,1)

LFMM_K1_ds  <- grilse_LFMM_K1 %>% filter(!SNP %in% c(Vgll3_SNPs_100Kwindo$SNP, six6_SNPs_100Kwindo$SNP, grilse_LFMM_OL_K1$SNP)) %>% 
  sample_n(20000)

LFMM_K1_ds_OLs  <- bind_rows(LFMM_K1_ds, Vgll3_SNPs_100Kwindo, six6_SNPs_100Kwindo, grilse_LFMM_OL_K1) %>% 
  distinct()

#plot

grilseplot_LFMMK1 <- ggman(LFMM_K1_ds_OLs , 
                        chrom = "CHROM", pvalue = "qval", snp = "SNP", bp="BP", 
                        pointSize = 1, title = "Proportion grilse per river",
                        xlabel = "Chromosome", ylabel = "-log10(qvalue)", 
                        logTransform = T,  sigLine = -log10(0.05)) + theme_classic()
ggmanHighlightGroup(grilseplot_LFMMK1, highlightDfm = SNPs_highlight, snp = "SNP", group = "gene")

grilse_LFMM_OL_K1$SNP[grilse_LFMM_OL_K1$SNP %in% PCA_qvals_map_OL$SNP]
grilse_LFMM_OL_K1$SNP[grilse_LFMM_OL_K1$SNP %in% RDA_grilse_OL_99]



#gene overlap
grilse_LFMM_OL_K1 <- inner_join(grilse_LFMM_OL_K1, ssa_chromconvert)
grilse_LFMM_OL_K1_bed <-  grilse_LFMM_OL_K1 %>%  select(CHROM_nam, BP, qval) %>% 
  mutate(BP_jr = BP + 1) %>%  select(CHROM_nam, BP, BP_jr, qval)

fwrite(grilse_LFMM_OL_K1_bed, "grilse_LFMM_OL_K1.bed",
       col.names = F, row.names = F, sep = "\t", quote = F )

system("bedtools intersect -a ../SSA_genes.bed -b grilse_LFMM_OL_K1.bed -wb > OL_grileK1_chipgrilse_geneoverlap")
K1_chipgrilse_geneoverlap <- fread("grileK1_chipgrilse_geneoverlap", header = F) %>% select(V1, V2, V3, V4, V8)
colnames(K1_chipgrilse_geneoverlap) <- c("CHROM", "POS", "BPtoo", "Symbol", "K1_log10p")




grilse_LFMM_K1$log10P <- -log10(grilse_LFMM_K1$pval)
grilse_LFMM_K1 <- inner_join(grilse_LFMM_K1, ssa_chromconvert)
grilse_LFMM_K1_bed <-  grilse_LFMM_K1 %>%  select(CHROM_nam, BP, log10P) %>% 
  mutate(BP_jr = BP + 1) %>%  select(CHROM_nam, BP, BP_jr, log10P)

fwrite(grilse_LFMM_K1_bed, "grilse_LFMM_K1.bed",
       col.names = F, row.names = F, sep = "\t", quote = F )

system("bedtools intersect -a ../SSA_genes.bed -b grilse_LFMM_K1.bed -wb > grileK1_chipgrilse_geneoverlap")
K1_chipgrilse_geneoverlap <- fread("grileK1_chipgrilse_geneoverlap", header = F) %>% select(V1, V2, V3, V4, V8)
colnames(K1_chipgrilse_geneoverlap) <- c("CHROM", "POS", "BPtoo", "Symbol", "K1_log10p")

#K2
## ## Fit an LFMM, i.e, compute B, U, V estimates:
mod.lfmm <- lfmm_ridge(Y = Grilse_geno_IMP, 
                       X = Grilse_prop_pheno_matrix$Grilse_prop, 
                       K = 2)
## performs association testing using the fitted model:
pv <- lfmm_test(Y = Grilse_geno_IMP, 
                X = Grilse_prop_pheno_matrix$Grilse_prop, 
                lfmm = mod.lfmm, 
                calibrate = "gif")

pvalues <- pv$calibrated.pvalue 
qvalues <- qvalue(pv$calibrated.pvalue)
grilse_LFMM_K2 <- data.frame(cbind(rownames(pvalues), pvalues, qvalues$qvalues))
colnames(grilse_LFMM_K2) <- c("SNP", "pval", "qval")
grilse_LFMM_K2$qval <- as.numeric(as.character((grilse_LFMM_K2$qval)))
grilse_LFMM_K2$pval <- as.numeric(as.character((grilse_LFMM_K2$pval)))
grilse_LFMM_K2  <- grilse_LFMM_K2 %>%  mutate(SNP = str_replace(SNP, "_.*", ""))
grilse_LFMM_K2 <- inner_join(grilse_LFMM_K2, MAP)

grilse_LFMM_OL_K2 <- grilse_LFMM_K2 %>% filter(qval < 0.05)
Vgll3_SNPs_100Kwindo<- grilse_LFMM_K2  %>% filter(CHROM %in% 25,BP > (28654947 - 50000)	, BP < (28659019 + 50000))
six6_SNPs_100Kwindo <- grilse_LFMM_K2  %>% filter(CHROM %in% 9, BP > (24902777 - 50000), BP <	(24905552 + 50000))

qqplot(rexp(length(pvalues), rate = log(10)),
       -log10(pvalues), xlab = "Expected quantile",
       pch = 19, cex = .4)
abline(0,1)

LFMM_K2_ds  <- grilse_LFMM_K2 %>% filter(!SNP %in% c(Vgll3_SNPs_100Kwindo$SNP, six6_SNPs_100Kwindo$SNP, grilse_LFMM_OL_K2$SNP)) %>% 
  sample_n(20000)

LFMM_K2_ds_OLs  <- bind_rows(LFMM_K2_ds, Vgll3_SNPs_100Kwindo, six6_SNPs_100Kwindo, grilse_LFMM_OL_K2) %>% 
  distinct()

#plot

grilseplot_LFMMK2 <- ggman(LFMM_K2_ds_OLs , 
                           chrom = "CHROM", pvalue = "qval", snp = "SNP", bp="BP", 
                           pointSize = 1, title = "Proportion grilse per river",
                           xlabel = "Chromosome", ylabel = "-log10(qvalue)", 
                           logTransform = T,  sigLine = -log10(0.05)) + theme_classic()
ggmanHighlightGroup(grilseplot_LFMMK2, highlightDfm = SNPs_highlight, snp = "SNP", group = "gene")

grilse_LFMM_OL_K2$SNP[grilse_LFMM_OL_K2$SNP %in% PCA_qvals_map_OL$SNP]
grilse_LFMM_OL_K2$SNP[grilse_LFMM_OL_K2$SNP %in% RDA_grilse_OL_99]


#gene overlap
#gene overlap
grilse_LFMM_OL_K2 <- inner_join(grilse_LFMM_OL_K2, ssa_chromconvert)
grilse_LFMM_OL_K2_bed <-  grilse_LFMM_OL_K2 %>%  select(CHROM_nam, BP, qval) %>% 
  mutate(BP_jr = BP + 1) %>%  select(CHROM_nam, BP, BP_jr, qval)

fwrite(grilse_LFMM_OL_K2_bed, "grilse_LFMM_OL_K2.bed",
       col.names = F, row.names = F, sep = "\t", quote = F )

system("bedtools intersect -a ../SSA_genes.bed -b grilse_LFMM_OL_K2.bed -wb > OL_grileK2_chipgrilse_geneoverlap")
OL_grileK2_chipgrilse_geneoverlap <- fread("OL_grileK2_chipgrilse_geneoverlap", header = F) %>% select(V1, V2, V3, V4, V8)
colnames(OL_grileK2_chipgrilse_geneoverlap) <- c("CHROM", "POS", "BPtoo", "Symbol", "qval")
arrange(OL_grileK2_chipgrilse_geneoverlap, qval)

OL_grileK2_chipgrilse_geneoverlap %>% select(Symbol) %>% distinct()

grilse_LFMM_K2$log10P <- -log10(grilse_LFMM_K2$pval)
grilse_LFMM_K2 <- inner_join(grilse_LFMM_K2, ssa_chromconvert)
grilse_LFMM_K2_bed <-  grilse_LFMM_K2 %>%  select(CHROM_nam, BP, log10P) %>% 
  mutate(BP_jr = BP + 1) %>%  select(CHROM_nam, BP, BP_jr, log10P)

fwrite(grilse_LFMM_K2_bed, "grilse_LFMM_K2.bed",
       col.names = F, row.names = F, sep = "\t", quote = F )

system("bedtools intersect -a ../SSA_genes.bed -b grilse_LFMM_K2.bed -wb > grileK2_chipgrilse_geneoverlap")
K2_chipgrilse_geneoverlap <- fread("grileK2_chipgrilse_geneoverlap", header = F) %>% select(V1, V2, V3, V4, V8)
colnames(K2_chipgrilse_geneoverlap) <- c("CHROM", "POS", "BPtoo", "Symbol", "K2_log10p")


K1K2_geneoverlap <- inner_join(K2_chipgrilse_geneoverlap, K1_chipgrilse_geneoverlap)
Allstat_overlap <- inner_join(PCA_chipgrilse_geneoverlap, K1K2_geneoverlap)
write.table(Allstat_overlap, "Chip_allstat_overlap.tsv", col.names = T, row.names = F, sep = "\t", quote = F)





#prep for polysel
Chip_topPC1_genes <- Allstat_overlap %>% 
  group_by(Symbol) %>%
  filter(PC1 == max(PC1)) 
Chip_topPC1_genes <- as.data.frame(unique(setDT(Chip_topPC1_genes), by = "Symbol"))

Chip_topK1_log10p_genes <- Allstat_overlap %>% 
  group_by(Symbol) %>%
  filter(K1_log10p == max(K1_log10p)) 
Chip_topK1_log10p_genes <- as.data.frame(unique(setDT(Chip_topK1_log10p_genes), by = "Symbol"))

Chip_topK2_log10p_genes <- Allstat_overlap %>% 
  group_by(Symbol) %>%
  filter(K2_log10p == max(K2_log10p)) 
Chip_topK2_log10p_genes <- as.data.frame(unique(setDT(Chip_topK2_log10p_genes), by = "Symbol"))


#get gene info
ssalar_genes <- fread("../ssalar_gene_result.txt")
Obj_Info <- inner_join(Chip_topK1_log10p_genes, ssalar_genes)

Obj_Info$SNPcount <- 1
colnames(Obj_Info)
Obj_Info$GeneLength <- Obj_Info$end_position_on_the_genomic_accession - Obj_Info$start_position_on_the_genomic_accession
#Make an object info file with charr positions and FST, but used humann gene IDs
Obj_Info <- data.frame(cbind(Obj_Info$GeneID, Obj_Info$K2_log10p, Obj_Info$Symbol, Obj_Info$GeneLength, Obj_Info$CHROM, Obj_Info$start_position_on_the_genomic_accession, Obj_Info$end_position_on_the_genomic_accession, Obj_Info$orientation), stringsAsFactors = F)
colnames(Obj_Info) <- c("objID",	"objStat",	"objName",	"GeneLength",	"chr",	"startpos",	"endpos",	"strand")
Obj_Info$strand <- NA

#SetInfo - TBD

SetInfo <- read.csv("../ssalar_kegg.csv")

colnames(SetInfo) <- c("setID", "setName", "setSource")


#SetObj
#turns into SetObj after filtering by GeneIDs in Salmon dataset
SetObj <- read.delim("../biosystems_gene", header=FALSE)
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

system("mkdir ~/Desktop/Software/polysel/data/salmon_chip_K2/")

write.table(Obj_Info, "~/Desktop/Software/polysel/data/salmon_chip_K2/ObjInfo.txt", col.names = T, row.names = F, sep = "\t", quote = F)
write.table(SetInfo, "~/Desktop/Software/polysel/data/salmon_chip_K2/SetInfo.txt", col.names = T, row.names = F, sep = "\t", quote = F)
write.table(SetObj, "~/Desktop/Software/polysel/data/salmon_chip_K2/SetObj.txt", col.names = T, row.names = F, sep = "\t", quote = F)




