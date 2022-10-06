#load libraries
#data manip and plot
library(data.table)  
library(tidyverse)
library(ggman)
library(marmap)
library(parallel)

#pop genomics
library(pcadapt)
library(LEA)
library(lfmm)
library(qvalue)
library(vegan)

#set a working dir
setwd("~/Desktop/Working/chip_2022_seaage/")

#### data required ####
##plink format files - .bed 
                      ##.ped in 12 format with --recode 12 option in plink, 
                      ###.raw format --recode A option in plink
                      ##.fam with individual info
                      ### .bim for marker information

##metadata with grilse proportion per river, latitude and longitude

##for polysel - ssalar_genes_result.txt gene info and ssalar_kegg.csv kegg pathways, biosystems_gene file

#read in metadata
grilse.meta <- fread("salmo_220K_site_group_meta_2022_NA_grilse_DFOregions.tsv")

#Map
bathydata <- getNOAA.bathy(-71,-43,40,57, res=3,keep=T)
plot(bathydata)
map=autoplot(bathydata, geom=c("r", "c"), colour="grey", size=0.1) +
  scale_fill_etopo(guide = FALSE) +
  geom_contour(aes(z=z), breaks=c(-100, -200, -500, -1000, -2000, -4000), colour="grey90", size=0.01) +
  xlab("Degrees Longitude") +
  ylab("Degrees Latitude") +
  theme(legend.position = "none")

map + geom_point(data = grilse.meta, aes(x = Lon, y =Lat, col = p1SW), size = 3, inherit.aes = F) + theme(legend.position="right") +
  scale_colour_gradient(low="blue", high="red") #+ geom_text(data = Pheno_grilse_prop, aes(x = Long, y = Lat, label =FID))

map + geom_point(data = grilse.meta, aes(x = Lon, y =Lat, col = Region), size = 3, inherit.aes = F) + theme(legend.position="right")


#match metadata to individual IDs in order of genotype file
grilse.meta.fam <- fread("Salmo_220K_merge2022_grilse.fam") %>%  mutate(SiteCode = V1, ID = V2) %>%  
  select(SiteCode, ID) %>% 
  inner_join(., grilse.meta)

#get SNP map
chrom_snp_map <- fread("Salmo_220K_merge2022_grilse.bim") %>%  
  mutate(CHROM =V1, SNP = V2, BP = V4) %>%  
  select(CHROM, BP, SNP)

#get alt chromosome IDs
chrom_snp_map <- fread("ssa_chromconvert.txt") %>% 
  inner_join(., chrom_snp_map)

#pop structure inference ##
#PCA functions
PCADAPT_import_PCA <- function(filename, Knum) {
  to_pc <- read.pcadapt(paste0(filename, ".bed"), type = "bed")
  PCA <- pcadapt(to_pc, K = Knum, min.maf = 0.05)
  return(PCA)
} # a function to import data and run a PCA

PCADAPT_scores_meta <- function(PCAobj, filename, Knum, meta) {
  FAM <- fread(paste0(filename, ".fam")) %>% 
    select(V1,V2) %>% 
    mutate(SiteCode = V1, ID = V2) %>%  
    select(-V1, -V2)
  meta_ordered <- inner_join(FAM, meta)
  PCAscores <- data.frame(PCAobj$scores)
  colnames(PCAscores) <- paste0("PC", rep(1:Knum))
  PCA_meta <- bind_cols(meta_ordered, PCAscores)
  return(PCA_meta)
} ## a function to match PCA scores to metadata

PCADAPT_var_per_axis <- function(PCAobj) {
  var_per_axis <- (PCAobj$singular.values^2)
  return(var_per_axis)} #a function to get PCA per axis variation

#run PCA
PC_K5 <- PCADAPT_import_PCA(filename = "Salmo_220K_merge2022_grilse", Knum = 5)
PC_K10 <- PCADAPT_import_PCA(filename = "Salmo_220K_merge2022_grilse", Knum = 10)

#screeplot
plot(PC_K10, option = "screeplot") + theme_classic()  #plateaus at ~K = 5

##Population structure analysis with PCA
PC_K5_meta <-  PCADAPT_scores_meta(PCAobj = PC_K5, filename =  "Salmo_220K_merge2022_grilse", Knum = 5, meta = grilse.meta.fam)
PCADAPT_var_per_axis(PC_K5)


ggplot() + 
  geom_point(data = PC_K5_meta, aes(x = PC1, y = PC2, colour = Region)) + 
  theme_classic()

ggplot() + 
  geom_point(data = PC_K5_meta, aes(x = PC1, y = PC2, colour = p1SW)) + 
  theme_classic() + scale_colour_gradient(low = "blue", high = "red")

##LEA snmf
###SNMF #ALL  ###

ped2lfmm(input.file = "Salmo_220K_merge2022_grilse.ped") #conversion only needs to be done once

#snmf for admixture proportions
grilse.lea = NULL
grilse.lea= snmf(input.file = "Salmo_220K_merge2022_grilse.lfmm", K = 1:15, entropy = TRUE, project = "new", CPU = 16 )

#plot CV scores
plot(grilse.lea, col = "blue", pch = 19, cex = 1.2)  #big declines until ~6? will use 5 for consistency, could use 10.. art not a science, it's just resolving rivers

snmf_Qscores_meta <- function(leaproj, Knum, meta) {
  leaQ <-Q(object = leaproj, K = Knum, run = 1)
  colnames(leaQ) <- paste0("Q", rep(1:Knum))
  leaQ_meta <- bind_cols(meta, leaQ)
  return(leaQ_meta) }  # a function to import Q scores from snmf and match to metadata

grilse.snmf.K5.meta<- snmf_Qscores_meta(leaproj = grilse.lea, Knum = 5, meta = grilse.meta.fam)

Make_admix_table<-  function(Admix_meta_object, Knum){Admix_table <- Admix_meta_object
rownames(Admix_table) <-Admix_table$ID
plot_data <-  Admix_table  %>% 
  gather('pop', 'prob', Q1:paste0("Q", Knum))
return(plot_data)} #a function to make a plottable table for snmf 

grilse.snmf.K5.admixtable <- Make_admix_table(grilse.snmf.K5.meta, Knum = 5)

ggplot(grilse.snmf.K5.admixtable, aes(ID, prob, fill = pop)) +
  geom_col() +
  facet_grid(~Region, scales = 'free', space = 'free')

#RDA - grilse proportion
#get geno matrix for RDA
genos.dose <- fread("Salmo_220K_merge2022_grilse.raw", sep = " ") %>%  select (-FID, -IID, -PAT, -MAT, -SEX, -PHENOTYPE)
genos.dose.imp<- apply(genos.dose, 2, function(x) replace(x, is.na(x), as.numeric(names(which.max(table(x))))))

##RDA
grilse.rda <- rda(genos.dose.imp ~ grilse.meta.fam$p1SW  ,scale=T )

#adjusted R2 and significance check 
RsquareAdj(grilse.rda)

# estimate variance partitioning between trait and pop structure PCs
grilse.pc.varpart <- varpart(grilse.meta.fam$p1SW, ~ PC_K5_meta$PC1, ~ PC_K5_meta$PC4, ~ PC_K5_meta$PC5) ##here PC1 explains lots of variance...
plot(grilse.pc.varpart)
VP_sigtest <-rda(grilse.meta.fam$p1SW ~ PC_K5_meta$PC1 + PC_K5_meta$PC4 + PC_K5_meta$PC5)
anova.cca(VP_sigtest, parallel= 16,by = "terms")
RsquareAdj(VP_sigtest)

####SNP associations with population structure and grilse proportion
##PCA
PCADAPT_scores_qval_map  <- function(PCAobj, mapobj, Knum) {
  qvalues <- qvalue(PCAobj$pvalues)
  PCAobj_map_pval_qval <- bind_cols(chrom_snp_map, PCAobj$pvalues, qvalues$qvalues)
  colnames(PCAobj_map_pval_qval)[5:6] <- c("pvalues", "qvalues")
  loadings <- data.frame(PCAobj$loadings)
  names(loadings) <- paste0("PC", rep(1:Knum), "_loading")
  PCAobj_map_pval_qval_loadings <- bind_cols(PCAobj_map_pval_qval, loadings)
  return( PCAobj_map_pval_qval_loadings)} #a function to get PCA pvalues, qvalues, loadings and SNP locations

PC_K5_scores <- PCADAPT_scores_qval_map(PCAobj = PC_K5, mapobj = chrom_snp_map, Knum = 5)

#make a downsampled object for plotting
PC_K5_scores_downsamp <-PC_K5_scores %>%   sample_n(10000)

#get PC outliers and combine with downsamp for plot comparison
PC_K5_qval05_OL <- PC_K5_scores %>%  filter(qvalues < 0.05)  

PC_K5_qval05_OL_downsamp<-  PC_K5_qval05_OL %>%  bind_rows(., PC_K5_scores_downsamp) %>%  
  distinct()

#get six6, vgll3, SV 1-23, location and make list for plotting
#vgll3 - 100 kbp due to low marker density
vgll3_100kbwindow_SNPs_list <- PC_K5_scores  %>% filter(CHROM %in% 25,BP > (28654947 - 50000)	, BP < (28659019 + 50000)) %>% 
  mutate(region = "vgll3") %>% 
  select(SNP, region)

#sic6 - 100 kbp due to low marker density
six6_100kbwindow_SNPs_list <- PC_K5_scores %>% filter(CHROM %in% 9, BP > (24902777 - 50000), BP <	(24905552 + 50000)) %>% 
  mutate(region = "six6") %>% 
  select(SNP, region)

#EU SV chroms 1-23
SV_1<- PC_K5_scores %>%  filter(CHROM %in% "1", BP > 44000000,  BP < 53000000)
SV_23<- PC_K5_scores %>%  filter(CHROM %in% "23", BP < 8000000)
SV1_23 <- bind_rows(SV_1, SV_23) %>% mutate(region = "SV1_23") %>% 
  select(SNP, region)


SNPs_highlight <- bind_rows(six6_100kbwindow_SNPs_list, 
                            vgll3_100kbwindow_SNPs_list,
                            SV1_23)

#plot q vals
grilseplot_PCA_K5qval <- ggman(PC_K5_qval05_OL_downsamp, 
                        chrom = "CHROM", pvalue = "qvalues", snp = "SNP", bp="BP", 
                        pointSize = 1, title = "PCA K5",
                        xlabel = "Chromosome", ylabel = "-log10(qvalue)", 
                        logTransform = T,  sigLine = 1.3) + theme_classic()
ggmanHighlightGroup(grilseplot_PCA_K5qval, highlightDfm = SNPs_highlight, snp = "SNP", group = "region")

#get proportion of PC outliers overlapping gene/SVs
SNPs_highlight %>%  filter(SNP %in% PC_K5_qval05_OL$SNP, region %in% "vgll3") #0
SNPs_highlight %>%  filter(SNP %in% PC_K5_qval05_OL$SNP, region %in% "six6") # 12
SNPs_highlight %>%  filter(SNP %in% PC_K5_qval05_OL$SNP, region %in% "SV1_23") #559

##get gene overlap

#make bed
PC_K5_qval05_OL_bed <- PC_K5_qval05_OL %>%  select(CHROM_nam, BP, qvalues) %>% 
  mutate(BP_jr = BP + 1) %>%  select(CHROM_nam, BP, BP_jr, qvalues)

fwrite(PC_K5_qval05_OL_bed , "PC_K5_qval05_OL.bed",
       col.names = F, row.names = F, sep = "\t", quote = F )

system("bedtools intersect -a SSA_genes.bed -b PC_K5_qval05_OL.bed -wb >  PC_K5_qval05_geneoverlap")

PC_K5_qval05_geneoverlap <- fread("PC_K5_qval05_geneoverlap") %>% 
  mutate (Chromosome = V1, Position = V2, gene = V4,  qvalues = V8) %>%  
  select (Chromosome, Position, gene ,  qvalues)

fwrite(PC_K5_qval05_geneoverlap, "Supplementary_Table_1_Axiom_pcadapt_genes.txt", sep = "\t", col.names = F, row.names = F, quote = F)
PC_K5_qval05_geneoverlap_unique <- PC_K5_qval05_geneoverlap %>%  arrange(gene) %>% 
  distinct(gene)

##GWAS LFMM 

##RIDGE REGRESSION LFMM
lfmmfunctionmap <- function(pheno, genomat, Knum, mapobj) {
## ## Fit an LFMM, i.e, compute B, U, V estimates:
mod.lfmm <- lfmm_ridge(Y = genomat, 
                       X = pheno, 
                       K = Knum)
## performs association testing using the fitted model:
pv <- lfmm_test(Y = genomat, 
                X = pheno, 
                lfmm = mod.lfmm, 
                calibrate = "gif") 

#get pvalues and qvalues
qvals <- qvalue(pv$calibrated.pvalue)
qvalues = qvals$qvalues[,1]
pvalues = pv$calibrated.pvalue[,1]
#join to map
mapped_lfmm <- bind_cols(mapobj, pvalues = pvalues, qvalues = qvalues)
return(mapped_lfmm)} # a function to fit and run a ridge regression model, FDR correction, matching to SNP map 

#across Ks
#K1
lfmm_grilse_K1 <- lfmmfunctionmap(pheno = grilse.meta.fam$p1SW, genomat = genos.dose.imp, Knum = 1, mapobj = chrom_snp_map)

#OL and downsamp for plot
lfmm_grilse_K1_OL <- lfmm_grilse_K1 %>%  filter(qvalues < 0.05) 
lfmm_grilse_K1_downsamp <- lfmm_grilse_K1  %>% sample_n(10000) %>%  bind_rows(., lfmm_grilse_K1_OL) %>%  distinct()

## qqplot
qqplot(rexp(length(lfmm_grilse_K1$pvalues), rate = log(10)),
       -log10(lfmm_grilse_K1$pvalues), xlab = "Expected quantile",
       pch = 19, cex = .4)
abline(0,1)


# gene overlap
lfmm_grilse_K1_OL.bed <- lfmm_grilse_K1_OL %>%  mutate(BP_jr = BP + 1) %>% 
  select(CHROM_nam, BP, BP_jr, qvalues) 

write.table(lfmm_grilse_K1_OL.bed , "lfmm_grilse_K1_OL.bed", col.names = F, row.names = F, sep = "\t", quote = F)
system("bedtools intersect -a SSA_genes.bed -b lfmm_grilse_K1_OL.bed -wb > lfmm_grilse_K1_OL_geneoverlap.tsv")

lfmm_grilse_K1_OL_geneoverlap <- fread("lfmm_grilse_K1_OL_geneoverlap.tsv") %>%   
  mutate (Chromosome = V1, Position = V2, gene = V4,  qvalues = V8) %>%  
  select (Chromosome, Position, gene ,  qvalues) %>% 
  mutate(Method = "LFMM K1")

#K5
lfmm_grilse_K5 <- lfmmfunctionmap(pheno = grilse.meta.fam$p1SW, genomat = genos.dose.imp, Knum = 5, mapobj = chrom_snp_map)

#OL and downsamp for plot
lfmm_grilse_K5_OL <-  lfmm_grilse_K5 %>%  filter(qvalues < 0.05) 
lfmm_grilse_K5_downsamp <- lfmm_grilse_K5  %>% sample_n(10000) %>%  bind_rows(., lfmm_grilse_K5_OL) %>%  distinct()


##qqplot
qqplot(rexp(length(lfmm_grilse_K5$pvalues), rate = log(10)),
       -log10(lfmm_grilse_K5$pvalues), xlab = "Expected quantile",
       pch = 19, cex = .4)
abline(0,1)

# gene overlap
lfmm_grilse_K5_OL.bed <- lfmm_grilse_K5_OL %>%  mutate(BP_jr = BP + 1) %>% 
  select(CHROM_nam, BP, BP_jr, qvalues) 

write.table(lfmm_grilse_K5_OL.bed , "lfmm_grilse_K5_OL.bed", col.names = F, row.names = F, sep = "\t", quote = F)
system("bedtools intersect -a SSA_genes.bed -b lfmm_grilse_K5_OL.bed -wb > lfmm_grilse_K5_OL_geneoverlap.tsv")

lfmm_grilse_K5_OL_geneoverlap <- fread("lfmm_grilse_K5_OL_geneoverlap.tsv") %>%   
  mutate (Chromosome = V1, Position = V2, gene = V4,  qvalues = V8) %>%  
  select (Chromosome, Position, gene ,  qvalues) %>% 
  mutate(Method = "LFMM K5")

lfmm_grilse_OL_geneoverlap <- bind_rows(lfmm_grilse_K1_OL_geneoverlap, lfmm_grilse_K5_OL_geneoverlap) %>%   arrange(qvalues)
write.table(lfmm_grilse_OL_geneoverlap, "Supplementary_Table_3_LFMMGrilse_OL_genes.tsv", col.names = T, row.names = F, sep = "\t", quote = F)


#check overlaps with six6, vgll3, 1-23 SV, between sets
SNPs_highlight %>% filter(SNP %in% lfmm_grilse_K1_OL$SNP)
SNPs_highlight %>% filter(SNP %in% lfmm_grilse_K5_OL$SNP)
lfmm_grilse_K1_OL %>%  filter(SNP %in% lfmm_grilse_K5_OL$SNP)
inner_join(lfmm_grilse_K1_OL_geneoverlap, lfmm_grilse_K5_OL_geneoverlap, by = gene)

lfmm_grilse_K1_OL_geneoverlap_unique <- lfmm_grilse_K1_OL_geneoverlap  %>%  arrange(gene) %>% 
  distinct(gene)

lfmm_grilse_K5_OL_geneoverlap_unique <- lfmm_grilse_K5_OL_geneoverlap  %>%  arrange(gene) %>% 
  distinct(gene)


##plot
grilseplot_lfmm_K1 <- ggman(lfmm_grilse_K1_downsamp, 
                                chrom = "CHROM", pvalue = "qvalues", snp = "SNP", bp="BP", 
                                pointSize = 1, title = "K1",
                                xlabel = "Chromosome", ylabel = "-log10(qvalue)", 
                                logTransform = T,  sigLine = 1.3) + theme_classic()
ggmanHighlightGroup(grilseplot_lfmm_K1 , highlightDfm = SNPs_highlight, snp = "SNP", group = "region")


grilseplot_lfmm_K5 <- ggman(lfmm_grilse_K5_downsamp, 
                            chrom = "CHROM", pvalue = "qvalues", snp = "SNP", bp="BP", 
                            pointSize = 1, title = "K5",
                            xlabel = "Chromosome", ylabel = "-log10(qvalue)", 
                            logTransform = T,  sigLine = 1.3) + theme_classic()
ggmanHighlightGroup(grilseplot_lfmm_K5 , highlightDfm = SNPs_highlight, snp = "SNP", group = "region")


###RDA scores 

RDA_scores_map  <- function(RDAobj, mapobj, Knum) {
  RDAscores <- data.frame(abs(RDAobj$CCA$v[,1:Knum])) #get absolute values of RDA scores for easy plotting
  colnames(RDAscores) <- paste0("absRDA", rep(1:Knum))
  RDAobj_map <- bind_cols(chrom_snp_map, RDAscores)
  return(RDAobj_map)}


grilse.rda.SNPscores <- RDA_scores_map(RDAobj = grilse.rda, mapobj = chrom_snp_map, Knum = 1)

#RDA outliers - top 1%
grilse.rda.SNPscores_OL99 <-grilse.rda.SNPscores %>%  
  filter(absRDA1 > quantile(absRDA1, 0.99))

#gene overlap
grilse.rda.SNPscores_OL99.bed   <- grilse.rda.SNPscores_OL99    %>%  select(CHROM_nam, BP, absRDA1) %>% 
  mutate(BP_jr = BP + 1) %>%  select(CHROM_nam, BP, BP_jr, absRDA1)

fwrite(grilse.rda.SNPscores_OL99.bed, "grilse.rda.SNPscores_OL99.bed",
       col.names = F, row.names = F, sep = "\t", quote = F )

system("bedtools intersect -a SSA_genes.bed -b grilse.rda.SNPscores_OL99.bed -wb > grilse.rda.SNPscores_OL99_geneoverlap.tsv")

grilse.rda.SNPscores_OL99_geneoverlap <- fread("grilse.rda.SNPscores_OL99_geneoverlap.tsv", header = F) %>% select(V1, V2, V4, V8)
colnames(grilse.rda.SNPscores_OL99_geneoverlap) <- c("Chromosome", "Position", "gene" ,  "absRDA1")
write.table(grilse.rda.SNPscores_OL99_geneoverlap, "Supplementary_Table_4_RDA_geneoverlap.tsv", col.names = T, row.names = F, quote = F, sep = "\t")

grilse.rda.SNPscores_OL99_geneoverlap_unique <- grilse.rda.SNPscores_OL99_geneoverlap %>%  arrange(desc(absRDA1)) %>%  distinct(gene)


# gene overlaps between analyses
grilse.rda.SNPscores_OL99_geneoverlap_unique %>%  filter(gene %in% PC_K10_qval05_geneoverlap_unique$gene)
grilse.rda.SNPscores_OL99_geneoverlap_unique %>%  filter(gene %in% lfmm_grilse_K1_OL_geneoverlap_unique$gene)
grilse.rda.SNPscores_OL99_geneoverlap_unique %>%  filter(gene %in% lfmm_grilse_K5_OL_geneoverlap_unique$gene)
lfmm_grilse_K1_OL_geneoverlap_unique  %>%  filter(gene %in% PC_K10_qval05_geneoverlap_unique$gene)
lfmm_grilse_K5_OL_geneoverlap_unique  %>%  filter(gene %in% PC_K10_qval05_geneoverlap_unique$gene)


# SNP overlaps between analyses
grilse.rda.SNPscores_OL99 %>%  filter(SNP %in% PC_K10_qval05_OL$SNP)
grilse.rda.SNPscores_OL99 %>%  filter(SNP %in% lfmm_grilse_K1_OL$SNP)
grilse.rda.SNPscores_OL99 %>%  filter(SNP %in% lfmm_grilse_K5_OL$SNP)
lfmm_grilse_K1_OL  %>%  filter(SNP %in% PC_K10_qval05_OL$SNP)
lfmm_grilse_K5_OL  %>%  filter(SNP %in% PC_K10_qval05_OL$SNP)
lfmm_grilse_K5_OL  %>%  filter(SNP %in% lfmm_grilse_K1_OL$SNP)

##polysel prep for gene set enrichment

#prep for polysel - SET info 
SetInfo <- read.csv("ssalar_kegg.csv")
colnames(SetInfo) <- c("setID", "setName", "setSource")

#SetObj
#turns into SetObj after filtering by GeneIDs in Salmon dataset
SetObj <- read.delim("biosystems_gene", header=FALSE)
colnames(SetObj) <- c("setID",	"objID")
SetObj <- SetObj[1:2]
SetObj <- SetObj[SetObj$setID %in% SetInfo$setID,]

#OBJ info
ssalar_genes <- fread("ssalar_gene_result.txt")
get.topgene.perstat.objinfo <- function(genescoreobj, statvar){
  topvar_genes <- genescoreobj %>% 
  group_by(Symbol) %>%
  filter({{statvar}} == max({{statvar}}))
  topvar_genes.uniq <- as.data.frame(unique(setDT(topvar_genes), by = "Symbol"))
  #get gene info
  
  Obj_Info <- inner_join(topvar_genes.uniq, ssalar_genes)
  Obj_Info$SNPcount <- 1
  Obj_Info$GeneLength <- Obj_Info$end_position_on_the_genomic_accession - Obj_Info$start_position_on_the_genomic_accession
  #Make an object info file with charr positions and FST, but used humann gene IDs
  Obj_Info <- Obj_Info %>%  select(GeneID, {{statvar}}, Symbol ,GeneLength, Chromosome, start_position_on_the_genomic_accession, 
                                   end_position_on_the_genomic_accession, orientation)
  colnames(Obj_Info) <- c("objID",	"objStat",	"objName",	"GeneLength",	"chr",	"startpos",	"endpos",	"strand")
  Obj_Info$strand <- NA
  return(Obj_Info)
} # a function to make an obj info file

#PC5
PC_K5.bed <- PC_K5_scores %>%  mutate(BP_jr = BP + 1) %>% 
  select(CHROM_nam, BP, BP_jr, PC1_loading)

fwrite(PC_K5.bed , "PC_K5.bed",
       col.names = F, row.names = F, sep = "\t", quote = F )

system("bedtools intersect -a SSA_genes.bed -b PC_K5.bed -wb >  PC_K5_geneoverlap")

PC_K5_geneoverlap <- fread("PC_K5_geneoverlap") %>% 
  mutate (Chromosome = V1, Position = V2, Symbol = V4, PC1 = V8) %>%  ##we use symbol instead of gene here for consistency with ssalar_gene_result.txt
  select (Chromosome, Position, Symbol ,  PC1)

#get ObjInfo
PC1_ObjInfo <- get.topgene.perstat.objinfo(genescoreobj = PC_K5_geneoverlap, statvar = PC1)

#ensure consistency of genes, pathways across comparisons
PC1_ObjInfo <- PC1_ObjInfo %>%  filter(objID  %in% SetObj$objID)
PC1_SetObj <- SetObj %>% filter(objID  %in% PC1_ObjInfo$objID) 
PC1_SetInfo <- SetInfo[SetInfo$setID %in% SetObj$setID,]

make.polysel.proj <- function(ObjInfo, SetObj, SetInfo, projname){
  mkdir_cmd <- paste0("mkdir ~/Desktop/Software/polysel/data/", projname,  "/")
  system(mkdir_cmd)
  write.table(ObjInfo, paste0("~/Desktop/Software/polysel/data/", projname, "/ObjInfo.txt"), col.names = T, row.names = F, sep = "\t", quote = F) 
  write.table(SetInfo, paste0("~/Desktop/Software/polysel/data/", projname, "/SetInfo.txt"), col.names = T, row.names = F, sep = "\t", quote = F)
  write.table(SetObj, paste0("~/Desktop/Software/polysel/data/", projname, "/SetObj.txt"), col.names = T, row.names = F, sep = "\t", quote = F)
}


make.polysel.proj(ObjInfo = PC1_ObjInfo,
                  SetObj = PC1_SetObj,
                  SetInfo = PC1_SetInfo,
                  projname = "grilse_PC1")


#Lfmm_K1
lfmm_grilse_K1.bed <- lfmm_grilse_K1 %>%  mutate(BP_jr = BP + 1) %>% 
  mutate(log10pvals = -log10(pvalues)) %>% 
  select(CHROM_nam, BP, BP_jr, log10pvals)

fwrite(lfmm_grilse_K1.bed, "lfmm_grilse_K1.bed",
       col.names = F, row.names = F, sep = "\t", quote = F )

system("bedtools intersect -a SSA_genes.bed -b lfmm_grilse_K1.bed -wb >  lfmm_grilse_K1_geneoverlap")

lfmm_grilse_K1_geneoverlap <- fread("lfmm_grilse_K1_geneoverlap") %>% 
  mutate (Chromosome = V1, Position = V2, Symbol = V4, log10pvals = V8) %>%  ##we use symbol instead of gene here for consistency with ssalar_gene_result.txt
  select (Chromosome, Position, Symbol ,  log10pvals)

#get ObjInfo
lfmm_K1_ObjInfo <- get.topgene.perstat.objinfo(genescoreobj = lfmm_grilse_K1_geneoverlap, statvar = log10pvals)

#ensure consistency of genes, pathways across comparisons
lfmm_K1_ObjInfo <- lfmm_K1_ObjInfo %>%  filter(objID  %in% SetObj$objID)
lfmm_K1_SetObj <- SetObj %>% filter(objID  %in% lfmm_K1_ObjInfo$objID) 
lfmm_K1_SetInfo <- SetInfo[SetInfo$setID %in% SetObj$setID,]

make.polysel.proj(ObjInfo = lfmm_K1_ObjInfo,
                  SetObj = lfmm_K1_SetObj,
                  SetInfo = lfmm_K1_SetInfo,
                  projname = "grilse_lfmm_K1")
##lfmm K5

lfmm_grilse_K5.bed <- lfmm_grilse_K5 %>%  mutate(BP_jr = BP + 1) %>% 
  mutate(log10pvals = -log10(pvalues)) %>% 
  select(CHROM_nam, BP, BP_jr, log10pvals)

fwrite(lfmm_grilse_K5.bed, "lfmm_grilse_K5.bed",
       col.names = F, row.names = F, sep = "\t", quote = F )

system("bedtools intersect -a SSA_genes.bed -b lfmm_grilse_K5.bed -wb >  lfmm_grilse_K5_geneoverlap")

lfmm_grilse_K5_geneoverlap <- fread("lfmm_grilse_K5_geneoverlap") %>% 
  mutate (Chromosome = V1, Position = V2, Symbol = V4, log10pvals = V8) %>%  ##we use symbol instead of gene here for consistency with ssalar_gene_result.txt
  select (Chromosome, Position, Symbol ,  log10pvals)

#get ObjInfo
lfmm_K5_ObjInfo <- get.topgene.perstat.objinfo(genescoreobj = lfmm_grilse_K5_geneoverlap, statvar = log10pvals)

#ensure consistency of genes, pathways across comparisons
lfmm_K5_ObjInfo <- lfmm_K5_ObjInfo %>%  filter(objID  %in% SetObj$objID)
lfmm_K5_SetObj <- SetObj %>% filter(objID  %in% lfmm_K5_ObjInfo$objID) 
lfmm_K5_SetInfo <- SetInfo[SetInfo$setID %in% SetObj$setID,]

make.polysel.proj(ObjInfo = lfmm_K5_ObjInfo,
                  SetObj = lfmm_K5_SetObj,
                  SetInfo = lfmm_K5_SetInfo,
                  projname = "grilse_lfmm_K5")


#go to analysis.3.Salmon_Seaage_polysel.R  
