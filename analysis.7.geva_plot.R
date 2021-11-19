setwd("/Users/ianbradbury/Desktop/Working/Sea_Age/geva")
#lsit files
XPnSL_files <- list.files(path = "./", pattern = "XPnsl.*.sites.txt")

# get chrom IDs
gene_IDs <- str_replace(XPnSL_files,  "XPnsl", "") %>% 
  str_replace("-", "") %>% 
  str_replace(".geva.sites.txt", "")

#read in files as list
All_XPnslgene_geva <- map(.x =XPnSL_files, .f = fread)

#change chrom names
names(All_XPnslgene_geva) <- gene_IDs
All_XPnslgene_geva$kmt2a
#start breaking down lists

All_XPnslgene_geva <- map_df(All_XPnslgene_geva, ~as.data.frame(.x), .id="gene") 

All_XPnslgene_geva <- All_XPnslgene_geva %>% filter(Filtered == 1, Clock %in% "J") %>% 
  mutate(MeanYears = 5 * PostMean) %>%  select(gene, MarkerID, MeanYears)

vgll3_geva <- All_XPnslgene_geva %>%  filter(gene %in% "vgll3")
six6_geva <- All_XPnslgene_geva %>%  filter(gene %in% "six6")

#lsit files
Random1K_files <- list.files(path = "./", pattern = "rand.*.sites.txt")

# get chrom IDs
chrom_IDs <- str_replace(Random1K_files,  "random1k", "") %>% 
  str_replace("-", "") %>% 
  str_replace(".geva.sites.txt", "")

#read in files as list
All_random_geva <- map(.x =Random1K_files , .f = fread)

#change chrom names
names(All_random_geva) <- chrom_IDs

#start breaking down lists

All_random_geva  <- map_df(All_random_geva , ~as.data.frame(.x), .id="chrom") 

All_random_geva <- All_random_geva %>% filter(Filtered == 1, Clock %in% "J") %>% 
  mutate(MeanYears = 5 * PostMean) %>%  select(chrom, MarkerID, MeanYears)


#plot em

ggplot() + geom_density(data = All_XPnslgene_geva, aes( x = MeanYears, fill = "XPnSL outlier", alpha = 0.5)) + geom_vline(xintercept = 25000) + geom_vline(xintercept = 100000) +
   geom_density(data =All_random_geva , aes( x = MeanYears, fill = "Random1k", alpha = 0.5)) 


ggplot() +  stat_ecdf(data =All_XPnslgene_geva, aes(x = MeanYears, colour = "XPnSL outlier"), geom = "step") +
  stat_ecdf(data =All_random_geva, aes(x = MeanYears, colour = "Random1K"), geom = "step")

mean(All_XPnslgene_geva$MeanYears)
mean(All_random_geva$MeanYears)
wilcox.test(All_XPnslgene_geva$MeanYears, All_random_geva$MeanYears)


