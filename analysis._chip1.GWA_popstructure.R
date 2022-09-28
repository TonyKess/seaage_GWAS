###seting up regional groups based on map position and PCA
library(data.table)
library(tidyverse)
library(pcadapt)

setwd("~/Desktop/Working/")

#functions
get_meta_merge_order <- function(filename, metafilename){
  FAM <- fread(paste0(filename, ".fam")) %>% 
    select(V1,V2) %>% 
    mutate(SiteCode = V1, ID = V2) %>%  
    select(-V1, -V2)
  meta <- fread(paste0(metafilename))
  meta_ordered <- inner_join(FAM, meta)
  return(meta_ordered)
} #get and merge metadata

salmo_220K_2022_meta <- get_meta_merge_order(filename = "Salmo_220K_220K_Merged2021_2022_Polyhigh_nodup_nohyb_maf01",
                     metafilename = "CIGENE_220K_Metadata_All_May_2022.tsv")

salmo_220K_site_group_meta_2022_NA <- salmo_220K_2022_meta %>%  
  select(SiteCode, Site, Group, Lat, Lon) %>%  
  filter(Group %in% "NorthAm") %>%  
  distinct()

write.table(salmo_220K_site_group_meta_2022_NA, "salmo_220K_site_group_meta_2022_NA.tsv", col.names = T, 
            row.names = F,
            sep = "\t", quote = F)

