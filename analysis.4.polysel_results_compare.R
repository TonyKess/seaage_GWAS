library(superheat)

SetInfo <- read.csv("ssalar_kegg.csv")

colnames(SetInfo)[1:3] <- c("setID", "setName", "setSource")

All_PC1_polysel <- fread("salmon_PC1_All_polysel.txt")
All_PC1_polysel_set <- inner_join( All_PC1_polysel, SetInfo, BY = "setID", "setName") %>% 
  mutate(AllPC1_q = setQ) %>% 
  select(setSource, AllPC1_q )

All_LRT_polysel <- fread("salmon_LRT_All_polysel.txt")
All_LRT_polysel_set <- inner_join( All_LRT_polysel, SetInfo, BY = "setID", "setName") %>% 
  mutate(All_LRT_q = setQ) %>% 
  select(setSource, All_LRT_q )
All_LRT_polysel_sig <- All_LRT_polysel_set$setSource[which(All_LRT_polysel_set$All_LRT_q < 0.05)]

MQN_LRT_polysel <- fread("salmon_MQN_LRT_polysel.txt")
MQN_LRT_polysel_set <- inner_join( MQN_LRT_polysel, SetInfo, BY = "setID", "setName") %>% 
  mutate(MQN_LRT_q = setQ) %>% 
  select(setSource, MQN_LRT_q )

MQN_LRT_polysel_sig <- MQN_LRT_polysel_set$setSource[which(MQN_LRT_polysel_set$MQN_LRT_q < 0.05)]

FQN_LRT_polysel <- fread("salmon_FQN_LRT_polysel.txt")
FQN_LRT_polysel_set <- inner_join( FQN_LRT_polysel, SetInfo, BY = "setID", "setName") %>% 
  mutate(FQN_LRT_q = setQ) %>% 
  select(setSource, FQN_LRT_q )

FQN_LRT_polysel_sig <- FQN_LRT_polysel_set$setSource[which(FQN_LRT_polysel_set$FQN_LRT_q < 0.05)]


FNL_LRT_polysel <- fread("salmon_FNL_LRT_polysel.txt")
FNL_LRT_polysel_set <- inner_join( FNL_LRT_polysel, SetInfo, BY = "setID", "setName") %>% 
  mutate(FNL_LRT_q = setQ) %>% 
  select(setSource, FNL_LRT_q )

FNL_LRT_polysel_sig <- FNL_LRT_polysel_set$setSource[which(FNL_LRT_polysel_set$FNL_LRT_q < 0.05)]

MNL_LRT_polysel <- fread("salmon_MNL_LRT_polysel.txt")
MNL_LRT_polysel_set <- inner_join( MNL_LRT_polysel, SetInfo, BY = "setID", "setName") %>% 
  mutate(MNL_LRT_q = setQ) %>% 
  select(setSource, MNL_LRT_q )

MNL_LRT_polysel_sig <- MNL_LRT_polysel_set$setSource[which(MNL_LRT_polysel_set$MNL_LRT_q < 0.05)]

All_sig <- unique(c(MNL_LRT_polysel_sig , FNL_LRT_polysel_sig, FQN_LRT_polysel_sig, MQN_LRT_polysel_sig, All_LRT_polysel_sig))

#many inner joins
Polysel_sig <- inner_join(All_PC1_polysel_set, All_LRT_polysel_set,) %>% 
  inner_join(., FNL_LRT_polysel_set) %>% 
  inner_join(., MQN_LRT_polysel_set) %>% 
  inner_join(., FQN_LRT_polysel_set) %>% 
  inner_join(.,  MNL_LRT_polysel_set)  %>% 
  filter(setSource %in% All_sig)

Polysel_sig[,2:7] <- -log10(Polysel_sig[,2:7])
library("wesanderson")
setscoremat <- Polysel_sig[,2:7]
colnames(setscoremat) <- colnames(Polysel_sig[,2:7]) 
rownames(setscoremat) <- Polysel_sig$setSource
superheat(X = setscoremat, pretty.order.rows = T, pretty.order.cols = T, row.dendrogram = T,col.dendrogram = T, scale = F, force.left.label = TRUE, heat.pal = wes_palette("Zissou1"))

