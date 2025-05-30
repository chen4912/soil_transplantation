getwd()
setwd("../soil_transplantation/data/")

##############################PACKAGES##############################
library(dplyr)
library(vegan)
library(tidyr)
library(openxlsx)
library(adespatial)
library(phyloseq)
library(picante)

set.seed(1234)
####The function of finding the mean and standard deviation of different treatments####
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE),
      se = sd(x[[col]]/sqrt(length(x[[col]]))))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
} 

####code of persistent and recruited communities#####
OTU_group = read.xlsx("field.xlsx", sheet = "group", rowNames = T, colNames = T)
OTU_group$Sample_ID = rownames(OTU_group)
head(OTU_group)
#### 
OTU_group$Site_code = ifelse(OTU_group$Site == "High", "H", ifelse(OTU_group$Site == "Mid", "M", "L"))
OTU_group$Popu_code = ifelse(OTU_group$Seedling_source == "High", "H", ifelse(OTU_group$Seedling_source == "Mid", "M", "L"))
OTU_group$Soil_code = ifelse(OTU_group$Soil == "High", "H", ifelse(OTU_group$Soil == "Mid", "M", "L"))
head(OTU_group)

#### type code
OTU_group$type = paste(OTU_group$Soil_code, OTU_group$Site_code, sep = "_")

#### rhizosphere
#bacteria_OTU_abun = read.xlsx("field.xlsx", sheet = "rhizos", rowNames = T, colNames = T)

bacteria_OTU_abun = read.xlsx("field.xlsx", sheet = "fungi", rowNames = T, colNames = T)###OTU data, where rows are OTUs, columns are samples

bacteria_OTU_abun[1:6,1:6]
#### Obtain persistent and non-persistent OTU tables: OTU classes present in all samples before and after transplantation
Species = unique(OTU_group$Species); Species
#i = "P. asiatica"

####bacteria_OTU_abun
group_otu_abun = bacteria_OTU_abun

###
final_results_list = list()

for (i in Species) {
  select_data = subset(OTU_group, Species == i & Time == "Aug")
  #### L_L_L(home)
  L_L_L_home = select_data[select_data$type %in% "L_L", ]$Sample_ID
  sample_low_OTU = group_otu_abun[, L_L_L_home]
  sample_low_OTU <- sample_low_OTU[rowSums(sample_low_OTU)>0,]
  #### M_M_M(home)
  M_M_M_home = select_data[select_data$type %in% "M_M", ]$Sample_ID
  sample_mid_OTU = group_otu_abun[, M_M_M_home]
  sample_mid_OTU <- sample_mid_OTU[rowSums(sample_mid_OTU)>0,]
  #### H_H_H(home)
  H_H_H_home = select_data[select_data$type %in% "H_H", ]$Sample_ID
  sample_high_OTU = group_otu_abun[, H_H_H_home]
  sample_high_OTU <- sample_high_OTU[rowSums(sample_high_OTU)>0,]
  ####
  #### 
  #### away sample( M -> H)
  MtoH_away = c((subset(select_data, type == "M_H"))$Sample_ID)
  MtoH_away_OTU = group_otu_abun[, MtoH_away]
  MtoH_away_OTU <- MtoH_away_OTU[rowSums(MtoH_away_OTU)>0,]
  #### away sample( L -> H)
  LtoH_away = c((subset(select_data, type == "L_H"))$Sample_ID)
  LtoH_away_OTU = group_otu_abun[, LtoH_away]
  LtoH_away_OTU <- LtoH_away_OTU[rowSums(LtoH_away_OTU)>0,]
  #### away sample( L -> M)
  LtoM_away = c((subset(select_data, type == "L_M"))$Sample_ID)
  LtoM_away_OTU = group_otu_abun[, LtoM_away]
  LtoM_away_OTU <- LtoM_away_OTU[rowSums(LtoM_away_OTU)>0,]
  #### 
  #### away sample( H -> M)
  HtoM_away = c((subset(select_data, type == "H_M"))$Sample_ID)
  HtoM_away_OTU = group_otu_abun[, HtoM_away]
  HtoM_away_OTU <- HtoM_away_OTU[rowSums(HtoM_away_OTU)>0,]
  #### away sample( H -> L)
  HtoL_away = c((subset(select_data, type == "H_L"))$Sample_ID)
  HtoL_away_OTU = group_otu_abun[, HtoL_away]
  HtoL_away_OTU <- HtoL_away_OTU[rowSums(HtoL_away_OTU)>0,]
  #### away sample( M -> L)
  MtoL_away = c((subset(select_data, type == "M_L"))$Sample_ID)
  MtoL_away_OTU = group_otu_abun[, MtoL_away]
  MtoL_away_OTU <- MtoL_away_OTU[rowSums(MtoL_away_OTU)>0,]
  #### 
  ### low
  sample_low = unique(c(L_L_L_home, LtoM_away, LtoH_away))
  sample_low_OTU = group_otu_abun[, sample_low]
  sample_low_nohome = unique(c(LtoM_away, LtoH_away))
  sample_low_nohome_OTU = group_otu_abun[, sample_low_nohome]
  
  # Calculate the number of samples for each group
  half_L_L_L_home <- length(L_L_L_home)
  half_LtoM_away <- length(LtoM_away)
  half_LtoH_away <- length(LtoH_away)
  # Keep the OTU that appears in all samples in L_L_L_home
  L_L_L_home_OTU1 <- rownames(sample_low_OTU)[rowSums(sample_low_OTU[, L_L_L_home] > 0) >= half_L_L_L_home]
  # Keep the OTU that appears in all samples in L to M_away
  LtoM_away_OTU1 <- rownames(group_otu_abun[, LtoM_away])[rowSums(group_otu_abun[, LtoM_away] > 0) >= half_LtoM_away]
  # Keep the OTU that appears in all samples in L to H_away
  LtoH_away_OTU1 <- rownames(group_otu_abun[, LtoH_away])[rowSums(group_otu_abun[, LtoH_away] > 0) >= half_LtoH_away]
  # Find the OTU shared by the three groups
  common_OTUs <- Reduce(intersect, list(L_L_L_home_OTU1, LtoM_away_OTU1, LtoH_away_OTU1))
  #common_OTUs_LM <- Reduce(intersect, list(L_L_L_home_OTU1, LtoM_away_OTU1))
  #common_OTUs_LH <- Reduce(intersect, list(L_L_L_home_OTU1, LtoH_away_OTU1))
  LtoM_home_OTU <- LtoM_away_OTU[common_OTUs, ]
  LtoH_home_OTU <- LtoH_away_OTU[common_OTUs, ]
  df1 <- tibble::rownames_to_column(LtoM_home_OTU, var = "ID")
  df2 <- tibble::rownames_to_column(LtoH_home_OTU, var = "ID")
  low_sample_home <- full_join(df1, df2, by = "ID")
  low_sample_home[is.na(low_sample_home)] <- 0
  
  ###Newly recruited OTU abundance matrix
  LtoM_away_OTU <- LtoM_away_OTU[setdiff(rownames(LtoM_away_OTU), rownames(LtoM_home_OTU)),]
  LtoH_away_OTU <- LtoH_away_OTU[setdiff(rownames(LtoH_away_OTU), rownames(LtoH_home_OTU)),]
  df3 <- tibble::rownames_to_column(LtoM_away_OTU, var = "ID")
  df4 <- tibble::rownames_to_column(LtoH_away_OTU, var = "ID")
  low_sample_away <- full_join(df3, df4, by = "ID")
  low_sample_away[is.na(low_sample_away)] <- 0
  
  ###########Same operation for medium altitude and high altitude operation
  ###middle
  sample_mid = unique(c(M_M_M_home, MtoH_away, MtoL_away))
  sample_mid_OTU = group_otu_abun[, sample_mid]
  sample_mid_nohome = unique(c(MtoH_away, MtoL_away))
  sample_mid_nohome_OTU = group_otu_abun[, sample_mid_nohome]
  
  # Calculate the number of samples for each group
  half_M_M_M_home <- length(M_M_M_home)
  half_MtoH_away <- length(MtoH_away)
  half_MtoL_away <- length(MtoL_away)
  # Keep the OTU that appears in all samples in M_M_M_home
  M_M_M_home_OTU1 <- rownames(sample_mid_OTU)[rowSums(sample_mid_OTU[, M_M_M_home] > 0) >= half_M_M_M_home]
  # Keep the OTU that appears in all samples in M to H_away 
  MtoH_away_OTU1 <- rownames(group_otu_abun[, MtoH_away])[rowSums(group_otu_abun[, MtoH_away] > 0) >= half_MtoH_away]
  # Keep the OTU that appears in all samples in M to L_away
  MtoL_away_OTU1 <- rownames(group_otu_abun[, MtoL_away])[rowSums(group_otu_abun[, MtoL_away] > 0) >= half_MtoL_away]
  # Find the OTU shared by the three groups
  common_OTUs <- Reduce(intersect, list(M_M_M_home_OTU1, MtoH_away_OTU1, MtoL_away_OTU1))
  #common_OTUs_MH <- Reduce(intersect, list(M_M_M_home_OTU1, MtoH_away_OTU1))
  #common_OTUs_ML <- Reduce(intersect, list(M_M_M_home_OTU1, MtoL_away_OTU1))
  MtoH_home_OTU <- MtoH_away_OTU[common_OTUs, ]
  MtoL_home_OTU <- MtoL_away_OTU[common_OTUs, ]
  df1 <- tibble::rownames_to_column(MtoH_home_OTU, var = "ID")
  df2 <- tibble::rownames_to_column(MtoL_home_OTU, var = "ID")
  mid_sample_home <- full_join(df1, df2, by = "ID")
  mid_sample_home[is.na(mid_sample_home)] <- 0
  
  ###Newly recruited OTU abundance matrix
  MtoH_away_OTU <- MtoH_away_OTU[setdiff(rownames(MtoH_away_OTU), rownames(MtoH_home_OTU)),]
  MtoL_away_OTU <- MtoL_away_OTU[setdiff(rownames(MtoL_away_OTU), rownames(MtoL_home_OTU)),]
  df3 <- tibble::rownames_to_column(MtoH_away_OTU, var = "ID")
  df4 <- tibble::rownames_to_column(MtoL_away_OTU, var = "ID")
  mid_sample_away <- full_join(df3, df4, by = "ID")
  mid_sample_away[is.na(mid_sample_away)] <- 0
  
  ### high
  sample_high = unique(c(H_H_H_home, HtoL_away, HtoM_away))
  sample_high_OTU = group_otu_abun[, sample_high]
  sample_high_nohome = unique(c(HtoL_away, HtoM_away))
  sample_high_nohome_OTU = group_otu_abun[, sample_high_nohome]
  
  # Calculate the number of samples for each group
  half_H_H_H_home <- length(M_M_M_home)
  half_HtoL_away <- length(HtoL_away)
  half_HtoM_away <- length(HtoM_away)
  # Keep the OTU that appears in all samples in H_H_H_home 
  H_H_H_home_OTU1 <- rownames(sample_high_OTU)[rowSums(sample_high_OTU[, H_H_H_home] > 0) >= half_H_H_H_home]
  #Keep the OTU that appears in all samples in H to L_away 
  HtoL_away_OTU1 <- rownames(group_otu_abun[, HtoL_away])[rowSums(group_otu_abun[, HtoL_away] > 0) >= half_HtoL_away]
  #Keep the OTU that appears in all samples in H to M_away
  HtoM_away_OTU1 <- rownames(group_otu_abun[, HtoM_away])[rowSums(group_otu_abun[, HtoM_away] > 0) >= half_HtoM_away]
  # Find the OTU shared by the three groups
  common_OTUs <- Reduce(intersect, list(H_H_H_home_OTU1, HtoL_away_OTU1, HtoM_away_OTU1))
  #common_OTUs_HM <- Reduce(intersect, list(H_H_H_home_OTU1, HtoM_away_OTU1))
  #common_OTUs_HL <- Reduce(intersect, list(H_H_H_home_OTU1, HtoL_away_OTU1))
  HtoM_home_OTU <- HtoM_away_OTU[common_OTUs, ]
  HtoL_home_OTU <- HtoL_away_OTU[common_OTUs, ]
  df1 <- tibble::rownames_to_column(HtoM_home_OTU, var = "ID")
  df2 <- tibble::rownames_to_column(HtoL_home_OTU, var = "ID")
  high_sample_home <- full_join(df1, df2, by = "ID")
  high_sample_home[is.na(high_sample_home)] <- 0
  
  ###Newly recruited OTU abundance matrix
  HtoM_away_OTU <- HtoM_away_OTU[setdiff(rownames(HtoM_away_OTU), rownames(HtoM_home_OTU)),]
  HtoL_away_OTU <- HtoL_away_OTU[setdiff(rownames(HtoL_away_OTU), rownames(HtoL_home_OTU)),]
  df3 <- tibble::rownames_to_column(HtoM_away_OTU, var = "ID")
  df4 <- tibble::rownames_to_column(HtoL_away_OTU, var = "ID")
  high_sample_away <- full_join(df3, df4, by = "ID")
  high_sample_away[is.na(high_sample_away)] <- 0
  
  home_OTU_low_mid <- full_join(low_sample_home, mid_sample_home, by = "ID")
  home_OTU <- full_join(home_OTU_low_mid, high_sample_home, by = "ID")
  home_OTU[is.na(home_OTU)] <- 0
  
  away_OTU_low_mid <- full_join(low_sample_away, mid_sample_away, by = "ID")
  away_OTU <- full_join(away_OTU_low_mid, high_sample_away, by = "ID")
  away_OTU[is.na(away_OTU)] <- 0
  
  final_results_list[[paste("home_OTU", i, sep = "_")]] <- home_OTU
  final_results_list[[paste("away_OTU", i, sep = "_")]] <- away_OTU
}
final_results_list
observed_species_fnil <- colSums(final_results_list$`home_OTU_F. nilgerrensis` > 0)
observed_species_pasi <- colSums(final_results_list$`home_OTU_P. asiatica` > 0)
observed_species_trep <- colSums(final_results_list$`home_OTU_T. repens` > 0)

home_Fnil <- final_results_list$`home_OTU_F. nilgerrensis`
home_Pasi <- final_results_list$`home_OTU_P. asiatica`
home_Trep <- final_results_list$`home_OTU_T. repens`
nohome_Fnil <- final_results_list$`away_OTU_F. nilgerrensis`
nohome_Pasi <- final_results_list$`away_OTU_P. asiatica`
nohome_Trep <- final_results_list$`away_OTU_T. repens`
###merge persisted OTU
merged_1_2 <- full_join(home_Fnil, home_Pasi, by = "ID")
home_result <- full_join(merged_1_2, home_Trep, by = "ID")
home_result[is.na(home_result)] <- 0
#write.table(home_result, file ="home_result_fungi.csv",sep =",", quote =FALSE)
###merge non-persisted OTU
merged_1_2 <- full_join(nohome_Fnil, nohome_Pasi, by = "ID")
nohome_result <- full_join(merged_1_2, nohome_Trep, by = "ID")
nohome_result[is.na(nohome_result)] <- 0
#write.table (nohome_result, file ="nohome_result_fungi.csv",sep =",", quote =FALSE)

####Get persistent and newly introduced pcoa
####persisted
tree = read.tree("fungi.nwk")
#tree = read.tree("bacterial.nwk")
rownames(home_result) <- home_result[,1]
home_result <- home_result[,-1]
fungi_hell <- decostand(data.frame(t(home_result)),method = "hellinger")
fungi_tree1 <- prune.sample(fungi_hell,tree)
fungi_unifrac <- phyloseq(otu_table(t(fungi_hell),taxa_are_rows = T),phy_tree(fungi_tree1))#create phyloeq index
fungi_unifrac1 <- distance(fungi_unifrac,method = "wunifrac")

#####PCoA
pcoa <- cmdscale(fungi_unifrac1, k = 2, eig = TRUE)
pcoa <- as.data.frame(pcoa$points)
names(pcoa)[1:2] <- c("PCoA1", "PCoA2")
pcoa

##end
