getwd()
setwd("../soil_transplantation/data/")

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

##############################PACKAGES##############################
library(openxlsx)
library(vegan)
library(magrittr)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(patchwork)
library(ggsci)
library(ggpmisc)
library(mgcv) 
library(reshape)
library(emmeans)
library(agricolae)
library(ggbreak)
library(multcomp)
set.seed(1234)
##############Tables S6############
group1 <- read.xlsx("field.xlsx", sheet = "group",colNames = T, rowNames = T)
group1 <- subset(group1, Time == "Aug" & Type != "CK")##Samples from August were kept and transplanted at elevation
group1$Site <- factor(group1$Site)
group1$Seedling_source <- factor(group1$Seedling_source)

species_list <- c("F. nilgerrensis", "P. asiatica", "T. repens")##Species
response_vars <- c("rich_phyllo", "rich_rhizos", "rich_fungi")

anova_results_list <- list()# Create an empty list to store the results

# Loop through each species and each response variable
for (species in species_list) {
  for (response_var in response_vars) {
    
    formula <- as.formula(paste("log(", response_var, ") ~ Site * Soil_type * Seedling_source"))
    model <- aov(formula, data = subset(group1, Species == species))
    
    anova_result <- car::Anova(model, type = 2)
    anova_result_df <- as.data.frame(anova_result)
    
    anova_result_df$`Pr(>F)` <- round(anova_result_df$`Pr(>F)`, 3)
    anova_result_df$`F value` <- round(anova_result_df$`F value`, 3)
    anova_result_df$`Sum Sq` <- round(anova_result_df$`Sum Sq`, 3)
    
    rownames(anova_result_df) <- c("Site", "Soil_type", "Seedling_source", 
                                   "Site x Soil_type", "Site x Seedling_source", 
                                   "Soil_type x Seedling_source", 
                                   "Site x Soil_type x Seedling_source", "Residuals")
    
    anova_results_list[[paste(species, response_var, sep = "&")]] <- anova_result_df
  }
}

anova_results_list

###############Tables S7############
###community composition
set.seed(12345678)
group1 <- read.xlsx("field.xlsx", sheet = "group", colNames = TRUE, rowNames = TRUE)
group1 <- subset(group1, Time == "Aug" & Type != "CK")# Subset group1 for the August time and exclude "CK" type

# Define the list of species and corresponding matrices
species_list <- c("F. nilgerrensis", "P. asiatica", "T. repens")
microbial_matrices <- list(
  phyllo = "phyllo_trans.rds",
  rhizos = "rhizos_trans.rds",
  fungi = "fungi_trans.rds"
)

# Create an empty list to store results
adonis_results_list <- list()

# Loop through each species and microbial matrix
for (species in species_list) {
  for (microbe in names(microbial_matrices)) {
    microbial_data <- readRDS(microbial_matrices[[microbe]])
    microbial_data <- as.matrix(microbial_data)
    
    group2 <- subset(group1, Species == species)
    
    # Filter the microbial data for the current species
    data <- microbial_data[rownames(microbial_data) %in% rownames(group2), ]
    microbial_filtered <- data[, colnames(data) %in% rownames(group2)]
    
    # Convert Site and Seedling_source to factors
    group2$Site <- factor(group2$Site)
    group2$Seedling_source <- factor(group2$Seedling_source)
    
    # Perform the adonis2 analysis
    model <- adonis2(microbial_filtered ~ Site * Soil_type * Seedling_source, data = group2, permutations = 999, by = "terms")
    
    result_df <- as.data.frame(model)
    result_df$SumOfSqs <- round(result_df$SumOfSqs, 3)
    result_df$R2 <- round(result_df$R2, 3)
    result_df$F <- round(result_df$F, 3)
    result_df$`Pr(>F)` <- round(result_df$`Pr(>F)`, 3)
    
    adonis_results_list[[paste(species, microbe, sep = "_")]] <- result_df
  }
}

# Print the results list
adonis_results_list

#######S × SS × ST in fungi richness of P. asiatica############
group1 <- read.xlsx("field.xlsx", sheet = "group",colNames = T, rowNames = T)

group1 <- subset(group1, Time == "Aug" & Type != "CK")
group <- subset(group1, Species == "P. asiatica")
group$rich_fungi <- log(group$rich_fungi)
group$Type <- paste0(group$Site,"_",group$Seedling_source)
unique(group$Type)
df2 <- data_summary(group, varname = "rich_fungi",groupnames = c("Seedling_source","Type","Soil","Soil_type"))
df2$Seedling_source <- factor(df2$Seedling_source, levels = c("Low","Mid","High"))
df2$Type <- factor(df2$Type, levels = c("Low_Mid","Low_High","Mid_Low","Mid_High","High_Low","High_Mid"))

df2$Soil_type <- factor(df2$Soil_type, levels = c("home","away"))
ggplot(data = df2) +
  geom_errorbar(mapping = aes(x = Type, y = rich_fungi, ymin = rich_fungi - se, ymax = rich_fungi + se, 
                              color = Seedling_source, group = interaction(Soil_type, Seedling_source)),
                position = position_dodge(width = 0.85), width = 0, size = 0.5, alpha = 1) +
  geom_bar(mapping = aes(x = Type, y = rich_fungi, fill = interaction(Soil_type, Seedling_source), group = interaction(Soil_type, Seedling_source), color = Seedling_source),
           position = position_dodge(width = 0.85), stat = "identity", size = 0.7, width = 0.7, alpha = 0.85) +
  scale_fill_manual(values = c("#2C91E0","NA", "#F0A73A","NA", "#3ABF99","NA")) +
  scale_color_manual(values = c("#2C91E0", "#F0A73A", "#3ABF99")) + 
  labs(y = "Effect size", title = NULL) +
  theme_bw() +
  geom_hline(aes(yintercept = 0), linetype = "dashed", size = 0.5) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "right",
        axis.line = element_line(linetype = 1, color = "black", size = 0.5),
        axis.ticks = element_line(color = "black", size = 0.5, lineend = "round"),
        axis.text = element_text(colour = 'black', size = 9),
        axis.title = element_text(colour = 'black', size = 10)) -> Fig1a;Fig1a

####Multiple comparisons
model <- aov(rich_fungi ~ Site * Soil_type * Seedling_source, data = group)
group$Site <- as.factor(group$Site)
group$Seedling_source <- as.factor(group$Seedling_source)
car::Anova(model, type =2)
model_means <- emmeans(model, specs = ~ Soil_type | Seedling_source * Site)
pairwise_comparisons <- pairs(model_means, adjust = "tukey")
model_means_cld <- cld(object = model_means,
                       Letters = letters,
                       alpha = 0.05,
                       decreasing = T)
model_means_cld

###################Fig 2a-c####################
###phyllosphere
phyllo <- readRDS("phyllo_trans_Jun.rds")
group1 <- read.xlsx("field.xlsx", sheet = "group",colNames = T, rowNames = T)
group1 <- subset(group1, Type != "CK")
phyllo <- as.matrix(phyllo)
data <- phyllo[rownames(phyllo) %in% rownames(group1),]
rhizos1 <- data[,colnames(data) %in% rownames(group1)]
pcoa <- cmdscale(rhizos1, k = 2, eig = TRUE)
poi <- pcoa$points
eigval <- pcoa$eig
pcoa_eig <- (pcoa$eig)[1:2]/sum(pcoa$eig)
pcoa_eig <- round(pcoa_eig, 3)
pcoa <- as.data.frame(poi)
names(pcoa)[1:2] <- c("PCoA1", "PCoA2")
pcoa$SampleID <- rownames(pcoa)
sample_site <- cbind(pcoa,group1[,c(1,2,3,4,5,6)])
sample_site$group1 <- paste0(sample_site$Site,"_",sample_site$Soil,"_",sample_site$Seedling_source)
sample_site_Jun <- subset(sample_site, Time == "Jun")
sample_site <- subset(sample_site, Time == "Aug")

###Organize the data set
df_trans <- aggregate(cbind(PCoA1, PCoA2) ~ group1+Species+Soil_type, sample_site, FUN = mean)
df2_trans <- df_trans %>%
  filter(sapply(strsplit(as.character(group1), "_"), function(x) x[1] == x[2]))
df3_trans <- df_trans %>%
  filter(!group1 %in% df2_trans$group1)
result_trans <- df3_trans %>%
  mutate(matched_value = NA, PCoA1_matched = NA, PCoA2_matched = NA)
for (i in 1:nrow(df3_trans)) {
  group1_values <- unlist(strsplit(as.character(df3_trans$group1[i]), "_"))
  matching_rows <- df2_trans %>%
    filter(sapply(strsplit(as.character(group1), "_"), function(x) {
      x[1] == group1_values[1] && x[3] == group1_values[3]
    }) & Species == df3_trans$Species[i])
  if (nrow(matching_rows) > 0) {
    result_trans$matched_value[i] <- matching_rows$group1[1] 
    result_trans$PCoA1_matched[i] <- matching_rows$PCoA1[1]
    result_trans$PCoA2_matched[i] <- matching_rows$PCoA2[1]
  }
}
df_trans1 <- df_trans %>%
  separate(group1, into = c("Site", "Soil", "Seedling_source"), sep = "_")
sample_site_Jun$Site <- factor(sample_site_Jun$Site)
sample_site$Site <- factor(sample_site$Site)
#df_situ1$group <- paste0(df_situ1$Site,"_",df_situ1$Soilsource1)
df_trans1$group <- paste0(df_trans1$Site,"_",df_trans1$Soil_type)
sample_site$group <- paste0(sample_site$Site,"_",sample_site$Soil_type)
ggplot() +
  #geom_segment(data = result_situ,mapping = aes(x = PCoA1, y = PCoA2, 
  #                          xend = PCoA1_matched, yend = PCoA2_matched),curvature = 0.3, size = 0.1, alpha = 0.7, color = "#FF8C00") +
  geom_segment(data = result_trans,mapping = aes(x = PCoA1_matched,  y = PCoA2_matched, 
                                                 xend = PCoA1, yend = PCoA2),curvature = 0.3, size = 0.1, alpha = 0.7, color = "#3D505A") +
  #geom_point(data = df_situ1,  mapping = aes(x = PCoA1, y = PCoA2, 
  #                        color = Species, shape = group, fill = Soilsource1),size = 2) +
  geom_point(data = df_trans1, mapping = aes(x = PCoA1, y = PCoA2, 
                                             color = Species, shape = group, fill = Soil_type),size = 2,alpha = 0.8) +
  geom_point(data = sample_site, mapping = aes(x = PCoA1, y = PCoA2, 
                                               color = Species, fill = Soil_type, shape = group),size = 1.2,alpha = 0.2) +
  geom_point(data = sample_site_Jun, mapping = aes(x = PCoA1, y = PCoA2, 
                                               color = Species,shape = Site),size = 2,stroke = 1.2,alpha = 0.8) +
  scale_color_manual(values = c("#974FA3", "#4CAB3B", "#B7352A")) +
  scale_shape_manual(values = c(Low_away=15,Low_home=22,Mid_away=16,Mid_home=21,High_away=17,High_home=24,Low=3,Mid=4,High=8)) +
  scale_fill_manual(values = c("white","white","white"))+
  labs(x=paste("PCoA1 (", format(100 * pcoa_eig[1], digits=3), "%)", sep=""),
       y=paste("PCoA2 (", format(100 * pcoa_eig[2], digits=3), "%)", sep=""))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(text=element_text(size=12))+
  theme(panel.background = element_rect(fill='white', colour='black'),
        panel.grid=element_blank(), 
        legend.position = "none",
        axis.text = element_text(color = "black"),
        legend.title=element_text(size = 12,face = "bold"),
        legend.text=element_text(size= 10))+
  geom_vline(aes(xintercept = 0),linetype="dashed",size = 0.5)+
  geom_hline(aes(yintercept = 0),linetype="dashed",size = 0.5)+
  scale_x_continuous(labels = scales::label_comma(accuracy =0.01))+
  scale_y_continuous(labels = scales::label_comma(accuracy =0.01)) -> Fig2a;Fig2a

####Rhizosphere
rhizos <- readRDS("rhizos_trans_Jun.rds")
group1 <- read.xlsx("field.xlsx", sheet = "group",colNames = T, rowNames = T)
group1 <- subset(group1, Type != "CK")
rhizos <- as.matrix(rhizos)
data <- rhizos[rownames(rhizos) %in% rownames(group1),]
rhizos1 <- data[,colnames(data) %in% rownames(group1)]
pcoa <- cmdscale(rhizos1, k = 2, eig = TRUE)
poi <- pcoa$points
eigval <- pcoa$eig
pcoa_eig <- (pcoa$eig)[1:2]/sum(pcoa$eig)
pcoa_eig <- round(pcoa_eig, 3)
pcoa <- as.data.frame(poi)
names(pcoa)[1:2] <- c("PCoA1", "PCoA2")
pcoa$SampleID <- rownames(pcoa)
sample_site <- cbind(pcoa,group1[,c(1,2,3,4,5,6)])
sample_site$group1 <- paste0(sample_site$Site,"_",sample_site$Soil,"_",sample_site$Seedling_source)
sample_site_Jun <- subset(sample_site, Time == "Jun")
sample_site <- subset(sample_site, Time == "Aug")

###Organize the data set
df_trans <- aggregate(cbind(PCoA1, PCoA2) ~ group1+Species+Soil_type, sample_site, FUN = mean)
df2_trans <- df_trans %>%
  filter(sapply(strsplit(as.character(group1), "_"), function(x) x[1] == x[2]))
df3_trans <- df_trans %>%
  filter(!group1 %in% df2_trans$group1)
result_trans <- df3_trans %>%
  mutate(matched_value = NA, PCoA1_matched = NA, PCoA2_matched = NA)
for (i in 1:nrow(df3_trans)) {
  group1_values <- unlist(strsplit(as.character(df3_trans$group1[i]), "_"))
  matching_rows <- df2_trans %>%
    filter(sapply(strsplit(as.character(group1), "_"), function(x) {
      x[1] == group1_values[1] && x[3] == group1_values[3]
    }) & Species == df3_trans$Species[i])
  if (nrow(matching_rows) > 0) {
    result_trans$matched_value[i] <- matching_rows$group1[1] 
    result_trans$PCoA1_matched[i] <- matching_rows$PCoA1[1]
    result_trans$PCoA2_matched[i] <- matching_rows$PCoA2[1]
  }
}
df_trans1 <- df_trans %>%
  separate(group1, into = c("Site", "Soil", "Seedling_source"), sep = "_")
sample_site_Jun$Site <- factor(sample_site_Jun$Site)
sample_site$Site <- factor(sample_site$Site)
#df_situ1$group <- paste0(df_situ1$Site,"_",df_situ1$Soilsource1)
df_trans1$group <- paste0(df_trans1$Site,"_",df_trans1$Soil_type)
sample_site$group <- paste0(sample_site$Site,"_",sample_site$Soil_type)
ggplot() +
  #geom_segment(data = result_situ,mapping = aes(x = PCoA1, y = PCoA2, 
  #                          xend = PCoA1_matched, yend = PCoA2_matched),curvature = 0.3, size = 0.1, alpha = 0.7, color = "#FF8C00") +
  geom_segment(data = result_trans,mapping = aes(x = PCoA1_matched,  y = PCoA2_matched, 
                                                 xend = PCoA1, yend = PCoA2),curvature = 0.3, size = 0.1, alpha = 0.7, color = "#3D505A") +
  #geom_point(data = df_situ1,  mapping = aes(x = PCoA1, y = PCoA2, 
  #                        color = Species, shape = group, fill = Soilsource1),size = 2) +
  geom_point(data = df_trans1, mapping = aes(x = PCoA1, y = PCoA2, 
                                             color = Species, shape = group, fill = Soil_type),size = 2,alpha = 0.8) +
  geom_point(data = sample_site, mapping = aes(x = PCoA1, y = PCoA2, 
                                               color = Species, fill = Soil_type, shape = group),size = 1.2,alpha = 0.2) +
  geom_point(data = sample_site_Jun, mapping = aes(x = PCoA1, y = PCoA2, 
                                                  color = Species,shape = Site),size = 2,stroke = 1.2,alpha = 0.8) +
  scale_color_manual(values = c("#974FA3", "#4CAB3B", "#B7352A")) +
  scale_shape_manual(values = c(Low_away=15,Low_home=22,Mid_away=16,Mid_home=21,High_away=17,High_home=24,Low=3,Mid=4,High=8)) +
  scale_fill_manual(values = c("white","white","white"))+
  labs(x=paste("PCoA1 (", format(100 * pcoa_eig[1], digits=3), "%)", sep=""),
       y=paste("PCoA2 (", format(100 * pcoa_eig[2], digits=3), "%)", sep=""))+
  #xlim(-0.06, 0.04)+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(text=element_text(size=12))+
  theme(panel.background = element_rect(fill='white', colour='black'),
        panel.grid=element_blank(), 
        legend.position = "none",
        axis.text = element_text(color = "black"),
        legend.title=element_text(size = 12,face = "bold"),
        legend.text=element_text(size= 10))+
  geom_vline(aes(xintercept = 0),linetype="dashed",size = 0.5)+
  geom_hline(aes(yintercept = 0),linetype="dashed",size = 0.5)+ 
  scale_x_continuous(labels = scales::label_comma(accuracy =0.01))+
  scale_y_continuous(labels = scales::label_comma(accuracy =0.01)) -> Fig2b;Fig2b

#####Fungi
fungi <- readRDS("fungi_trans_Jun.rds")
group1 <- read.xlsx("field.xlsx", sheet = "group",colNames = T, rowNames = T)
group1 <- subset(group1, Type != "CK")
fungi <- as.matrix(fungi)
data <- fungi[rownames(fungi) %in% rownames(group1),]
fungi1 <- data[,colnames(data) %in% rownames(group1)]
pcoa <- cmdscale(fungi1, k = 2, eig = TRUE)
poi <- pcoa$points
eigval <- pcoa$eig
pcoa_eig <- (pcoa$eig)[1:2]/sum(pcoa$eig)
pcoa_eig <- round(pcoa_eig, 3)
pcoa <- as.data.frame(poi)
names(pcoa)[1:2] <- c("PCoA1", "PCoA2")
pcoa$SampleID <- rownames(pcoa)
sample_site <- cbind(pcoa,group1[,c(1,2,3,4,5,6)])
sample_site$group1 <- paste0(sample_site$Site,"_",sample_site$Soil,"_",sample_site$Seedling_source)
sample_site_Jun <- subset(sample_site, Time == "Jun")
sample_site <- subset(sample_site, Time == "Aug")

###Organize the data set
df_trans <- aggregate(cbind(PCoA1, PCoA2) ~ group1+Species+Soil_type, sample_site, FUN = mean)
df2_trans <- df_trans %>%
  filter(sapply(strsplit(as.character(group1), "_"), function(x) x[1] == x[2]))
df3_trans <- df_trans %>%
  filter(!group1 %in% df2_trans$group1)
result_trans <- df3_trans %>%
  mutate(matched_value = NA, PCoA1_matched = NA, PCoA2_matched = NA)
for (i in 1:nrow(df3_trans)) {
  group1_values <- unlist(strsplit(as.character(df3_trans$group1[i]), "_"))
  matching_rows <- df2_trans %>%
    filter(sapply(strsplit(as.character(group1), "_"), function(x) {
      x[1] == group1_values[1] && x[3] == group1_values[3]
    }) & Species == df3_trans$Species[i])
  if (nrow(matching_rows) > 0) {
    result_trans$matched_value[i] <- matching_rows$group1[1] 
    result_trans$PCoA1_matched[i] <- matching_rows$PCoA1[1]
    result_trans$PCoA2_matched[i] <- matching_rows$PCoA2[1]
  }
}
df_trans1 <- df_trans %>%
  separate(group1, into = c("Site", "Soil", "Seedling_source"), sep = "_")
sample_site_Jun$Site <- factor(sample_site_Jun$Site)
sample_site$Site <- factor(sample_site$Site)
#df_situ1$group <- paste0(df_situ1$Site,"_",df_situ1$Soilsource1)
df_trans1$group <- paste0(df_trans1$Site,"_",df_trans1$Soil_type)
sample_site$group <- paste0(sample_site$Site,"_",sample_site$Soil_type)
ggplot() +
  #geom_segment(data = result_situ,mapping = aes(x = PCoA1, y = PCoA2, 
  #                          xend = PCoA1_matched, yend = PCoA2_matched),curvature = 0.3, size = 0.1, alpha = 0.7, color = "#FF8C00") +
  geom_segment(data = result_trans,mapping = aes(x = PCoA1_matched,  y = PCoA2_matched, 
                                                 xend = PCoA1, yend = PCoA2),curvature = 0.3, size = 0.1, alpha = 0.7, color = "#3D505A") +
  #geom_point(data = df_situ1,  mapping = aes(x = PCoA1, y = PCoA2, 
  #                        color = Species, shape = group, fill = Soilsource1),size = 2) +
  geom_point(data = df_trans1, mapping = aes(x = PCoA1, y = PCoA2, 
                                             color = Species, shape = group, fill = Soil_type),size = 2,alpha = 0.8) +
  geom_point(data = sample_site, mapping = aes(x = PCoA1, y = PCoA2, 
                                               color = Species, fill = Soil_type, shape = group),size = 1.2,alpha = 0.2) +
  geom_point(data = sample_site_Jun, mapping = aes(x = PCoA1, y = PCoA2, 
                                                   color = Species,shape = Site),size = 2,stroke = 1.2,alpha = 0.8) +
  scale_color_manual(values = c("#974FA3", "#4CAB3B", "#B7352A")) +
  scale_shape_manual(values = c(Low_away=15,Low_home=22,Mid_away=16,Mid_home=21,High_away=17,High_home=24,Low=3,Mid=4,High=8)) +
  scale_fill_manual(values = c("white","white","white"))+
  labs(x=paste("PCoA1 (", format(100 * pcoa_eig[1], digits=3), "%)", sep=""),
       y=paste("PCoA2 (", format(100 * pcoa_eig[2], digits=3), "%)", sep=""))+
  #xlim(-0.05, 0.13)+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(text=element_text(size=12))+
  theme(panel.background = element_rect(fill='white', colour='black'),
        panel.grid=element_blank(), 
        legend.position = "none",
        axis.text = element_text(color = "black"),
        legend.title=element_text(size = 12,face = "bold"),
        legend.text=element_text(size= 10))+
  geom_vline(aes(xintercept = 0),linetype="dashed",size = 0.5)+
  geom_hline(aes(yintercept = 0),linetype="dashed",size = 0.5)+
  scale_x_continuous(labels = scales::label_comma(accuracy =0.01))+
  scale_y_continuous(labels = scales::label_comma(accuracy =0.01)) -> Fig2c;Fig2c
Fig2b + Fig2c + Fig2a 

#################Fig 2d & S3############
trans_group <- read.xlsx("field.xlsx",
                         sheet = "feast_result", colNames = T, rowNames = T)

trans <- melt(trans_group, measure.vars = c("initial", "local",
                                            "loca_leaf", "Unknown"),
              variable_name = "Proportion")
df <- data_summary(trans, varname="value",
                   groupnames=c("Species","type","Proportion"))
df
Species <- unique(df$Species)
Type <- unique(df$type)
plots <- list()
for (i in Species) {
  for (j in Type) {
    filtered_df <- df[df$Species == i & df$type == j, ]
    p <- ggplot(filtered_df, aes(x="", y= value, fill=Proportion)) +
      geom_bar(width=1, stat="identity", color = "black",size = 0.5) +
      theme_classic(base_size=12) +
      theme(axis.line=element_blank(),
            axis.text=element_blank(),
            legend.position="none",
            axis.ticks=element_blank(),
            panel.grid=element_blank(),
            plot.title=element_text(hjust=0.5)) +
      coord_polar(theta="y", start=0) +
      scale_fill_manual(values=c("initial"="#E09513", 
                                 "loca_leaf"="#BDB985", 
                                 "local"="#C68267", 
                                 "Unknown"="grey"))
    plots[[paste(i, j, sep="_")]] <- p
  }
}
FigS3 <- wrap_plots(plots, ncol = 3)
FigS3

##non speciation
filtered_df <- data_summary(trans, varname="value",
                   groupnames=c("type","Proportion"))
filtered_df
Type <- unique(filtered_df$type)
plots <- list()

for (j in Type) {

  filtered_df1 <- filtered_df[filtered_df$type == j, ]
  
  p <- ggplot(filtered_df1, aes(x="", y= value, fill=Proportion)) +
    geom_bar(width=1, stat="identity", color = "black",size = 0.5) +
    theme_classic(base_size=12) +
    theme(axis.line=element_blank(),
          axis.text=element_blank(),
          legend.position="none",
          axis.ticks=element_blank(),
          panel.grid=element_blank(),
          plot.title=element_text(hjust=0.5)) +
    coord_polar(theta="y", start=0) +
    scale_fill_manual(values=c("initial"="#E09513", 
                               "loca_leaf"="#BDB985", 
                               "local"="#C68267", 
                               "Unknown"="grey"))
  
  plots[[paste(i, j, sep="_")]] <- p
}

Fig2d <- wrap_plots(plots, ncol = 3)
Fig2d

#############Fig 2e-f  Relationship between similarity distance and elevation#####
###rhizos_bacteria
rhizos <- readRDS("rhizos_trans.rds")# Load rhizosphere microbial data
rhizos <- as.matrix(rhizos)# Convert data to matrix format
group <- read.xlsx("field.xlsx", sheet = "group",colNames = T, rowNames = T)# Load grouping information
group <- subset(group, Time == "Aug")
group$group <- paste0(group$Site,"_",group$Soil)
unique(group$group)
#i = "F. nilgerrensis"
Species = unique(group$Species)
results = data.frame()
# Loop through each species
for (i in Species) {
  group_trans <- rownames(subset(group, Site == "Low" & Soil == "Low" & Species == i))
  group_initial1 <-rownames(subset(group, Site == "Mid" & Soil == "Low" & Species == i))
  group_initial2 <- rownames(subset(group, Site == "High" & Soil == "Low" & Species == i))
  # Calculate community differences (distance)
  distance <- melt(rhizos[group_trans, group_initial1, drop = FALSE])
  distance$Site <- "Mid"
  distance$Soil <- "Low"
  distance$Species <- i
  
  distance1 <- melt(rhizos[group_trans, group_initial2, drop = FALSE])
  distance1$Site <- "High"
  distance1$Soil <- "Low"
  distance1$Species <- i
  # Combine results
  distance2 <- rbind(distance, distance1)
  results = rbind(distance2, results) 
}
results
results1 = data.frame()
for (i in Species) {
  group_trans <- rownames(subset(group, Site == "Mid" & Soil == "Mid" & Species == i))
  group_initial1 <-rownames(subset(group, Site == "Low" & Soil == "Mid" & Species == i))
  group_initial2 <- rownames(subset(group, Site == "High" & Soil == "Mid" & Species == i))
  distance <- melt(rhizos[group_trans, group_initial1, drop = FALSE])
  distance$Site <- "Low"
  distance$Soil <- "Mid"
  distance$Species <- i
  
  distance1 <- melt(rhizos[group_trans, group_initial2, drop = FALSE])
  distance1$Site <- "High"
  distance1$Soil <- "Mid"
  distance1$Species <- i
  
  distance2 <- rbind(distance, distance1)
  results1 = rbind(distance2, results1) 
}
results1 
results2 = data.frame()
for (i in Species) {
  group_trans <- rownames(subset(group, Site == "High" & Soil == "High" & Species == i))
  group_initial1 <-rownames(subset(group, Site == "Low" & Soil == "High" & Species == i))
  group_initial2 <- rownames(subset(group, Site == "Mid" & Soil == "High" & Species == i))
  distance <- melt(rhizos[group_trans, group_initial1, drop = FALSE])
  distance$Site <- "Low"
  distance$Soil <- "High"
  distance$Species <- i
  
  distance1 <- melt(rhizos[group_trans, group_initial2, drop = FALSE])
  distance1$Site <- "Mid"
  distance1$Soil <- "High"
  distance1$Species <- i
  
  distance2 <- rbind(distance, distance1)
  results2 = rbind(distance2, results2) 
}
results2
results3 <- rbind(results, results1,results2) # Combine results, results1, and results2
# Convert Site and Soil to numeric values
results3 <- results3 %>%
  mutate(
    Site = case_when(
      Site == "Low" ~ 1590,
      Site == "Mid" ~ 1946,
      Site == "High" ~ 2641,
      TRUE ~ as.numeric(Site)
    ),
    Soil = case_when(
      Soil == "Low" ~ 1590,
      Soil == "Mid" ~ 1946,
      Soil == "High" ~ 2641,
      TRUE ~ as.numeric(Soil)
    )
  )

results3$elevation <- as.numeric(results3$Site)-as.numeric(results3$Soil)

fit_leaf <- gam(value ~ s(elevation, bs = 'cr', k = 3), data = subset(results3, Species == "F. nilgerrensis"))
fit_leaf <- gam(value ~ s(elevation, bs = 'cr', k = 3), data = subset(results3, Species == "P. asiatica"))
fit_leaf <- gam(value ~ s(elevation, bs = 'cr', k = 3), data = subset(results3, Species == "T. repens"))
summary(fit_leaf)
# Calculate elevational offset based on species
results3 <- results3 %>%
  mutate(Elevational_offset = case_when(
    Species == "F. nilgerrensis" ~ elevation - 30, 
    Species == "P. asiatica" ~ elevation,       
    Species == "T. repens" ~ elevation + 30   
  ))

# Plot scatter points and smooth curves
ggplot()+
  geom_point(results3,mapping = aes(Elevational_offset, value,color = Species), shape= 21)+
  labs(x = "Elevation differrence", y = "Community dissimilarity")+
  scale_color_manual(values = c("#974FA3", "#4CAB3B", "#B7352A"))+
  geom_smooth(data = results3,
              aes(x = Elevational_offset, y = value,color = Species),
              method = "gam", formula = y ~ s(x, k = 3),se = F, linewidth = 1) +
  theme(panel.background = element_rect(fill='white', colour='black'),
        panel.grid=element_blank(), 
        axis.title.x=element_text(colour='black', size=12),
        axis.title.y=element_text(colour='black', size=12),
        axis.text=element_text(colour='black',size=10),
        legend.title=element_text(size = 10,face = "bold"),
        legend.text=element_text(size= 10),
        legend.position = "top") -> Fig2e;Fig2e

###FUNGI
rhizos <- readRDS("fungi_trans.rds")
rhizos <- as.matrix(rhizos)
group <- read.xlsx("field.xlsx", sheet = "group",colNames = T, rowNames = T)
group <- subset(group, Time == "Aug")
group$group <- paste0(group$Site,"_",group$Soil)
unique(group$group)

Species = unique(group$Species)
results = data.frame()
for (i in Species) {
  group_trans <- rownames(subset(group, Site == "Low" & Soil == "Low" & Species == i))
  group_initial1 <-rownames(subset(group, Site == "Mid" & Soil == "Low" & Species == i))
  group_initial2 <- rownames(subset(group, Site == "High" & Soil == "Low" & Species == i))
  distance <- melt(rhizos[group_trans, group_initial1, drop = FALSE])
  distance$Site <- "Mid"
  distance$Soil <- "Low"
  distance$Species <- i
  
  distance1 <- melt(rhizos[group_trans, group_initial2, drop = FALSE])
  distance1$Site <- "High"
  distance1$Soil <- "Low"
  distance1$Species <- i
  
  distance2 <- rbind(distance, distance1)
  results = rbind(distance2, results) 
}
results
results1 = data.frame()
for (i in Species) {
  group_trans <- rownames(subset(group, Site == "Mid" & Soil == "Mid" & Species == i))
  group_initial1 <-rownames(subset(group, Site == "Low" & Soil == "Mid" & Species == i))
  group_initial2 <- rownames(subset(group, Site == "High" & Soil == "Mid" & Species == i))
  distance <- melt(rhizos[group_trans, group_initial1, drop = FALSE])
  distance$Site <- "Low"
  distance$Soil <- "Mid"
  distance$Species <- i
  
  distance1 <- melt(rhizos[group_trans, group_initial2, drop = FALSE])
  distance1$Site <- "High"
  distance1$Soil <- "Mid"
  distance1$Species <- i
  
  distance2 <- rbind(distance, distance1)
  results1 = rbind(distance2, results1) 
}
results1 
results2 = data.frame()
for (i in Species) {
  group_trans <- rownames(subset(group, Site == "High" & Soil == "High" & Species == i))
  group_initial1 <-rownames(subset(group, Site == "Low" & Soil == "High" & Species == i))
  group_initial2 <- rownames(subset(group, Site == "Mid" & Soil == "High" & Species == i))
  distance <- melt(rhizos[group_trans, group_initial1, drop = FALSE])
  distance$Site <- "Low"
  distance$Soil <- "High"
  distance$Species <- i
  
  distance1 <- melt(rhizos[group_trans, group_initial2, drop = FALSE])
  distance1$Site <- "Mid"
  distance1$Soil <- "High"
  distance1$Species <- i
  
  distance2 <- rbind(distance, distance1)
  results2 = rbind(distance2, results2) 
}
results2
results3 <- rbind(results, results1,results2) 
results3 <- results3 %>%
  mutate(
    Site = case_when(
      Site == "Low" ~ 1590,
      Site == "Mid" ~ 1946,
      Site == "High" ~ 2641,
      TRUE ~ as.numeric(Site)
    ),
    Soil = case_when(
      Soil == "Low" ~ 1590,
      Soil == "Mid" ~ 1946,
      Soil == "High" ~ 2641,
      TRUE ~ as.numeric(Soil)
    )
  )
results3$elevation <- as.numeric(results3$Site)-as.numeric(results3$Soil)

fit_leaf <- gam(value ~ s(elevation, bs = 'cr', k = 3), data = subset(results3, Species == "F. nilgerrensis"))
fit_leaf <- gam(value ~ s(elevation, bs = 'cr', k = 3), data = subset(results3, Species == "P. asiatica"))
fit_leaf <- gam(value ~ s(elevation, bs = 'cr', k = 3), data = subset(results3, Species == "T. repens"))
summary(fit_leaf)

results3 <- results3 %>%
  mutate(Elevational_offset = case_when(
    Species == "F. nilgerrensis" ~ elevation - 30, 
    Species == "P. asiatica" ~ elevation,       
    Species == "T. repens" ~ elevation + 30  
  ))

ggplot()+
  geom_point(results3,mapping = aes(Elevational_offset, value,color = Species), shape= 21)+
  labs(x = "Elevation differrence", y = "Community dissimilarity")+
  scale_color_manual(values = c("#974FA3", "#4CAB3B", "#B7352A"))+
  geom_smooth(data = results3,
              aes(x = Elevational_offset, y = value,color = Species),
              method = "gam", formula = y ~ s(x, k = 3),se = F, linewidth = 1) +
  theme(panel.background = element_rect(fill='white', colour='black'),
        panel.grid=element_blank(), 
        axis.title.x=element_text(colour='black', size=12),
        axis.title.y=element_text(colour='black', size=12),
        axis.text=element_text(colour='black',size=10),
        legend.title=element_text(size = 10,face = "bold"),
        legend.text=element_text(size= 10),
        legend.position = "top") -> Fig2f;Fig2f
Fig2e+Fig2f

##########Fig S1########
rm(list=ls())
library(openxlsx)
library(vegan)
library(adespatial)
library(ggplot2)
library(phyloseq)
library(picante)
library(ggpmisc)
############community dissimilarity & elevation##############
PATH <- readRDS("PATH.rds")
AMF <- readRDS("AMF.rds")
group <- read.xlsx("path_AMF.xlsx", sheet = "group",colNames = T, rowNames = T)
#i = "F"
Species <- unique(group$Species)
data2 <- data.frame()
for (i in Species) {
  library(NST)
  group1 = subset(group, Species == i)
  Euce <- vegdist(group1[,2], method = "euclidean")
  Euce1 <- dist.3col(Euce)
  
  PATH <- as.matrix(PATH)
  data <- PATH[rownames(PATH) %in% rownames(group1),]
  PATH1 <- data[,colnames(data) %in% rownames(group1)]
  PATH2 <- dist.3col(PATH1)
  
  AMF <- as.matrix(AMF)
  data <- AMF[rownames(AMF) %in% rownames(group1),]
  AMF1 <- data[,colnames(data) %in% rownames(group1)]
  AMF2 <- dist.3col(AMF1)
  
  data <- cbind(Euce1,PATH2[,3],AMF2[,3])
  data$Species <- i
  data2 <- rbind(data2, data)
}
data2
colnames(data2) <- c("name1","name2","elevation","path","AMF","Species")

data3 <- subset(data2,elevation != "0")

ggplot(data3, aes(x = elevation, y = path, color = Species, group = Species)) +
  geom_point(size = 1.5, shape = 21) +  
  geom_smooth(method = "loess", se = F, size = 0.6) + 
  labs(y = "Community dissimilarity \n weighted Unifrac",
       x = "Elevation difference (m)", title = NULL) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.line = element_line(linetype = 1, color = "black", size = 0.5),
        axis.ticks = element_line(color = "black", size = 0.5, lineend = 2),
        axis.text = element_text(colour = 'black', size = 9),
        axis.title = element_text(colour = 'black', size = 10)) +
  scale_fill_manual(values = c("#974FA3", "#4CAB3B", "#B7352A")) +
  scale_color_manual(values = c("#974FA3", "#4CAB3B", "#B7352A")) -> Fig_S1a;Fig_S1a

ggplot(group, aes(x = elevation, y = AMF, color = Species, group = Species)) +
  geom_point(size = 1.5, shape = 21) + 
  geom_smooth(method = "loess", se = F, size = 0.6) +  
  labs(y = "Community dissimilarity \n weighted Unifrac",
       x = "Elevation difference (m)", title = NULL) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.line = element_line(linetype = 1, color = "black", size = 0.5),
        axis.ticks = element_line(color = "black", size = 0.5, lineend = 2),
        axis.text = element_text(colour = 'black', size = 9),
        axis.title = element_text(colour = 'black', size = 10)) +
  scale_fill_manual(values = c("#974FA3", "#4CAB3B", "#B7352A")) +
  scale_color_manual(values = c("#974FA3", "#4CAB3B", "#B7352A")) -> Fig_S1b;Fig_S1b

ggplot(group,mapping = aes(x = Elevation, y = ra_path,color = Species,group = Species)) +
  geom_point(group,mapping = aes(x = Elevation, y = ra_path,color = Species,group = Species),
             position = position_dodge(width = 0.7), size = 1.5, shape = 21) + 
  geom_smooth(method = "loess", se = F, size = 0.6) + 
  labs(y = "Relative abundance of pathogen (%)",x = "Elevation (m)",title = NULL) +
  theme_bw()+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        legend.position = "none",
        axis.line=element_line(linetype=1,color="black",size=0.5),
        axis.ticks=element_line(color="black",size=0.5,lineend = 2),
        axis.text=element_text(colour='black',size=9),
        axis.title=element_text(colour='black', size=10))+
  scale_fill_manual(values = c("#974FA3", "#4CAB3B", "#B7352A"))+
  scale_color_manual(values = c("#974FA3", "#4CAB3B", "#B7352A"))+
  scale_x_continuous(breaks = c(1590,1946,2641))  -> Fig_S1c;Fig_S1c

ggplot(group,mapping = aes(x = Elevation, y = ra_AMF,color = Species,group = Species)) +
  geom_point(group,mapping = aes(x = Elevation, y = ra_AMF,color = Species,group = Species),
             position = position_dodge(width = 0.7), size = 1.5, shape = 21) + 
  geom_smooth(method = "loess", se = F, size = 0.6) + 
  labs(y = "Relative abundance of AMF (%)",x = "Elevation (m)",title = NULL) +
  theme_bw()+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        legend.position = "none",
        axis.line=element_line(linetype=1,color="black",size=0.5),
        axis.ticks=element_line(color="black",size=0.5,lineend = 2),
        axis.text=element_text(colour='black',size=9),
        axis.title=element_text(colour='black', size=10))+
  scale_fill_manual(values = c("#974FA3", "#4CAB3B", "#B7352A"))+
  scale_color_manual(values = c("#974FA3", "#4CAB3B", "#B7352A"))+
  scale_x_continuous(breaks = c(1590,1946,2641))  -> Fig_S1d;Fig_S1d

(pFig_S1a+Fig_S1b) / (Fig_S1c+Fig_S1d)


####end
