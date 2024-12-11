##############################Fig 2###################
library(openxlsx)
library(vegan)
library(magrittr)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(patchwork)
library(ggsci)
library(ggpmisc)
library(reshape)
library(patchwork)
library(emmeans)
library(agricolae)
library(multcomp)
set.seed(1234)

#######PCoA analysis
###phyllosphere
phyllo <- readRDS("phyllo_trans.rds")
group1 <- read.xlsx("field.xlsx", sheet = "group",colNames = T, rowNames = T)
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
sample_site <- cbind(pcoa,group1[,c(2,3,4,5,6)])
sample_site$group1 <- paste0(sample_site$Site,"_",sample_site$Soil,"_",sample_site$Popu)

###Organize the data set
df_trans <- aggregate(cbind(PCoA1, PCoA2) ~ group1+Species+Soilsource, sample_site, FUN = mean)
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
  separate(group1, into = c("Site", "Soil", "Popu"), sep = "_")

sample_site$Site <- factor(sample_site$Site)
#df_situ1$group <- paste0(df_situ1$Site,"_",df_situ1$Soilsource1)
df_trans1$group <- paste0(df_trans1$Site,"_",df_trans1$Soilsource)
sample_site$group <- paste0(sample_site$Site,"_",sample_site$Soilsource)
ggplot() +
  #geom_segment(data = result_situ,mapping = aes(x = PCoA1, y = PCoA2, 
  #                          xend = PCoA1_matched, yend = PCoA2_matched),curvature = 0.3, size = 0.1, alpha = 0.7, color = "#FF8C00") +
  geom_segment(data = result_trans,mapping = aes(x = PCoA1_matched,  y = PCoA2_matched, 
                                                 xend = PCoA1, yend = PCoA2),curvature = 0.3, size = 0.1, alpha = 0.7, color = "#3D505A") +
  #geom_point(data = df_situ1,  mapping = aes(x = PCoA1, y = PCoA2, 
  #                        color = Species, shape = group, fill = Soilsource1),size = 2) +
  geom_point(data = df_trans1, mapping = aes(x = PCoA1, y = PCoA2, 
                                             color = Species, shape = group, fill = Soilsource),size = 2) +
  geom_point(data = sample_site, mapping = aes(x = PCoA1, y = PCoA2, 
                                               color = Species, fill = Soilsource, shape = group),size = 1.2,alpha = 0.2) +
  scale_color_manual(values = c("#974FA3", "#4CAB3B", "#3C72AF")) +
  scale_shape_manual(values = c(15,22,16,21,17,24)) +
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
rhizos <- readRDS("rhizos_trans.rds")
group <- read.xlsx("field.xlsx", sheet = "group",colNames = T, rowNames = T)
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
sample_site <- cbind(pcoa,group1[,c(2,3,4,5,6)])
sample_site$group1 <- paste0(sample_site$Site,"_",sample_site$Soil,"_",sample_site$Popu)

###Organize the data set
df_trans <- aggregate(cbind(PCoA1, PCoA2) ~ group1+Species+Soilsource, sample_site, FUN = mean)
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
  separate(group1, into = c("Site", "Soil", "Popu"), sep = "_")

sample_site$Site <- factor(sample_site$Site)
df_trans1$group <- paste0(df_trans1$Site,"_",df_trans1$Soilsource)
sample_site$group <- paste0(sample_site$Site,"_",sample_site$Soilsource)
ggplot() +
  #geom_segment(data = result_situ,mapping = aes(x = PCoA1, y = PCoA2, 
  #                          xend = PCoA1_matched, yend = PCoA2_matched),curvature = 0.3, size = 0.1, alpha = 0.7, color = "#FF8C00") +
  geom_segment(data = result_trans,mapping = aes(x = PCoA1_matched,  y = PCoA2_matched, 
                                                 xend = PCoA1, yend = PCoA2),curvature = 0.3, size = 0.1, alpha = 0.7, color = "#3D505A") +
  #geom_point(data = df_situ1,  mapping = aes(x = PCoA1, y = PCoA2, 
  #                        color = Species, shape = group, fill = Soilsource1),size = 2) +
  geom_point(data = df_trans1, mapping = aes(x = PCoA1, y = PCoA2, 
                                             color = Species, shape = group, fill = Soilsource),size = 2) +
  geom_point(data = sample_site, mapping = aes(x = PCoA1, y = PCoA2, 
                                               color = Species, fill = Soilsource, shape = group),size = 1.2,alpha = 0.2) +
  scale_color_manual(values = c("#974FA3", "#4CAB3B", "#3C72AF")) +
  scale_shape_manual(values = c(15,22,16,21,17,24)) +
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
  scale_y_continuous(labels = scales::label_comma(accuracy =0.01)) -> Fig2b;Fig2b

#####Fungi
fungi <- readRDS("fungi_trans.rds")
group <- read.xlsx("field.xlsx", sheet = "group",colNames = T, rowNames = T)
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
sample_site <- cbind(pcoa,group1[,c(2,3,4,5,6)])
sample_site$group1 <- paste0(sample_site$Site,"_",sample_site$Soil,"_",sample_site$Popu)
###Organize the data set
df_trans <- aggregate(cbind(PCoA1, PCoA2) ~ group1+Species+Soilsource, sample_site, FUN = mean)
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
  separate(group1, into = c("Site", "Soil", "Popu"), sep = "_")

sample_site$Site <- factor(sample_site$Site)
#df_situ1$group <- paste0(df_situ1$Site,"_",df_situ1$Soilsource1)
df_trans1$group <- paste0(df_trans1$Site,"_",df_trans1$Soilsource)
sample_site$group <- paste0(sample_site$Site,"_",sample_site$Soilsource)
ggplot() +
  #geom_segment(data = result_situ,mapping = aes(x = PCoA1, y = PCoA2, 
  #                          xend = PCoA1_matched, yend = PCoA2_matched),curvature = 0.3, size = 0.1, alpha = 0.7, color = "#FF8C00") +
  geom_segment(data = result_trans,mapping = aes(x = PCoA1_matched,  y = PCoA2_matched, 
                                                 xend = PCoA1, yend = PCoA2),curvature = 0.3, size = 0.1, alpha = 0.7, color = "#3D505A") +
  #geom_point(data = df_situ1,  mapping = aes(x = PCoA1, y = PCoA2, 
  #                        color = Species, shape = group, fill = Soilsource1),size = 2) +
  geom_point(data = df_trans1, mapping = aes(x = PCoA1, y = PCoA2, 
                                             color = Species, shape = group, fill = Soilsource),size = 2) +
  geom_point(data = sample_site, mapping = aes(x = PCoA1, y = PCoA2, 
                                               color = Species, fill = Soilsource, shape = group),size = 1.2,alpha = 0.2) +
  scale_color_manual(values = c("#974FA3", "#4CAB3B", "#3C72AF")) +
  scale_shape_manual(values = c(15,22,16,21,17,24)) +
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
  scale_y_continuous(labels = scales::label_comma(accuracy =0.01)) -> Fig2c;Fig2c

#######Feast_result
library(ggplot2)
library(ggsci)
library(magrittr)
library(patchwork)
library(ggpmisc)
library(reshape)
library(openxlsx)
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
FigS2 <- wrap_plots(plots, ncol = 3)
FigS2

##non speciation
Type <- unique(df$type)
plots <- list()

for (j in Type) {

  filtered_df <- df[df$type == j, ]
  
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

Fig2d <- wrap_plots(plots, ncol = 3)
Fig2d


###Relationship between similarity distance and altitude
###rhizos_bacteria
rhizos <- readRDS("rhizos_trans.rds")
rhizos <- as.matrix(rhizos)
group <- read.xlsx("field.xlsx", sheet = "group",colNames = T, rowNames = T)
group$group <- paste0(group$Site,"_",group$Soil)
unique(group$group)
i = "F. nilgerrensis"
Species = unique(group$Species)
results = data.frame()
for (i in Species) {
  group_trans <- rownames(subset(group, Site == "1590" & Soil == "1590" & Species == i))
  group_initial1 <-rownames(subset(group, Site == "1946" & Soil == "1590" & Species == i))
  group_initial2 <- rownames(subset(group, Site == "2641" & Soil == "1590" & Species == i))
  distance <- melt(rhizos[group_trans, group_initial1, drop = FALSE])
  distance$Site <- "1946"
  distance$Soil <- "1590"
  distance$Species <- i
  
  distance1 <- melt(rhizos[group_trans, group_initial2, drop = FALSE])
  distance1$Site <- "2641"
  distance1$Soil <- "1590"
  distance1$Species <- i
  
  distance2 <- rbind(distance, distance1)
  results = rbind(distance2, results) 
}
results
results1 = data.frame()
for (i in Species) {
  group_trans <- rownames(subset(group, Site == "1946" & Soil == "1946" & Species == i))
  group_initial1 <-rownames(subset(group, Site == "1590" & Soil == "1946" & Species == i))
  group_initial2 <- rownames(subset(group, Site == "2641" & Soil == "1946" & Species == i))
  distance <- melt(rhizos[group_trans, group_initial1, drop = FALSE])
  distance$Site <- "1590"
  distance$Soil <- "1946"
  distance$Species <- i
  
  distance1 <- melt(rhizos[group_trans, group_initial2, drop = FALSE])
  distance1$Site <- "2641"
  distance1$Soil <- "1946"
  distance1$Species <- i
  
  distance2 <- rbind(distance, distance1)
  results1 = rbind(distance2, results1) 
}
results1 
results2 = data.frame()
for (i in Species) {
  group_trans <- rownames(subset(group, Site == "2641" & Soil == "2641" & Species == i))
  group_initial1 <-rownames(subset(group, Site == "1590" & Soil == "2641" & Species == i))
  group_initial2 <- rownames(subset(group, Site == "1946" & Soil == "2641" & Species == i))
  distance <- melt(rhizos[group_trans, group_initial1, drop = FALSE])
  distance$Site <- "1590"
  distance$Soil <- "2641"
  distance$Species <- i
  
  distance1 <- melt(rhizos[group_trans, group_initial2, drop = FALSE])
  distance1$Site <- "1946"
  distance1$Soil <- "2641"
  distance1$Species <- i
  
  distance2 <- rbind(distance, distance1)
  results2 = rbind(distance2, results2) 
}
results2
results3 <- rbind(results, results1,results2) 
results3$elevation <- as.numeric(results3$Site)-as.numeric(results3$Soil)

ggplot()+
  geom_point(results3,mapping = aes(elevation, value,color = Species), shape= 21)+
  labs(x = "Elevation differrence", y = "Community dissimilarity")+
  scale_color_manual(values = c("#974FA3", "#4CAB3B", "#3C72AF"))+
  stat_smooth(results3,mapping = aes(elevation, value),
              method = "lm",
              formula = y ~ I(x^2),
              se = T,
              alpha = 0.5)+
  stat_poly_eq(results3,mapping = aes(elevation, value, label =  paste(..rr.label..,..p.value.label..,  sep = "~~~~")),
               formula = y~I(x^2), parse = TRUE)+
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
group$group <- paste0(group$Site,"_",group$Soil)
unique(group$group)
i = "F. nilgerrensis"
Species = unique(group$Species)
results = data.frame()
for (i in Species) {
  group_trans <- rownames(subset(group, Site == "1590" & Soil == "1590" & Species == i))
  group_initial1 <-rownames(subset(group, Site == "1946" & Soil == "1590" & Species == i))
  group_initial2 <- rownames(subset(group, Site == "2641" & Soil == "1590" & Species == i))
  distance <- melt(rhizos[group_trans, group_initial1, drop = FALSE])
  distance$Site <- "1946"
  distance$Soil <- "1590"
  distance$Species <- i
  
  distance1 <- melt(rhizos[group_trans, group_initial2, drop = FALSE])
  distance1$Site <- "2641"
  distance1$Soil <- "1590"
  distance1$Species <- i
  
  distance2 <- rbind(distance, distance1)
  results = rbind(distance2, results) 
}
results
results1 = data.frame()
for (i in Species) {
  group_trans <- rownames(subset(group, Site == "1946" & Soil == "1946" & Species == i))
  group_initial1 <-rownames(subset(group, Site == "1590" & Soil == "1946" & Species == i))
  group_initial2 <- rownames(subset(group, Site == "2641" & Soil == "1946" & Species == i))
  distance <- melt(rhizos[group_trans, group_initial1, drop = FALSE])
  distance$Site <- "1590"
  distance$Soil <- "1946"
  distance$Species <- i
  
  distance1 <- melt(rhizos[group_trans, group_initial2, drop = FALSE])
  distance1$Site <- "2641"
  distance1$Soil <- "1946"
  distance1$Species <- i
  
  distance2 <- rbind(distance, distance1)
  results1 = rbind(distance2, results1) 
}
results1 
results2 = data.frame()
for (i in Species) {
  group_trans <- rownames(subset(group, Site == "2641" & Soil == "2641" & Species == i))
  group_initial1 <-rownames(subset(group, Site == "1590" & Soil == "2641" & Species == i))
  group_initial2 <- rownames(subset(group, Site == "1946" & Soil == "2641" & Species == i))
  distance <- melt(rhizos[group_trans, group_initial1, drop = FALSE])
  distance$Site <- "1590"
  distance$Soil <- "2641"
  distance$Species <- i
  
  distance1 <- melt(rhizos[group_trans, group_initial2, drop = FALSE])
  distance1$Site <- "1946"
  distance1$Soil <- "2641"
  distance1$Species <- i
  
  distance2 <- rbind(distance, distance1)
  results2 = rbind(distance2, results2) 
}
results2
results3 <- rbind(results, results1,results2) 
results3$elevation <- as.numeric(results3$Site)-as.numeric(results3$Soil)

ggplot()+
  geom_point(results3,mapping = aes(elevation, value,color = Species), shape= 21)+
  labs(x = "Elevation differrence", y = "Community dissimilarity")+
  scale_color_manual(values = c("#974FA3", "#4CAB3B", "#3C72AF"))+
  stat_smooth(results3,mapping = aes(elevation, value),
              method = "lm",
              formula = y ~ I(x^2),
              se = T,
              alpha = 0.5)+
  stat_poly_eq(results3,mapping = aes(elevation, value, label =  paste(..rr.label..,..p.value.label..,  sep = "~~~~")),
               formula = y~I(x^2), parse = TRUE)+
  theme(panel.background = element_rect(fill='white', colour='black'),
        panel.grid=element_blank(), 
        axis.title.x=element_text(colour='black', size=12),
        axis.title.y=element_text(colour='black', size=12),
        axis.text=element_text(colour='black',size=10),
        legend.title=element_text(size = 10,face = "bold"),
        legend.text=element_text(size= 10),
        legend.position = "top") -> Fig2f;Fig2f

####putative fungal pathogens
###不同海拔的病原菌
group <- read.xlsx("field.xlsx", sheet = "group",colNames = T, rowNames = T)
group$Site <- factor(group$Site, levels = c("1590","1946","2641"))
group$Popu <- factor(group$Popu, levels = c("1590","1946","2641"))
group$Soil <- factor(group$Soil, levels = c("1590","1946","2641"))

###rich_path
group <- subset(group, Type == "CK")
group$rich_path <- log(group$rich_path)

df2 <- data_summary(group, varname = "rich_path",groupnames = c("Site","Species"))
df2$Site <- factor(df2$Site, levels = c("1590","1946","2641"))

ggplot(data = df2) +
  geom_errorbar(df2, mapping = aes(x= Site, y = rich_path, ymin = rich_path-se, ymax = rich_path+se, color = Species),
                position=position_dodge(width = 0.8), width = 0, size = 0.5, alpha = 1) +
  geom_bar(df2,mapping = aes(x = Site, y = rich_path, fill = Species),
           position = position_dodge(width = 0.8), stat = "identity", size = 0.7, width = 0.7, alpha = 0.8) +
  scale_fill_manual(values = c("#974FA3", "#4CAB3B", "#3C72AF")) +
  scale_color_manual(values = c("#974FA3", "#4CAB3B", "#3C72AF")) +
  #scale_shape_manual(values = c(16,21)) +
  labs(y = "Richness of pathogen (log)", title = NULL) +
  theme_bw() +
  #facet_wrap( ~ Popu, scales = "free") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.line = element_line(linetype = 1, color = "black", size = 0.5),
        axis.ticks = element_line(color = "black", size = 0.5, lineend = 2),
        axis.text = element_text(colour = 'black', size = 9),
        axis.title = element_text(colour = 'black', size = 10)) -> FigS3a;FigS3a

###ra_path
df2 <- data_summary(group, varname = "ra_path",groupnames = c("Site","Species"))
df2$Site <- factor(df2$Site, levels = c("1590","1946","2641"))

ggplot(data = df2) +
  geom_errorbar(df2, mapping = aes(x= Site, y = ra_path, ymin = ra_path-se, ymax = ra_path+se, color = Species),
                position=position_dodge(width = 0.8), width = 0, size = 0.5, alpha = 1) +
  geom_bar(df2,mapping = aes(x = Site, y = ra_path, fill = Species),
           position = position_dodge(width = 0.8), stat = "identity", size = 0.7, width = 0.7, alpha = 0.8) +
  scale_fill_manual(values = c("#974FA3", "#4CAB3B", "#3C72AF")) +
  scale_color_manual(values = c("#974FA3", "#4CAB3B", "#3C72AF")) +
  #scale_shape_manual(values = c(16,21)) +
  labs(y = "Abundance of pathogen (%)", title = NULL) +
  theme_bw() +
  #facet_wrap( ~ Popu, scales = "free") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "right",
        axis.line = element_line(linetype = 1, color = "black", size = 0.5),
        axis.ticks = element_line(color = "black", size = 0.5, lineend = 2),
        axis.text = element_text(colour = 'black', size = 9),
        axis.title = element_text(colour = 'black', size = 10)) -> FigS3b;FigS3b

####end
