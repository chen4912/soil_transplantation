####################Fig 3######################
####plant growth in common garden###
library(openxlsx)
library(ggplot2)
library(ggh4x)
library(agricolae)
library(multcomp)
library(patchwork)
library(emmeans)
library(multcomp)
library(lmerTest)
library(sjPlot)
library(tidyverse)
library(glmm.hp)
library(MASS)
library(rfPermute)
library(rfUtilities)
library(MuMIn)
library(car)
library(reshape)
library(Hmisc)
group <- read.xlsx("greenhouse.xlsx", sheet = "greenhouse_PGR",colNames = T, rowNames = T)
group$plant <- factor(group$plant, levels = c("1590","1946","2641"))
group$soil <- factor(group$soil, levels = c("1590","1946","2641"))

df <- data_summary(group, varname = "PGR",groupnames = c("plant","soil","species"))
df$plant <- factor(df$plant, levels = c("1590","1946","2641"))
df$soil <- factor(df$soil, levels = c("1590","1946","2641"))
ggplot(data = df) +
  geom_errorbar(mapping = aes(x = plant, y = PGR, ymin = PGR - se, ymax = PGR + se, color = soil),
                position = position_dodge(width = 0.7), width = 0, size = 0.5, alpha = 1) +
  geom_bar(mapping = aes(x = plant, y = PGR, fill = soil),
           position = position_dodge(width = 0.7),width = 0.6, stat = "identity", size = 3, alpha = 0.8) +
  scale_fill_manual(values = c("#974FA3","#4CAB3B","#3C72AF")) +
  scale_color_manual(values = c("#974FA3","#4CAB3B","#3C72AF")) +
  scale_shape_manual(values = c(16,21)) +
  labs(y = "Plant growth rate (cm/day)", title = NULL) +
  theme_bw() +
  facet_wrap( ~ species, scales = "free") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "right",
        axis.line = element_line(linetype = 1, color = "black", size = 0.5),
        axis.ticks = element_line(color = "black", size = 0.5, lineend = 2),
        axis.text = element_text(colour = 'black', size = 9),
        axis.title = element_text(colour = 'black', size = 10)) -> FigS1;FigS1

###Multiple comparison
group$Type <- paste0(group$soil,"_",group$plant)
type <- unique(group$Type)
species <- unique(group$species)
result <- data.frame()
#j = "T. repens"
for (j in species) {
  low <- subset(group, species == j)
  model <- aov(PGR ~ Type, low)
  summary(model)
  model_means <- emmeans(object = model,
                         specs = "Type")
  model_means_cld <- cld(object = model_means,
                         adjust = "none",
                         Letters = letters,
                         alpha = 0.05,
                         decreasing = T)
  model <- data.frame(model_means_cld)
  model$species <- j
  result <- rbind(result, model)
}
result

#####deming
group <- read.xlsx("greenhouse.xlsx", sheet = "deming",colNames = T, rowNames = T)
Pasi <- subset(group, Species == "P. asiatica")
cor.test(Pasi$field,Pasi$house)
###deming regression analysis of P. asiatica
library(deming)
fit <- deming(field ~ house, data=Pasi, xstd=se1, ystd=se)
print(fit)

ggplot(Pasi, aes(x = house, y = field, color = species)) +
  geom_errorbar(aes(ymin = field - se, ymax = field + se), 
                width = 0, position = position_dodge(0.5)) +
  geom_errorbarh(aes(xmin = house - se1, xmax = house + se1, y = field), 
                 height = 0, position = position_dodge(0.5)) +
  geom_point(size = 2.5, position = position_dodge(0.5)) +
  scale_color_manual(values = c("#4CAB3B")) +
  scale_shape_manual(values = c(16,21)) +
  labs(y = "Filed", x = "Greenhouse", title = NULL) +
  xlim(0.85,1.15)+
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  geom_abline(slope = 0.94, intercept = -0.08, linetype = "dashed", color = "black") +
  theme_bw() +
  #facet_wrap( ~ species, scales = "free") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.line = element_line(linetype = 1, color = "black", size = 0.5),
        axis.ticks = element_line(color = "black", size = 0.5, lineend = 2),
        axis.text = element_text(colour = 'black', size = 9),
        axis.title = element_text(colour = 'black', size = 10)) -> Fig3a;Fig3a

###deming regression analysis of T. repens
Trep <- subset(group, Species == "T. repens")
cor.test(Trep$field,Trep$house)

fit <- deming(field ~ house, data=Trep, xstd=se1, ystd=se)
print(fit)
ggplot(Trep, aes(x = house, y = field, color = species)) +
  geom_errorbar(aes(ymin = field - se, ymax = field + se), 
                width = 0, position = position_dodge(0.5)) +
  geom_errorbarh(aes(xmin = house - se1, xmax = house + se1, y = field), 
                 height = 0, position = position_dodge(0.5)) +
  geom_point(size = 2.5, position = position_dodge(0.5)) +
  scale_color_manual(values = c("#3C72AF")) +
  scale_shape_manual(values = c(16,21)) +
  labs(y = "Filed", x = "Greenhouse", title = NULL) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  geom_abline(slope = 1.70, intercept = -0.84, linetype = "dashed", color = "black") +
  theme_bw() +
  #facet_wrap( ~ species, scales = "free") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.line = element_line(linetype = 1, color = "black", size = 0.5),
        axis.ticks = element_line(color = "black", size = 0.5, lineend = 2),
        axis.text = element_text(colour = 'black', size = 9),
        axis.title = element_text(colour = 'black', size = 10)) -> Fig3b;Fig3b

#####Multiple regression of LMM
options(na.action ="na.fail")
trans = read.xlsx("field.xlsx", sheet = "Multiple_regression", rowNames = T, colNames = T)
colnames(trans)
###强相关筛选
#trans[,c(7:27)] <- scale(trans[,c(7:27)])
soilData <- trans[,c(5,10,12,14,16,18,24:27)] 

cormat <- round(cor(soilData, method = "spearman"),2)

reorder_cormat <- function(cormat){
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  cormat <-cormat[hc$order, hc$order]
}

# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

reordered_cormat <- reorder_cormat(cormat)

upper_tri <- get_upper_tri(reordered_cormat)

# Melt the correlation matrix
melted_cormat <- melt(upper_tri, na.rm = TRUE)

#write.table(upper_tri, file ="cor_field.csv",sep =",", quote =FALSE)

p_heatmap<- ggplot(melted_cormat, aes(X1, X2, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Spearman\nCorrelation") +
  theme_bw()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1,
                                   size = 14, hjust = 1),
        axis.text.y = element_text(size=14),
        legend.text = element_text(size=14),
        legend.title = element_text(size=14),
        axis.title = element_blank())+
  coord_fixed();p_heatmap

###PGR
trans = read.xlsx("trans.xlsx", sheet = "MIRM", rowNames = T, colNames = T)
colnames(trans)
trans[,c(10,12,16,18,25,27)]  <- scale(trans[,c(10,12,16,18,25,27)])
fit<-lmer(RGR ~ Phyllo_PCoA1+Home_rhizos_PCoA1+Away_rhizos_PCoA1+Away_fungi_PCoA1+away_fungi+
            away_rhizos+(1|Species),data=trans,REML=FALSE)###full model
vif(fit)
fit_select <- dredge(fit)
head(fit_select)
best_est<-model.avg(fit_select,subset= delta <2, revised.var = TRUE)
summary(best_est)
sw(best_est)
fit<-lmer(RGR ~ Home_rhizos_PCoA1+Away_rhizos_PCoA1+
            Phyllo_PCoA1+(1|Species),data=trans)###final_model
summary(fit)
r.squaredGLMM(fit)
glmm.hp::glmm.hp(fit,type = "R2")

anova_results <- car::Anova(fit,type=2)
p_values <- anova_results$`Pr(>Chisq)`
adjusted_p_values <- p.adjust(p_values, method = "BH")

trans = read.xlsx("trans.xlsx", sheet = "MIRM", rowNames = T, colNames = T)
pd_attributes_variable1 <- attributes(scale(trans[c("Home_rhizos_PCoA1")]))
pd_attributes_variable2 <- attributes(scale(trans[c("Away_rhizos_PCoA1")]))
pd_attributes_variable3 <- attributes(scale(trans[c("Phyllo_PCoA1")]))
trans[,c(10,12,14,16,18,25:27)]  <- scale(trans[,c(10,12,14,16,18,25:27)])
fit<-lmer(RGR ~ Home_rhizos_PCoA1+Away_rhizos_PCoA1+
            Phyllo_PCoA1+(1|Species),data=trans)###final_model

###Home_rhizos_PCoA1
get.data = get_model_data(fit, type = "pred", terms = "Home_rhizos_PCoA1[all]")
get.data = as.data.frame(get.data); colnames(get.data)[1] = "Home_rhizos_PCoA1"
##### rotation
get.data["Home_rhizos_PCoA1_row"] <- (pd_attributes_variable1$`scaled:center`["Home_rhizos_PCoA1"] + pd_attributes_variable1$`scaled:scale`["Home_rhizos_PCoA1"]*get.data["Home_rhizos_PCoA1"])
trans["Home_rhizos_PCoA1_row"] <- (pd_attributes_variable1$`scaled:center`["Home_rhizos_PCoA1"] + pd_attributes_variable1$`scaled:scale`["Home_rhizos_PCoA1"]*trans["Home_rhizos_PCoA1"])

ggplot(get.data, aes(x = Home_rhizos_PCoA1_row, y = predicted)) +
  geom_ribbon(aes(ymin=conf.low,ymax=conf.high),fill="#6D6D67", alpha = I(0.3)) +
  geom_line(size=1, linetype = 1, color = "black") +
  geom_point(trans, mapping = aes(x = Home_rhizos_PCoA1_row, y = RGR, color = Species),size = 2.5,shape = 21) +
  theme_bw()+
  scale_color_manual(values = c("#974FA3", "#4CAB3B", "#3C72AF"))+
  labs(x = "PCoA1 of rhizosphere bacteria persisted \n from origin site",y = "Effect of plant growth (cm/day)", title = NULL) +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        legend.position = "none",
        axis.line=element_line(linetype=1,color="black",size=0.5),
        axis.ticks=element_line(color="black",size=0.5,lineend = 2),
        axis.text=element_text(colour='black',size=9),
        axis.title=element_text(colour='black', size=10)) -> Fig3a; Fig3a

###Phyllo_PCoA1
get.data = get_model_data(fit, type = "pred", terms = "Phyllo_PCoA1[all]")
get.data = as.data.frame(get.data); colnames(get.data)[1] = "Phyllo_PCoA1"
##### rotation
get.data["Phyllo_PCoA1_row"] <- (pd_attributes_variable3$`scaled:center`["Phyllo_PCoA1"] + pd_attributes_variable3$`scaled:scale`["Phyllo_PCoA1"]*get.data["Phyllo_PCoA1"])
trans["Phyllo_PCoA1_row"] <- (pd_attributes_variable3$`scaled:center`["Phyllo_PCoA1"] + pd_attributes_variable3$`scaled:scale`["Phyllo_PCoA1"]*trans["Phyllo_PCoA1"])

ggplot(get.data, aes(x = Phyllo_PCoA1_row, y = predicted)) +
  geom_ribbon(aes(ymin=conf.low,ymax=conf.high),fill="#6D6D67", alpha = I(0.3)) +
  geom_line(size=1, linetype = 1, color = "black") +
  geom_point(trans, mapping = aes(x = Phyllo_PCoA1_row, y = RGR, color = Species),size = 2.5,shape = 21) +
  theme_bw()+
  scale_color_manual(values = c("#974FA3", "#4CAB3B", "#3C72AF"))+
  labs(x = "PCoA1 of phyllosphere bacteria",y = "Effect of plant growth (cm/day)", title = NULL) +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        legend.position = "none",
        axis.line=element_line(linetype=1,color="black",size=0.5),
        axis.ticks=element_line(color="black",size=0.5,lineend = 2),
        axis.text=element_text(colour='black',size=9),
        axis.title=element_text(colour='black', size=10)) -> Fig3b; Fig3b

###Away_rhizos_PCoA1
get.data = get_model_data(fit, type = "pred", terms = "Away_rhizos_PCoA1[all]")
get.data = as.data.frame(get.data); colnames(get.data)[1] = "Away_rhizos_PCoA1"
##### rotation
get.data["Away_rhizos_PCoA1_row"] <- (pd_attributes_variable2$`scaled:center`["Away_rhizos_PCoA1"] + pd_attributes_variable2$`scaled:scale`["Away_rhizos_PCoA1"]*get.data["Away_rhizos_PCoA1"])
trans["Away_rhizos_PCoA1_row"] <- (pd_attributes_variable2$`scaled:center`["Away_rhizos_PCoA1"] + pd_attributes_variable2$`scaled:scale`["Away_rhizos_PCoA1"]*trans["Away_rhizos_PCoA1"])

ggplot(get.data, aes(x = Away_rhizos_PCoA1_row, y = predicted)) +
  geom_ribbon(aes(ymin=conf.low,ymax=conf.high),fill="#6D6D67", alpha = I(0.3)) +
  geom_line(size=1, linetype = 1, color = "black") +
  geom_point(trans, mapping = aes(x = Away_rhizos_PCoA1_row, y = RGR, color = Species),size = 2.5,shape = 21) +
  theme_bw()+
  scale_color_manual(values = c("#974FA3", "#4CAB3B", "#3C72AF"))+
  labs(x = "PCoA1 of rhizosphere bacteria taxa \n from transplantation site (%)",y = "Effect of plant growth (cm/day)", title = NULL) +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        legend.position = "none",
        axis.line=element_line(linetype=1,color="black",size=0.5),
        axis.ticks=element_line(color="black",size=0.5,lineend = 2),
        axis.text=element_text(colour='black',size=9),
        axis.title=element_text(colour='black', size=10)) -> Fig3c; Fig3c

####end
