getwd()
setwd("../soil_transplantation/data/")
##############################PACKAGES##############################
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
library(MuMIn)
library(car)
library(reshape)
library(Hmisc)

#####Multiple regression of LMM########
options(na.action ="na.fail")
trans = read.xlsx("field.xlsx", sheet = "Multiple_regression", rowNames = T, colNames = T)
colnames(trans)

###########correlation analysis
#trans[,c(7:27)] <- scale(trans[,c(7:27)])
soilData <- trans[,c(5:14)] 

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
  coord_fixed();p_heatmap###we selected parameter:Phyllo_PCoA1/Home_rhizos_PCoA1/Away_rhizos_PCoA1/Away_fungi_PCoA1/away_fungi/away_rhizos
away_rhizos

######Filter the optimal model
trans = read.xlsx("field.xlsx", sheet = "Multiple_regression", rowNames = T, colNames = T)
colnames(trans)# View column names to verify the data
trans[,c(6,7,9,10,12,14)]  <- scale(trans[,c(6,7,9,10,12,14)])# Scale specific columns (standardize to mean = 0, sd = 1)

# Fit a full linear mixed-effects model with PGR as the response variable
# Fixed effects: Phyllo_PCoA1, Home_rhizos_PCoA1, Away_rhizos_PCoA1, Away_fungi_PCoA1, away_fungi, away_rhizos
# Random effect: Species (random intercept)
# REML = FALSE is used for model comparison and selection
fit<-lmer(PGR ~ Phyllo_PCoA1+Home_rhizos_PCoA1+Away_rhizos_PCoA1+Away_fungi_PCoA1+away_fungi+
            away_rhizos+(1|Species),data=trans,REML=FALSE)###full model
vif(fit)# Calculate Variance Inflation Factor (VIF) to check for multicollinearity among predictors
# Perform model selection using the dredge function to explore all possible subsets of the full model
fit_select <- dredge(fit,rank = "AIC")
head(fit_select)

get.models(fit_select,1)[[1]]

# Fit the final model based on the selected predictors
# Fixed effects: Home_rhizos_PCoA1, Away_rhizos_PCoA1, Phyllo_PCoA1
# Random effect: Species (random intercept)
fit<-lmer(PGR ~ Home_rhizos_PCoA1+Away_rhizos_PCoA1+
            Phyllo_PCoA1+(1|Species),data=trans)###final_model
#AIC(fit)# Calculate AIC for the final model
summary(fit)# Calculate AIC for the final model
r.squaredGLMM(fit)# Calculate R-squared values (marginal and conditional) for the final model
glmm.hp::glmm.hp(fit,type = "R2")# Perform hierarchical partitioning of variance for the final model

anova_results <- car::Anova(fit,type=2)# Perform Type II ANOVA to test the significance of fixed effects
p_values <- anova_results$`Pr(>Chisq)`# Extract p-values from the ANOVA results
adjusted_p_values <- p.adjust(p_values, method = "BH")# Adjust p-values using the Benjamini-Hochberg (BH) method to control for multiple testing

#########Fig 3###########
trans = read.xlsx("field.xlsx", sheet = "Multiple_regression", rowNames = T, colNames = T)
pd_attributes_variable1 <- attributes(scale(trans[c("Home_rhizos_PCoA1")]))
pd_attributes_variable2 <- attributes(scale(trans[c("Away_rhizos_PCoA1")]))
pd_attributes_variable3 <- attributes(scale(trans[c("Phyllo_PCoA1")]))
trans[,c(6,7,9,10,12,14)]  <- scale(trans[,c(6,7,9,10,12,14)])
fit<-lmer(PGR ~ Home_rhizos_PCoA1+Away_rhizos_PCoA1+
            Phyllo_PCoA1+(1|Species),data=trans)###final_model

###Home_rhizos_PCoA1
get.data = get_model_data(fit, type = "pred", terms = "Home_rhizos_PCoA1[all]")
get.data = as.data.frame(get.data); colnames(get.data)[1] = "Home_rhizos_PCoA1"
##### rotation the scale data
get.data["Home_rhizos_PCoA1_row"] <- (pd_attributes_variable1$`scaled:center`["Home_rhizos_PCoA1"] + pd_attributes_variable1$`scaled:scale`["Home_rhizos_PCoA1"]*get.data["Home_rhizos_PCoA1"])
trans["Home_rhizos_PCoA1_row"] <- (pd_attributes_variable1$`scaled:center`["Home_rhizos_PCoA1"] + pd_attributes_variable1$`scaled:scale`["Home_rhizos_PCoA1"]*trans["Home_rhizos_PCoA1"])

ggplot(get.data, aes(x = Home_rhizos_PCoA1_row, y = predicted)) +
  geom_ribbon(aes(ymin=conf.low,ymax=conf.high),fill="#6D6D67", alpha = I(0.3)) +
  geom_line(size=1, linetype = 1, color = "black") +
  geom_point(trans, mapping = aes(x = Home_rhizos_PCoA1_row, y = PGR, color = Species),size = 2.5,shape = 21) +
  theme_bw()+
  scale_color_manual(values = c("#974FA3", "#4CAB3B", "#B7352A"))+
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
##### rotation the scale data
get.data["Phyllo_PCoA1_row"] <- (pd_attributes_variable3$`scaled:center`["Phyllo_PCoA1"] + pd_attributes_variable3$`scaled:scale`["Phyllo_PCoA1"]*get.data["Phyllo_PCoA1"])
trans["Phyllo_PCoA1_row"] <- (pd_attributes_variable3$`scaled:center`["Phyllo_PCoA1"] + pd_attributes_variable3$`scaled:scale`["Phyllo_PCoA1"]*trans["Phyllo_PCoA1"])

ggplot(get.data, aes(x = Phyllo_PCoA1_row, y = predicted)) +
  geom_ribbon(aes(ymin=conf.low,ymax=conf.high),fill="#6D6D67", alpha = I(0.3)) +
  geom_line(size=1, linetype = 1, color = "black") +
  geom_point(trans, mapping = aes(x = Phyllo_PCoA1_row, y = PGR, color = Species),size = 2.5,shape = 21) +
  theme_bw()+
  scale_color_manual(values = c("#974FA3", "#4CAB3B", "#B7352A"))+
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
##### rotation the scale data
get.data["Away_rhizos_PCoA1_row"] <- (pd_attributes_variable2$`scaled:center`["Away_rhizos_PCoA1"] + pd_attributes_variable2$`scaled:scale`["Away_rhizos_PCoA1"]*get.data["Away_rhizos_PCoA1"])
trans["Away_rhizos_PCoA1_row"] <- (pd_attributes_variable2$`scaled:center`["Away_rhizos_PCoA1"] + pd_attributes_variable2$`scaled:scale`["Away_rhizos_PCoA1"]*trans["Away_rhizos_PCoA1"])

ggplot(get.data, aes(x = Away_rhizos_PCoA1_row, y = predicted)) +
  geom_ribbon(aes(ymin=conf.low,ymax=conf.high),fill="#6D6D67", alpha = I(0.3)) +
  geom_line(size=1, linetype = 1, color = "black") +
  geom_point(trans, mapping = aes(x = Away_rhizos_PCoA1_row, y = PGR, color = Species),size = 2.5,shape = 21) +
  theme_bw()+
  scale_color_manual(values = c("#974FA3", "#4CAB3B", "#B7352A"))+
  labs(x = "PCoA1 of rhizosphere bacteria taxa \n from transplantation site",y = "Effect of plant growth (cm/day)", title = NULL) +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        legend.position = "none",
        axis.line=element_line(linetype=1,color="black",size=0.5),
        axis.ticks=element_line(color="black",size=0.5,lineend = 2),
        axis.text=element_text(colour='black',size=9),
        axis.title=element_text(colour='black', size=10)) -> Fig3c; Fig3c
Fig3a +Fig3b +Fig3c

######end
trans = read.xlsx("D:/桌面/soil_transplantation250508/data/greenhouse.xlsx", sheet = "greenhouse_PGR", rowNames = T, colNames = T)
trans$PGR <- log(trans$PGR)
df2 <- data_summary(trans, varname = "PGR", groupnames = c("Soil_source", "Seedling_source", "species"))
write.csv(df2, file = "D:/桌面/greenhouse_PGR.csv", row.names = FALSE)
