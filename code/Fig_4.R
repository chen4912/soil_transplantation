getwd()
setwd("D:/桌面/soil_transplantation250508/data/")
##############################PACKAGES##############################
library(openxlsx)
library(deming)
library(ggplot2)
library(dplyr)

##############tables S9############
########   PGR
####T. repens
group <- read.xlsx("greenhouse.xlsx", sheet = "greenhouse_PGR",colNames = T, rowNames = T)
group$plant <- factor(group$Seedling_source, levels = c("Low","Mid","High"))
group$Soil_source <- factor(group$Soil_source, levels = c("Low","Mid","High"))
greenhouse_mod1 <- aov(PGR ~ Seedling_source * Soil_source, data = subset(group, species == "T. repens"))

greenhouse_mod1 <- as.data.frame(car::Anova(greenhouse_mod1, type =2))
greenhouse_mod1$`Sum Sq` <- round(greenhouse_mod1$`Sum Sq`,3)
greenhouse_mod1$`F value` <- round(greenhouse_mod1$`F value`,3)
greenhouse_mod1$`Pr(>F)` <- round(greenhouse_mod1$`Pr(>F)`,3)
print(greenhouse_mod1)

####P. asiatica
greenhouse_mod2 <- aov(PGR ~ Seedling_source * Soil_source, data = subset(group,species == "P. asiatica"))

greenhouse_mod2 <- as.data.frame(car::Anova(greenhouse_mod2, type =2))
greenhouse_mod2$`Sum Sq` <- round(greenhouse_mod2$`Sum Sq`,3)
greenhouse_mod2$`F value` <- round(greenhouse_mod2$`F value`,3)
greenhouse_mod2$`Pr(>F)` <- round(greenhouse_mod2$`Pr(>F)`,3)
print(greenhouse_mod2)

########   ABOVE BIOMASS
####T. repens
greenhouse_mod1 <- aov(above_biomass ~ Seedling_source * Soil_source, data = subset(group,species == "T. repens"))

greenhouse_mod1 <- as.data.frame(car::Anova(greenhouse_mod1, type =2))
greenhouse_mod1$`Sum Sq` <- round(greenhouse_mod1$`Sum Sq`,3)
greenhouse_mod1$`F value` <- round(greenhouse_mod1$`F value`,3)
greenhouse_mod1$`Pr(>F)` <- round(greenhouse_mod1$`Pr(>F)`,3)
print(greenhouse_mod1)

####P. asiatica
greenhouse_mod2 <- aov(above_biomass ~ Seedling_source * Soil_source, data = subset(group,species == "P. asiatica"))

greenhouse_mod2 <- as.data.frame(car::Anova(greenhouse_mod2, type =2))
greenhouse_mod2$`Sum Sq` <- round(greenhouse_mod2$`Sum Sq`,3)
greenhouse_mod2$`F value` <- round(greenhouse_mod2$`F value`,3)
greenhouse_mod2$`Pr(>F)` <- round(greenhouse_mod2$`Pr(>F)`,3)
print(greenhouse_mod2)

###############Fig S4############
####   PGR 
greenhouse_mod1 <- aov(PGR ~ Seedling_source * Soil_source, data = subset(group, species == "T. repens"))
emm_results <- emmeans(greenhouse_mod1, ~ Seedling_source * Soil_source)
df2 <- summary(emm_results)

ggplot(data = df2) +
  geom_errorbar(mapping = aes(x = Seedling_source, y = emmean, ymin = emmean - SE, ymax = emmean + SE, color = Soil_source),
                position = position_dodge(width = 0.7), width = 0, size = 0.5, alpha = 1) +
  geom_bar(mapping = aes(x = Seedling_source, y = emmean, fill = Soil_source),
           position = position_dodge(width = 0.7),width = 0.6, stat = "identity", size = 3, alpha = 0.8) +
  scale_fill_manual(values = c("#974FA3","#4CAB3B","#B7352A")) +
  scale_color_manual(values = c("#974FA3","#4CAB3B","#B7352A")) +
  labs(y = "Plant growth rate (cm/day)", title = NULL) +
  theme_bw() +
  scale_y_continuous(expand=c(0,0))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.line = element_line(linetype = 1, color = "black", size = 0.5),
        axis.ticks = element_line(color = "black", size = 0.5, lineend = 2),
        axis.text = element_text(colour = 'black', size = 9),
        axis.title = element_text(colour = 'black', size = 10)) -> FigS1a;FigS1a

model_means <- emmeans(greenhouse_mod1, specs = ~ Seedling_source * Soil_source)
pairwise_comparisons <- pairs(model_means, adjust = "none")
model_means_cld <- cld(object = model_means,
                       Letters = letters,
                       alpha = 0.05,
                       decreasing = T)
model_means_cld

####P. asiatica
greenhouse_mod2 <- aov(PGR ~ Seedling_source * Soil_source, data = subset(group,species == "P. asiatica"))
emm_results <- emmeans(greenhouse_mod2, ~ Seedling_source * Soil_source)
df2 <- summary(emm_results)

ggplot(data = df2) +
  geom_errorbar(mapping = aes(x = Seedling_source, y = emmean, ymin = emmean - SE, ymax = emmean + SE, color = Soil_source),
                position = position_dodge(width = 0.7), width = 0, size = 0.5, alpha = 1) +
  geom_bar(mapping = aes(x = Seedling_source, y = emmean, fill = Soil_source),
           position = position_dodge(width = 0.7),width = 0.6, stat = "identity", size = 3, alpha = 0.8) +
  scale_fill_manual(values = c("#974FA3","#4CAB3B","#B7352A")) +
  scale_color_manual(values = c("#974FA3","#4CAB3B","#B7352A")) +
  labs(y = "Plant growth rate (cm/day)", title = NULL) +
  theme_bw() +
  scale_y_continuous(expand=c(0,0))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.line = element_line(linetype = 1, color = "black", size = 0.5),
        axis.ticks = element_line(color = "black", size = 0.5, lineend = 2),
        axis.text = element_text(colour = 'black', size = 9),
        axis.title = element_text(colour = 'black', size = 10)) -> FigS1b;FigS1b

model_means <- emmeans(greenhouse_mod2, specs = ~ Seedling_source * Soil_source)
pairwise_comparisons <- pairs(model_means, adjust = "none")
model_means_cld <- cld(object = model_means,
                       Letters = letters,
                       alpha = 0.05,
                       decreasing = T)
model_means_cld

####   biomass   ##
####T. repens
greenhouse_mod1 <- aov(above_biomass ~ Seedling_source * Soil_source, data = subset(group,species == "T. repens"))
emm_results <- emmeans(greenhouse_mod1, ~ Seedling_source * Soil_source)
df2 <- summary(emm_results)

ggplot(data = df2) +
  geom_errorbar(mapping = aes(x = Seedling_source, y = emmean, ymin = emmean - SE, ymax = emmean + SE, color = Soil_source),
                position = position_dodge(width = 0.7), width = 0, size = 0.5, alpha = 1) +
  geom_bar(mapping = aes(x = Seedling_source, y = emmean, fill = Soil_source),
           position = position_dodge(width = 0.7),width = 0.6, stat = "identity", size = 3, alpha = 0.8) +
  scale_fill_manual(values = c("#974FA3","#4CAB3B","#B7352A")) +
  scale_color_manual(values = c("#974FA3","#4CAB3B","#B7352A")) +
  labs(y = "Plant growth rate (cm/day)", title = NULL) +
  theme_bw() +
  scale_y_continuous(expand=c(0,0))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.line = element_line(linetype = 1, color = "black", size = 0.5),
        axis.ticks = element_line(color = "black", size = 0.5, lineend = 2),
        axis.text = element_text(colour = 'black', size = 9),
        axis.title = element_text(colour = 'black', size = 10)) -> FigS1c;FigS1c

model_means <- emmeans(greenhouse_mod1, specs = ~ Seedling_source * Soil_source)
pairwise_comparisons <- pairs(model_means, adjust = "none")
model_means_cld <- cld(object = model_means,
                       Letters = letters,
                       alpha = 0.05,
                       decreasing = T)
model_means_cld

####P. asiatica
greenhouse_mod2 <- aov(above_biomass ~ Seedling_source * Soil_source, data = subset(group,species == "P. asiatica"))
emm_results <- emmeans(greenhouse_mod2, ~ Seedling_source * Soil_source)
df2 <- summary(emm_results)

ggplot(data = df2) +
  geom_errorbar(mapping = aes(x = Seedling_source, y = emmean, ymin = emmean - SE, ymax = emmean + SE, color = Soil_source),
                position = position_dodge(width = 0.7), width = 0, size = 0.5, alpha = 1) +
  geom_bar(mapping = aes(x = Seedling_source, y = emmean, fill = Soil_source),
           position = position_dodge(width = 0.7),width = 0.6, stat = "identity", size = 3, alpha = 0.8) +
  scale_fill_manual(values = c("#974FA3","#4CAB3B","#B7352A")) +
  scale_color_manual(values = c("#974FA3","#4CAB3B","#B7352A")) +
  labs(y = "Plant growth rate (cm/day)", title = NULL) +
  theme_bw() +
  scale_y_continuous(expand=c(0,0))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.line = element_line(linetype = 1, color = "black", size = 0.5),
        axis.ticks = element_line(color = "black", size = 0.5, lineend = 2),
        axis.text = element_text(colour = 'black', size = 9),
        axis.title = element_text(colour = 'black', size = 10)) -> FigS1d;FigS1d

model_means <- emmeans(greenhouse_mod2, specs = ~ Seedling_source * Soil_source)
pairwise_comparisons <- pairs(model_means, adjust = "none")
model_means_cld <- cld(object = model_means,
                       Letters = letters,
                       alpha = 0.05,
                       decreasing = T)
model_means_cld

(FigS1a + FigS1b) / (FigS1c + FigS1d)

#############Fig 4#############
group <- read.xlsx("greenhouse.xlsx", sheet = "deming",colNames = T, rowNames = T)
Pasi <- subset(group, species == "P. asiatica")
cor.test(Pasi$field,Pasi$house)
###deming regression analysis of P. asiatica
fit <- deming(field ~ house, data=Pasi, xstd=se1, ystd=se)
print(fit)

ggplot(Pasi, aes(x = house, y = field)) +
  geom_errorbar(aes(ymin = field - se, ymax = field + se), 
                width = 0, position = position_dodge(0.5), color= "#4CAB3B") +
  geom_errorbarh(aes(xmin = house - se1, xmax = house + se1, y = field), 
                 height = 0, position = position_dodge(0.5), color= "#4CAB3B") +
  geom_point(size = 2.5, position = position_dodge(0.5), color= "#4CAB3B",shape = 16) +
  labs(y = "Filed", x = "Greenhouse", title = NULL) +
  geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "red") +
  geom_abline(slope = -14.38, intercept = 0.07, linetype = "dashed", color = "black") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.line = element_line(linetype = 1, color = "black", size = 0.5),
        axis.ticks = element_line(color = "black", size = 0.5, lineend = 2),
        axis.text = element_text(colour = 'black', size = 9),
        axis.title = element_text(colour = 'black', size = 10)) -> Fig4a;Fig4a

###deming regression analysis of T. repens
Trep <- subset(group, species == "T. repens")
cor.test(Trep$field,Trep$house)

fit <- deming(field ~ house, data=Trep, xstd=se1, ystd=se)
print(fit)

ggplot(Trep, aes(x = house, y = field, color = species)) +
  geom_errorbar(aes(ymin = field - se, ymax = field + se), 
                width = 0, position = position_dodge(0.5), color= "#B7352A") +
  geom_errorbarh(aes(xmin = house - se1, xmax = house + se1, y = field), 
                 height = 0, position = position_dodge(0.5), color= "#B7352A") +
  geom_point(size = 2.5, position = position_dodge(0.5), color= "#B7352A",shape = 16) +
  labs(y = "Filed", x = "Greenhouse", title = NULL) +
  geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "red") +
  geom_abline(slope = -2.00, intercept = -0.195, linetype = "dashed", color = "black") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.line = element_line(linetype = 1, color = "black", size = 0.5),
        axis.ticks = element_line(color = "black", size = 0.5, lineend = 2),
        axis.text = element_text(colour = 'black', size = 9),
        axis.title = element_text(colour = 'black', size = 10)) -> Fig4b;Fig4b

Fig4a+Fig4b
####end

