getwd()
setwd("D:/桌面/soil_transplantation250508/data/")

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
library(ggplot2)
library(ggh4x)
library(patchwork)
library(emmeans)
library(agricolae)
library(multcomp)
library(sjPlot)
library(lme4)
library(mgcv)
library(tidyverse)
library(dplyr)
library(purrr)
library(tidyr)
library(car)
library(effectsize)
library(metap)

##############Table S4################
group <- read.xlsx("field.xlsx", sheet = "field_PGR",colNames = T, rowNames = T)
mod1 <- aov(GR ~ Site * Seedling_source * Soil_type, data = subset(group, Species == "F. nilgerrensis"))
car::Anova(mod1, type = 2)

mod1 <- as.data.frame(car::Anova(mod1, type =2))
mod1$`Sum Sq` <- round(mod1$`Sum Sq`,3)
mod1$`F value` <- round(mod1$`F value`,3)
mod1$`Pr(>F)` <- round(mod1$`Pr(>F)`,3)
print(mod1)


mod2 <- aov(GR ~ Seedling_source * Soil_type, data = subset(group, Species == "P. asiatica"))
car::Anova(mod2, type = 2)
mod2 <- as.data.frame(car::Anova(mod2, type =2))
mod2$`Sum Sq` <- round(mod2$`Sum Sq`,3)
mod2$`F value` <- round(mod2$`F value`,3)
mod2$`Pr(>F)` <- round(mod2$`Pr(>F)`,3)
print(mod2)

mod3 <- aov(GR ~ Seedling_source * Soil_type, data = subset(group, Species == "T. repens"))
car::Anova(mod3, type = 2)
mod3 <- as.data.frame(car::Anova(mod3, type =2))
mod3$`Sum Sq` <- round(mod3$`Sum Sq`,3)
mod3$`F value` <- round(mod3$`F value`,3)
mod3$`Pr(>F)` <- round(mod3$`Pr(>F)`,3)
print(mod2)

#############Table S5##############
set.seed(123456)

# 1. ANOVA and compare models for main and interactive effects and main effects only
#####F. nilgerrensis
mod1 <- aov(GR ~ Site * Seedling_source * Soil_type, 
            data = subset(group, Species == "F. nilgerrensis"))###models for main and interactive effects
car::Anova(mod1, type = 2)

mod1_mian_only <- aov(GR ~ Site + Seedling_source + Soil_type, 
            data = subset(group, Species == "F. nilgerrensis"))###models for main effects only
car::Anova(mod1_mian_only, type = 2)

# Model comparison
anova(mod1, mod1_mian_only)
AIC(mod1, mod1_mian_only)

#####P. asiatica
# ANOVA and compare models for main and interactive effects and main effects only
mod2 <- aov(GR ~ Site * Seedling_source * Soil_type, 
            data = subset(group, Species == "P. asiatica"))###models for main and interactive effects
car::Anova(mod2, type = 2)

mod2_mian_only <- aov(GR ~ Site + Seedling_source + Soil_type, 
                     data = subset(group, Species == "P. asiatica"))###models for main effects only
car::Anova(mod2_mian_only, type = 2)

# Model comparison
anova(mod2, mod2_mian_only)
AIC(mod2, mod2_mian_only)

#####T. repens
# ANOVA and compare models for main and interactive effects and main effects only
mod3 <- aov(GR ~ Site * Seedling_source * Soil_type, 
            data = subset(group, Species == "T. repens"))###models for main and interactive effects
car::Anova(mod3, type = 2)

mod3_mian_only <- aov(GR ~ Site + Seedling_source + Soil_type, 
                     data = subset(group, Species == "T. repens"))###models for main effects only
car::Anova(mod3_mian_only, type = 2)

# Model comparison
anova(mod3, mod3_mian_only)
AIC(mod3, mod3_mian_only)

#p<0.05, the AIC of the mix_model is smaller, and the logLik is larger 
#Complex models showing support for retaining interaction terms

### 2.Calculate the statistical effectiveness of each item
# Calculateη² and Cohen's f
############F. nilgerrensis
mod1 <- aov(PGR ~ Site * Seedling_source * Soil_type, 
            data = subset(group, Species == "F. nilgerrensis"))

eta_results <- eta_squared(mod1, partial = T)
f_results <- cohens_f(mod1, partial = T)

effect1_size_table <- cbind(eta_results, Cohen_f = f_results$Cohens_f)
effect1_size_table$Eta2_partial <- round(effect1_size_table$Eta2_partial, 3)
effect1_size_table$CI_low <- round(effect1_size_table$CI_low, 3)
effect1_size_table$Cohen_f <- round(effect1_size_table$Cohen_f, 3)
print(effect1_size_table)

############P. asiatica
mod2 <- aov(PGR ~ Site * Seedling_source * Soil_type, 
            data = subset(group, Species == "P. asiatica"))

eta_results <- eta_squared(mod2, partial = T)
f_results <- cohens_f(mod2, partial = T)

effect2_size_table <- cbind(eta_results, Cohen_f = f_results$Cohens_f)
effect2_size_table$Eta2_partial <- round(effect2_size_table$Eta2_partial, 3)
effect2_size_table$CI_low <- round(effect2_size_table$CI_low, 3)
effect2_size_table$Cohen_f <- round(effect2_size_table$Cohen_f, 3)
print(effect2_size_table)

############T. repens
mod3 <- aov(PGR ~ Site * Seedling_source * Soil_type, 
            data = subset(group, Species == "T. repens"))

eta_results <- eta_squared(mod3, partial = T)
f_results <- cohens_f(mod3, partial = T)

effect3_size_table <- cbind(eta_results, Cohen_f = f_results$Cohens_f)
effect3_size_table$Eta2_partial <- round(effect3_size_table$Eta2_partial, 3)
effect3_size_table$CI_low <- round(effect3_size_table$CI_low, 3)
effect3_size_table$Cohen_f <- round(effect3_size_table$Cohen_f, 3)
print(effect3_size_table)
#Eta2: <0.06, Weak effect; 0.06-0.14, Medium effect; ≥0.14, Strong effect 
#Cohens_f: <0.25, Weak effect; 0.25-0.40, Medium effect; ≥ 0.40, Strong effect

###############resampling analysis###############
set.seed(123456)
group <- read.xlsx("field.xlsx", sheet = "field_PGR",colNames = T, rowNames = T)
species_list <- unique(group$Species)  # Get all unique species

# Define parameters
n_boot <- 999  # Number of resamples
# Initialize the independent storage structure
anova_storage <- list()
# Cycling for each species
for (species in species_list) {
  data_subset <- subset(group, Species == species)
  # Get the original sample size for each process combination
  original_sample_sizes <- data_subset %>%
    group_by(Site, Soil_type, Seedling_source) %>%
    summarise(n = n(), .groups = "drop")
  for (i in 1:n_boot) {
    if (i %% 50 == 0) message("Processing iteration ", i, " for Species: ", species)
    # Stratified sampling by raw sample size
    boot_data <- data_subset %>%
      group_by(Site, Soil_type, Seedling_source) %>%
      group_modify(~ {
        # Get the original sample size of the current group
        current_size <- original_sample_sizes %>%
          filter(Site == .y$Site, 
                 Soil_type == .y$Soil_type, 
                 Seedling_source == .y$Seedling_source) %>%
          pull(n)
        # Return sampling by original sample size
        sample_n(.x, size = current_size, replace = TRUE)
      }) %>%
      ungroup()
    tryCatch({
      # Model fitting
      mix_model <- aov(GR ~ Site * Soil_type * Seedling_source, data = boot_data)
      # Store ANOVA results
      anova_res <- car::Anova(mix_model, type = 2) %>%
        as.data.frame() %>%
        tibble::rownames_to_column("Effect") %>%
        mutate(
          BootstrapID = i,
          Species = species
        )
      
      anova_storage[[length(anova_storage) + 1]] <- anova_res
    }, error = function(e) {
      message("Error in iteration ", i, " for Species ", species, ": ", e$message)
    })
  }
}

anova_results_df <- bind_rows(anova_storage)
anova_results_cleaned <- anova_results_df %>% filter(!is.na(`Pr(>F)`))

####Calculation of significant percentages
# p-value adjust
anova_results_cleaned <- anova_results_df %>%
  filter(!is.na(`Pr(>F)`))

final_results1 <- anova_results_cleaned %>%
  group_by(Species) %>%
  mutate(Adjusted_p = p.adjust(`Pr(>F)`, method = "BH"))

p_value_summary <- final_results1 %>%
  group_by(Species, Effect) %>%
  summarise(
    Count_less_than_0.05 = sum(Adjusted_p < 0.05, na.rm = TRUE),
    Count_greater_equal_0.05 = sum(Adjusted_p >= 0.05, na.rm = TRUE),
    Total_Count = n(),
    Proportion_less_than_0.05 = Count_less_than_0.05 / Total_Count,
    Proportion_greater_equal_0.05 = Count_greater_equal_0.05 / Total_Count,
    .groups = 'drop'
  )

p_value_long <- p_value_summary %>%
  pivot_longer(cols = starts_with("Proportion"),
               names_to = "Significance",
               values_to = "Proportion") %>%
  mutate(Significance = ifelse(Significance == "Proportion_less_than_0.05", "< 0.05", ">= 0.05"))

p_value_long <- p_value_long %>%
  mutate(Effect = as.character(Effect),
         Effect = ifelse(Effect == "Site", "Site", Effect),
         Effect = ifelse(Effect == "Seedling_source", "Seedling_source", Effect),
         Effect = ifelse(Effect == "Soil_type", "Soil_type", Effect),
         Effect = ifelse(Effect == "Site:Seedling_source", "Site x Seedling_source", Effect),
         Effect = ifelse(Effect == "Site:Soil_type", "Site x Soil_type", Effect),
         Effect = ifelse(Effect == "Soil_type:Seedling_source", "Soil_type x Seedling_source", Effect),
         Effect = ifelse(Effect == "Site:Soil_type:Seedling_source", "Site x Soil_type x Seedling_source", Effect) %>% 
           factor(levels = c("Site", "Seedling_source", "Soil_type", 
                             "Site x Seedling_source", "Site x Soil_type", 
                             "Soil_type x Seedling_source", 
                             "Site x Soil_type x Seedling_source")))

p_value_long$Proportion <- round(p_value_long$Proportion*100,1)
p_value_long$Significance <- factor(p_value_long$Significance, levels = c(">= 0.05", "< 0.05"))

ggplot(p_value_long, aes(x = Effect, y = Proportion, fill = Significance)) +
  geom_bar(stat = "identity", position = "stack") +  
  scale_fill_manual(values = c("#6888F5", "#D77071")) + 
  facet_wrap(~ Species) + 
  labs(
    y = "Proportion of total models (%)"
  ) +
  theme_bw() +
  facet_grid(. ~ Species) +
  geom_hline(aes(yintercept = 50), linetype = "dashed", size = 1) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "right",
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.line = element_line(linetype = 1, color = "black", size = 0.5),
        axis.ticks = element_line(color = "black", size = 0.5, lineend = "round"),
        axis.text = element_text(colour = 'black', size = 9),
        axis.title = element_text(colour = 'black', size = 10)) -> Fig_S3; Fig_S3

###############Fig 1a##############
group <- read.xlsx("field.xlsx", sheet = "raw_plant_data",colNames = T, rowNames = T)

df2 <- data_summary(group, varname = "GR", groupnames = c("Site", "Seedling_source", "Soil_type", "Species", "Type"))
#write.xlsx(df2, file = "field_PGR_summary1.xlsx", sheetName = "Sheet1", row.names = TRUE) 

df2$Seedling_source <- factor(df2$Seedling_source, levels = c("Low", "Mid", "High"))
df2$Type <- factor(paste0(df2$Site, "_", df2$Seedling_source), 
                   levels = c("Low_Low","Mid_Mid","High_High","Low_Mid", "Low_High", "Mid_Low", 
                              "Mid_High", "High_Low", "High_Mid"))

ggplot(data = df2) +
  geom_errorbar(aes(x = Type, y = GR, ymin = GR - se, 
                    ymax = GR + se, color = Seedling_source, 
                    group = interaction(Soil_type, Seedling_source)),
                position = position_dodge(width = 0.85), width = 0, 
                size = 0.5, alpha = 1) +
  geom_bar(aes(x = Type, y = GR, 
               fill = interaction(Soil_type, Seedling_source), 
               group = interaction(Soil_type, Seedling_source), 
               color = Seedling_source),
           position = position_dodge(width = 0.85), stat = "identity", 
           size = 0.7, width = 0.7, alpha = 0.85) +
  scale_fill_manual(values = c("#2C91E0", "NA", "#F0A73A", "NA", "#3ABF99", "NA")) +
  scale_color_manual(values = c("#2C91E0", "#F0A73A", "#3ABF99")) + 
  labs(y = "Plant growth rate (cm/day)", title = NULL) +
  theme_bw() +
  scale_y_continuous(expand=c(0,0))+
  facet_wrap(~ Species, scales = "free_y") + 
  geom_hline(aes(yintercept = 0), linetype = "dashed", size = 0.5) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.line = element_line(linetype = 1, color = "black", size = 0.5),
        axis.ticks = element_line(color = "black", size = 0.5, lineend = "round"),
        axis.text = element_text(colour = 'black', size = 9),
        axis.title = element_text(colour = 'black', size = 10)) -> Fig_1a; Fig_1a

###############Fig 1b##############
#group <- read.xlsx("field.xlsx", sheet = "raw_plant_data",colNames = T, rowNames = T)
# Define a function to process data and plots
#group$GR <- log(group$GR)
#df2 <- data_summary(group, varname = "GR", groupnames = c("Site", "Seedling_source", "Soil_type", "Type","Species"))
#write.xlsx(df2, file = "field_PGR_summary1.xlsx", sheetName = "Sheet1", rowNames = TRUE) 

###Relationship between relative plant growth rate and transplantation distance
group <- read.xlsx("field.xlsx", sheet = "effect1",colNames = T, rowNames = T)

ggplot(data = group, aes(Elevational, Effect1, group = Species)) +
  labs(x = "Elevational",y = "ln(RG in home soil/RG in CK soil)", title = NULL) +
  scale_fill_manual(values = c("#974FA3", "#4CAB3B", "#B7352A")) +
  scale_color_manual(values = c("#974FA3", "#4CAB3B", "#B7352A")) +
  geom_smooth(data = group, aes(Elevational, Effect1, color = Species),
              method = "loess", se = FALSE, span = 0.4, position = position_dodge(width = 100)) +
  geom_errorbar(data = group, aes(x = Elevational, y = Effect1, ymin = Effect1 - CI, 
                                  ymax = Effect1 + CI),color = "#999595", 
                width = 0, size = 0.7, alpha = 1, position = position_dodge(width = 100)) +
  geom_errorbar(data = group, aes(x = Elevational, y = Effect1, ymin = Effect1 - se2, 
                                  ymax = Effect1 + se2),color = "black", 
                width = 0, size = 1.2, alpha = 1, position = position_dodge(width = 100)) +
  geom_point(data = group, 
             mapping = aes(x = Elevational, y = Effect1, fill = Species), 
             shape = 21, size = 3, position = position_dodge(width = 100)) +
  theme(panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid = element_blank(),
        axis.title.x = element_text(colour = 'black', size = 12),
        axis.title.y = element_text(colour = 'black', size = 12),
        axis.text = element_text(colour = 'black', size = 10),
        legend.title = element_text(size = 10, face = "bold"),
        legend.text = element_text(size = 10),
        legend.position = "top") -> Fig1b; Fig1b

###############Fig 1c##############
group <- read.xlsx("field.xlsx", sheet = "effect2",colNames = T, rowNames = T)

ggplot(data = group, aes(Elevational, Effect2, group = Species)) +
  labs(x = "Elevational",y = "ln(RG in home soil/RG in away soil)", title = NULL) +
  scale_fill_manual(values = c("#974FA3", "#4CAB3B", "#B7352A")) +
  scale_color_manual(values = c("#974FA3", "#4CAB3B", "#B7352A")) +
  geom_smooth(data = group, aes(Elevational, Effect2, color = Species),
              method = "loess", se = FALSE, span = 0.4, position = position_dodge(width = 100)) +
  geom_errorbar(data = group, aes(x = Elevational, y = Effect2, ymin = Effect2 - CI, 
                                  ymax = Effect2 + CI),color = "#999595", 
                width = 0, size = 0.7, alpha = 1, position = position_dodge(width = 100)) +
  geom_errorbar(data = group, aes(x = Elevational, y = Effect2, ymin = Effect2 - se2, 
                                  ymax = Effect2 + se2),color = "black", 
                width = 0, size = 1.2, alpha = 1, position = position_dodge(width = 100)) +
  geom_point(data = group, 
             mapping = aes(x = Elevational, y = Effect2, fill = Species), 
             shape = 21, size = 3, position = position_dodge(width = 100)) +
  theme(panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid = element_blank(),
        axis.title.x = element_text(colour = 'black', size = 12),
        axis.title.y = element_text(colour = 'black', size = 12),
        axis.text = element_text(colour = 'black', size = 10),
        legend.title = element_text(size = 10, face = "bold"),
        legend.text = element_text(size = 10),
        legend.position = "top") -> Fig1c; Fig1c

######relationship with plant height and biomass in greenhouse experiment#####
group <- read.xlsx("greenhouse.xlsx", sheet = "greenhouse_PGR",colNames = T, rowNames = T)
group$above_biomass <- log(group$above_biomass)
group$height_2 <- log(group$height_2)
summary(lm(above_biomass ~ height_2, data = subset(group, species == "P. asiatica")))
summary(lm(above_biomass ~ height_2, data = subset(group, species == "T. repens")))
####P. asiatica
ggplot() +
  geom_smooth(subset(group, species == "P. asiatica"), mapping = aes(x = height_2, y = above_biomass), color = "#4CAB3B",method = lm,size = 1) +
  geom_point(subset(group, species == "P. asiatica"), mapping = aes(x = height_2, y = above_biomass), color = "#4CAB3B",size = 2,shape = 21) +
  theme_bw()+
  labs(x = "Plant height (cm)",y = "Aboveground biomass (g)", title = NULL) +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        legend.position = "none",
        axis.line=element_line(linetype=1,color="black",size=0.5),
        axis.ticks=element_line(color="black",size=0.5,lineend = 2),
        axis.text=element_text(colour='black',size=9),
        axis.title=element_text(colour='black', size=10))
####T. repens
ggplot() +
  geom_smooth(subset(group, species == "T. repens"), mapping = aes(x = height_2, y = above_biomass), color = "#B7352A",method = lm,size = 1) +
  geom_point(subset(group, species == "T. repens"), mapping = aes(x = height_2, y = above_biomass), color = "#B7352A",size = 2,shape = 21) +
  theme_bw()+
  labs(x = "Plant height (cm)",y = "Aboveground biomass (g)", title = NULL) +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        legend.position = "none",
        axis.line=element_line(linetype=1,color="black",size=0.5),
        axis.ticks=element_line(color="black",size=0.5,lineend = 2),
        axis.text=element_text(colour='black',size=9),
        axis.title=element_text(colour='black', size=10))

####end
