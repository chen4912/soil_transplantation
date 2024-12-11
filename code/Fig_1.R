###################Fig 1###########
####plant growth in field###
library(openxlsx)
library(ggplot2)
library(ggh4x)
library(patchwork)
library(emmeans)
library(agricolae)
library(multcomp)
library(sjPlot)
group <- read.xlsx("field.xlsx", sheet = "trans_PGR",colNames = T, rowNames = T)
#group <- read.xlsx("trans.xlsx", sheet = "trans_CK",colNames = T, rowNames = T)
group$Site <- factor(group$Site, levels = c("1590","1946","2641"))
group$Popu <- factor(group$Popu, levels = c("1590","1946","2641"))

group$Type <- paste0(group$Popu,"_",group$Site)
df2 <- data_summary(group, varname = "effect1",groupnames = c("Popu","Type","Soil","Species","Soilsource"))
df2$Popu <- factor(df2$Popu, levels = c("1590","1946","2641"))
df2$Type <- factor(df2$Type, levels = c("1590_1590","1590_1946","1590_2641","1946_1946","1946_1590","1946_2641","2641_2641","2641_1590","2641_1946"))
df2$CI <- 1.96*df2$se

ggplot(data = df2) +
  geom_errorbar(mapping = aes(x = Type, y = effect1, ymin = effect1 - CI, ymax = effect1 + CI, 
                              color = Species, group = interaction(Soilsource, Species)),
                position = position_dodge(width = 1), width = 0, size = 0.5, alpha = 1) +
  geom_bar(mapping = aes(x = Type, y = effect1, fill = interaction(Soilsource, Species), group = interaction(Soilsource, Species), color = Species),
           position = position_dodge(width = 1), stat = "identity", size = 0.7, width = 0.7, alpha = 0.8) +
  scale_fill_manual(values = c("#974FA3","NA", "#4CAB3B","NA", "#3C72AF","NA")) +
  scale_color_manual(values = c("#974FA3", "#4CAB3B", "#3C72AF")) + 
  labs(y = "Effect size", title = NULL) +
  theme_bw() +
  geom_hline(aes(yintercept = 0), linetype = "dashed", size = 0.5) +
  facet_wrap(~ Popu, scales = "free") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "right",
        axis.line = element_line(linetype = 1, color = "black", size = 0.5),
        axis.ticks = element_line(color = "black", size = 0.5, lineend = "round"),
        axis.text = element_text(colour = 'black', size = 9),
        axis.title = element_text(colour = 'black', size = 10)) -> Fig1a;Fig1a


###Relationship between relative plant growth rate and transplantation distance
group1 <- subset(group, Soilsource == "home")

############Elevational~~effect(home/CK)
fit<-lmer(effect1 ~ Elevational+(1|Species),data=group1)###final_model
summary(fit)
MuMIn::r.squaredGLMM(fit)

get.data = get_model_data(fit, type = "pred", terms = "Elevational[all]")
get.data = as.data.frame(get.data); colnames(get.data)[1] = "Elevational"

ggplot(get.data, aes(x = Elevational, y = predicted)) +
  geom_ribbon(aes(ymin=conf.low,ymax=conf.high),fill="#6D6D67", alpha = I(0.2)) +
  geom_line(size=1, linetype = 1, color = "black") +
  geom_point(group1, mapping = aes(x = Elevational, y = effect1, color = Species),size = 2,shape = 21) +
  theme_bw()+
  scale_color_manual(values = c("#974FA3", "#4CAB3B", "#3C72AF"))+
  labs(x = "Elevational",y = "Effect(Home/CK)", title = NULL) +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        legend.position = "none",
        axis.line=element_line(linetype=1,color="black",size=0.5),
        axis.ticks=element_line(color="black",size=0.5,lineend = 2),
        axis.text=element_text(colour='black',size=9),
        axis.title=element_text(colour='black', size=10)) -> Fig1b; Fig1b

####Elevational~~effect(home/away)
fit<-lmer(effect2 ~ Elevational+(1|Species),data=group1)###final_model
summary(fit)
MuMIn::r.squaredGLMM(fit)

get.data = get_model_data(fit, type = "pred", terms = "Elevational[all]")
get.data = as.data.frame(get.data); colnames(get.data)[1] = "Elevational"

ggplot(get.data, aes(x = Elevational, y = predicted)) +
  geom_ribbon(aes(ymin=conf.low,ymax=conf.high),fill="#6D6D67", alpha = I(0.2)) +
  geom_line(size=1, linetype = 1, color = "#6D6D67") +
  geom_point(group1, mapping = aes(x = Elevational, y = effect2, color = Species),size = 2,shape = 21) +
  theme_bw()+
  scale_color_manual(values = c("#974FA3", "#4CAB3B", "#3C72AF"))+
  labs(x = "Elevational",y = "Effect(Home/away)", title = NULL) +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        legend.position = "none",
        axis.line=element_line(linetype=1,color="black",size=0.5),
        axis.ticks=element_line(color="black",size=0.5,lineend = 2),
        axis.text=element_text(colour='black',size=9),
        axis.title=element_text(colour='black', size=10)) -> Fig1c; Fig1c

###Multiple comparison
group$Type <- paste0(group$Popu,"_",group$Site)
type <- unique(group$Type)
species <- unique(group$Species)
result <- data.frame()
#i="2641_1590"
#j="F. nilgerrensis"

for (i in type) {
  for (j in species) {
    low <- subset(group, Type == i & Species == j)
    low$group <- paste0(low$Site,"_",low$Soil)
    model <- aov(effect1 ~ Soilsource, low)
    summary(model)
    model_means <- emmeans(object = model,
                           specs = "Soilsource")
    model_means_cld <- cld(object = model_means,
                           adjust = "none",
                           Letters = letters,
                           alpha = 0.05,
                           decreasing = T)
    model <- data.frame(model_means_cld)
    model$Type <- i
    model$Species <- j
    result <- rbind(result, model)
  }
}
result1

####t test compared to 0
group$Type <- paste0(group$Popu,"_",group$Site,"_",group$Soilsource)
type <- unique(group$Type)
species <- unique(group$Species)
result <- data.frame()
#i="2641_1590"
#j="F. nilgerrensis"

for (i in type) {
  for (j in species) {
    low <- subset(group, Type == i & Species == j)
    t_test_result <- t.test(low$effect1, mu = 0)
    model <- data.frame(t_test_result$p.value)
    model$Type <- i
    model$Species <- j
    result <- rbind(result, model)
  }
}
result2
