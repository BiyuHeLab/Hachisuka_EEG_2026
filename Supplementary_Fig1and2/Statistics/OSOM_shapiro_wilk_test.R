library(R.matlab)
library(tidyverse)
library(rstatix)

library(afex)
library(emmeans)
library(lme4)
library(lmerTest)
library(ARTool)

setwd("C:/Users/ahach/OneDrive/Desktop/C2F_manuscript/Behavior_statsv2/OSOM_Exp1/")

data <- read.csv("Exp1_SuperBasicHR.csv")

df <- data %>%
  pivot_longer(
    cols = c(Normal_Super, OSM_Super, OSOM_Super, Normal_Basic, OSM_Basic, OSOM_Basic),
    names_to = c("Super", "Basic"),
    names_pattern = "(Normal|OSM|OSOM)_(Super|Basic)",
    values_to = "PercCorr"
  )

names(df)[names(df) == "Super"] <- "Condition"
names(df)[names(df) == "Basic"] <- "Category"

df$Subject <- factor(rep(1:(nrow(df)/6), each = 6))
df$Category <- factor(df$Category, levels = c("Super", "Basic"))
df$Condition <- factor(df$Condition, levels = c("Normal", "OSM", "OSOM"))


anova_model <- aov_ez(data = df, dv = "PercCorr", id = "Subject", 
                      within = c("Category", "Condition"))
summary(anova_model)

#Test of normality for 2-way RM
res <- residuals(anova_model$lm) #residuals
shapiro.test(res) # Shapiro-Wilk test
qqnorm(res)
qqline(res)

############################################################
data <- read.csv("Exp1_PAS.csv")

df <- data %>%
  pivot_longer(
    cols = c(Normal, OSM, OSOM),
    values_to = "PAS"
  )

names(df)[names(df) == "name"] <- "Condition"

df$Subject <- factor(rep(1:(nrow(df)/3), each = 3))
df$Condition <- factor(df$Condition, levels = c("Normal", "OSM", "OSOM"))


anova_model <- aov_ez(data = df, dv = "PAS", id = "Subject", 
                      within = c("Condition"))
summary(anova_model)

#Test of normality for 1-way RM
anova_model <- aov_ez(data = df, dv = "PAS", id = "Subject", 
                      within = c("Condition"))
res <- residuals(anova_model$lm) #residuals
shapiro.test(res) # Shapiro-Wilk test
qqnorm(res)
qqline(res)

#############################################################
setwd("C:/Users/ahach/OneDrive/Desktop/C2F_manuscript/Behavior_statsv2/OSOM_Exp2/")
data <- read.csv("Exp2_SuperBasicHR.csv")

df <- data %>%
  pivot_longer(
    cols = c(OSM_Super, OSOMstable_Super, OSOMmoved_Super, OSM_Basic, OSOMstable_Basic, OSOMmoved_Basic),
    names_to = c("Super", "Basic"),
    names_pattern = "(OSM|OSOMstable|OSOMmoved)_(Super|Basic)",
    values_to = "PercCorr"
  )

names(df)[names(df) == "Super"] <- "Condition"
names(df)[names(df) == "Basic"] <- "Category"

df$Subject <- factor(rep(1:(nrow(df)/6), each = 6))
df$Category <- factor(df$Category, levels = c("Super", "Basic"))
df$Condition <- factor(df$Condition, levels = c("OSM", "OSOMstable", "OSOMmoved"))


anova_model <- aov_ez(data = df, dv = "PercCorr", id = "Subject", 
                      within = c("Category", "Condition"))
summary(anova_model)

#Test of normality for 2-way RM
res <- residuals(anova_model$lm) #residuals
shapiro.test(res) # Shapiro-Wilk test
qqnorm(res)
qqline(res)

############################################################

data <- read.csv("Exp2_PAS.csv")

df <- data %>%
  pivot_longer(
    cols = c(OSM, OSOMstable, OSOMmoved),
    values_to = "PAS"
  )

names(df)[names(df) == "name"] <- "Condition"

df$Subject <- factor(rep(1:(nrow(df)/3), each = 3))
df$Condition <- factor(df$Condition, levels = c("OSM", "OSOMstable", "OSOMmoved"))


anova_model <- aov_ez(data = df, dv = "PAS", id = "Subject", 
                      within = c("Condition"))
summary(anova_model)

#Test of normality for 1-way RM
anova_model <- aov_ez(data = df, dv = "PAS", id = "Subject", 
                      within = c("Condition"))
res <- residuals(anova_model$lm) #residuals
shapiro.test(res) # Shapiro-Wilk test
qqnorm(res)
qqline(res)

friedman.test(PAS ~ Condition | Subject, data = df)
pairwise.wilcox.test(
  x = df$PAS,
  g = df$Condition,
  paired = TRUE,           
  p.adjust.method = "bonferroni"
)
