library(R.matlab)
library(tidyverse)
library(rstatix)

library(afex)
library(emmeans)
library(lme4)
library(lmerTest)
library(permuco)

setwd("C:/Users/ahach/OneDrive/Desktop/C2F_manuscript/Behavior_statsv2/C2F/")

mat_data <- readMat("corrbyPAS.mat")

super <- mat_data$corrbyPAS[[1]][[1]] #cell 1 is for "superordinate"
basic <- mat_data$corrbyPAS[[2]][[1]] #cell 2 is for "basic"

condition_labels <- c("Normal","Masked","LSF")
reshape_scores <- function(score_array, Category_level) {
  as.data.frame.table(score_array, responseName = "PercCorr") %>%
    rename(Subject = Var1, Condition = Var2, PAS = Var3) %>%
    mutate(
      Subject = as.integer(Subject),
      Condition = factor(condition_labels[as.integer(Condition)],
                         levels = condition_labels),
      PAS = as.numeric(PAS),
      Category_level = factor(Category_level)
    )
}

df_A <- reshape_scores(super, "super")
df_B <- reshape_scores(basic, "basic")

# Combine both into a single table:
df <- bind_rows(df_A, df_B)
df$PAS <- as.factor(df$PAS)
df <- na.omit(df)  # removes any rows with NA factors

# 3-way RM anova, with effect sizess via partial eta squared:
anova_model <- anova_test(data = df, dv = PercCorr, wid = Subject, 
                          within = c("Category_level", "Condition", "PAS"),
                          effect.size = "pes")
get_anova_table(anova_model)

#another ANOVA package for post-hoc comparison (main effects are the same)
anova_model <- aov_ez(data = df, dv = "PercCorr", id = "Subject", 
                      within = c("Category_level", "Condition", "PAS"))
summary(anova_model)
emmeans(anova_model, pairwise ~ Category_level | PAS * Condition,adjust = "bonferroni")

###################################################
#Test of normality for 1-way RM
anova_model <- aov_ez(data = df, dv = "PercCorr", id = "Subject", 
                      within = c("PAS"))
res <- residuals(anova_model$lm) #residuals
shapiro.test(res) # Shapiro-Wilk test
qqnorm(res)
qqline(res)

#Test of normality for 2-way RM
anova_model <- aov_ez(data = df, dv = "PercCorr", id = "Subject", 
                      within = c("Category_level", "Condition"))
res <- residuals(anova_model$lm) #residuals
shapiro.test(res) # Shapiro-Wilk test
qqnorm(res)
qqline(res)

#Test of normality for 3-way RM
anova_model <- aov_ez(data = df, dv = "PercCorr", id = "Subject", 
                      within = c("Category_level", "Condition", "PAS"))
res <- residuals(anova_model$lm) #residuals
shapiro.test(res) # Shapiro-Wilk test
qqnorm(res)
qqline(res)

# Non-parametric test instead of ANOVA for the 3-way RM:
# Permutation approach w/ permuco
df$Subject <- factor(df$Subject)
fit <- aovperm(
  formula =PercCorr ~ Category_level * Condition * PAS + Error(Subject/(Category_level * Condition * PAS)),
  data = df,
  np = 10000, #10,000 MC permutation
)
summary(fit)

# Non-parametric post-hoc with Wilcoxon signed-rank test
wilcox_A_by_BC <- df %>%
  group_by(PAS, Condition) %>%
  wilcox_test(
    PercCorr ~ Category_level,
    paired = TRUE,
    subject = "Subject"
  ) %>%
  adjust_pvalue(method = "bonferroni")


###################################

# Wilcoxon-signed rank test against chance:
wilcox_result <- full_df %>%
  group_by(Category_level, Condition, PAS) %>%
  summarize(
    p = wilcox.test(PercCorr, mu = 0.5, alternative = "greater")$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    p_fdr = p.adjust(p, method = "fdr"),
    sig_label = case_when(
      p_fdr < 0.001 ~ "***",
      p_fdr < 0.01  ~ "**",
      p_fdr < 0.05  ~ "*",
      TRUE          ~ ""
    )
  )

print(wilcox_result,n = Inf)


