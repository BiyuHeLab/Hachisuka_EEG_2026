library(lme4)
library(MASS)
library(ordinal)
library(dplyr)
library(ggplot2)
library(broom)
library(tidyr)

edges <- seq(1, 241, by = 1) 
nWindows <- length(edges);
setwd("/isilon/LFMI/VMdrive/Ayaka/EEG/Distance2Bound_LMMstats/Signed_dist/SuperCorrect/")

######################## Load files:
csv_file <- "lmm_superV2.csv"

df <- read.csv(csv_file)
df$Subject <- as.integer(df$Subject)
#df$Dist <- as.numeric(scale(df$Dist))
super_df <- df

csv_file <- "lmm_basicV2.csv"

df <- read.csv(csv_file)
df$Subject <- as.integer(df$Subject)
#df$Dist <- as.numeric(scale(df$Dist))
basic_df <- df

  super_df <- super_df %>%
    group_by(Subject) %>%
    mutate(Dist = scale(Dist))%>%
    ungroup()

  basic_df <- basic_df %>%
    group_by(Subject) %>%
    mutate(Dist = scale(Dist))%>%
    ungroup()

# Average sig. timepoints for Basic, 62ms-112ms (92:103) and 338-378 (148:156), 62ms (93) is peak, 368 (154) is second peak.
dat_basic <- basic_df %>% filter(timebin %in% 153:153) %>%
  group_by(Trial,Subject,PAS, Cond) %>%
  summarize(avg_dist_basic = mean(Dist))

# Average sig. timepoints for Superordinate, 348 ms - 388ms (150:158), 363ms (153) is peak; 88 is peak for early (36 ms)
dat_super <-super_df %>% filter(timebin %in% 153:153) %>%
  group_by(Trial,Subject) %>%
  summarize(avg_dist_super = mean(Dist))

# Join everything into a single table.
dat <- dat_basic %>%
  left_join(dat_super, by = c("Trial","Subject"))

##################################
# store coefficients per subject
subjects <- unique(dat$Subject)
results <- data.frame(Subject = subjects, a = NA, b = NA, indirect = NA)

for (i in seq_along(subjects)) {
  subdat <- filter(dat, Subject == subjects[i])
  
  catB <- subdat$avg_dist_basic
  catA <- subdat$avg_dist_super
  
  # path a:
  model_a <- lm(catA ~ catB + as.factor(Cond), data = subdat)
  a_i <- coef(model_a)["catB"]
  
  # path b:
  #model_b <- clm(ordered(PAS) ~ catB + catA, data = subdat)
  model_b <- clm(ordered(PAS) ~ 1 + as.factor(Cond) + catB + catA, data = subdat)
  b_i <- coef(model_b)["catA"]
  cprime_i <- coef(model_b)["catB"]
  
  # path c:
 model_c <- clm(ordered(PAS) ~ 1 + as.factor(Cond) + catB, data = subdat)
 # model_c <- clm(ordered(PAS) ~ catB, data = subdat)
  c_i <- coef(model_c)["catB"]
  
  results[i, "a"] <- a_i
  results[i, "b"] <- b_i
  results[i, "direct"] <- cprime_i
  results[i, "indirect"] <- a_i * b_i
  results[i, "total"] <- c_i #cprime_i + a_i * b_i#
}

vars <- c("a", "b", "direct", "indirect", "total")

# Run t-tests and store p-values
pvals <- sapply(vars, function(v) wilcox.test(results[[v]], mu = 0)$p.value)

# Create a data frame with significance labels
sig_labels <- data.frame(
  variable = vars,
  p_value = pvals,
  label = ifelse(pvals < 0.001, "***",
                 ifelse(pvals < 0.01, "**",
                        ifelse(pvals < 0.05, "*",
                               ifelse(pvals < 0.08, "~", ""))))
)

sig_labels

df <- data.frame(
  subject = 1:length(subjects) ,
  a = results$a,
  b = results$b,
  indirect = results$indirect,
  direct = results$direct,
  total = results$total
)

df_long <- df %>%
  pivot_longer(cols = a:total, names_to = "variable", values_to = "value")

summary_df <- df_long %>%
  group_by(variable) %>%
  summarise(
    mean = mean(value),
    se = sd(value)/sqrt(n())
  ) %>%
  left_join(sig_labels, by = "variable")

ggplot(summary_df, aes(x = variable, y = mean, fill = variable)) +
  scale_fill_okabe_ito() +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2) +
  geom_text(aes(y = mean + se + 0.01, label = label), size = 12) +  # add asterisks above bars
  theme_minimal(base_size =22) +
  labs(y = "Value", x = "Variable")

