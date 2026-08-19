library(dplyr)
library(ordinal)
library(lme4)
library(ggplot2)
library(doParallel)
library(foreach)
library(R.matlab)

setwd("/isilon/LFMI/VMdrive/Ayaka/EEG/Revision_round1/Fig5/AllCond/PASrating1/")

edges <- seq(1, 241, by = 1) 
nWindows <- length(edges);
######################## Load files:
#csv_file <- "lmm_super.csv"
csv_file <- "lmm_basic.csv"

df <- read.csv(csv_file)

df$Subject <- as.integer(df$Subject)

df <- df %>%
  group_by(Subject) %>%
  mutate(Dist = scale(Dist))%>%
  ungroup()

betavals <- numeric(nWindows)
pvals <- numeric(nWindows)

# Setup parallel backend
n_cores <- 10 #as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
if (is.na(n_cores) || n_cores < 1) n_cores <- 1

cl <- makeCluster(n_cores)
registerDoParallel(cl)

on.exit({
  try(stopCluster(cl), silent = TRUE)
}, add = TRUE)

results <- foreach(w = 2:nWindows-1, .packages = c("lme4")) %dopar% {
  timeStart = edges[w] - 1
  timeEnd = edges[w] + 1
  thisGroup <- df[df$timebin >= timeStart & df$timebin <= timeEnd, ]
  # General:
  model <- glmer(BasicCorr ~ 1 + Dist + as.factor(Cond) + (1 | Subject), data = thisGroup, family = binomial(link="probit"))

  coefs <- summary(model)$coefficients
  
  coef_x1 <- coefs["Dist", "Estimate"]
  z_x1    <- coefs["Dist", "z value"]
  p_x1    <- coefs["Dist", "Pr(>|z|)"]

  c(beta1 = coef_x1, pval1 = p_x1, tval_dist = z_x2)
}

stopCluster(cl)
results_df <- do.call(rbind, results)
results_df <- as.data.frame(results_df)
results_df$time <- seq(-0.4, 0.8, length.out = 240)

results_df$p_value_adj <- results_df$pval1
ggplot(results_df, aes(x = time, y = beta1)) +
  geom_line(color = "blue") +  # Line for beta values
  geom_point(data = subset(results_df, p_value_adj < 0.05), color = "red", size = 2) +
  #geom_hline(yintercept = 0, linetype = "dotted", color = "black") +  # Horizontal line at y = 0
  #geom_vline(xintercept = 0, linetype = "dotted", color = "black") +  # Vertical line at x = 0
  labs(x = "Time", y = "Beta Values") +
  theme_minimal() +
  ggtitle("Basic ALL, cross-pos trials: Signed Dist (Dist & Cond fixed, 3 timebin)")

write.csv(results_df, "Basic_ALL_3timebin_signeddist.csv", row.names = FALSE)
