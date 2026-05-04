library(dplyr)
library(ordinal)
library(ggplot2)
library(doParallel)
library(foreach)
library(R.matlab)

setwd("/isilon/LFMI/VMdrive/Ayaka/EEG/Distance2Bound_LMMstats/Signed_dist/SuperCorrect/")

edges <- seq(1, 241, by = 1) 
nWindows <- length(edges);
######################## Load files:
#csv_file <- "lmm_superV2.csv"
csv_file <- "lmm_basicV2.csv"

df <- read.csv(csv_file)

df$Subject <- as.integer(df$Subject)

df <- df %>%
  group_by(Subject) %>%
  mutate(Dist = scale(Dist))%>%
  ungroup()

betavals <- numeric(nWindows)
pvals <- numeric(nWindows)

n_cores <- parallel::detectCores() - 2
cl <- makeCluster(n_cores)
registerDoParallel(cl)
results <- foreach(w = 2:nWindows-1, .packages = c("ordinal")) %dopar% {
  timeStart = edges[w] - 1
  timeEnd = edges[w] + 1
  thisGroup <- df[ df$timebin >= timeStart & df$timebin <= timeEnd, ]
  # General:
  model <- clmm(as.ordered(PAS) ~ 1 + Dist + Cond + (1 | Subject), data = thisGroup)

  vc <- vcov(model)
  
  coef_x1 <- coef(model)["Dist"]
  se1 <- sqrt(vc["Dist", "Dist"])
  t_value1 <- coef_x1 / se1
  
  pval_x1 <- summary(model)$coefficients["Dist", "Pr(>|z|)"]
  
  coef_x2 <- coef(model)["Cond"]
  se2 <- sqrt(vc["Cond", "Cond"])
  t_value2 <- coef_x2 / se2
  
  pval_x2 <- summary(model)$coefficients["Cond", "Pr(>|z|)"]
  c(beta1 = coef_x1, pval1 = pval_x1, beta2 = coef_x2, pval2 = pval_x2, tval_dist = t_value1, tval_cond = t_value2)
}

stopCluster(cl)
results_df <- do.call(rbind, results)
results_df <- as.data.frame(results_df)
results_df$time <- seq(-0.4, 0.8, length.out = 240)

results_df$p_value_adj <- results_df$pval1
ggplot(results_df, aes(x = time, y = beta1.Dist)) +
  geom_line(color = "blue") +  # Line for beta values
  geom_point(data = subset(results_df, p_value_adj < 0.05), color = "red", size = 2) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "black") +  # Horizontal line at y = 0
  geom_vline(xintercept = 0, linetype = "dotted", color = "black") +  # Vertical line at x = 0
  labs(x = "Time", y = "Beta Values") +
  theme_minimal() +
  ggtitle("Basic SuperCorrect trials: Signed Dist (Dist & Cond fixed, 3 timebin)")

# Sanity check plot:
means_V1 <- tapply(df$Dist, df$PAS, mean)  # mean dist per PAS
bar_centers <- barplot(means_V1,
                       names.arg = sort(unique(subject_means$Y)),
                       col = "lightblue",
                       #ylim = range(c(subject_means$V1, means_V1)),
                       ylab = "Basic distances",
                       xlab = "PAS",
                       main = "Basic vs PAS")

# add individual subject points
for (i in 1:nrow(df)) {
  y_level <- df$PAS[i]
  x_pos <- bar_centers[which(sort(unique(df$PAS)) == y_level)]
  points(x_pos, df$Dist[i], pch = 16, col = rgb(0,0,0,0.5))
}

# Example time vector and p-values
realtime<-seq(-0.4,0.8,by = 0.005)
time <- seq(0, 240, by = 1)  # e.g., 0 ms to 1000 ms in 10 ms bins
# Threshold for significance
alpha <- 0.05
p = results_df$pval1
# Find indices of significant time bins
sig_indices <- which(p < alpha)

# Print the corresponding time bins
sig_timebins <- time[sig_indices]
print(sig_timebins)

write.csv(results_df, "Basic_SuperCorrect_3timebin_signeddist_zscore.csv", row.names = FALSE)
