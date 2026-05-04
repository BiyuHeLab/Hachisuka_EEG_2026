#install.packages("dplyr")
#install.packages("ordinal")
#install.packages("ggplot2")
#install.packages("doParallel")
#install.packages("foreach")
#install.packages("dplyr")
#install.packages("tibble")

library(dplyr)
library(ordinal)
library(ggplot2)
library(doParallel)
library(foreach)
library(dplyr)
library(tibble)

setwd("/isilon/LFMI/VMdrive/Ayaka/EEG/CosmoMVPA_results/Final_Resubmission_results_April16/Signed_dist/CrossPosition/SuperCorrect/AllCond/")

# Read .csv
csv_file  <- "lmm_basic.csv"

# It is advised that you make multiple copies of this script, and run the permutation in chunks (by 250 perms, for example)
start_perm <- 1
end_perm   <- 1000

edges <- seq(1, 241, by = 1) 
nWindows <- length(edges);
  
df <- read.csv(csv_file)

# Setup parallel backend
n_cores <- parallel::detectCores() - 2
cl <- makeCluster(n_cores)
registerDoParallel(cl)

nPerms <- end_perm  # number of permutations
for (i in start_perm:nPerms) {
  print(sprintf("Running permutation %d...", i))
  #Step 0: z-score the distances, per subject.
  df <- df %>%
    group_by(Subject) %>%
    mutate(Dist = scale(Dist)) %>%
    ungroup()
  
  # Step 1: Assign unique TrialID per Subject
  df <- df %>%
    group_by(Subject, Trial) %>%
    mutate(TrialID = cur_group_id()) %>%
    ungroup()
  
  # Step 2: Shuffle trial blocks WITHIN each subject
  df_shuffled <- df %>%
    group_by(Subject) %>%
    group_split() %>%
    lapply(function(subject_df) {
      # Split into trial blocks within subject
      trial_blocks <- split(subject_df, subject_df$TrialID)
      # Shuffle these blocks
      shuffled_blocks <- sample(trial_blocks)
      # Replace Dist only
      subject_df$Dist <- unlist(lapply(shuffled_blocks, function(block) block$Dist))
      return(subject_df)
    }) %>%
    bind_rows() %>%
    select(-TrialID)
  
  results <- foreach(w = 2:nWindows-1, .packages = c("ordinal")) %dopar% {
    timeStart = edges[w] - 1
    timeEnd = edges[w] + 1
    thisGroup <- df_shuffled[ df_shuffled$timebin >= timeStart & df_shuffled$timebin <= timeEnd, ]
    
    # General:
    model <- clmm(as.ordered(PAS) ~ 1 + Dist + as.factor(Cond) + (1 | Subject), data = thisGroup)
    vc <- vcov(model)
    
    coef_x1 <- coef(model)["Dist"]
    se1 <- sqrt(vc["Dist", "Dist"])
    t_value1 <- coef_x1 / se1
    
    pval_x1 <- summary(model)$coefficients["Dist", "Pr(>|z|)"]
    c(beta1 = coef_x1, pval1 = pval_x1, tval_dist = t_value1)
    
    }
  
  results_df <- do.call(rbind, results)
  results_df <- as.data.frame(results_df)
  results_df$time <- seq(-0.4, 0.8, length.out = 240)
  write.csv(results_df, paste0("perms/Basic/Basic_SuperCorr_perm_", i, ".csv"), row.names = TRUE)
}

stopCluster(cl)

}
