# simulation.R (for Subtype Differentiation)
#
# Description:
# This script simulates a scenario to test if a model can differentiate between
# two disease subtypes (A and B). The model is trained on data containing only
# Subtype A vs. Healthy Controls, but leverages an external "source" GWAS
# from a general disease (a mix of A and B).
#
# The key comparison is between P+T, standard Elastic Net, and ePRS to see
# which method can better learn to distinguish A from B.

# --- 1. Load Libraries ---
library(glmnet)
library(MASS)
library(pROC)
library(ggplot2)
library(dplyr)
library(tidyr)

# --- 2. Simulation Parameters ---
sim_count <- 40           # Number of full simulation runs
iter_count <- 15          # Number of r_g steps per run
n_row <- 400              # Base sample size
n_col <- 1000             # Number of SNPs

# Define the SNP correlation structure (AR1 with rho=0)
ar1_cor <- function(n, rho) {
  exponent <- abs(matrix(1:n - 1, nrow = n, ncol = n, byrow = TRUE) - (1:n - 1))
  rho^exponent
}
corr_matrix <- ar1_cor(n_col, 0)

# --- 3. Initialize Result Matrices ---
# We will compare three main methods: P+T, Elastic Net, and ePRS.
r_g_subtype <- r_g_eprs <- p_t <- en <- eprs <- matrix(nrow = sim_count, ncol = iter_count)


# --- 4. Main Simulation Loop ---
for (sim in 1:sim_count) {
  # --- Simulate Genotype Data ---
  df_source <- mvrnorm(n = n_row * 8, mu = rep(0, n_col), Sigma = corr_matrix)
  df_target <- mvrnorm(n = n_row, mu = rep(0, n_col), Sigma = corr_matrix)
  df_test <- mvrnorm(n = n_row, mu = rep(0, n_col), Sigma = corr_matrix)
  
  for (iter in 1:iter_count) {
    # --- Define True Effect Sizes (Betas) ---
    # beta_a: Effects for Subtype A
    beta_a <- c(rep(5, 15), rep(0, n_col - 15))
    # beta_b: Effects for Subtype B, with varying overlap with Subtype A
    beta_b <- c(rep(0, iter), rep(5, 15), rep(0, n_col - 15 - iter))
    # beta_source: Effects for the general disease (sum of A and B)
    beta_source <- beta_a + beta_b
    
    
    # Genetic correlation between subtypes
    r_g_subtype[sim, iter] = sum(beta_a*beta_b) / sqrt((sum(beta_a^2) * sum(beta_b^2)))
    # Genetic correlation between Subtype A and the general source phenotype.
    r_g_eprs[sim, iter] <- sum(beta_a * beta_source) / sqrt(sum(beta_a^2) * sum(beta_source^2))
    
    # --- Simulate Genetic Risk & Phenotype Status ---
    # Calculate underlying genetic risk scores for each subtype
    source_a_risk <- df_source %*% beta_a + rnorm(n_row * 8, 0, 0)
    source_b_risk <- df_source %*% beta_b + rnorm(n_row * 8, 0, 0)
    target_a_risk <- df_target %*% beta_a + rnorm(n_row, 0, 10)
    target_b_risk <- df_target %*% beta_b + rnorm(n_row, 0, 10)
    test_a_risk <- df_test %*% beta_a + rnorm(n_row, 0, 10)
    test_b_risk <- df_test %*% beta_b + rnorm(n_row, 0, 10)
    
    # An individual is a "case" of general phenotype if the combined risk is high
    # General Phenotype
    source_g_risk <- source_a_risk + source_b_risk + rnorm(n_row * 8, 0, 30)
    target_g_risk <- target_a_risk + target_b_risk + rnorm(n_row, 0, 30)
    test_g_risk <- test_a_risk + test_b_risk + rnorm(n_row, 0, 30)
    
    # Convert risks to binary status 
    y_source <- ifelse(source_g_risk > median(source_g_risk), 1, 0)
    target_status = ifelse(target_g_risk > median(target_g_risk), 1, 0)
    
    # --- Construct the Final Training and Test Datasets ---
    
    ##-- Create the Test Set: Subtype A vs. Subtype B --##
    # Identify which subtype has higher genetic risk for each person in the test set
    diff_test <- test_a_risk - test_b_risk
    # Classify as A if risk_A >> risk_B, and as B if risk_B >> risk_A
    subtype_test <- ifelse(diff_test > 0.2*sd(diff_test), 1,  # Subtype A
                           ifelse(diff_test < -0.2*sd(diff_test), 0, -1)) # Subtype B, -1 is ambiguous
    
    test_set_indices <- which(subtype_test != -1 & test_g_risk > median(test_g_risk))
    df_test_final <- df_test[test_set_indices, ]
    y_test_final <- subtype_test[test_set_indices]
    
    ##-- Create the Target Training/Validation Set: Subtype A vs. Healthy Controls --##
    diff_target <- target_a_risk - target_b_risk
    subtype_target <- ifelse(diff_target > 0.2*sd(diff_target), 1,
                             ifelse(diff_target < -0.2*sd(diff_target), 0, -1))   # Is it more likely A or B?
    
    # We can only confidently identify "Subtype A" cases and "Healthy" controls
    indices_subtype_A <- which(target_status == 1 & subtype_target == 1)
    indices_healthy <- which(target_status == 0)
    
    df_target_A_vs_healthy <- rbind(df_target[indices_subtype_A, ], df_target[indices_healthy, ])
    y_target_A_vs_healthy <- c(rep(1, length(indices_subtype_A)), rep(0, length(indices_healthy)))
    
    # Split this target data for hyperparameter tuning
    n_target_final <- nrow(df_target_A_vs_healthy)
    
    train_size <- floor(0.7 * n_target_final)
    train_idx <- sample(1:n_target_final, train_size)
    valid_idx <- setdiff(1:n_target_final, train_idx)
    
    # --- Run "External" Source GWAS ---
    # Generates p-values from the general disease GWAS
    source_pvals <- sapply(1:n_col, function(i) {
      model <- summary(glm(y_source ~ df_source[, i], family = "binomial"))
      return(model$coef[2, 4])
    })
    
    # Generates effect sizes from the target GWAS
    target_effects <- sapply(1:n_col, function(i) {
      model <- summary(glm(y_target_A_vs_healthy[train_idx] ~ df_target_A_vs_healthy[train_idx, i]))
      return(model$coef[2, 1])
    })
    
    # --- Run and Evaluate Models ---
    
    # Define hyperparameter grid for tuning
    alpha_thr <- c(0, 0.2, 0.4, 0.6, 0.8, 1)
    
    ##-- Method 1: P+T --##
    p_threshold <- c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5, 1)
    valid_aucs_pt <- sapply(p_threshold, function(threshold) {
      selected_snps <- which(source_pvals < threshold)
      pred <- df_target_A_vs_healthy[valid_idx, selected_snps, drop = FALSE] %*% as.matrix(target_effects[selected_snps])
      roc(y_target_A_vs_healthy[valid_idx], as.vector(pred), direction = "<", quiet = TRUE)$auc
    })
    best_p_threshold <- p_threshold[which.max(valid_aucs_pt)]
    
    selected_snps <- which(source_pvals < best_p_threshold)
    pred_pt <- df_test_final[, selected_snps, drop = FALSE] %*% as.matrix(target_effects[selected_snps])
    p_t[sim, iter] <- roc(y_test_final, as.vector(pred_pt), direction = "<", quiet = TRUE)$auc
    
    ##-- Method 2: Standard Elastic Net (EN) --##
    valid_aucs_en <- sapply(alpha_thr, function(alpha) {
      m <- cv.glmnet(df_target_A_vs_healthy[train_idx,], y_target_A_vs_healthy[train_idx],
                     family = "binomial", alpha = alpha)
      pred <- predict(m, df_target_A_vs_healthy[valid_idx,], s = "lambda.min", type = "response")
      roc(y_target_A_vs_healthy[valid_idx], as.vector(pred), direction = "<", quiet=TRUE)$auc
    })
    best_alpha_en <- alpha_thr[which.max(valid_aucs_en)]
    
    # Train final EN model on all available target data
    final_en_model <- cv.glmnet(df_target_A_vs_healthy, y_target_A_vs_healthy, family = "binomial", alpha = best_alpha_en)
    pred_en <- predict(final_en_model, df_test_final, s = "lambda.min", type = "response")
    en[sim, iter] <- roc(y_test_final, as.vector(pred_en), direction = "<", quiet=TRUE)$auc
    
    ##-- Method 3: ePRS with Elastic Net --##
    # E_j <- 1 / -log10(source_pvals)
    E_j <- 10*(1 - (1-source_pvals)^8)
    E_j[is.infinite(E_j)] <- 1e6 # Handle p=1 case
    pf <- E_j * r_g_eprs[sim, iter] + (1 - r_g_eprs[sim, iter])
    
    valid_aucs_eprs <- sapply(alpha_thr, function(alpha) {
      m <- cv.glmnet(df_target_A_vs_healthy[train_idx,], y_target_A_vs_healthy[train_idx],
                     family = "binomial", alpha = alpha, penalty.factor = pf)
      pred <- predict(m, df_target_A_vs_healthy[valid_idx,], s = "lambda.min", type = "response")
      roc(y_target_A_vs_healthy[valid_idx], as.vector(pred), direction = "<", quiet=TRUE)$auc
    })
    best_alpha_eprs <- alpha_thr[which.max(valid_aucs_eprs)]
    
    # Train final ePRS model
    final_eprs_model <- cv.glmnet(df_target_A_vs_healthy, y_target_A_vs_healthy, family = "binomial",
                                  alpha = best_alpha_eprs, penalty.factor = pf)
    pred_eprs <- predict(final_eprs_model, df_test_final, s = "lambda.min", type = "response")
    eprs[sim, iter] <- roc(y_test_final, as.vector(pred_eprs), direction = "<", quiet=TRUE)$auc
    
    print(paste("sim:", sim, "| Iter:", iter, "| r_g:", round(r_g_subtype[sim, iter], 2),
                "| P+T AUC:", round(p_t[sim, iter], 3),
                "| EN AUC:", round(en[sim, iter], 3),
                "| ePRS AUC:", round(eprs[sim, iter], 3)))
  }
}

# --- 5. Plot Results ---
# Convert wide matrices to a long (tidy) dataframe for ggplot
df_results <- data.frame(
  auc = c(as.vector(p_t), as.vector(en), as.vector(eprs)),
  rg = c(as.vector(r_g_subtype), as.vector(r_g_subtype), as.vector(r_g_subtype)),
  method = rep(c("P+T", "Elastic Net", "ePRS"), each = length(en))
)

# Create the plot
ggplot(df_results, aes(x = rg, y = auc, color = method)) +
  geom_smooth(se = TRUE, alpha = 0.2, span = 0.5) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "red") +
  scale_color_manual(values = c("P+T" = "#E69F00", "Elastic Net" = "#0072B2", "ePRS" = "#D55E00")) +
  labs(
    title = "ePRS Improves Subtype Differentiation",
    subtitle = "Model trained on 'Subtype A vs. Healthy', tested on 'Subtype A vs. Subtype B'",
    x = "Genetic Correlation between Subtypes",
    y = "AUC",
    color = "Method"
  ) +
  coord_cartesian(ylim = c(0.45, 0.9)) +
  theme_bw(base_size = 14)
