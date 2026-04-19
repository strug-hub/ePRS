# run_simulation.R
#
# Description:
# This script simulates genotype-phenotype data to compare the predictive
# performance of different methods. The core comparison is between
# standard Elastic Net and ePRS, which leverages p-values from an external
# source GWAS to improve prediction in a smaller target cohort.
# The simulation varies the genetic correlation (r_g) between the source and
# target phenotypes to assess robustness.

# --- 1. Load Libraries---
library(glmnet)
library(MASS)
library(ggplot2)

# --- 2. Simulation Parameters ---
set.seed(1000)
sim_count <- 40           # Number of full simulation runs
iter_count <- 40          # Number of different r_g values per run
n_row <- 300              # Sample size of target cohort
n_col <- 1000             # Number of simulated SNPs

# Train-test splits
train_prop <- 0.55
valid_prop <- 0.15

# Define the SNP correlation structure (AR1 with rho=0)
ar1_cor <- function(n, rho) {
  exponent <- abs(matrix(1:n - 1, nrow = n, ncol = n, byrow = TRUE) - (1:n - 1))
  rho^exponent
}
corr_matrix <- ar1_cor(n_col, 0)

# --- 3. Initialize Result Matrices ---
# We will compare four main methods:
# - source: External-source benchmark
# - p_t: P+T using target data for beta estimation
# - en: Standard Elastic Net
# - eprs: ePRS with Elastic Net

r_g_values <- source <- p_t <- en <- eprs <- matrix(nrow = sim_count, ncol = iter_count)

# --- 4. Main Simulation Loop ---
for (sim in 1:sim_count) {

  # Source data (larger N)
  df_source <- mvrnorm(n = n_row * 5, mu = rep(0, n_col), Sigma = corr_matrix)
  # Target data (smaller N)
  df_target <- mvrnorm(n = n_row, mu = rep(0, n_col), Sigma = corr_matrix)

  # True effect sizes for source phenotype
  beta_source <- c(rep(1, iter_count), rep(0, n_col - iter_count))

  for (iter in 1:iter_count) {
    # --- Define Target Effect Sizes and Genetic Correlation ---
    # Slide the block of causal variants to control overlap with beta_source
    beta_target <- c(rep(0, iter), rep(1, iter_count), rep(0, n_col - iter_count - iter))
    r_g_values[sim, iter] <- sum(beta_source * beta_target) / sqrt((sum(beta_source^2) * sum(beta_target^2)))

    # --- Simulate Phenotypes ---
    target_noise_sd <- 0    # Figure 1 main-text setting, h2 = 1
    # target_noise_sd <- 4  # residual-noise comparison, h2 ≈ 0.71
    y_source <- df_source %*% beta_source + rnorm(nrow(df_source), 0, 0)   
    y_target <- df_target %*% beta_target + rnorm(n_row, 0, target_noise_sd)              

    # --- Run "External" Source GWAS ---
    # The corresponding p-values are used as external evidence for ePRS.
    source_gwas <- vapply(1:n_col, function(i) {
      coef(summary(lm(y_source ~ df_source[, i])))[2, c(1, 4)]
    }, numeric(2))
    beta_source_hat <- source_gwas[1, ]
    pval <- source_gwas[2, ]

    # --- Split Target Data into Training/Validation/Test Sets ---
    # For tuning alpha and final evaluation
    n_train <- floor(train_prop * n_row)
    n_valid <- floor(valid_prop * n_row)
    n_test  <- n_row - n_train - n_valid

    train_idx <- seq_len(n_train)
    valid_idx <- seq_len(n_valid) + n_train
    test_idx  <- seq_len(n_test) + n_train + n_valid

    y_train <- y_target[train_idx]
    y_valid <- y_target[valid_idx]
    y_test <- y_target[test_idx]

    # --- Method 1: Source PRS benchmark ---
    # External-source benchmark: apply oracle source effect sizes directly
    # to the target cohort without refitting in the target data.
    y_pred_source <- as.vector(df_target[test_idx, ] %*% beta_source)
    source[sim, iter] <- cor(y_test, y_pred_source)

    # --- Method 2: P+T (Clumping and Thresholding) using Target Betas ---
    effect_size = pval_target = array(dim=n_col)
    for(i in 1:n_col){
      model = summary(lm(y_train ~ df_target[train_idx,i]))$coef
      effect_size[i] = model[2,1]; pval_target[i] = model[2,4]
    }

    p_threshold = c(1, 0.5, 0.4, 0.3, 0.2, 0.1, 0.05, 0.01, 0.001); to_compare = NULL
    for(pthr in p_threshold){
      selected_snps = which(pval < pthr)
      y_pred = df_target[valid_idx, selected_snps] %*% as.matrix(effect_size[selected_snps])
      to_compare = append(to_compare, (cor(y_valid, y_pred)))
    }
    pthr = p_threshold[which.max(to_compare)]
    selected_snps = which(pval < pthr)
    y_pred = df_target[test_idx, selected_snps] %*% as.matrix(effect_size[selected_snps])
    p_t[sim, iter] = cor(y_test, y_pred)

    # --- Method 3: Standard Elastic Net (EN) ---
    # Tune alpha on a validation set
    alpha_thr <- c(0, 0.2, 0.4, 0.6, 0.8, 1)
    validation_cors <- sapply(alpha_thr, function(alpha_tune) {
      m <- cv.glmnet(df_target[train_idx, ], y_train, nfolds = 3, alpha = alpha_tune)
      y_pred_valid <- predict(m, df_target[valid_idx, ], s = "lambda.min")
      cor(y_valid, y_pred_valid)
    })
    # Fall back to ridge regression if no predictors are selected during cross-validation
    best_alpha_en <- alpha_thr[which.max(validation_cors)]
    best_alpha_en <- ifelse(sum(is.na(alpha_thr[validation_cors]))==length(alpha_thr), 0, alpha_thr[which.max(validation_cors)])
    
    # Train final model on combined train+valid set and evaluate on test set
    final_en_model <- cv.glmnet(df_target[c(train_idx, valid_idx), ], y_target[c(train_idx, valid_idx)], nfolds = 3, alpha = best_alpha_en)
    y_pred_en <- predict(final_en_model, df_target[test_idx, ], s = "lambda.min")
    en[sim, iter] <- cor(y_test, y_pred_en)

    # --- Method 4: ePRS with Elastic Net ---
    # Define the evidence term E_j from source p-values.
    # We use a transformation that is more sensitive to p-value changes than -log10(p).
    E_j <- 1/-log10(pval)
    # E_j <- 10 * (1 - (1 - pval)^8)

    # The penalty factor combines external evidence (E_j) with genetic correlation (r_g)
    pf <- E_j * r_g_values[sim, iter] + (1 - r_g_values[sim, iter])

    # Tune alpha for ePRS
    validation_cors_eprs <- sapply(alpha_thr, function(alpha_tune) {
      m <- cv.glmnet(df_target[train_idx, ], y_train, nfolds = 3, alpha = alpha_tune, penalty.factor = pf)
      y_pred_valid <- predict(m, df_target[valid_idx, ], s = "lambda.min")
      cor(y_valid, y_pred_valid)
    })
    best_alpha_eprs <- alpha_thr[which.max(validation_cors_eprs)]

    # Train final ePRS model and evaluate on test set
    final_eprs_model <- cv.glmnet(df_target[c(train_idx, valid_idx), ], y_target[c(train_idx, valid_idx)], nfolds = 3, alpha = best_alpha_eprs, penalty.factor = pf)
    y_pred_eprs <- predict(final_eprs_model, df_target[test_idx, ], s = "lambda.min")
    eprs[sim, iter] <- cor(y_test, y_pred_eprs)
    
    print(paste("sim:", sim, "| Iter:", iter, "| r_g:", round(r_g_values[sim, iter], 2)))
  }
}

# --- 5. Plotting Results ---
df = data.frame(r_sq = c(source, p_t, en, eprs),
                rg = rep(r_g_values, 4),
                method = rep(c("Source", "P+T", "EN", "ePRS"), each = sim_count * iter_count))

ggplot(df, aes(x = rg, y = r_sq)) + 
  geom_point(aes(col = method)) +
  geom_smooth(col = "black", lty = "dashed") +
  facet_grid(method ~ .)

ggplot(df, aes(x = rg, y = r_sq)) +
  geom_smooth(aes(col = method), lty = "dashed", se = FALSE)
