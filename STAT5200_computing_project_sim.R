library(glmnet)
library(ggplot2)

##########################################################################
# function: generate_data
generate_data <- function(n, p, s) {
  
  # observed predictors
  X <- matrix(rnorm(n * p), n, p)
  
  # ---- Create true beta coefficients (length s) ----
  betas_true <- runif(s, 1, 3) * sample(c(-1, 1), s, replace = TRUE)
  
 
  if (s <= p) {
    betas_obs <- c(betas_true[1:s], rep(0, p - s))
    y <- X %*% betas_obs + rnorm(n, 0, 1)
    
    # ---- Case 2: s > p ----
    # All p observed betas are nonzero
    # Remaining (s - p) signals affect y but are unobserved latent variables
  } else {
    betas_obs <- betas_true[1:p]              # p nonzero betas
    X_latent <- matrix(rnorm(n * (s - p)), n, s - p)
    y <- X %*% betas_obs + X_latent %*% betas_true[(p+1):s] + rnorm(n, 0, 1)
  }
  y <- as.vector(y)
 
  # final data frame
  sim_data <- data.frame(y = y, X = X)
  colnames(sim_data) <- c("y", paste0("X", 1:p))
  
  return(list(sim_data = sim_data, betas_obs = betas_obs, betas_true = betas_true))
}




#########################################################################
# function: generate_data 2, makes a dense model with a set signal to noise ratio
generate_data2 <- function(n, p, SNR = 5, sigma = 1) {
  
  # X ~ N(0,1)
  X <- matrix(rnorm(n * p), n, p)
  
  # Draw random direction for beta
  beta_raw <- rnorm(p)
  norm_raw <- sqrt(sum(beta_raw^2))
  
  # Set ||beta||^2 / sigma^2 = SNR
  r <- sqrt(SNR) * sigma
  beta <- beta_raw * (r / norm_raw)
  
  # Generate y = X beta + noise
  y <- as.vector(X %*% beta + rnorm(n, 0, sigma))
  
  sim_data <- data.frame(y = y, X)
  colnames(sim_data) <- c("y", paste0("X", 1:p))
  
  return(list(sim_data = sim_data, beta = beta))
}


##########################################################################
# Ridgeless OLS using pseudoinverse
ridgeless_fit <- function(X, y, tol = 1e-8, lambda = 0) {
  # Compute SVD
  svd_X <- svd(X)
  
  # Compute D_plus, handling small singular values
  D_plus <- diag(ifelse(svd_X$d > tol, 1 / svd_X$d, 0))
  
  # Optionally include tiny ridge regularization
  if (lambda > 0) {
    D_plus <- diag(ifelse(svd_X$d > tol, svd_X$d / (svd_X$d^2 + lambda), 0))
  }
  
  # Compute pseudoinverse
  X_pinv <- svd_X$v %*% D_plus %*% t(svd_X$u)
  
  # Compute coefficients
  beta_hat <- X_pinv %*% y
  
  return(beta_hat)
}


##########################################################################
# Fit models and compute test MSE
model_fits <- function(train_X, train_y, test_X, test_y, train) {
  
  
  # Lasso
  fit_lasso <- cv.glmnet(train_X, train_y, alpha = 1, nfolds = 5, intercept = FALSE,
                         standardize = FALSE)
  lasso_pred <- as.numeric(predict(fit_lasso, test_X, s = "lambda.min"))
  lasso_mse <- mean((test_y - lasso_pred)^2)
  
  # Ridge
  fit_ridge <- cv.glmnet(train_X, train_y, alpha = 0, nfolds = 5,
                         intercept = FALSE, standardize = FALSE)
  ridge_pred <- as.numeric(predict(fit_ridge, test_X, s = "lambda.min"))
  ridge_mse <- mean((test_y - ridge_pred)^2)
  
 
  
  # Ridgeless OLS
  beta_rl <- ridgeless_fit(train_X, train_y)
  rl_pred <- as.numeric(test_X %*% beta_rl)
  ridgeless_mse <- mean((test_y - rl_pred)^2)
  
  if (train == 1) {
    train_lasso_pred <- as.numeric(predict(fit_lasso, train_X, s = "lambda.min"))
    train_lasso_mse <- mean((train_y - train_lasso_pred)^2)
    
    train_ridge_pred <- as.numeric(predict(fit_ridge, train_X, s = "lambda.min"))
    train_ridge_mse <- mean((train_y - train_ridge_pred)^2)
    
    
    # to get a nice graph showing interpolation when p >= n
    train_rl_pred <- as.numeric(train_X %*% beta_rl)
    train_ridgeless_mse <- mean((train_y - train_rl_pred)^2)
    
    
    
    return(list(lasso_mse = lasso_mse, ridge_mse = ridge_mse, 
                ridgeless_mse = ridgeless_mse,
                train_lasso_mse = train_lasso_mse,
                train_ridge_mse = train_ridge_mse,
                train_ridgeless_mse = train_ridgeless_mse))
  }
  
  
  return(list(lasso_mse = lasso_mse, ridge_mse = ridge_mse, 
               ridgeless_mse = ridgeless_mse,
              train_ridgless_mse = train_ridgless_mse))
  
}
##########################################################################
# Compute single simulation
# train variable controls whether we include training mse also to get another graph
compute_simulation <- function(n, p, s, seed, train = 0, SNR) {
  set.seed(seed)
  
  if (SNR == 0) {
    dat <- generate_data(n, p, s)
  } else {
    dat <- generate_data2(n, p, SNR)
  }
  
  sim_df <- dat$sim_data
  X <- as.matrix(sim_df[, -1])
  y <- sim_df$y
  
  ###############################################
  # train/test split
  ###########################################
  n_train <- floor(0.8 * n)
  train_idx <- sample(seq_len(n), n_train)
  test_idx <- setdiff(seq_len(n), train_idx)
  
  train_X <- X[train_idx, , drop = FALSE]
  test_X <- X[test_idx, , drop = FALSE]
  train_y <- y[train_idx]
  test_y <- y[test_idx]

  
  mse_vec <- model_fits(train_X, train_y, test_X, test_y, train)
  return(mse_vec)
}

# Simulation study over varying p
sim_study <- function(n, J, p_seq, s_fixed = 10, train, SNR) {
  num_p <- length(p_seq)
  
  lasso_mse <-  matrix(NA, nrow = J, ncol= num_p)
  ridge_mse <-  matrix(NA, nrow = J, ncol= num_p)
  ridgeless_mse <-  matrix(NA, nrow = J, ncol= num_p)
  
  if (train == 1) {
    train_lasso_mse <-  matrix(NA, nrow = J, ncol= num_p)
    train_ridge_mse <-  matrix(NA, nrow = J, ncol= num_p)
    train_ridgeless_mse <-  matrix(NA, nrow = J, ncol= num_p)
  }

  
  for (i in 1:num_p) {
    cat("Running p =", p_seq[i], "of", p_seq[num_p], "\n")
    for(j in 1:J) {
      fit <- compute_simulation(n = n, p = p_seq[i], s = s_fixed, seed = 100 + j, train, SNR)
      lasso_mse[j, i] <- fit$lasso_mse
      ridge_mse[j, i] <- fit$ridge_mse
      ridgeless_mse[j, i] <- fit$ridgeless_mse
      if (train == 1) {
        train_lasso_mse[j, i] <- fit$train_lasso_mse
        train_ridge_mse[j, i] <- fit$train_ridge_mse
        train_ridgeless_mse[j, i] <- fit$train_ridgeless_mse
      }
      
    }
  }
  
  if(train == 1) {
    return(list(
      lasso_mse = lasso_mse, ridge_mse = ridge_mse, ridgeless_mse = ridgeless_mse,
      train_lasso_mse = train_lasso_mse, train_ridge_mse = train_ridge_mse,
      train_ridgeless_mse = train_ridgeless_mse
    ))
  }
  
  
  return(list(lasso_mse = lasso_mse, ridge_mse = ridge_mse,
              ridgeless_mse = ridgeless_mse))
}

##########################################################################
# calls the sim study and returns plot of test mse, optiol to include plot of train mse
complete_run <- function(n, J, p_seq, s_fixed, train, SNR) {
  # Run simulation
  results <- sim_study(n, J, p_seq, s_fixed, train, SNR)
  
  # Compute means and SDs
  lasso_mean <- colMeans(results$lasso_mse)
  ridge_mean <- colMeans(results$ridge_mse)
  ridgeless_mean <- colMeans(results$ridgeless_mse)
  
  lasso_sd <- apply(results$lasso_mse, 2, sd)
  ridge_sd <- apply(results$ridge_mse, 2, sd)
  ridgeless_sd <- apply(results$ridgeless_mse, 2, sd)
  
  test_mse_df <- data.frame(
    p = rep(p_seq, times = 3),
    Method = rep(c("Lasso", "Ridge", "Ridgeless"), each = length(p_seq)),
    MSE = c(lasso_mean, ridge_mean, ridgeless_mean),
    SD = c(lasso_sd, ridge_sd,  ridgeless_sd)
  )
  
  # Plot

  
  if (train == 1) {
    # Compute means and SDs
    train_lasso_mean <- colMeans(results$train_lasso_mse)
    train_ridge_mean <- colMeans(results$train_ridge_mse)
    train_ridgeless_mean <- colMeans(results$train_ridgeless_mse)
    
    train_lasso_sd <- apply(results$train_lasso_mse, 2, sd)
    train_ridge_sd <- apply(results$train_ridge_mse, 2, sd)
    train_ridgeless_sd <- apply(results$train_ridgeless_mse, 2, sd)
    
    train_mse_df <- data.frame(
      p = rep(p_seq, times = 3),
      Method = rep(c("Lasso", "Ridge", "Ridgeless"), each = length(p_seq)),
      MSE = c(train_lasso_mean, train_ridge_mean, train_ridgeless_mean),
      SD = c(train_lasso_sd, train_ridge_sd,  train_ridgeless_sd)
    )
    
    
    return(list(test_mse_df = test_mse_df, train_mse_df = train_mse_df))
  }
  return(list(test_mse_df= test_mse_df))
}



####################################################################################
# *********************************************************************************
#####################################################################################

# end of functions

# Parameters
n <- 250
J <- 30
p_seq <- seq(10, 600, by = 10)
s_fixed <-  50  # small sparsity to show double descent

plot_dfs <- complete_run(n, J, p_seq, s_fixed, train = 1, SNR = 5)

test_mse_df <- plot_dfs$test_mse_df

# potentially wont exist
train_mse_df <- plot_dfs$train_mse_df

# plot
test_plot <- ggplot(test_mse_df, aes(x = p, y = MSE, color = Method)) +
  geom_line(size = 1.2) +
  geom_point(size = 1.5) +
  geom_errorbar(aes(ymin = MSE - SD, ymax = MSE + SD), width = 5, alpha = 0.3) +
  labs(
    title = "Double Descent Simulation",
    subtitle = paste0("n = 250, J = 30, SNR = 5"),
    x = "Number of Predictors (p)",
    y = "Test MSE"
  ) +
  coord_cartesian(ylim = c(0, 10))  +
  theme_minimal(base_size = 14) +
  theme(legend.position = "top", plot.title = element_text(face = "bold"))

test_plot

# Plot
train_plot <- ggplot(train_mse_df, aes(x = p, y = MSE, color = Method)) +
  geom_line(size = 1.2) +
  geom_point(size = 1.5) +
  geom_errorbar(aes(ymin = MSE - SD, ymax = MSE + SD), width = 5, alpha = 0.3) +
  labs(
    title = "Double Descent Simulation",
    subtitle = paste0("n = 250, J = 30, SNR = 1"),
    x = "Number of Predictors (p)",
    y = "Train MSE"
  ) +
  coord_cartesian(ylim = c(0, 2))  +
  theme_minimal(base_size = 14) +
  theme(legend.position = "top", plot.title = element_text(face = "bold"))

train_plot
