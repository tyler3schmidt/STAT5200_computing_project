library(glmnet)
library(ggplot2)

##########################################################################
# function: generate_data
generate_data <- function(n, p, s){
  
  # covariates
  X <- matrix(rnorm(n * p, mean = 0, sd = 1), nrow = n, ncol = p)
  
  # select s important features
  important_idx <- sample(seq_len(p), size = min(s, p), replace = FALSE)
  
  # beta coefficients
  betas <- rep(0, p)
  betas[important_idx] <- runif(length(important_idx), min = 1, max = 3) * sample(c(-1,1), length(important_idx), replace = TRUE)
  
  # outcome with Gaussian noise
  y <- X %*% betas + rnorm(n, mean = 0, sd = 0.1)
  
  sim_data <- data.frame(y = as.vector(y), X)
  colnames(sim_data) <- c("y", paste0("X", 1:p))
  
  return(list(sim_data = sim_data, betas = betas))
}

##########################################################################
# Ridgeless OLS using pseudoinverse
ridgeless_fit <- function(X, y, tol = 1e-8) {
  n <- nrow(X)
  p <- ncol(X)
  
  if (p < n) {
    # Underparameterized: use standard OLS formula
    beta_hat <- solve(t(X) %*% X, t(X) %*% y)
  } else {
    # Overparameterized or p >= n: use Moore-Penrose pseudoinverse
    svd_X <- svd(X)
    D_plus <- diag(ifelse(svd_X$d > tol, 1 / svd_X$d, 0))
    X_pinv <- svd_X$v %*% D_plus %*% t(svd_X$u)
    beta_hat <- X_pinv %*% y
  }
  
  return(beta_hat)
}

##########################################################################
# Fit models and compute test MSE
model_fits <- function(train_X, train_y, test_X, test_y) {
  # Lasso
  fit_lasso <- cv.glmnet(train_X, train_y, alpha = 1, nfolds = 5)
  lasso_pred <- as.numeric(predict(fit_lasso, test_X, s = "lambda.min"))
  lasso_mse <- mean((test_y - lasso_pred)^2)
  
  # Ridge
  fit_ridge <- cv.glmnet(train_X, train_y, alpha = 0, nfolds = 5)
  ridge_pred <- as.numeric(predict(fit_ridge, test_X, s = "lambda.min"))
  ridge_mse <- mean((test_y - ridge_pred)^2)
  
  # Elastic Net (alpha = 0.5)
  fit_en <- cv.glmnet(train_X, train_y, alpha = 0.5, nfolds = 5)
  en_pred <- as.numeric(predict(fit_en, test_X, s = "lambda.min"))
  en_mse <- mean((test_y - en_pred)^2)
  
  # Ridgeless OLS
  beta_rl <- ridgeless_fit(train_X, train_y)
  rl_pred <- as.numeric(test_X %*% beta_rl)
  ridgeless_mse <- mean((test_y - rl_pred)^2)
  
  return(list(lasso_mse = lasso_mse, ridge_mse = ridge_mse, 
              elasticnet_mse = en_mse, ridgeless_mse = ridgeless_mse))
  
}
##########################################################################
# Compute single simulation
compute_simulation <- function(n, p, s, seed) {
  set.seed(seed)
  dat <- generate_data(n, p, s)
  sim_df <- dat$sim_data
  X <- as.matrix(sim_df[, -1])
  y <- sim_df$y
  
  # train/test split
  n_train <- floor(0.8 * n)
  train_idx <- sample(seq_len(n), n_train)
  test_idx <- setdiff(seq_len(n), train_idx)
  
  train_X <- X[train_idx, , drop = FALSE]
  test_X <- X[test_idx, , drop = FALSE]
  train_y <- y[train_idx]
  test_y <- y[test_idx]
  
  mse_vec <- model_fits(train_X, train_y, test_X, test_y)
  return(mse_vec)
}

# Simulation study over varying p
sim_study <- function(n, J, p_seq, s_fixed = 10) {
  num_p <- length(p_seq)
  
  lasso_mse <-  matrix(NA, nrow = J, ncol= num_p)
  ridge_mse <-  matrix(NA, nrow = J, ncol= num_p)
  en_mse <- matrix(NA, nrow = J, ncol = num_p)
  ridgeless_mse <-  matrix(NA, nrow = J, ncol= num_p)
  
  for (i in 1:num_p) {
    cat("Running p =", p_seq[i], "of", num_p, "\n")
    for(j in 1:J) {
      fit <- compute_simulation(n = n, p = p_seq[i], s = s_fixed, seed = 100 + j)
      lasso_mse[j, i] <- fit$lasso_mse
      ridge_mse[j, i] <- fit$ridge_mse
      en_mse[j, i] <- fit$elasticnet_mse
      ridgeless_mse[j, i] <- fit$ridgeless_mse
    }
  }
  
  return(list(lasso_mse = lasso_mse, ridge_mse = ridge_mse,
              elasticnet_mse = en_mse, ridgeless_mse = ridgeless_mse))
}

##########################################################################
# Parameters
n <- 125
J <- 50
p_seq <- seq(10, 500, by = 5)
s_fixed <- 10  # small sparsity to show double descent

# Run simulation
results <- sim_study(n, J, p_seq, s_fixed)

# Compute means and SDs
lasso_mean <- colMeans(results$lasso_mse)
ridge_mean <- colMeans(results$ridge_mse)
en_mean <- colMeans(results$elasticnet_mse)
ridgeless_mean <- colMeans(results$ridgeless_mse)

lasso_sd <- apply(results$lasso_mse, 2, sd)
ridge_sd <- apply(results$ridge_mse, 2, sd)
en_sd <- apply(results$elasticnet_mse, 2, sd)
ridgeless_sd <- apply(results$ridgeless_mse, 2, sd)

mse_df <- data.frame(
  p = rep(p_seq, times = 4),
  Method = rep(c("Lasso", "Ridge", "ElasticNet", "Ridgeless"), each = length(p_seq)),
  MSE = c(lasso_mean, ridge_mean, en_mean, ridgeless_mean),
  SD = c(lasso_sd, ridge_sd, en_sd, ridgeless_sd)
)

# Plot
ggplot(mse_df, aes(x = p, y = MSE, color = Method)) +
  geom_line(size = 1.2) +
  geom_point(size = 1.5) +
  geom_errorbar(aes(ymin = MSE - SD, ymax = MSE + SD), width = 5, alpha = 0.3) +
  labs(
    title = "Double Descent Simulation",
    subtitle = paste0("n = ", n, ", ", J, " replications per setting"),
    x = "Number of Predictors (p)",
    y = "Test MSE"
  ) +
  coord_cartesian(ylim = c(0, 0.5))  +
  theme_minimal(base_size = 14) +
  theme(legend.position = "top", plot.title = element_text(face = "bold"))

