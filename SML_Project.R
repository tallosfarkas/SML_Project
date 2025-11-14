#### SML Project ####

### Packages ###################################################################
required_packages <- c(
  "xts","zoo","frenchdata","mistr","ggplot2","quantmod","dplyr","purrr",
  "lubridate","stringr","tidyr","stats","gt","PerformanceAnalytics","kableExtra",
  "knitr","broom","tibble","tidyquant","tidymodels","glmnet","ranger","doParallel",
  "GGally","readxl","pROC","caret"
)

installed <- rownames(installed.packages())
missing <- setdiff(required_packages, installed)
if (length(missing) > 0) install.packages(missing, dependencies = TRUE)
invisible(lapply(required_packages, library, character.only = TRUE))

### Helper: safe FRED getter ###################################################
get_fred_safe <- function(sym, from) {
  tryCatch(
    quantmod::getSymbols(sym, src = "FRED", from = from, auto.assign = FALSE),
    error = function(e) stop(sprintf("FRED download failed for %s: %s", sym, e$message))
  )
}

### Parameters #################################################################
start <- as.Date("1950-01-01")
end   <- Sys.Date()


### Market data (S&P 500) ######################################################
sp500_xts <- quantmod::getSymbols("^GSPC", from = start, to = end,
                                  src = "yahoo", auto.assign = FALSE)

# first trading day of each month
sp500_first <- do.call(rbind, lapply(split(sp500_xts, f = "months"), head, 1))

sp500_df <- data.frame(date = index(sp500_first), coredata(sp500_first)) %>%
  transmute(
    year_month = floor_date(date, "month"),
    price      = GSPC.Adjusted,
    volume     = GSPC.Volume
  ) %>%
  arrange(year_month) %>%
  mutate(
    price_lag1  = dplyr::lag(price, 1),
    return      = price / price_lag1 - 1,
    lag1_return = dplyr::lag(return, 1),
    lag2_return = dplyr::lag(return, 2),
    lag3_return = dplyr::lag(return, 3),
    lag4_return = dplyr::lag(return, 4),
    lag5_return = dplyr::lag(return, 5),
    volume_change = volume / dplyr::lag(volume, 1) - 1,
    volume_change_lag = dplyr::lag(volume_change,1),
    UP_DOWN     = if_else(price >= price_lag1, "Up", "Down")
  )

### Macro data (FRED + NBER) ###################################################
cpi_xts <- get_fred_safe("CPIAUCSL", from = start)
fed_xts <- get_fred_safe("FEDFUNDS", from = start)

nber_tbl <- tq_get("USREC", get = "economic.data", from = start)
nber_xts <- xts(nber_tbl[["price"]], order.by = nber_tbl$date)

macro_xts <- merge(cpi_xts, fed_xts, nber_xts)
colnames(macro_xts) <- c("CPI","FedFundsRate","NBER")

macro_df <- data.frame(date = index(macro_xts), coredata(macro_xts)) %>%
  mutate(
    CPI_lag          = dplyr::lag(CPI, 1),
    FedFundsRate_lag = dplyr::lag(FedFundsRate, 1),
    NBER_lag         = dplyr::lag(NBER, 1),
    year_month       = floor_date(date, "month")
  ) %>%
  select(-date) %>%
  tidyr::drop_na()

macro_df <- macro_df %>%
  select(year_month, CPI_lag, FedFundsRate_lag, NBER_lag)


### VIX and Daily News Sentiment ##########################################

# 1. Volatility Index (VIX)
vix_xts <- getSymbols("^VIX", from = start, to = end,
                      src = "yahoo", auto.assign = FALSE)

vix_df <- data.frame(date = index(vix_xts), coredata(vix_xts)) %>%
  transmute(
    year_month = floor_date(date, "month"),
    VIX = VIX.Adjusted
  ) %>%
  group_by(year_month) %>%
  summarise(VIX = mean(VIX, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(
    VIX_lag = dplyr::lag(VIX, 1),                   # one-month lag
    VIX_change_lag = (VIX_lag / dplyr::lag(VIX_lag, 1)) - 1  # change from t-2 to t-1
  ) %>%
  select(year_month, VIX_change_lag) %>%
  drop_na()

# 2. Daily News Sentiment Index (FRB San Francisco)

dnsi_raw <- read_excel("news_sentiment_data.xlsx")


dnsi_df <- dnsi_raw %>%
  rename(date = 1, DNSI = 2) %>%
  mutate(
    date = as.Date(date),
    year_month = floor_date(date, "month")
  ) %>%
  group_by(year_month) %>%
  summarise(DNSI = mean(DNSI, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(
    DNSI_lag = dplyr::lag(DNSI, 1),
    DNSI_change_lag = DNSI_lag - dplyr::lag(DNSI_lag, 1)
  ) %>%
  select(year_month, DNSI_change_lag) %>%
  drop_na()



################################################################################
# --- Merge everything for final data set  -------------------------------------
################################################################################

final_data_df <- macro_df %>%
  left_join(sp500_df, by = "year_month") %>%
  left_join(vix_df,   by = "year_month") %>%
  left_join(dnsi_df,  by = "year_month") %>%
  select(-price, -price_lag1, -volume,-volume_change, -return) %>%
  tidyr::drop_na()

# Ensure response factor
final_data_df <- final_data_df %>%
  mutate(UP_DOWN = factor(UP_DOWN, levels = c("Down", "Up")))

final_data_df$Y <- ifelse(final_data_df$UP_DOWN == "Up", 1, 0)

# Create the plot and store it
pair_plot <- GGally::ggpairs(
  final_data_df %>% select(-year_month, -UP_DOWN)
)

# Show the plot (optional)
pair_plot

# Save as PNG
ggsave(
  filename = "pair_plot.png",
  plot = pair_plot,
  width = 12,
  height = 12,
  dpi = 300
)


################################################################################
# --- Analyse Imbalance in the underlying data set  ----------------------------
################################################################################

imbalance_df <- final_data_df %>%
  mutate(UP_DOWN = factor(UP_DOWN, levels = c("Down","Up")))

# Create imbalance plot and store it
imb_plot <- ggplot(imbalance_df, aes(UP_DOWN)) +
  geom_bar(fill = "steelblue", alpha = 0.8) +
  geom_text(
    stat = "count",
    aes(label = after_stat(count)),
    vjust = -0.5,
    size = 5
  ) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Class Imbalance in Monthly S&P 500 Direction",
    x = "Market Direction",
    y = "Frequency"
  )

# Save plot as PNG
ggsave("imbalance_plot.png",
       plot = imb_plot,
       width = 6,
       height = 6,
       dpi = 300)

################################################################################
################################################################################

### --- Logistic Regression with Elastic Net --- ###################

# --- Prepare data -------------------------------------------------------------

# We remove non-predictor variables and create interaction terms between the key
# macro and news variables

# Create all lagged + interaction features BEFORE the split
data_model_enet <- final_data_df %>%
  select(-UP_DOWN, -year_month) %>%
  mutate(
    # Interaction terms
    DNSI_VIX_lag       = DNSI_change_lag * VIX_change_lag,
    DNSI_FedFunds_lag  = DNSI_change_lag * FedFundsRate_lag,
    VIX_CPI_lag        = VIX_change_lag * CPI_lag,
    DNSI_NBER_lag      = DNSI_change_lag * NBER_lag
  ) %>%
  tidyr::drop_na()

# --- Chronological 70/30 split -----------------------------------------------

# We decided to use the first 30% of our time series data set for the calibrtion of the
# hyperparamter and coefficient estimation.
# Since it is a time series, we cannot randomize.

n <- nrow(data_model_enet)
split_point <- floor(0.7 * n)

train_df_enet <- data_model_enet[1:split_point, ]
test_df_enet  <- data_model_enet[(split_point + 1):n, ]

# --- Matrices for glmnet ------------------------------------------------------
x_train <- model.matrix(Y ~ . - 1, data = train_df_enet)
y_train <- train_df_enet$Y

x_test  <- model.matrix(Y ~ . - 1, data = test_df_enet)
y_test  <- test_df_enet$Y


################################################################################
# --- Time-Series Cross-Validation for Alpha and Lambda ------------------------
################################################################################

# For the CV, we test alphas (mixture of lasso and ridge) between [0,1], with using cv, for each
# alpha to define for each alphas the optimal lambda.
# createTimeSlices() creates rolling time windows

### --- Time-Series CV Folds ---------------------------------------------------

window_length <- 60   # 5 years of monthly data
horizon       <- 12   # 1-year validation

folds <- createTimeSlices(
  1:nrow(x_train),
  initialWindow = window_length,
  horizon       = horizon,
  fixedWindow   = TRUE
)

### --- Grid Search ------------------------------------------------------------
# We'll perform a small grid search for alpha and lambda

# Lambda grid (log-spaced)
lambda_grid <- exp(seq(log(1e-4), log(10), length.out = 40))

# Alpha grid: ridge -> lasso spectrum
alpha_grid  <- seq(0, 1, by = 0.1)

### --- Storage ----------------------------------------------------------------

cv_results <- data.frame(alpha = numeric(), lambda = numeric(),
                         Accuracy = numeric(), AUC = numeric())


### --- Time-Series CV ---------------------------------------------------------

# The loop performes "nested" time-series cross validation

# For each alpha, run several expanding window folds and within each fold find
# the best lambda

# Then fit the model and record the Accuracy and AUC on the validation slice

# Average the fold metrics and store them in cv_results.
# This gives us a grid of (alpha, lambda, performance) combinations.

for (a in alpha_grid) {
  cat("Testing alpha =", a, "\n")

  auc_mat <- matrix(NA, nrow = length(folds$train), ncol = length(lambda_grid))
  acc_mat <- matrix(NA, nrow = length(folds$train), ncol = length(lambda_grid))

  for (i in seq_along(folds$train)) {

    tr  <- folds$train[[i]]
    val <- folds$test[[i]]

    # Skip folds with only one class
    if (length(unique(y_train[val])) < 2) next

    x_tr <- x_train[tr, ]
    y_tr <- y_train[tr]

    x_val <- x_train[val, ]
    y_val <- y_train[val]

    # Fit path for all lambdas (TS-safe)
    fit <- glmnet(
      x = x_tr,
      y = y_tr,
      family = "binomial",
      alpha = a,
      lambda = lambda_grid
    )

    # Predict whole lambda path at once (efficient)
    pred_prob_mat <- predict(fit, newx = x_val, type = "response")

    for (j in seq_along(lambda_grid)) {
      p <- pred_prob_mat[, j]
      pred_class <- ifelse(p > 0.5, 1, 0)

      acc_mat[i, j] <- mean(pred_class == y_val)

      auc_mat[i, j] <- tryCatch({
        as.numeric(auc(roc(y_val, p, quiet = TRUE)))
      }, error = function(e) NA)
    }
  }

  mean_auc <- colMeans(auc_mat, na.rm = TRUE)
  best_idx <- which.max(mean_auc)

  cv_results <- rbind(cv_results, data.frame(
    alpha = a,
    lambda = lambda_grid[best_idx],
    AUC = mean_auc[best_idx],
    Accuracy = colMeans(acc_mat, na.rm = TRUE)[best_idx]
  ))
}


# --- Select best alpha & lambda -----------------------------------------------

# We then pick the alpha & lambda combination that gave the the highest mean AUC
# (Area under the curve) i.e. the area under the ROC curve.

best_row     <- cv_results[which.max(cv_results$AUC), ]
best_alpha   <- best_row$alpha
best_lambda  <- best_row$lambda

cat("\nBest alpha:", best_alpha, "\nBest lambda:", best_lambda, "\n")

################################################################################
# --- Refit on full training data ---------------------------------------------
################################################################################

# The last step of the model fitting is the training of the final model on the
# entire training sample to receive the final coefficients using the chosen alpha
# and lambda

final_model <- glmnet(
  x = x_train,
  y = y_train,
  family = "binomial",
  alpha = best_alpha,
  lambda = best_lambda
)

# --- Threshold Tuning due to imbalances in the underlying data classes --------

# Train-set probability predictions
train_prob <- predict(final_model, newx = x_train, type = "response")

# ROC
roc_train <- roc(y_train, as.numeric(train_prob), quiet = TRUE)

opt_thresh <- coords(roc_train, "best", ret = "threshold")

opt_thresh_enet <- as.numeric(opt_thresh)

cat("Optimal threshold (training):", round(opt_thresh_enet, 4), "\n")

################################################################################
# --- Out-of-sample evaluation -------------------------------------------------
################################################################################

# For the OoS prediction, we use our final_model to predict the probabilities on the
# unseen test data. Finally, we just convert probabilities to a binary class predictions


pred_prob <- predict(final_model, newx = x_test, type = "response")
pred_class <- ifelse(pred_prob > opt_thresh_enet, 1, 0)

conf_mat <- table(Predicted = pred_class, Actual = y_test)
accuracy <- mean(pred_class == y_test)

# We can evaluate the classifications with a confusion matrix, accuracy and
# ROC curve

cat("\nConfusion Matrix:\n")
print(conf_mat)
cat("\nOut-of-sample Accuracy:", round(accuracy, 3), "\n")

roc_obj_enet <- roc(y_test, as.numeric(pred_prob))
cat("Out-of-sample AUC:", round(auc(roc_obj_enet), 3), "\n")

png("roc_enet.png", width = 800, height = 600, res = 120)
plot(roc_obj_enet, col = "blue", main = "ROC Curve - Out-of-sample (Elastic Net)")
dev.off()

# --- Extract coefficients at best lambda and alpha ----------------------------

# Because coefficients are shrunk and selected jointly through the penalty parameter,
# classical inference is not valid. Therefore, significance testing is replaced by
# OoS predictive performance metrics.


# Extract sparse coefficient matrix
coef_mat <- coef(final_model)

# Convert to a tidy data frame manually
coef_df <- data.frame(
  term = rownames(coef_mat),
  estimate = as.numeric(coef_mat)
)
coef_df$odds_ratio <- exp(coef_df$estimate)

coef_df <- coef_df %>%
  filter(term != "(Intercept)") %>%
  arrange(desc(abs(estimate)))

print(coef_df)

ggplot(coef_df, aes(x = reorder(term, odds_ratio), y = odds_ratio)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  theme_minimal(base_size = 12) +
  labs(
    title = "Elastic Net Coefficients",
    x = "Predictor",
    y = "Coefficient (Chnange in the Odds Ratio)"
  )

ggsave(
  filename = "enet_coeffs.png",
  plot = last_plot(),
  width = 7, height = 6, dpi = 300
)


### --- Random Forest Classification --- ###################

# What to tune:

# mtry: Features tried per split
# min.node.size: Leaf Size --> Controls tree depth/overfit
# Max.depth
# Sample Fraction

# What is fixed:

# num.tree = 1000
# splitrule = gini

# --- Create model matrices for R. Forrest and Gradient Boosting ---------------

# --- Prepare data (no interactions) -------------------------------------------
final_data_df$Y <- ifelse(final_data_df$UP_DOWN == "Up", 1, 0)

# Remove unused columns but keep the original lagged predictors
data_model_rf <- final_data_df %>%
  select(-UP_DOWN, -year_month) %>%
  tidyr::drop_na()

# --- Chronological 70/30 split -----------------------------------------------
n <- nrow(data_model_rf)
split_point <- floor(0.7 * n)

rf_train_df <- data_model_rf[1:split_point, ]
rf_test_df  <- data_model_rf[(split_point + 1):n, ]

# --- Factorize target for classification -------------------------------------
rf_train_df$Y <- factor(rf_train_df$Y, levels = c(0, 1))
rf_test_df$Y  <- factor(rf_test_df$Y, levels = c(0, 1))

################################################################################
# --- Time-Series Cross-Validation for Hyperparameter --------------------------
################################################################################

# we can use the folds from the elastic net tuning here again
# First create the hyperparameter tuning grid:

p <- ncol(rf_train_df) - 1 # not counting y
param_grid <- expand.grid(
  mtry            = unique(pmax(1, round(c(sqrt(p), p/3, p/2, p/1.5, p)))),
  min.node.size   = c(2, 5, 10, 20),
  max.depth       = c(NA, 5, 10, 15, 20),
  sample.fraction = c(0.6, 0.7, 0.8, 0.9)
)
nrow(param_grid)

# --- Cross-validation loop as before ------------------------------------------

rf_results <- data.frame(
  mtry = numeric(),
  min.node.size = numeric(),
  max.depth = numeric(),
  sample.fraction = numeric(),
  Accuracy = numeric(),
  AUC = numeric()
)

for (r in 1:nrow(param_grid)) {
  pset <- param_grid[r, ]
  cat("\nTesting combination:", r, "of", nrow(param_grid), "\n")
  print(pset)

  accuracies <- c()
  aucs <- c()

  for (i in seq_along(folds$train)) {
    train_idx <- folds$train[[i]]
    val_idx   <- folds$test[[i]]

    if (length(unique(y_train[val_idx])) < 2) next  # skip degenerate folds

    # skip folds with only one class
    if (length(unique(rf_train_df$Y[val_idx])) < 2) next

    df_train <- rf_train_df[train_idx, ]
    df_val   <- rf_train_df[val_idx, ]

    # Conditional argument for max.depth
    if (is.na(pset$max.depth)) {
      rf_model <- ranger(
        Y ~ .,
        data = df_train,
        num.trees = 1000,
        mtry = pset$mtry,
        min.node.size = pset$min.node.size,
        sample.fraction = pset$sample.fraction,
        splitrule = "gini",
        probability = TRUE,
        importance = "impurity",
        seed = 123
      )
    } else {
      rf_model <- ranger(
        Y ~ .,
        data = df_train,
        num.trees = 1000,
        mtry = pset$mtry,
        min.node.size = pset$min.node.size,
        max.depth = pset$max.depth,
        sample.fraction = pset$sample.fraction,
        splitrule = "gini",
        probability = TRUE,
        importance = "impurity",
        seed = 123
      )
    }

    # Predict probabilities and classes
    pred_prob <- predict(rf_model, data = df_val)$predictions[, 2]
    pred_class <- ifelse(pred_prob > 0.5, 1, 0)

    acc <- mean(pred_class == df_val$Y)
    auc_val <- tryCatch({
      roc_obj <- roc(df_val$Y, pred_prob, quiet = TRUE)
      as.numeric(auc(roc_obj))
    }, error = function(e) NA)

    accuracies <- c(accuracies, acc)
    aucs <- c(aucs, auc_val)
  }

  rf_results <- rbind(rf_results, data.frame(
    mtry = pset$mtry,
    min.node.size = pset$min.node.size,
    max.depth = pset$max.depth,
    sample.fraction = pset$sample.fraction,
    Accuracy = mean(accuracies, na.rm = TRUE),
    AUC = mean(aucs, na.rm = TRUE)
  ))
}

# --- Select best hyperparameter combo ----------------------------------------
best_row <- rf_results[which.max(rf_results$AUC), ]
best_params <- best_row[1:4]

cat("\nBest Parameters:\n")
print(best_params)

# If max.depth = NA it means, we should not restrict the tree depth in our model

################################################################################
# --- Refit on full training data ---------------------------------------------
################################################################################

if (is.na(best_params$max.depth)) {
  final_rf <- ranger(
    Y ~ .,
    data = rf_train_df,
    num.trees = 1000,
    mtry = best_params$mtry,
    min.node.size = best_params$min.node.size,
    sample.fraction = best_params$sample.fraction,
    splitrule = "gini",
    probability = TRUE,
    importance = "impurity",
    seed = 123
  )
} else {
  final_rf <- ranger(
    Y ~ .,
    data = rf_train_df,
    num.trees = 1000,
    mtry = best_params$mtry,
    min.node.size = best_params$min.node.size,
    max.depth = best_params$max.depth,
    sample.fraction = best_params$sample.fraction,
    splitrule = "gini",
    probability = TRUE,
    importance = "impurity",
    seed = 123
  )
}



roc_rf <- roc(rf_train_df$Y, predict(final_rf, data = rf_train_df)$predictions[, 2], quiet = TRUE)
opt_thresh_rf <- coords(roc_rf, "best", ret = "threshold")
opt_thresh_rf <- as.numeric(opt_thresh_rf)

################################################################################
# --- Out-of-Sample Evaluation -------------------------------------------------
################################################################################

pred_prob_rf <- predict(final_rf, data = rf_test_df)$predictions[, "1"]
pred_class_rf <- ifelse(pred_prob_rf > opt_thresh_rf, 1, 0)

conf_mat_rf <- table(Predicted = pred_class_rf, Actual = rf_test_df$Y)
accuracy_rf <- mean(pred_class_rf == as.numeric(as.character(rf_test_df$Y)))


cat("\nConfusion Matrix (RF):\n")
print(conf_mat_rf)
cat("\nOut-of-sample Accuracy (RF):", round(accuracy_rf, 3), "\n")

# --- ROC/AUC plot -------------------------------------------------------------
if (length(unique(rf_test_df$Y)) == 2) {
  roc_obj_rf_tuned <- roc(rf_test_df$Y, pred_prob_rf)
  cat("\nTest AUC:", round(auc(roc_obj_rf_tuned), 3), "\n")
  png("roc_rf.png", width = 800, height = 600, res = 120)
  plot(roc_obj_rf_tuned, col = "darkgreen",
       main = "ROC Curve - Tuned Random Forest")
  dev.off()
} else {
  cat("\nROC skipped: test set has only one class.\n")
}


# --- Feature Importance -------------------------------------------------------
importance_df <- data.frame(
  Variable = names(final_rf$variable.importance),
  Importance = final_rf$variable.importance
) %>%
  arrange(desc(Importance))

importance_df %>%
  top_n(20, Importance) %>%
  ggplot(aes(x = reorder(Variable, Importance), y = Importance)) +
  geom_col() +
  coord_flip() +
  labs(
    title = "Top 20 Features — Impurity Importance",
    x = "Feature",
    y = "Importance"
  ) +
  theme_minimal()

# --- Interaction strength

install.packages("randomForestExplainer")
library(randomForestExplainer)


# Prepare random forest explainer data
rf_expl <- explain_forest(final_rf, interactions = TRUE)


################################################################################
################################################################################

# We need the 'gbm' package
if (!"gbm" %in% installed) install.packages("gbm")
library(gbm)

# --- Parameter grid to tune for GBM ---
param_grid_gbm <- expand.grid(
  n.trees = c(100, 200, 300),
  interaction.depth = c(1, 2, 3), # As in slides07 (1=additive, >1=interactions)
  shrinkage = c(0.01, 0.1)     # Learning rate
)

# --- Rolling-window tuning (for gbm) ---
# Use the same window parameters
n_train <- nrow(rf_train_df)
window_size <- floor(0.7 * n_train)
horizon <- 12

results_gbm <- data.frame()
eps <- 1e-9 # for log-loss stability

print("Starting GBM time-series tuning...")
set.seed(123)
for (i in seq(window_size, n_train - horizon, by = horizon)) {

  # Get window data (dataframes)
  train_window <- rf_train_df[1:i, ]
  test_window <- rf_train_df[(i + 1):(i + horizon), ]

  # --- Prepare data for gbm ---
  # gbm for "bernoulli" requires a 0/1 numeric target
  train_window$Y_num <- as.numeric(as.character(train_window$Y))
  test_window$Y_num <- as.numeric(as.character(test_window$Y))

  # Skip invalid test windows
  if (nrow(test_window) == 0 || length(unique(test_window$Y_num)) < 2) next

  # Create a formula that uses Y_num and excludes the factor Y
  predictors <- names(train_window)[!names(train_window) %in% c("Y", "Y_num")]
  gbm_formula <- as.formula(paste("Y_num ~", paste(predictors, collapse = " + ")))

  for (j in 1:nrow(param_grid_gbm)) {
    p <- param_grid_gbm[j, ]

    # Fit the gbm model
    gbm_model <- gbm(
      gbm_formula,
      data = train_window,
      distribution = "bernoulli", # for 0/1 classification
      n.trees = p$n.trees,
      interaction.depth = p$interaction.depth,
      shrinkage = p$shrinkage,
      n.minobsinnode = 10 # As in slides07, a good default
    )

    # Predict on the validation window
    preds <- predict(gbm_model,
                     newdata = test_window,
                     n.trees = p$n.trees,
                     type = "response")
    y_true_window <- test_window$Y_num

    # --- Log-loss (cross-entropy) ---
    logloss <- -mean(y_true_window * log(preds + eps) + (1 - y_true_window) * log(1 - preds + eps))
    auc_val <- tryCatch(as.numeric(auc(y_true_window, preds, quiet = TRUE)), error = function(e) NA)

    results_gbm <- rbind(
      results_gbm,
      data.frame(
        n.trees = p$n.trees,
        interaction.depth = p$interaction.depth,
        shrinkage = p$shrinkage,
        LogLoss = logloss,
        AUC = auc_val
      )
    )
  }
}
print("GBM tuning complete.")

# --- Aggregate across rolling folds ---
tuning_summary_gbm <- results_gbm %>%
  group_by(n.trees, interaction.depth, shrinkage) %>%
  summarise(
    mean_LogLoss = mean(LogLoss, na.rm = TRUE),
    mean_AUC = mean(AUC, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(mean_LogLoss) # lower log-loss is better

cat("Top GBM tuning results (lowest log-loss):\n")
print(head(tuning_summary_gbm, 5))

# --- Select best parameters ---
best_params_gbm <- tuning_summary_gbm[1, ]
cat("\nSelected GBM parameters:\n")
print(best_params_gbm)

# --- Fit final model on full training data ---
# Prepare full rf_train_df for gbm
train_df_gbm <- rf_train_df
train_df_gbm$Y_num <- as.numeric(as.character(train_df_gbm$Y))
predictors <- names(train_df_gbm)[!names(train_df_gbm) %in% c("Y", "Y_num")]
gbm_formula <- as.formula(paste("Y_num ~", paste(predictors, collapse = " + ")))

set.seed(123)
gbm_final <- gbm(
  gbm_formula,
  data = train_df_gbm,
  distribution = "bernoulli",
  n.trees = best_params_gbm$n.trees,
  interaction.depth = best_params_gbm$interaction.depth,
  shrinkage = best_params_gbm$shrinkage,
  n.minobsinnode = 10
)

# --- Evaluate on hold-out test set ---
# Prepare test_df for gbm
test_df_gbm <- rf_test_df
test_df_gbm$Y_num <- as.numeric(as.character(test_df_gbm$Y))
y_true_test_gbm <- test_df_gbm$Y_num # Use the numeric 0/1 Y

gbm_pred_prob <- predict(gbm_final,
                         newdata = test_df_gbm,
                         n.trees = best_params_gbm$n.trees,
                         type = "response")
gbm_pred_class <- ifelse(gbm_pred_prob > 0.5, 1, 0)

# --- Metrics ---
gbm_conf_mat <- table(Predicted = gbm_pred_class, Actual = y_true_test_gbm)
gbm_accuracy <- mean(gbm_pred_class == y_true_test_gbm)
gbm_logloss_test <- -mean(y_true_test_gbm * log(gbm_pred_prob + eps) + (1 - y_true_test_gbm) * log(1 - gbm_pred_prob + eps))

cat("\n--- GBM Final Test Performance ---\n")
cat("\nConfusion Matrix:\n"); print(gbm_conf_mat)
cat("\nTest Accuracy:", round(gbm_accuracy, 3))
cat("\nTest Log-Loss:", round(gbm_logloss_test, 4))

# --- ROC/AUC plot ---
gbm_roc_obj <- roc(y_true_test_gbm, gbm_pred_prob, quiet = TRUE)
cat("\nTest AUC:", round(auc(gbm_roc_obj), 3), "\n")
plot(gbm_roc_obj, col = "darkorange", main = "ROC Curve - Tuned GBM")

# --- Feature Importance ---
cat("\nGBM Feature Importance:\n")
print(summary(gbm_final, plotit = FALSE))

###############################################################################
###############################################################################


################ PLOTS & FINAL SUMMARY #############################

# --- 1. Create Final Objects for Comparison ---

# --- FIX: Define y_true_test (the 0/1 numeric test target) ---
# This ensures the true test set target is available
y_true_test <- as.numeric(as.character(test_df$Y))

# --- Re-calculate all metrics from final models ---
# This makes the script runnable in chunks, as it only needs the final models
print("Re-calculating final metrics...")

# Elastic Net (enet_final was named 'final_model' in your script)
enet_pred_prob <- predict(final_model, newx = x_test, type = "response")
enet_accuracy  <- mean(ifelse(enet_pred_prob > opt_thresh_enet, 1, 0) == y_true_test)
enet_logloss_test <- -mean(y_true_test * log(enet_pred_prob + eps) + (1 - y_true_test) * log(1 - enet_pred_prob + eps))
# We already have 'enet_roc_obj' from the Elastic Net block, but we re-create it here
# just in case the block wasn't run.
enet_roc_obj   <- roc(y_true_test, as.numeric(enet_pred_prob), quiet = TRUE)

# Random Forest (final_rf)
rf_pred_prob_tuned <- predict(final_rf, data = test_df)$predictions[, 2]
rf_accuracy_tuned  <- mean(ifelse(rf_pred_prob_tuned > opt_thresh_rf, 1, 0) == y_true_test)
rf_logloss_test_tuned <- -mean(y_true_test * log(rf_pred_prob_tuned + eps) + (1 - y_true_test) * log(1 - rf_pred_prob_tuned + eps))
# We already have 'roc_obj_rf_tuned', but re-create for safety.
rf_roc_obj_tuned   <- roc(y_true_test, rf_pred_prob_tuned, quiet = TRUE)


length(rf_pred_prob_tuned)
length(y_true_test)

# Gradient Boosting (gbm_final)
# (Must re-create test_df_gbm as it's not saved globally)
test_df_gbm <- test_df
test_df_gbm$Y_num <- as.numeric(as.character(test_df_gbm$Y))

gbm_pred_prob <- predict(gbm_final,
                         newdata = test_df_gbm,
                         n.trees = best_params_gbm$n.trees,
                         type = "response")
gbm_accuracy  <- mean(ifelse(gbm_pred_prob > 0.5, 1, 0) == y_true_test)
gbm_logloss_test <- -mean(y_true_test * log(gbm_pred_prob + eps) + (1 - y_true_test) * log(1 - gbm_pred_prob + eps))
# We already have 'gbm_roc_obj', but re-create for safety.
gbm_roc_obj   <- roc(y_true_test, gbm_pred_prob, quiet = TRUE)


# --- 2. Save EDA Plot ---
print("Saving EDA plot...")
png("ggpairs_plot.png", width = 1200, height = 1000, res = 100)
# Create the subset for a cleaner plot
eda_subset <- data_model %>%
  mutate(UP_DOWN = factor(Y, levels = c(0, 1), labels = c("Down", "Up"))) %>%
  select(
    UP_DOWN,
    lag1_return,
    FedFundsRate_lag,
    VIX_change_lag,
    DNSI_change_lag,
    DNSI_VIX_lag
  )
print(
  GGally::ggpairs(
    eda_subset,
    aes(color = UP_DOWN, alpha = 0.6),
    title = "EDA: Predictor Relationships by Market Direction (Subset)"
  ) +
    theme(text = element_text(size = 8)) +
    scale_color_manual(values = c("Down" = "darkred", "Up" = "darkgreen")) +
    scale_fill_manual(values = c("Down" = "darkred", "Up" = "darkgreen"))
)
dev.off()


# --- 3. Save Tuned Elastic Net ROC Plot ---
print("Saving Elastic Net ROC plot...")
png("enet_roc_plot.png", width = 800, height = 600, res = 100)
plot(enet_roc_obj, col = "blue", main = "ROC Curve - Tuned Elastic Net")
dev.off()

# --- 4. Save Tuned Random Forest ROC & VarImp Plots ---
print("Saving Random Forest plots...")
png("rf_tuned_roc_plot.png", width = 800, height = 600, res = 100)
plot(rf_roc_obj_tuned, col = "darkgreen", main = "ROC Curve - Tuned Random Forest")
dev.off()

png("rf_var_imp_plot.png", width = 800, height = 700, res = 100)
rf_importance_df <- data.frame(
  Variable = names(final_rf$variable.importance),
  Importance = final_rf$variable.importance
) %>% arrange(desc(Importance))
print(
  rf_importance_df %>%
    top_n(10, Importance) %>%
    ggplot(aes(x = reorder(Variable, Importance), y = Importance, fill = Importance)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    labs(
      title = "Top 10 Most Important Predictors (Random Forest)",
      x = "Predictor",
      y = "Mean Decrease in Gini Impurity"
    ) +
    scale_fill_gradient(low = "lightblue", high = "darkblue") +
    theme_minimal() +
    theme(legend.position = "none")
)
dev.off()


# --- 5. Save Tuned GBM ROC Plot ---
print("Saving GBM ROC plot...")
png("gbm_roc_plot.png", width = 800, height = 600, res = 100)
plot(gbm_roc_obj, col = "darkorange", main = "ROC Curve - Tuned GBM")
dev.off()

# --- 6. Save Combined ROC Plot ---
print("Saving combined ROC plot...")
png("combined_roc_plot.png", width = 800, height = 700, res = 100)
plot(enet_roc_obj, col = "blue", main = "Final Model ROC Comparison (Test Set)")
plot(rf_roc_obj_tuned, add = TRUE, col = "darkgreen")
plot(gbm_roc_obj, add = TRUE, col = "darkorange")
graphics::legend("bottomright",
                 legend = c(
                   paste0("Elastic Net (AUC: ", round(auc(enet_roc_obj), 3), ")"),
                   paste0("Random Forest (AUC: ", round(auc(rf_roc_obj_tuned), 3), ")"),
                   paste0("GBM (AUC: ", round(auc(gbm_roc_obj), 3), ")")
                 ),
                 col = c("blue", "darkgreen", "darkorange"),
                 lwd = 2,
                 cex = 0.9
)
dev.off()

print("All plots saved as PNG files in your working directory.")

# --- 7. Create and Save Final Summary Table (FIXED) ---
final_summary <- data.frame(
  Model = c("Elastic Net (Tuned)", "Random Forest (Tuned)", "Gradient Boosting (Tuned)"),
  Test_AUC = c(
    auc(enet_roc_obj),
    auc(rf_roc_obj_tuned),
    auc(gbm_roc_obj)
  ),
  Test_Accuracy = c(
    enet_accuracy,
    rf_accuracy_tuned,
    gbm_accuracy
  ),
  Test_LogLoss = c(
    enet_logloss_test,
    rf_logloss_test_tuned,
    gbm_logloss_test
  )
)

# Print the table to the console
print(
  final_summary %>%
    arrange(desc(Test_AUC)) %>%
    mutate_if(is.numeric, round, 4) %>%
    kable(caption = "Final Model Performance on Hold-Out Test Set") %>%
    kable_styling(bootstrap_options = "striped", full_width = FALSE)
)

# Save the table to a file
library(gt)
final_summary %>%
  arrange(desc(Test_AUC)) %>%
  mutate_if(is.numeric, round, 4) %>%
  gt() %>%
  tab_header(title = "Final Model Performance on Hold-Out Test Set") %>%
  gtsave("final_summary_table.png")

print("Final summary table saved as final_summary_table.png")
