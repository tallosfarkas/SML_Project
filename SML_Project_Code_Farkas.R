#### SML Project ####

### Packages ###################################################################
required_packages <- c(
  "xts","zoo","frenchdata","mistr","ggplot2","quantmod","dplyr","purrr",
  "lubridate","stringr","tidyr","stats","gt","PerformanceAnalytics","kableExtra",
  "knitr","broom","tibble","tidyquant","tidymodels","glmnet","ranger","doParallel","GGally","readxl",
  "pROC"
)
installed <- rownames(installed.packages())
missing <- setdiff(required_packages, installed)
if (length(missing) > 0) install.packages(missing, dependencies = TRUE)
invisible(lapply(required_packages, library, character.only = TRUE))
tidymodels_prefer()

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



### Merge + features ###########################################################

final_data_df <- macro_df %>%
  left_join(sp500_df, by = "year_month") %>%
  left_join(vix_df,   by = "year_month") %>%
  left_join(dnsi_df,  by = "year_month") %>%
  select(-price, -price_lag1, -volume,-volume_change,-return) %>%
  tidyr::drop_na()

# Ensure response factor
final_data_df <- final_data_df %>%
  mutate(UP_DOWN = factor(UP_DOWN, levels = c("Down", "Up")))

GGally::ggpairs(final_data_df %>% select(-year_month,-UP_DOWN))

# --- Prepare data -------------------------------------------------------------
final_data_df$Y <- ifelse(final_data_df$UP_DOWN == "Up", 1, 0)

data_model <- final_data_df %>%
  select(-UP_DOWN, -year_month) %>%
  mutate(
    # Example interaction effects
    DNSI_VIX_lag       = DNSI_change_lag * VIX_change_lag,
    DNSI_FedFunds_lag  = DNSI_change_lag * FedFundsRate_lag,
    VIX_CPI_lag        = VIX_change_lag * CPI_lag,
    DNSI_NBER_lag      = DNSI_change_lag * NBER_lag
  ) %>%
  tidyr::drop_na()

# --- Chronological 70/30 split -----------------------------------------------
n <- nrow(data_model)
split_point <- floor(0.7 * n)

train_df <- data_model[1:split_point, ]
test_df  <- data_model[(split_point + 1):n, ]

# --- Matrices for glmnet ------------------------------------------------------
x_train <- model.matrix(Y ~ . - 1, data = train_df)
y_train <- train_df$Y

x_test  <- model.matrix(Y ~ . - 1, data = test_df)
y_test  <- test_df$Y

# --- Fit Elastic Net Logistic Regression -------------------------------------
# --- Parameter grid to tune for Elastic Net ---
# We tune alpha (the L1/L2 mix) and lambda (the penalty)
param_grid_enet <- expand.grid(
  alpha = c(0.0, 0.25, 0.5, 0.75, 1.0), # 0=Ridge, 1=Lasso
  lambda = 10^seq(-4, -1, length.out = 10)
)

# --- Rolling-window tuning within training data (for glmnet) ---
# Use the same window parameters as your RF loop
n_train <- nrow(train_df)
window_size <- floor(0.7 * n_train)
horizon <- 12

results_enet <- data.frame()
eps <- 1e-9 # for log-loss stability

print("Starting Elastic Net time-series tuning...")
set.seed(123)
for (i in seq(window_size, n_train - horizon, by = horizon)) {
  
  # Define the row indices for this fold
  train_indices <- 1:i
  test_indices <- (i + 1):(i + horizon)
  
  # --- Get the window data from the pre-split matrices ---
  # This is much faster as model.matrix isn't re-run
  x_train_window <- x_train[train_indices, ]
  y_train_window <- y_train[train_indices]
  
  x_test_window <- x_train[test_indices, ] # Use validation part of x_train
  y_test_window <- y_train[test_indices] # Use validation part of y_train
  
  # Skip invalid test windows
  if (nrow(x_test_window) == 0 || length(unique(y_test_window)) < 2) next
  
  for (j in 1:nrow(param_grid_enet)) {
    p <- param_grid_enet[j, ]
    
    # Fit the glmnet model
    enet_model <- glmnet(
      x_train_window,
      y_train_window,
      family = "binomial",
      alpha = p$alpha,
      lambda = p$lambda # Note: glmnet can tune lambda itself, but this is more explicit
    )
    
    # Predict on the validation window
    preds <- predict(enet_model, newx = x_test_window, type = "response")
    
    # --- Log-loss (cross-entropy) ---
    logloss <- -mean(y_test_window * log(preds + eps) + (1 - y_test_window) * log(1 - preds + eps))
    auc_val <- tryCatch(as.numeric(auc(y_test_window, as.numeric(preds))), error = function(e) NA)
    
    results_enet <- rbind(
      results_enet,
      data.frame(
        alpha = p$alpha,
        lambda = p$lambda,
        LogLoss = logloss,
        AUC = auc_val
      )
    )
  }
}
print("Elastic Net tuning complete.")

# --- Aggregate across rolling folds ---
tuning_summary_enet <- results_enet %>%
  group_by(alpha, lambda) %>%
  summarise(
    mean_LogLoss = mean(LogLoss, na.rm = TRUE),
    mean_AUC = mean(AUC, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(mean_LogLoss) # lower log-loss is better

cat("Top Elastic Net tuning results (lowest log-loss):\n")
print(head(tuning_summary_enet, 5))

# --- Select best parameters ---
best_params_enet <- tuning_summary_enet[1, ]
cat("\nSelected Elastic Net parameters:\n")
print(best_params_enet)

# --- Fit final model on full training data ---
set.seed(123)
enet_final <- glmnet(
  x_train,
  y_train,
  family = "binomial",
  alpha = best_params_enet$alpha,
  lambda = best_params_enet$lambda
)

# --- Evaluate on hold-out test set ---
enet_pred_prob <- predict(enet_final, newx = x_test, type = "response")
enet_pred_class <- ifelse(enet_pred_prob > 0.5, 1, 0)
y_true_test <- y_test # Use the numeric y_test

# --- Metrics ---
enet_conf_mat <- table(Predicted = enet_pred_class, Actual = y_true_test)
enet_accuracy <- mean(enet_pred_class == y_true_test)
enet_logloss_test <- -mean(y_true_test * log(enet_pred_prob + eps) + (1 - y_true_test) * log(1 - enet_pred_prob + eps))

cat("\n--- Elastic Net Final Test Performance ---\n")
cat("\nConfusion Matrix:\n"); print(enet_conf_mat)
cat("\nTest Accuracy:", round(enet_accuracy, 3))
cat("\nTest Log-Loss:", round(enet_logloss_test, 4))

# --- ROC/AUC plot ---
enet_roc_obj <- roc(y_true_test, as.numeric(enet_pred_prob))
cat("\nTest AUC:", round(auc(enet_roc_obj), 3), "\n")
plot(enet_roc_obj, col = "blue", main = "ROC Curve - Tuned Elastic Net")

# --- Coefficients ---
cat("\nFinal Model Coefficients:\n")
print(coef(enet_final))



### --- Random Forest Classification (time-series safe) --- ###################


# --- Prepare data -------------------------------------------------------------
# (You already created train_df and test_df above)

# Convert UP/DOWN target to factor for ranger
train_df$Y <- factor(train_df$Y, levels = c(0, 1))
test_df$Y  <- factor(test_df$Y, levels = c(0, 1))

# --- Fit Random Forest --------------------------------------------------------
set.seed(123)

rf_model <- ranger(
  formula = Y ~ .,
  data = train_df,
  num.trees = 500,
  mtry = floor(sqrt(ncol(train_df) - 1)),  # typical heuristic
  importance = "impurity",
  probability = TRUE
)

# --- Predict on test data -----------------------------------------------------
rf_pred <- predict(rf_model, data = test_df)

# Extract predicted probabilities and classes
pred_prob <- rf_pred$predictions[, "1"]
pred_class <- ifelse(pred_prob > 0.5, 1, 0)

# --- Evaluate performance -----------------------------------------------------
conf_mat <- table(Predicted = pred_class, Actual = test_df$Y)
accuracy <- mean(pred_class == as.numeric(as.character(test_df$Y)))

cat("\nConfusion Matrix:\n")
print(conf_mat)
cat("\nClassification Accuracy:", round(accuracy, 3), "\n")

# --- AUC ----------------------------------------------------------------------
roc_obj <- roc(as.numeric(as.character(test_df$Y)), pred_prob)
cat("AUC:", round(auc(roc_obj), 3), "\n")

plot(roc_obj, col = "forestgreen", main = "ROC Curve - Random Forest")

# --- Feature Importance -------------------------------------------------------
importance_df <- data.frame(
  Variable = names(rf_model$variable.importance),
  Importance = rf_model$variable.importance
) %>%
  arrange(desc(Importance))

cat("\nTop 10 Most Important Features:\n")
print(head(importance_df, 10))


### --- Random Forest (Time-Series Safe, Log-Loss Tuned) --- ###################

# --- Parameter grid to tune ---------------------------------------------------
param_grid <- expand.grid(
  mtry = c(2, 4, 6, 8, 10),
  min.node.size = c(1, 5, 10),
  sample.fraction = c(0.6, 0.8, 1.0)
)

# --- Rolling-window tuning within training data ------------------------------
n_train <- nrow(train_df)
window_size <- floor(0.7 * n_train)   # first rolling window (70%)
horizon <- 12                          # predict next month

results <- data.frame()


print(window_size)
print(n_train)
print(seq(window_size, n_train - horizon, by = horizon))
set.seed(123)
for (i in seq(window_size, n_train - horizon, by = horizon)) {
  
  train_window <- train_df[1:i, ]
  test_window  <- train_df[(i + 1):(i + horizon), ]
  
  # Skip invalid test windows
  if (nrow(test_window) == 0 || length(unique(test_window$Y)) < 2) next
  
  for (j in 1:nrow(param_grid)) {
    p <- param_grid[j, ]
    
    rf_model <- ranger(
      Y ~ .,
      data = train_window,
      num.trees = 500,
      mtry = p$mtry,
      min.node.size = p$min.node.size,
      sample.fraction = p$sample.fraction,
      probability = TRUE,
      seed = 123
    )
    
    preds <- predict(rf_model, data = test_window)$predictions[, "1"]
    y_true <- as.numeric(as.character(test_window$Y))
    
    # --- Log-loss (cross-entropy) --------------------------------------------
    eps <- 1e-9
    logloss <- -mean(y_true * log(preds + eps) + (1 - y_true) * log(1 - preds + eps))
    
    # --- AUC (for reference) -------------------------------------------------
    auc_val <- tryCatch(as.numeric(auc(y_true, preds)), error = function(e) NA)
    
    results <- rbind(
      results,
      data.frame(
        mtry = p$mtry,
        min.node.size = p$min.node.size,
        sample.fraction = p$sample.fraction,
        LogLoss = logloss,
        AUC = auc_val
      )
    )
  }
}

# --- Aggregate across rolling folds ------------------------------------------
tuning_summary <- results %>%
  group_by(mtry, min.node.size, sample.fraction) %>%
  summarise(
    mean_LogLoss = mean(LogLoss, na.rm = TRUE),
    mean_AUC = mean(AUC, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(mean_LogLoss)   # lower log-loss is better

cat("Top tuning results (lowest log-loss):\n")
print(head(tuning_summary, 5))

# --- Select best parameters ---------------------------------------------------
best_params <- tuning_summary[1, ]
cat("\nSelected parameters:\n")
print(best_params)

# --- Fit final model on full training data -----------------------------------
set.seed(123)
rf_final <- ranger(
  Y ~ .,
  data = train_df,
  num.trees = 500,
  mtry = best_params$mtry,
  min.node.size = best_params$min.node.size,
  sample.fraction = best_params$sample.fraction,
  probability = TRUE,
  importance = "impurity"
)

# --- Evaluate on hold-out test set -------------------------------------------
rf_pred <- predict(rf_final, data = test_df)
pred_prob <- rf_pred$predictions[, "1"]
pred_class <- ifelse(pred_prob > 0.5, 1, 0)
y_true <- as.numeric(as.character(test_df$Y))

# --- Metrics -----------------------------------------------------------------
conf_mat <- table(Predicted = pred_class, Actual = y_true)
accuracy <- mean(pred_class == y_true)

# Log-loss
eps <- 1e-9
logloss_test <- -mean(y_true * log(pred_prob + eps) + (1 - y_true) * log(1 - pred_prob + eps))

cat("\nConfusion Matrix:\n"); print(conf_mat)
cat("\nTest Accuracy:", round(accuracy, 3))
cat("\nTest Log-Loss:", round(logloss_test, 4))

# --- ROC/AUC plot -------------------------------------------------------------
if (length(unique(y_true)) == 2) {
  roc_obj <- roc(y_true, pred_prob)
  cat("\nTest AUC:", round(auc(roc_obj), 3), "\n")
  plot(roc_obj, col = "darkgreen", main = "ROC Curve - Tuned Random Forest (Lagged Features)")
} else {
  cat("\nROC skipped: test set has only one class.\n")
}

# --- Feature Importance -------------------------------------------------------
importance_df <- data.frame(
  Variable = names(rf_final$variable.importance),
  Importance = rf_final$variable.importance
) %>%
  arrange(desc(Importance))

cat("\nTop 10 Most Important Features:\n")
print(head(importance_df, 10))





# We need the 'gbm' package
if (!"gbm" %in% installed) install.packages("gbm")
library(gbm)

# --- Parameter grid to tune for GBM ---
param_grid_gbm <- expand.grid(
  n.trees = c(100, 200, 300),
  interaction.depth = c(1, 2, 3), # As in slides07 (1=additive, >1=interactions)
  shrinkage = c(0.01, 0.1)    # Learning rate
)

# --- Rolling-window tuning (for gbm) ---
# Use the same window parameters
n_train <- nrow(train_df)
window_size <- floor(0.7 * n_train)
horizon <- 12

results_gbm <- data.frame()
eps <- 1e-9 # for log-loss stability

print("Starting GBM time-series tuning...")
set.seed(123)
for (i in seq(window_size, n_train - horizon, by = horizon)) {
  
  # Get window data (dataframes)
  train_window <- train_df[1:i, ]
  test_window <- train_df[(i + 1):(i + horizon), ]
  
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
    auc_val <- tryCatch(as.numeric(auc(y_true_window, preds)), error = function(e) NA)
    
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
# Prepare full train_df for gbm
train_df_gbm <- train_df
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
test_df_gbm <- test_df
test_df_gbm$Y_num <- as.numeric(as.character(test_df_gbm$Y))
y_true_test <- test_df_gbm$Y_num # Use the numeric 0/1 Y

gbm_pred_prob <- predict(gbm_final,
                         newdata = test_df_gbm,
                         n.trees = best_params_gbm$n.trees,
                         type = "response")
gbm_pred_class <- ifelse(gbm_pred_prob > 0.5, 1, 0)

# --- Metrics ---
gbm_conf_mat <- table(Predicted = gbm_pred_class, Actual = y_true_test)
gbm_accuracy <- mean(gbm_pred_class == y_true_test)
gbm_logloss_test <- -mean(y_true_test * log(gbm_pred_prob + eps) + (1 - y_true_test) * log(1 - gbm_pred_prob + eps))

cat("\n--- GBM Final Test Performance ---\n")
cat("\nConfusion Matrix:\n"); print(gbm_conf_mat)
cat("\nTest Accuracy:", round(gbm_accuracy, 3))
cat("\nTest Log-Loss:", round(gbm_logloss_test, 4))

# --- ROC/AUC plot ---
gbm_roc_obj <- roc(y_true_test, gbm_pred_prob)
cat("\nTest AUC:", round(auc(gbm_roc_obj), 3), "\n")
plot(gbm_roc_obj, col = "darkorange", main = "ROC Curve - Tuned GBM")

# --- Feature Importance ---
cat("\nGBM Feature Importance:\n")
print(summary(gbm_final, plotit = FALSE))






# --- 1. Plot all ROC curves together ---
plot(enet_roc_obj, col = "blue", main = "Final Model ROC Comparison (Test Set)")
plot(roc_obj, add = TRUE, col = "darkgreen")
plot(gbm_roc_obj, add = TRUE, col = "darkorange")

graphics::legend("bottomright",
       legend = c(
         paste0("Elastic Net (AUC: ", round(auc(enet_roc_obj), 3), ")"),
         paste0("Random Forest (AUC: ", round(auc(roc_obj), 3), ")"),
         paste0("GBM (AUC: ", round(auc(gbm_roc_obj), 3), ")")
       ),
       col = c("blue", "darkgreen", "darkorange"),
       lwd = 2
)

# --- 2. Create a final comparison table ---
# Get probabilities from your RF block (you named it 'pred_prob')
rf_pred_prob <- predict(rf_final, data = test_df)$predictions[, "1"]
rf_accuracy <- mean(ifelse(rf_pred_prob > 0.5, 1, 0) == y_true_test)
rf_logloss_test <- -mean(y_true_test * log(rf_pred_prob + eps) + (1 - y_true_test) * log(1 - rf_pred_prob + eps))

final_summary <- data.frame(
  Model = c("Elastic Net", "Random Forest", "Gradient Boosting"),
  Test_AUC = c(
    auc(enet_roc_obj),
    auc(roc_obj),
    auc(gbm_roc_obj)
  ),
  Test_Accuracy = c(
    enet_accuracy,
    rf_accuracy,
    gbm_accuracy
  ),
  Test_LogLoss = c(
    enet_logloss_test,
    rf_logloss_test,
    ----------------------------------------------------------------------
      gbm_logloss_test
  )
)

# Print a nice table
final_summary %>%
  arrange(desc(Test_AUC)) %>%
  mutate_if(is.numeric, round, 4) %>%
  kable(caption = "Final Model Performance on Hold-Out Test Set") %>%
  kable_styling(bootstrap_options = "striped", full_width = FALSE)






################ PLOTS #############################



# --- 1. Save EDA Plot ---
# This one is large and complex, so we set a larger size
print("Saving EDA plot...")
png("ggpairs_plot.png", width = 1200, height = 1000, res = 100)
print(
  GGally::ggpairs(final_data_df %>% select(-year_month,-UP_DOWN)))
dev.off()

# --- 2. Save Tuned Elastic Net ROC Plot ---
print("Saving Elastic Net ROC plot...")
png("enet_roc_plot.png", width = 800, height = 600, res = 100)
plot(enet_roc_obj, col = "blue", main = "ROC Curve - Tuned Elastic Net")
dev.off()

# --- 3. Save Tuned Random Forest ROC & VarImp Plots ---
# First, create the correct objects for the TUNED Random Forest
rf_pred_prob_tuned <- predict(rf_final, data = test_df)$predictions[, "1"]
rf_roc_obj_tuned <- roc(y_true_test, rf_pred_prob_tuned, quiet = TRUE)
rf_importance_df <- data.frame(
  Variable = names(rf_final$variable.importance),
  Importance = rf_final$variable.importance
) %>% arrange(desc(Importance))

print("Saving Random Forest plots...")
png("rf_tuned_roc_plot.png", width = 800, height = 600, res = 100)
plot(rf_roc_obj_tuned, col = "darkgreen", main = "ROC Curve - Tuned Random Forest")
dev.off()

png("rf_var_imp_plot.png", width = 800, height = 700, res = 100)
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


# --- 4. Save Tuned GBM ROC Plot ---
print("Saving GBM ROC plot...")
png("gbm_roc_plot.png", width = 800, height = 600, res = 100)
plot(gbm_roc_obj, col = "darkorange", main = "ROC Curve - Tuned GBM")
dev.off()

# --- 5. Save Combined ROC Plot ---
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

# --- 6. Create and Save Final Summary Table ---
# This fixes the errors in your script
rf_accuracy_tuned <- mean(ifelse(rf_pred_prob_tuned > 0.5, 1, 0) == y_true_test)
rf_logloss_test_tuned <- -mean(y_true_test * log(rf_pred_prob_tuned + eps) + (1 - y_true_test) * log(1 - rf_pred_prob_tuned + eps))

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

# Print the table so you can screenshot it
print(
  final_summary %>%
    arrange(desc(Test_AUC)) %>%
    mutate_if(is.numeric, round, 4) %>%
    kable(caption = "Final Model Performance on Hold-Out Test Set") %>%
    kable_styling(bootstrap_options = "striped", full_width = FALSE)
)

# You can also save this table to a file
library(gt)
final_summary %>%
  arrange(desc(Test_AUC)) %>%
  mutate_if(is.numeric, round, 4) %>%
  gt() %>%
  gtsave("final_summary_table.png")

print("Final summary table saved as final_summary_table.png")

