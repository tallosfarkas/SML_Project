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
set.seed(123)
cvfit <- cv.glmnet(
  x = x_train,
  y = y_train,
  family = "binomial",
  alpha = 0.5,         # Elastic Net: 0 = ridge, 1 = lasso
  nfolds = 10,
  type.measure = "class"
)

best_lambda <- cvfit$lambda.min
cat("Optimal lambda:", best_lambda, "\n")

# --- Refit on training data using best lambda ---------------------------------
model_enet <- glmnet(
  x = x_train,
  y = y_train,
  family = "binomial",
  alpha = 0.5,
  lambda = best_lambda
)

# --- Predict on test set ------------------------------------------------------
pred_prob <- predict(model_enet, newx = x_test, type = "response")
pred_class <- ifelse(pred_prob > 0.5, 1, 0)

# --- Evaluate performance -----------------------------------------------------
conf_mat <- table(Predicted = pred_class, Actual = y_test)
accuracy <- mean(pred_class == y_test)

cat("\nConfusion Matrix:\n")
print(conf_mat)
cat("\nClassification Accuracy:", round(accuracy, 3), "\n")

# Optional ROC/AUC
roc_obj <- roc(y_test, as.numeric(pred_prob))
cat("AUC:", round(auc(roc_obj), 3), "\n")

plot(roc_obj, col = "blue", main = "ROC Curve - Elastic Net Logistic Regression")


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
