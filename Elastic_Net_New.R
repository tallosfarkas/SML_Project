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



### Merge + features ###########################################################

final_data_df <- macro_df %>%
  left_join(sp500_df, by = "year_month") %>%
  left_join(vix_df,   by = "year_month") %>%
  left_join(dnsi_df,  by = "year_month") %>%
  select(-price, -price_lag1, -volume,-volume_change, return) %>%
  tidyr::drop_na()

# Ensure response factor
final_data_df <- final_data_df %>%
  mutate(UP_DOWN = factor(UP_DOWN, levels = c("Down", "Up")))

final_data_df$Y <- ifelse(final_data_df$UP_DOWN == "Up", 1, 0)

GGally::ggpairs(final_data_df %>% select(-year_month,-UP_DOWN,-return))

################################################################################
################################################################################

# --- Prepare data -------------------------------------------------------------

# We remove non-predictor variables and create interaction terms between the key
# macro and news variables

data_model <- final_data_df %>%
  select(-UP_DOWN, -year_month, -return) %>%
  mutate(
    DNSI_VIX_lag      = DNSI_change_lag * VIX_change_lag,
    DNSI_FedFunds_lag = DNSI_change_lag * FedFundsRate_lag,
    VIX_CPI_lag       = VIX_change_lag * CPI_lag,
    DNSI_NBER_lag     = DNSI_change_lag * NBER_lag
  ) %>%
  tidyr::drop_na()

# --- Chronological 30/70 split ------------------------------------------------

# We decided to use the first 30% of our time series data set for the calibrtion of the
# hyperparamter and coefficient estimation.
# Since it is a time series, we cannot randomize.


n <- nrow(data_model)
split_point <- floor(0.3 * n)

train_df <- data_model[1:split_point, ]
test_df  <- data_model[(split_point + 1):n, ]

x_train <- model.matrix(Y ~ . - 1, data = train_df)
y_train <- train_df$Y

x_test  <- model.matrix(Y ~ . - 1, data = test_df)
y_test  <- test_df$Y

################################################################################
# --- Time-Series Cross-Validation for Alpha and Lambda ------------------------
################################################################################

# For the CV, we test alphas (mixture of lasso and ridge) betwee [0,1], with using cv, for each
# alpha to define for each alphas the optimal lambda.
# createTimeSlices() creates expanding time windows



# We'll perform a small grid search for alpha and use time-series CV for lambda
alpha_grid <- seq(0, 1, by = 0.25)  # ridge -> lasso spectrum
n_folds <- 5

# Custom time-series folds
folds <- createTimeSlices(1:nrow(x_train),
                          initialWindow = floor(0.5 * nrow(x_train)),
                          horizon = floor((0.5 * nrow(x_train)) / n_folds),
                          fixedWindow = FALSE)

cv_results <- data.frame(alpha = numeric(), lambda = numeric(),
                         Accuracy = numeric(), AUC = numeric())


# The loop performes "nested" time-series cross validation

# For each alpha, run several expanding window folds and within each fold use cv.glmnet
# to find the best lambda

# Then fit the model and record the Accuracy and AUC on the validation slice

# Average the fold metrics and store them in cv_results.
# This gives us a grid of (alpha, lambda, performance) combinations.

for (a in alpha_grid) {
  cat("Testing alpha =", a, "\n")

  # Define time-series cross-validation manually
  accuracies <- c()
  aucs <- c()

  for (i in seq_along(folds$train)) {
    train_idx <- folds$train[[i]]
    val_idx   <- folds$test[[i]]

    x_tr <- x_train[train_idx, ]
    y_tr <- y_train[train_idx]
    x_val <- x_train[val_idx, ]
    y_val <- y_train[val_idx]

    cvfit <- cv.glmnet(
      x = x_tr,
      y = y_tr,
      family = "binomial",
      alpha = a,
      type.measure = "class"
    )

    best_lambda <- cvfit$lambda.min
    model <- glmnet(x_tr, y_tr, family = "binomial",
                    alpha = a, lambda = best_lambda)

    pred_prob <- predict(model, newx = x_val, type = "response")
    pred_class <- ifelse(pred_prob > 0.5, 1, 0)

    acc <- mean(pred_class == y_val)
    auc_val <- tryCatch({
      roc_obj <- roc(y_val, as.numeric(pred_prob), quiet = TRUE)
      as.numeric(auc(roc_obj))
    }, error = function(e) NA)

    accuracies <- c(accuracies, acc)
    aucs <- c(aucs, auc_val)
  }

  cv_results <- rbind(cv_results, data.frame(
    alpha = a,
    lambda = best_lambda,
    Accuracy = mean(accuracies, na.rm = TRUE),
    AUC = mean(aucs, na.rm = TRUE)
  ))
}

# --- Select best alpha & lambda ----------------------------------------------

# We then pick the alpha & lambda combination that gave the the highest mean AUC
# (Area under the curve) i.e. the area under the ROC curve.

best_row <- cv_results[which.max(cv_results$AUC), ]
best_alpha <- best_row$alpha
best_lambda <- best_row$lambda

cat("\nBest Alpha:", best_alpha, "\nBest Lambda:", best_lambda, "\n")

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

################################################################################
# --- Out-of-sample evaluation -------------------------------------------------
################################################################################

# For the OoS prediction, we use our final_model to predict the probabilities on the
# unseen test data. Finally, we just convert probabilities to a binary class predictions


pred_prob <- predict(final_model, newx = x_test, type = "response")
pred_class <- ifelse(pred_prob > 0.5, 1, 0)

conf_mat <- table(Predicted = pred_class, Actual = y_test)
accuracy <- mean(pred_class == y_test)

# We can evaluate the classifications with a confusion matrix, accuracy and
# ROC curve

cat("\nConfusion Matrix:\n")
print(conf_mat)
cat("\nOut-of-sample Accuracy:", round(accuracy, 3), "\n")

roc_obj <- roc(y_test, as.numeric(pred_prob))
cat("Out-of-sample AUC:", round(auc(roc_obj), 3), "\n")

plot(roc_obj, col = "blue", main = "ROC Curve - Out-of-sample (Elastic Net)")
