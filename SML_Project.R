
#### SML Project ####

### Needed Packages ############################################################

#We download any not already installed R package to smoothly run the code.

required_packages <- c(
  "xts", "zoo", "frenchdata", "mistr", "ggplot2",
  "quantmod", "dplyr", "purrr", "lubridate",
  "stringr", "tidyr", "stats", "gt", "PerformanceAnalytics","kableExtra","knitr","broom", "tibble","tidyquant")
installed <- rownames(installed.packages())
missing <- setdiff(required_packages, installed)

if (length(missing) > 0) {
  install.packages(missing, dependencies = TRUE)
}
# Load all required packages
invisible(lapply(required_packages, library, character.only = TRUE))



### Data Download ##############################################################


# 1. Download daily S&P 500 data
start <- as.Date("1950-01-01")
end <- Sys.Date()

sp500_xts <- getSymbols("^GSPC", from = start, to = end,
                        src = "yahoo", auto.assign = FALSE)

sp500_first <- do.call(rbind, lapply(split(sp500_xts, f = "months"), head, 1))


# 6. Combine into final data frame
sp500_df <- data.frame(
  date = index(sp500_first), coredata(sp500_first)
)

sp500_df <- sp500_df %>%
  mutate(
    year_month = floor_date(date, "month"))

# 7. Add lagged variables

sp500_df <- sp500_df %>%
  dplyr::select(year_month,GSPC.Adjusted, GSPC.Volume) %>%
  rename(price = GSPC.Adjusted,
         volume = GSPC.Volume) %>%
  mutate(
    price_lag1 = lag(price,1),
    price_lag2 = lag(price,2),
    price_lag3 = lag(price,3),
    price_lag4 = lag(price,4),
    price_lag5 = lag(price,5),
  ) %>%
  drop_na


### Macro Data ################################################################


# Download macro data from FRED --> GDP, Consumer Price Index, Fed Funds Rate
getSymbols(c("GDP", "CPIAUCSL", "FEDFUNDS"), src = "FRED", from = start)

# Load NBER recession dates
nber_dates <- tq_get("USREC", get = "economic.data", from = start)
nber_xts <- xts(nber_dates[["price"]], order.by = nber_dates$date)

# Merge data (monthly basis)
macro_data <- merge(
  GDP,
  CPIAUCSL,
  FEDFUNDS,
  nber_xts
)
colnames(macro_data) <- c("GDP", "CPI", "FedFundsRate","NBER")

macro_data_df <- data.frame(date = index(macro_data), coredata(macro_data))

# Apply lags for information availability
macro_data_df <- macro_data_df %>%
  mutate(
    GDP = lag(GDP, 2),
    CPI = lag(CPI, 1),
  )

# Forward-fill for monthly completeness
macro_data_df[ , 2:5] <- na.locf(macro_data_df[ , 2:5], na.rm = FALSE)

# Caluclate YoY Growth
macro_data_df <- macro_data_df %>%
  mutate(
    GDP_YoY = ROC(GDP, n = 12, type = "discrete"),
    CPI_YoY = ROC(CPI, n = 12, type = "discrete")
  )



## Inflation Lagged and Surprise

macro_data_df <- macro_data_df %>%
  mutate(
    CPI_YoY_lag = lag(CPI_YoY, 12),
    inflation_surprise = CPI_YoY - CPI_YoY_lag
    )

## Fed Rate

macro_data_df <- macro_data_df %>%
  mutate(
    FedChange = FedFundsRate - lag(FedFundsRate),
    year_month = floor_date(date, "month")) %>%
  dplyr::select(-date)

macro_data_df <- na.omit(macro_data_df)

### Merge Everything ##########################################################

final_data_df <-macro_data_df  %>%
  left_join(sp500_df, by = "year_month")
