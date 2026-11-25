# Validation of the Python DLNM with the R package dlnm (Gasparrini 2011)
# Same data, same model specification as notebook 04.



if (!requireNamespace("dlnm", quietly = TRUE)) install.packages("dlnm")
library(dlnm)
library(splines)

df <- read.csv("data/london_merged_daily.csv")
df$date <- as.Date(df$date)
df$time <- seq_len(nrow(df))
df$dow  <- factor(weekdays(df$date))
n_years <- length(unique(format(df$date, "%Y")))

# Lag knots equally spaced on the log scale - identical to the Python code
lag_knots <- exp(seq(0, log(22), length.out = 5))[2:4] - 1

fit_dlnm <- function(temp_col) {
  x <- df[[temp_col]]
  cb <- crossbasis(
    x, lag = 21,
    argvar = list(fun = "ns", knots = quantile(x, c(0.10, 0.75, 0.90)), Boundary.knots = range(x)),
    arglag = list(fun = "ns", knots = lag_knots, Boundary.knots = c(0, 21))
  )
  model <- glm(numdeaths ~ cb + ns(time, df = 7 * n_years) + dow,
               family = quasipoisson, data = df)
  
  # Minimum-mortality temperature between the 1st and 99th percentiles
  grid <- round(seq(quantile(x, 0.01), quantile(x, 0.99), by = 0.1), 1)
  pred_grid <- crosspred(cb, model, at = grid, cen = median(x))
  mmt <- grid[which.min(pred_grid$allRRfit)]
  
  # Relative risks at the 1st and 99th percentiles versus the MMT
  pct <- quantile(x, c(0.01, 0.99))
  pred <- crosspred(cb, model, at = pct, cen = mmt)
  
  list(cb = cb, model = model, mmt = mmt, pct = pct, pred = pred,
       dispersion = summary(model)$dispersion)
}

fmt <- function(p, i) sprintf("%.2f (%.2f-%.2f)", p$allRRfit[i], p$allRRlow[i], p$allRRhigh[i])

for (col in c("temperature", "era5_tmean")) {
  r <- fit_dlnm(col)
  cat(sprintf("%-12s overdispersion %.2f | MMT %.1f C | heat (99th, %.1f C) RR %s | cold (1st, %.1f C) RR %s\n",
              col, r$dispersion, r$mmt, r$pct[2], fmt(r$pred, 2), r$pct[1], fmt(r$pred, 1)))
}