# DLNM effect plots: station vs ERA5-Land temperature, London 2002-2006



for (p in c("dlnm", "ggplot2")) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
library(dlnm)
library(splines)
library(ggplot2)

df <- read.csv("data/london_merged_daily.csv")
df$date <- as.Date(df$date)
df$time <- seq_len(nrow(df))
df$dow  <- factor(weekdays(df$date))
n_years <- length(unique(format(df$date, "%Y")))
lag_knots <- exp(seq(0, log(22), length.out = 5))[2:4] - 1

sources <- c("Station" = "temperature", "ERA5-Land" = "era5_tmean")
colours <- c("Station" = "black", "ERA5-Land" = "#C0392B")

overall <- list(); lagres <- list(); fits <- list()

for (src in names(sources)) {
  x  <- df[[sources[[src]]]]
  cb <- crossbasis(x, lag = 21,
                   argvar = list(fun = "ns", knots = quantile(x, c(0.10, 0.75, 0.90)), Boundary.knots = range(x)),
                   arglag = list(fun = "ns", knots = lag_knots, Boundary.knots = c(0, 21)))
  model <- glm(numdeaths ~ cb + ns(time, df = 7 * n_years) + dow, family = quasipoisson, data = df)
  
  # minimum-mortality temperature (1st-99th percentile range)
  grid <- round(seq(quantile(x, 0.01), quantile(x, 0.99), by = 0.1), 1)
  mmt  <- grid[which.min(crosspred(cb, model, at = grid, cen = median(x))$allRRfit)]
  
  # cumulative exposure-response over lags 0-21
  at <- round(seq(quantile(x, 0.005), quantile(x, 0.995), by = 0.1), 1)
  cp <- crosspred(cb, model, at = at, cen = mmt)
  overall[[src]] <- data.frame(source = src, temp = cp$predvar, rr = cp$allRRfit,
                               low = cp$allRRlow, high = cp$allRRhigh, mmt = mmt)
  
  # lag-response at the 99th (heat) and 1st (cold) percentiles
  pct <- quantile(x, c(0.01, 0.99))
  cpl <- crosspred(cb, model, at = pct, cen = mmt, bylag = 0.2)
  for (i in 1:2) {
    lagres[[paste(src, i)]] <- data.frame(
      source = src, effect = c("Cold (1st percentile)", "Heat (99th percentile)")[i],
      lag = seq(0, 21, by = 0.2), rr = cpl$matRRfit[i, ], low = cpl$matRRlow[i, ], high = cpl$matRRhigh[i, ])
  }
  fits[[src]] <- list(cb = cb, model = model, mmt = mmt)
}
overall <- do.call(rbind, overall)
lagres  <- do.call(rbind, lagres)
mmts    <- unique(overall[, c("source", "mmt")])

theme_set(theme_minimal(base_size = 12) + theme(panel.grid.minor = element_blank(), legend.position = "top"))
dir.create("figures", showWarnings = FALSE)

# 1. Cumulative exposure-response, both temperature sources
p1 <- ggplot(overall, aes(temp, rr, colour = source, fill = source)) +
  geom_hline(yintercept = 1, colour = "grey50") +
  geom_ribbon(aes(ymin = low, ymax = high), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 1) +
  geom_vline(data = mmts, aes(xintercept = mmt, colour = source), linetype = "dashed", linewidth = 0.5, show.legend = FALSE) +
  scale_colour_manual(values = colours, name = NULL) +
  scale_fill_manual(values = colours, name = NULL) +
  labs(x = "Daily mean temperature (\u00B0C)", y = "Cumulative relative risk (lag 0\u201321)",
       title = "Temperature\u2013mortality association in London, 2002\u20132006",
       subtitle = "Dashed lines: minimum-mortality temperature")
ggsave("figures/R_exposure_response.png", p1, width = 7.5, height = 5, dpi = 150, bg = "white")

# 2. Lag-response for heat and cold, both temperature sources
p2 <- ggplot(lagres, aes(lag, rr, colour = source, fill = source)) +
  geom_hline(yintercept = 1, colour = "grey50") +
  geom_ribbon(aes(ymin = low, ymax = high), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 1) +
  facet_wrap(~ effect, scales = "free_y") +
  scale_colour_manual(values = colours, name = NULL) +
  scale_fill_manual(values = colours, name = NULL) +
  labs(x = "Lag (days)", y = "Relative risk vs minimum-mortality temperature",
       title = "Lag structure of heat and cold effects")
ggsave("figures/R_lag_response.png", p2, width = 10, height = 4.5, dpi = 150, bg = "white")

# 3. 3D surface of the full exposure-lag-response (station temperature)
cb    <- fits[["Station"]]$cb      # crosspred matches coefficients by the name "cb"
model <- fits[["Station"]]$model
cp3 <- crosspred(cb, model, cen = fits[["Station"]]$mmt, by = 0.5)
png("figures/R_3d_station.png", width = 1000, height = 800, res = 150)
plot(cp3, xlab = "Temperature (\u00B0C)", zlab = "RR", theta = 220, phi = 30, lphi = 30,
     main = "Exposure\u2013lag\u2013response, station temperature")
dev.off()

p1
p2