# LOAD REQUIRED PACKAGES AND FUNCTIONS -----------------------------------------
if (!require("pacman")) install.packages("pacman")
pkgs = c("surveillance", "dplyr", "sp", "sf", "hhh4addon",
         "ggplot2", "hhh4addon") # package names
pacman::p_load(pkgs, character.only = T)

source("R/functions.R")

# LOAD DATA --------------------------------------------------------------------

# Cases at country level
counts <- readRDS("data/processed/daily_cases.rds")

# Shapefile for Africa
africa <- st_read("data/processed/geodata/africa.gpkg") 

# Weather
weather_clean <- readr::read_csv("data/original/AfricaCountries_2020-12-08_ALLEXTRACTEDDATA.csv")

# DATA PREPARATION -------------------------------------------------------------

# Policy variables
# Sindex is divided by 10 to show a 10% increase per unit increase
sindex <- readRDS("data/processed/stringency.rds") / 10
testing <- readRDS("data/processed/testing.rds")

# Weather data 
extrac_var_climate <- function(var, data) {
  data_clean <- data %>% 
    select(Date, name, var) %>% 
    filter(name %in% colnames(counts)) %>% 
    tidyr::spread(key = name, value = var) %>% 
    select(-Date) %>% 
    as.matrix()
  rownames(data_clean) <- as.character(unique(data$Date))
  return(data_clean)
}




rain_mean <- extrac_var_climate(var = "rain_mean", data = weather_clean)
temp_mean <- extrac_var_climate(var = "temp_mean", data = weather_clean)
sh_mean <- extrac_var_climate(var = "sh_mean", data = weather_clean) 

# Standardise climatic variables
rain_mean <- (rain_mean - mean(rain_mean)) / sd(rain_mean)
temp_mean <- (temp_mean - mean(temp_mean)) / sd(temp_mean)
sh_mean <- (sh_mean - mean(sh_mean)) / sd(sh_mean)


# See what the common dates are for the time varying datasets and censor
# accordingly 
final_dates <- Reduce(intersect, list(as.character(rownames(rain_mean)), rownames(counts), 
                                      rownames(sindex)))

counts <- counts[rownames(counts) %in% final_dates, ]
sindex <- sindex[rownames(sindex) %in% final_dates, ]
testing <- testing[rownames(testing) %in% final_dates, ]

# Check that the order of cases and countries in the shapefile are the same
all(colnames(counts) == africa$name)

map <- as(africa, "Spatial")
row.names(map) <- as.character(africa$name)

# Create adj mat and neighbours order
africa_adjmat <- poly2adjmat(map)
africa_nbOrder <- nbOrder(africa_adjmat, maxlag = Inf)

epi_sts <- sts(observed = counts,
               start = c(2020, 23),
               frequency = 365,
               population = africa$Pop2020 / sum(africa$Pop2020),
               neighbourhood = africa_nbOrder,
               map = map)

# Create covariates 

pop <- population(epi_sts)


# HDI by category
HDI_cat <- as.numeric(africa$HDI_Level)
# HDImedium <- ifelse(HDI_cat == 1, 1, 0) 
# HDIhigh <- ifelse(HDI_cat == 2, 1, 0) 
# HDI_medium <- matrix(HDImedium, ncol = ncol(epi_sts), nrow = nrow(epi_sts),
#                      byrow = T)
# 
# HDI_high <- matrix(HDIhigh, ncol = ncol(epi_sts), nrow = nrow(epi_sts),
#                    byrow = T)
HDI_cat <- matrix(HDI_cat, ncol = ncol(epi_sts), nrow = nrow(epi_sts),
                  byrow = T)

# Median age
mage <- matrix(africa$Median_age, ncol = ncol(epi_sts), nrow = nrow(epi_sts),
               byrow = T)

SSA <- matrix(1 - africa$North_Afri, 
              ncol = ncol(epi_sts), nrow = nrow(epi_sts),
              byrow = T)

LL <- matrix(africa$landlock, 
             ncol = ncol(epi_sts), nrow = nrow(epi_sts),
             byrow = T)
# MODEL ------------------------------------------------------------------------

k <- 7
sindex_lag <- xts::lag.xts(sindex, k = 7)
testing_lag <- xts::lag.xts(testing, k = 7)
temp_mean_lag <- xts::lag.xts(temp_mean, k = 7)
rain_mean_lag <- xts::lag.xts(rain_mean, k = 7)
sh_mean_lag <- xts::lag.xts(sh_mean, k = 7)

start_day <- "2020-03-28"
end_day <- "2020-08-06"

fit_start <- which(rownames(counts) == start_day) 
fit_end <- which(rownames(counts) == end_day) 


# MODEL ------------------------------------------------------------------------

# Best AR1 model with no RE 
f_end <- ~ 1 
f_ar <- ~ 1 + log(pop) + HDI_cat + LL + sindex_lag + testing_lag + 
  rain_mean_lag + temp_mean_lag + sh_mean_lag 
f_ne <- ~ 1 + log(pop) + HDI_cat + LL + sindex_lag + testing_lag 


model_basic <- list(
  end = list(f = f_end, offset = population(epi_sts)),
  ar = list(f = f_ar),
  ne = list(f = f_ne, weights = W_powerlaw(maxlag = 9)),
  optimizer = list(stop = list(iter.max = 50)),
  family = "NegBin1",
  subset = fit_start:fit_end)

fit_basic <- hhh4(epi_sts, control = model_basic)

# Lagged version of best AR1 model 
f_end <- ~ 1 
f_ar <- ~ 1 + log(pop) + HDI_cat + LL + sindex_lag + testing_lag + 
  rain_mean_lag + temp_mean_lag + sh_mean_lag 
f_ne <- ~ 1 + log(pop) + HDI_cat + LL + sindex_lag + testing_lag 
lags <- 8

AIC_poisson <- numeric(lags - 1)
for (i in 1:(lags - 1)) {
  model_lag <- list(
    end = list(f = f_end, offset = population(epi_sts)),
    ar = list(f = f_ar),
    ne = list(f = f_ne, weights = W_powerlaw(maxlag = 9)),
    optimizer = list(stop = list(iter.max = 50)),
    family = "NegBin1",
    subset = fit_start:fit_end,
    funct_lag = poisson_lag, 
    max_lag = i + 1)
  
  # Note that funct_lag = geometric_lag and max_lag = 5 are the defaults in hhh4lag 
  # and would not need to be specified explicitly.
  fit_lag_pois <- profile_par_lag(epi_sts, model_lag) # now use hhh4lag
  AIC_poisson[i] <- AIC(fit_lag_pois)
  print(i)
}

AIC_geom <- numeric(lags - 1)
for (i in 1:(lags - 1)) {
  model_lag <- list(
    end = list(f = f_end, offset = population(epi_sts)),
    ar = list(f = f_ar),
    ne = list(f = f_ne, weights = W_powerlaw(maxlag = 9)),
    optimizer = list(stop = list(iter.max = 50)),
    family = "NegBin1",
    subset = fit_start:fit_end,
    funct_lag = geometric_lag, 
    max_lag = i + 1)
  
  # Note that funct_lag = geometric_lag and max_lag = 5 are the defaults in hhh4lag 
  # and would not need to be specified explicitly.
  fit_lag_geom <- profile_par_lag(epi_sts, model_lag) # now use hhh4lag
  AIC_geom[i] <- AIC(fit_lag_geom)
  print(i)
}


# AIC table
tibble(p = 2:lags, Geometric = AIC_geom, Poisson = AIC_poisson, 
       aic_baseline = AIC(fit_basic)) %>% 
  tidyr::gather(key = "dist", value = "AIC", -p, -aic_baseline) %>% 
  mutate(diff = AIC - aic_baseline) %>% 
  ggplot(aes(x = p, y = diff, col = dist)) +
  geom_line() +
  geom_point() +
  labs(y = "Improvement in AIC", x = "D") +
  scale_color_brewer("", type = "q", palette = 6) +
  scale_x_continuous(breaks = 2:8) +
  theme_gray(base_size = 13) +
  theme(legend.position = "top") 

ggsave("figs/paper/AIC_ud.pdf", width = 7, height = 5)

# FIT BEST POISSON AND GEOMETRIC MODEL WITH OPTIMAL LAG AND PLOT WEIGHTS
# Geometric
model_lag <- list(
  end = list(f = f_end, offset = population(epi_sts)),
  ar = list(f = f_ar),
  ne = list(f = f_ne, weights = W_powerlaw(maxlag = 9)),
  optimizer = list(stop = list(iter.max = 50)),
  family = "NegBin1",
  subset = fit_start:fit_end,
  funct_lag = geometric_lag, 
  max_lag = 7)

fit_lag <- profile_par_lag(epi_sts, model_lag)
summary(fit_lag, idx2Exp = T)
confint(fit_lag, parm = "overdisp")

AIC(fit_lag)

wgeom <- fit_lag$distr_lag

model_lag <- list(
  end = list(f = f_end, offset = population(epi_sts)),
  ar = list(f = f_ar),
  ne = list(f = f_ne, weights = W_powerlaw(maxlag = 9)),
  optimizer = list(stop = list(iter.max = 50)),
  family = "NegBin1",
  subset = fit_start:fit_end,
  funct_lag = poisson_lag, 
  max_lag = 7)

fit_lag <- profile_par_lag(epi_sts, model_lag)
summary(fit_lag, idx2Exp = T)
confint(fit_lag, parm = "overdisp")

AIC(fit_lag)

wpois <- fit_lag$distr_lag

tibble(dist = rep(c("Geometric", "Poisson"), each = length(wgeom)),
       weight = c(wgeom, wpois), lag = rep(1:7, 2)) %>% 
  ggplot(aes(x = lag, y = weight, fill = dist)) +
  geom_bar(stat = "identity", position = position_dodge(), col = "black") +
  labs(x = "d", y = expression(w[d])) + 
  scale_x_continuous(breaks = 1:7) +
  scale_fill_brewer("", type = "q", palette = 6) +
  theme_gray(base_size = 13) +
  theme(legend.position = "top") 

ggsave("figs/paper/weights.pdf", width = 7, height = 5)


# LAGGED POISSON WITH AND WITHOUT RE WITH OPTIMAL NUMBER OF LAGS  
lag_optimal <- 7
f_end <- ~ 1 
f_ar <- ~ 1 + log(pop) + HDI_cat + LL + sindex_lag + testing_lag + 
  rain_mean_lag + temp_mean_lag + sh_mean_lag 
f_ne <- ~ 1 + log(pop) + HDI_cat + LL + sindex_lag + testing_lag 
  
  
model_lag <- list(
  end = list(f = f_end, offset = population(epi_sts)),
  ar = list(f = f_ar),
  ne = list(f = f_ne, weights = W_powerlaw(maxlag = 9)),
  optimizer = list(stop = list(iter.max = 50)),
  family = "NegBin1",
  data = list(pop = pop, HDI_cat = HDI_cat,
              SSA = SSA, LL = LL, sindex_lag = sindex_lag,
              rain_mean_lag = rain_mean_lag, temp_mean_lag = temp_mean_lag, 
              testing_lag = testing_lag,
              sh_mean_lag = sh_mean_lag),
  subset = fit_start:fit_end,
  funct_lag = poisson_lag, 
  max_lag = lag_optimal)
  
fit_lag <- profile_par_lag(epi_sts, model_lag)
summary(fit_lag, idx2Exp = T)
confint(fit_lag, parm = "overdisp")
wpois <-   fit_lag$distr_lag
  
saveRDS(fit_lag, paste0("output/fitted_model_LAG", lag_optimal, ".rds"))
  
fit <- fit_lag
nterms <- terms(fit)$nGroups + 2
coefs <- exp(coef(fit)[1:nterms])
CIs <- exp(confint(fit)[1:nterms, ])
id_log <-  c(grep("over", names(coefs)), grep("neweights.d", names(coefs)))
coefs[id_log] <- log(coefs[id_log])
CIs[id_log, ] <- log(CIs[id_log, ])
tab <- round(cbind(coefs, CIs), 3)
  
# Calculate scores of predictive performance
tp <- c(fit_end, nrow(epi_sts) - 1)
forecast <- oneStepAhead_hhh4lag(fit, tp = tp, type = "final")
fitScores <- colMeans(scores(forecast))
# Produce final summary table
tab <- rbind(cbind(Params = rownames(tab), tab), 
             c("", "", "", ""),
             c("AIC", round(AIC(fit), 2), "", ""),
             c("", "", "", ""),
             names(fitScores),
             round(as.numeric(fitScores), 3))
write.csv(tab, paste0("output/tab_params_LAG", lag_optimal, ".csv"), row.names = F)
  
# INTRODUCE REs
f_end <- ~ 1 
f_ar <- ~ -1 + log(pop) + HDI_cat + LL + sindex_lag + testing_lag + 
  rain_mean_lag + temp_mean_lag + sh_mean_lag + ri()
f_ne <- ~ -1 + log(pop) + HDI_cat + LL + sindex_lag + testing_lag + ri() 
  
model_lag <- list(
  end = list(f = f_end, offset = population(epi_sts)),
  ar = list(f = f_ar),
  ne = list(f = f_ne, weights = W_powerlaw(maxlag = 9)),
  optimizer = list(stop = list(iter.max = 10)),
  family = "NegBin1",
  data = list(pop = pop, HDI_cat = HDI_cat,
              SSA = SSA, LL = LL, sindex_lag = sindex_lag,
              rain_mean_lag = rain_mean_lag, temp_mean_lag = temp_mean_lag, 
              testing_lag = testing_lag,
              sh_mean_lag = sh_mean_lag),
  subset = fit_start:fit_end,
  funct_lag = poisson_lag, 
  par_lag = fit_lag$par_lag,
  max_lag = lag_optimal)
  
fit_lag <- hhh4_lag(epi_sts, model_lag)
saveRDS(fit_lag, paste0("output/fitted_model_LAG", lag_optimal, "_RE.rds"))
  
fit <- fit_lag
nterms <- terms(fit)$nGroups + 2
coefs <- exp(coef(fit)[1:nterms])
CIs <- exp(confint(fit)[1:nterms, ])
id_log <-  c(grep("over", names(coefs)), grep("neweights.d", names(coefs)))
coefs[id_log] <- log(coefs[id_log])
CIs[id_log, ] <- log(CIs[id_log, ])
tab <- round(cbind(coefs, CIs), 3)

# Calculate scores of predictive performance
tp <- c(fit_end, nrow(epi_sts) - 1)
forecast <- oneStepAhead_hhh4lag(fit, tp = tp, type = "final")
fitScores <- colMeans(scores(forecast))
# Produce final summary table
tab <- rbind(cbind(Params = rownames(tab), tab), 
             c("", "", "", ""),
             c("AIC", round(AIC(fit), 2), "", ""),
             c("", "", "", ""),
             names(fitScores),
             round(as.numeric(fitScores), 3))
write.csv(tab, paste0("output/tab_params_LAG", lag_optimal, "_RE.csv"), row.names = F)

# # BACKCASTING SIMULATION -------------------------------------------------------
# 
# nrows <- sum(seq(fit_start, fit_end, by = 7) <= fit_end - 14)
# 
# out <- matrix(NA, nrow = nrows, ncol = ncol(epi_sts))
# 
# days <- seq(fit_start, fit_end, by = 7)[1:nrows]
# 
# for(i in 1:length(days)) {
#   try({
#     f_end <- ~ 1 
#     f_ar <- ~ 1 + log(pop) + HDI_cat + LL + sindex_lag + testing_lag +
#       rain_mean_lag + temp_mean_lag + sh_mean_lag
#     f_ne <- ~ 1 + log(pop) + HDI_cat + LL + sindex_lag + testing_lag
# 
#     message("Fitting model without REs")
#     
#     model_lag <- list(
#       end = list(f = f_end, offset = population(epi_sts)),
#       ar = list(f = f_ar),
#       ne = list(f = f_ne, weights = W_powerlaw(maxlag = 9)),
#       optimizer = list(stop = list(iter.max = 50),
#                        variance = list(method = "Nelder-Mead")),
#       family = "NegBin1",
#       data = list(pop = pop, HDI_cat = HDI_cat,
#                   SSA = SSA, LL = LL, sindex_lag = sindex_lag,
#                   rain_mean_lag = rain_mean_lag, temp_mean_lag = temp_mean_lag, 
#                   testing_lag = testing_lag,
#                   sh_mean_lag = sh_mean_lag),
#       subset = days[i]:fit_end,
#       funct_lag = poisson_lag, 
#       max_lag = lag_optimal)
#     
#     fit_lag <- profile_par_lag(epi_sts, model_lag)
#     wpois <-   fit_lag$distr_lag
#     
#     message("Fitting model with REs")
#     
#     # INTRODUCE REs
#     f_end <- ~ 1 
#     f_ar <- ~ -1 + log(pop) + HDI_cat + LL + sindex_lag + testing_lag +
#       rain_mean_lag + temp_mean_lag + sh_mean_lag + ri()
#     f_ne <- ~ 1 + log(pop) + HDI_cat + LL + sindex_lag + testing_lag 
#     
#     model_lag <- list(
#       end = list(f = f_end, offset = population(epi_sts)),
#       ar = list(f = f_ar),
#       ne = list(f = f_ne, weights = W_powerlaw(maxlag = 9)),
#       optimizer = list(stop = list(iter.max = 50)),
#       family = "NegBin1",
#       data = list(pop = pop, HDI_cat = HDI_cat,
#                   SSA = SSA, LL = LL, sindex_lag = sindex_lag,
#                   rain_mean_lag = rain_mean_lag, temp_mean_lag = temp_mean_lag, 
#                   testing_lag = testing_lag,
#                   sh_mean_lag = sh_mean_lag),
#       subset = days[i]:fit_end,
#       funct_lag = poisson_lag, 
#       par_lag = fit_lag$par_lag,
#       max_lag = lag_optimal)
#     
#     fit_lag <- hhh4_lag(epi_sts, model_lag)
#     
#     fit <- fit_lag
#     tp <- c(fit_end, nrow(epi_sts) - 1)
#     forecast <- oneStepAhead_hhh4lag(fit, tp = tp, type = "final")
#     fitScores <- colMeans(scores(forecast, which = "logs", individual = T))
#     out[i, ] <- as.numeric(fitScores)
#   })
#   print(i)
# }
# 
# colnames(out) <- names(fitScores)
# out <- bind_cols(start_day = rownames(counts)[days], as_tibble(out))
# readr::write_csv(out, path = "output/backcasting.csv")
