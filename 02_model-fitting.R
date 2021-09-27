# LOAD REQUIRED PACKAGES AND FUNCTIONS -----------------------------------------
if (!require("pacman")) install.packages("pacman")
pkgs = c("surveillance", "dplyr", "sp", "sf", "hhh4addon",
         "ggplot2", "hhh4addon") # package names
pacman::p_load(pkgs, character.only = T)

source("R/functions.R")


# CLEAR OUTPUT FOLDER
file.remove(list.files("output/models",full.names=TRUE))

# LOAD DATA --------------------------------------------------------------------

# Cases at country level
counts <- readRDS("data/processed/daily_cases.rds")

# Shapefile for Africa
africa <- st_read("data/processed/geodata/africa.gpkg") 

# Weather
# weather_clean <- readr::read_csv("data/original/AfricaCountries_2020-12-08_ALLEXTRACTEDDATA.csv")

# DATA PREPARATION -------------------------------------------------------------

# Policy variables
# Sindex is divided by 10 to show a 10% increase per unit increase
sindex <- readRDS("data/processed/stringency.rds") / 10
testing <- readRDS("data/processed/testing.rds")

# Weather data 
# rain_mean <- extrac_var_climate(var = "rain_mean", data = weather_clean)
# temp_mean <- extrac_var_climate(var = "temp_mean", data = weather_clean)
# sh_mean <- extrac_var_climate(var = "sh_mean", data = weather_clean) 


# See what the common dates are for the time varying datasets and censor
# accordingly 
# final_dates <- Reduce(intersect, list(rownames(rain_mean), rownames(counts), 
#                                       rownames(sindex)))
final_dates <- Reduce(intersect, list(rownames(counts),
                                      rownames(sindex)))

counts <- counts[rownames(counts) %in% final_dates, ]
sindex <- sindex[rownames(sindex) %in% final_dates, ]
testing <- testing[rownames(testing) %in% final_dates, ]


# rain_mean <- rain_mean[rownames(rain_mean) %in% final_dates,]
# temp_mean <- temp_mean[rownames(temp_mean) %in% final_dates,]
# sh_mean <- sh_mean[rownames(sh_mean) %in% final_dates,]


end_day = final_dates[length(final_dates)]
shinydata = as.data.frame(matrix(nrow=as.numeric(as.Date(end_day)-as.Date(final_dates[1])+1)*dim(counts)[2],ncol=0))
shinydata$time = rep(seq(from=as.Date(final_dates[1]),to=as.Date(end_day),by=1),dim(counts)[2])
shinydata$country = as.vector(t(matrix(colnames(counts),dim(counts)[2],as.numeric(as.Date(end_day)-as.Date(final_dates[1])+1))))
shinydata$observed =  as.vector(counts[rownames(counts)<=as.Date(end_day), ])
# shinydata$rain = as.vector(rain_mean)
# shinydata$humidity = as.vector(sh_mean)
# shinydata$temperature = as.vector(temp_mean)
shinydata$stringency = as.vector(sindex[as.Date(rownames(sindex)) <= end_day,])
shinydata$testing = as.vector(testing[as.Date(rownames(testing)) <= end_day,])


# Standardise climatic variables
# rain_mean <- (rain_mean - mean(rain_mean)) / sd(rain_mean)
# temp_mean <- (temp_mean - mean(temp_mean)) / sd(temp_mean)
# sh_mean <- (sh_mean - mean(sh_mean)) / sd(sh_mean)


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

# Population
pop <- population(epi_sts)


# HDI by category
HDI_cat <- as.numeric(africa$HDI_Level)
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

# Lag time varying variables
k <- 7
sindex_lag <- xts::lag.xts(sindex, k = 7)
testing_lag <- xts::lag.xts(testing, k = 7)
# temp_mean_lag <- xts::lag.xts(temp_mean, k = 7)
# rain_mean_lag <- xts::lag.xts(rain_mean, k = 7)
# sh_mean_lag <- xts::lag.xts(sh_mean, k = 7)

# MODEL ------------------------------------------------------------------------

# Days used to conduct inference
days_for_inferece <- 90
days_to_predict <- 7
fit_start <- nrow(epi_sts) - days_to_predict - days_for_inferece
fit_end <- nrow(epi_sts) - days_to_predict


start_day <- rownames(counts)[fit_start] 
saveRDS(list(shinydata, africa, start_day, end_day), 
        "output/models/shinydata.rds")


# Best AR1 model with no RE 
f_end <- ~ 1 
# f_ar <- ~ 1 + log(pop) + HDI_cat + LL + sindex_lag + testing_lag + 
#   rain_mean_lag + temp_mean_lag + sh_mean_lag 
# f_ne <- ~ 1 + log(pop) + HDI_cat + LL + sindex_lag + testing_lag 

f_ar <- ~ 1 + log(pop) + HDI_cat + LL + sindex_lag + testing_lag  
f_ne <- ~ 1 + HDI_cat + LL + sindex_lag + testing_lag 

model_basic <- list(
  end = list(f = f_end, offset = population(epi_sts)),
  ar = list(f = f_ar),
  ne = list(f = f_ne, weights = W_powerlaw(maxlag = 9)),
  optimizer = list(stop = list(iter.max = 50)),
  family = "NegBin1",
  subset = fit_start:fit_end)

fit_basic <- hhh4(epi_sts, control = model_basic)
summary(fit_basic)

# Lagged version of previous model

lags <- 14


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
  scale_x_continuous(breaks = 2:lags) +
  theme_gray(base_size = 13) +
  theme(legend.position = "top")

ggsave("figs/figureS10A.pdf", width = 7, height = 5)

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

ggsave("figs/figureS10B.pdf", width = 7, height = 5)


# LAGGED POISSON WITH AND WITHOUT RE WITH OPTIMAL NUMBER OF LAGS  
lag_optimal <- 7
f_end <- ~ 1 
# f_ar <- ~ 1 + log(pop) + HDI_cat + LL + sindex_lag + testing_lag +
#   rain_mean_lag + temp_mean_lag + sh_mean_lag 
f_ar <- ~ 1 + log(pop) + HDI_cat + LL + sindex_lag + testing_lag 
f_ne <- ~ 1 + HDI_cat + LL + sindex_lag + testing_lag 
  
  
model_lag <- list(
  end = list(f = f_end, offset = population(epi_sts)),
  ar = list(f = f_ar),
  ne = list(f = f_ne, weights = W_powerlaw(maxlag = 9)),
  optimizer = list(stop = list(tol = 1e-5, niter = 100), # stop rules
                   regression = list(method = "nlminb"), # for penLogLik
                   variance = list(method = "nlminb")), # for marLogLik
  family = "NegBin1",
  data = list(pop = pop, HDI_cat = HDI_cat,
              SSA = SSA, LL = LL, sindex_lag = sindex_lag,
              # rain_mean_lag = rain_mean_lag, temp_mean_lag = temp_mean_lag, 
              testing_lag = testing_lag #sh_mean_lag = sh_mean_lag
              ),
  subset = fit_start:fit_end,
  funct_lag = poisson_lag, 
  max_lag = lag_optimal)
  
fit_lag <- profile_par_lag(epi_sts, model_lag)
summary(fit_lag, idx2Exp = T)
confint(fit_lag, parm = "overdisp")
wpois <-   fit_lag$distr_lag
  
saveRDS(fit_lag, paste0("output/models/fitted_model_LAG", lag_optimal, ".rds"))
  
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
write.csv(tab, paste0("output/tables/tab_params_LAG", lag_optimal, ".csv"),
          row.names = F)
  
# INTRODUCE REs
f_end <- ~ 1 
# f_ar <- ~ -1 + log(pop) + HDI_cat + LL + sindex_lag + testing_lag +
#   rain_mean_lag + temp_mean_lag + sh_mean_lag + ri()
f_ar <- ~ -1 + log(pop) + HDI_cat + LL + sindex_lag + testing_lag + ri()
f_ne <- ~ 1 + HDI_cat + LL + sindex_lag + testing_lag 
  
model_lag <- list(
  end = list(f = f_end, offset = population(epi_sts)),
  ar = list(f = f_ar),
  ne = list(f = f_ne, weights = W_powerlaw(maxlag = 9)),
  optimizer = list(stop = list(tol = 1e-5, niter = 100), # stop rules
                   regression = list(method = "nlminb"), # for penLogLik
                   variance = list(method = "nlminb")),
  family = "NegBin1",
  data = list(pop = pop, HDI_cat = HDI_cat,
              SSA = SSA, LL = LL, sindex_lag = sindex_lag,
              # rain_mean_lag = rain_mean_lag, temp_mean_lag = temp_mean_lag, 
              # sh_mean_lag = sh_mean_lag,
              testing_lag = testing_lag),
  subset = fit_start:fit_end,
  funct_lag = poisson_lag, 
  par_lag = fit_lag$par_lag,
  max_lag = lag_optimal)
  
fit_lag <- hhh4_lag(epi_sts, model_lag)
saveRDS(fit_lag, paste0("output/models/fitted_model_LAG", lag_optimal, "_RE.rds"))
  
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
write.csv(tab, paste0("output/tables/tab_params_LAG", lag_optimal, "_RE.csv"),
          row.names = F)

# P - VALUES -------------------------------------------------------------------
beta_hat <- fit$coefficients[1:16]
sd_hat <- fit$se[1:16]
all.equal(names(beta_hat), names(sd_hat))

zscores <- beta_hat / sd_hat
pvalues <- 2 * pnorm(abs(zscores), lower.tail = F)
cbind(round(pvalues, 3), round(pvalues, 4))
