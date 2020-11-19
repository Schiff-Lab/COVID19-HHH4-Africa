# LOAD REQUIRED PACKAGES AND FUNCTIONS -----------------------------------------
if (!require("pacman")) install.packages("pacman")
pkgs = c("surveillance", "sf", "tmap", "ggplot2", "dplyr", "hhh4addon") # package names
pacman::p_load(pkgs, character.only = T)

source("R/functions.R")

# LOAD DATA --------------------------------------------------------------------

# Cases at country level
counts <- readRDS("data/processed/daily_cases_plot.rds")

# Shapefile for Africa
africa <- st_read("data/processed/geodata/africa_plot.gpkg") 

# Weather
weather_clean <- readr::read_csv("data/original/AfricaCountries_2020-12-08_ALLEXTRACTEDDATA.csv")

# DATA PREPARATION -------------------------------------------------------------

# Policy variables
sindex <- readRDS("data/processed/stringency_plot.rds")
testing <- readRDS("data/processed/testing_plot.rds")

# Weather variables
rain_mean <- extrac_var_climate(var = "rain_mean", data = weather_clean)
temp_mean <- extrac_var_climate(var = "temp_mean", data = weather_clean)
sh_mean <- extrac_var_climate(var = "sh_mean", data = weather_clean) * 1000

# See what the common dates are for the time varying datasets and censor
# accordingly 
final_dates <- Reduce(intersect, list(rownames(rain_mean), rownames(counts), rownames(sindex)))

counts <- counts[rownames(counts) %in% final_dates, ]
sindex <- sindex[rownames(sindex) %in% final_dates, ]
testing <- testing[rownames(testing) %in% final_dates, ]


# Reshape data for plotting
all_input <- reshape_df(counts, "observed") %>% 
  inner_join(reshape_df(sindex, "sindex")) %>% 
  inner_join(reshape_df(testing, "testing")) %>% 
  inner_join(reshape_df(temp_mean, "temp")) %>%
  inner_join(reshape_df(sh_mean, "sh")) %>% 
  inner_join(reshape_df(rain_mean, "rain")) %>% 
  inner_join(africa[c("name", "Pop2020")], by = c("COUNTRY" = "name")) 

all_input$rp100k <- (all_input$observed / all_input$Pop2020) * 100000

first_day <- as.Date("2020-03-28")
last_day <- as.Date("2020-08-05")

df <- all_input %>% 
  filter(COUNTRY %in% c("Senegal", "Egypt", "Uganda", "South Africa")) %>% 
  group_by(COUNTRY) %>% 
  summarise(time = min(time), value = mean(rp100k))

options(scipen=999)

all_input %>% 
  select(COUNTRY, time, sindex:rain, rp100k) %>% 
  filter(COUNTRY %in% c("Senegal", "Egypt", "Uganda", "South Africa")) %>% 
  tidyr::gather(key = "var", value = "value", -COUNTRY, -time) %>% 
  ggplot(aes(x = time, y = value)) +
  geom_rect(data = df, aes(xmin = first_day, xmax = last_day,
                           ymin = -Inf, ymax = Inf), fill = "grey", alpha = 0.1) +
  geom_rect(data = df, aes(xmin = last_day, xmax = last_day + 7,
                           ymin = -Inf, ymax = Inf), fill = "orange", alpha = 0.1) +
  geom_line(aes(col = COUNTRY)) +
  facet_grid(rows = vars(factor(var, 
                                levels = c("rp100k", "sindex", "testing",
                                           "temp", "rain", "sh"),
                                labels = c(expression("Cases per"~100000),
                                           "Stringency~Index", 
                                           "Testing",
                                           expression("Temperature"~(C^o)),
                                           "Rain~(mm)", 
                                           expression("Humidity"~(g/kg))))),
             scales = "free_y", labeller = label_parsed) +
  labs(x = "Day", y = "", col = "") +
  ggsci::scale_color_nejm() +
  scale_x_date(date_breaks = "2 week", date_labels = "%b %d", expand = c(0, 0)) +
  theme_bw(base_size = 12) +
  theme(legend.position = "top")

ggsave("figs/figure1A.pdf", width = 7, height = 9)

# Time series of time varying variables

# - Daily cases per 100,000
ggplot(all_input, aes(x = time, y = rp100k)) +
  geom_line() +
  facet_wrap(~ COUNTRY, scales = "free_y") + 
  labs(x = "", y = "Daily cases per 100,000") +
  scale_x_date(date_labels = "%b", expand = c(0, 0)) +
  theme_bw(base_size = 13) 

ggsave("figs/figureS1.pdf", width = 18, height = 10)

# - Sindex
ggplot(all_input, aes(x = time, y = sindex)) +
  geom_line() +
  facet_wrap(~ COUNTRY) + 
  labs(x = "", y = "Stringency index") +
  scale_x_date(date_labels = "%b", expand = c(0, 0)) +
  theme_bw(base_size = 13) 

ggsave("figs/figureS2.pdf", width = 18, height = 10)

# - Testing
ggplot(all_input, aes(x = time, y = testing)) +
  geom_line() +
  facet_wrap(~ COUNTRY) + 
  labs(x = "", y = "Testing regime") +
  scale_x_date(date_labels = "%b", expand = c(0, 0)) +
  scale_y_continuous(breaks = 0:3) +
  theme_bw(base_size = 13) 

ggsave("figs/figureS3.pdf", width = 18, height = 10)

# - Temperature
ggplot(all_input, aes(x = time, y = temp)) +
  geom_line() +
  facet_wrap(~ COUNTRY) + 
  labs(x = "", y = expression("Temperature"~(C^o))) +
  scale_x_date(date_labels = "%b", expand = c(0, 0)) +
  theme_bw(base_size = 13) 

ggsave("figs/figureS4.pdf", width = 18, height = 10)

# - Rain
ggplot(all_input, aes(x = time, y = rain)) +
  geom_line() +
  facet_wrap(~ COUNTRY) + 
  labs(x = "", y = "Rain (mm)") +
  scale_x_date(date_labels = "%b", expand = c(0, 0)) +
  theme_bw(base_size = 13) 

ggsave("figs/figureS5.pdf", width = 18, height = 10)

# Humidity
ggplot(all_input, aes(x = time, y = sh)) +
  geom_line() +
  facet_wrap(~ COUNTRY) + 
  labs(x = "", y = "Specific Humidity (g/kg)") +
  scale_x_date(date_labels = "%b", expand = c(0, 0)) +
  theme_bw(base_size = 13) 

ggsave("figs/figureS6.pdf", width = 18, height = 10)

# Map with total cases
cases_total <- reshape_df(counts, "observed") %>% 
  group_by(COUNTRY) %>% 
  summarise(total = sum(observed)) 


africa <- africa %>% 
  inner_join(cases_total, by = c("name" = "COUNTRY"))

africa$rpk100 <- (africa$total / africa$Pop2020) * 100000

rpk <- map(x = "rpk100", shape = africa, palette = "-RdYlBu",
           n = 10, digits = 0, legend_title = "", 
           panel_title = "Total cases reported per 100,000 up to 08-13-2020")

tmap_save(rpk, "figs/figure1B.pdf", width = 7, height = 8)

# Map with predictive scores
fit <- readRDS("output/fitted_model_LAG7_RE.rds") # final model
fit_end <- fit$control$subset[length(fit$control$subset)]
tp <- c(fit_end, nrow(fit$stsObj) - 1)
forecast <- oneStepAhead_hhh4lag(fit, tp = tp, type = "final")
fitScores <- colMeans(scores(forecast, which = "logs", individual = T))
logs <- tibble(name = names(fitScores), score = as.numeric(fitScores))

fit$control$ar$lag <- 1
fit$control$ne$lag <- 1
calibrationTest(fit, which = "logs")

calib_test <- tibble(name = colnames(forecast$observed),
                     z = NA,
                     pvalue = NA)
size <- as.numeric(unique(exp(forecast$psi)))
for (i in 1:ncol(forecast$observed)) {
  test <- x <- calibrationTest(x = forecast$observed[, i], 
                               mu = forecast$pred[, i], 
                               size = size, which = "logs")
  calib_test$z[i] <- test$statistic
  calib_test$pvalue[i] <- test$p.value
}

calibrationTest(x = forecast$observed, 
                mu = forecast$pred, 
                size = size, which = "logs")

africa <- africa %>% 
  left_join(calib_test)

breaks <- c(0, 0.01, 0.05, 0.10, 1)
pv_map <- map(x = "pvalue", shape = africa, 
              palette = "RdBu", 
              labs = c("< 0.01", "0.01 - 0.05", "0.05 - 0.10", "> 0.10"),
              breaks = breaks, digits = 2, legend_title = "P-value", 
              panel_title = "Calibration Test")
pv_map$tm_fill$palette <- c("#8C0C25", "#C7433F", "#EB9072", "#F7F7F7")

tmap_save(pv_map, "figs/figureS13.pdf", width = 7, height = 8)

# Forest plot with random effects 
res <- exp(ranef(fit))
ci_res <- exp(confint(fit))
out_res <-tibble(country = factor(substr(names(res), 12, 100)), 
                 type = toupper(substr(names(res), 1, 2)),
                 estimate = res,
                 low = ci_res[rownames(ci_res) %in% names(res), 1],
                 up = ci_res[rownames(ci_res) %in% names(res), 2])

levels(out_res$country)[levels(out_res$country) == "Democratic Republic of the Congo"] <- "DRC"

out_res$type <- ifelse(out_res$type == "AR", "Within country", "Between country") %>% 
  factor(levels = c("Within country", "Between country"))

ggplot(out_res) +
  geom_vline(xintercept = 1, linetype = 2, col = "blue", size = .6) +
  geom_pointrange(aes(x = estimate, y = reorder(country, desc(country)), xmin = low, xmax = up),
                  size = .2) +
  facet_wrap(~ type, scales = "free_x") +
  labs(x = "Relative risk of random effects", y = "Country") +
  theme_bw() 

ggsave("figs/figure3CD.pdf", width = 7 * 2, height = 8)  

# Term contributions plot ------------------------------------------------------
nterms <- terms(fit)$nGroups
coefs <- coef(fit)[1:nterms]

# AR contributions 
pop_ar <- (coefs["ar.log(pop)"] * log(fit$control$data$pop)) %>% 
  tidymat(type = "Time constant", var = "Population")
HDI_ar <- (coefs["ar.HDI_cat"] * fit$control$data$HDI_cat) %>% 
  tidymat(type = "Time constant", var = "HDI")
LL_ar <- (coefs["ar.LL"] * fit$control$data$LL) %>% 
  tidymat(type = "Time constant", var = "LL")
Sindex_ar <- (coefs["ar.sindex_lag"] * fit$control$data$sindex_lag) %>% 
  tidymat(type = "Time varying", var = "Sindex")
Testing_ar <- (coefs["ar.testing_lag"] * fit$control$data$testing_lag) %>% 
  tidymat(type = "Time varying", var = "Testing")
Rain_ar <- (coefs["ar.rain_mean_lag"] * fit$control$data$rain_mean_lag) %>% 
  tidymat(type = "Time varying", var = "Rain")
Temp_ar <- (coefs["ar.temp_mean_lag"] * fit$control$data$temp_mean_lag) %>% 
  tidymat(type = "Time varying", var = "Temperature")
Sh_ar <- (coefs["ar.sh_mean_lag"] * fit$control$data$sh_mean_lag) %>% 
  tidymat(type = "Time varying", var = "Humidity")

# NE contributions 
pop_ne <- (coefs["ne.log(pop)"] * log(fit$control$data$pop)) %>% 
  tidymat(type = "Time constant", var = "Population")
HDI_ne <- (coefs["ne.HDI_cat"] * fit$control$data$HDI_cat) %>% 
  tidymat(type = "Time constant", var = "HDI")
LL_ne <- (coefs["ne.LL"] * fit$control$data$LL) %>% 
  tidymat(type = "Time constant", var = "LL")
Sindex_ne <- (coefs["ne.sindex_lag"] * fit$control$data$sindex_lag) %>% 
  tidymat(type = "Time varying", var = "Sindex")
Testing_ne <- (coefs["ne.testing_lag"] * fit$control$data$testing_lag) %>% 
  tidymat(type = "Time varying", var = "Testing")

# Put everything together
contrib_ar <- bind_rows(pop_ar, HDI_ar, LL_ar, Sindex_ar, Testing_ar,
                        Rain_ar, Temp_ar, Sh_ar)

contrib_ne <- bind_rows(pop_ne, HDI_ne, LL_ne, Sindex_ne, Testing_ne)

forlegend1 <- expand.grid(time = as.Date("2020-07-15"), 
                          type = "Time constant",
                          var = unique(contrib_ar$var[contrib_ar$type == "Time varying"]),
                          country = unique(contrib_ar$country),
                          contrib = NA) %>% 
  as_tibble()

forlegend2 <- expand.grid(time = as.Date("2020-07-15"), 
                          type = "Time varying",
                          var = unique(contrib_ar$var[contrib_ar$type == "Time constant"]),
                          country = unique(contrib_ar$country),
                          contrib = NA) %>% 
  as_tibble()

A1 <- contrib_ar %>% 
  na.omit() %>% 
  bind_rows(forlegend1) %>% 
  filter(country %in% c("Egypt", "Senegal", "South Africa", "Uganda"),
         type == "Time constant", time >= "2020-07-12") %>% 
  ggplot(aes(x = time, y = contrib)) +
  geom_line(aes(col = var)) +
  facet_grid(country ~ type, scales = "free_x") +
  labs(col = "", x = "", 
       y = expression("Contribution to"~log(lambda[it]))) +
  ylim(c(-0.6, 0.6)) +
  colorblindr::scale_color_OkabeIto(use_black = T) +
  scale_x_date(expand = c(0, 0), 
               breaks = as.Date(c("2020-07-30")),
               labels = c("March to August")) +
  theme_bw() +
  theme(legend.position = "bottom",
        strip.background.y = element_blank(),
        strip.text.y = element_blank())

A2 <- contrib_ar %>% 
  na.omit() %>% 
  bind_rows(forlegend2) %>% 
  filter(country %in% c("Egypt", "Senegal", "South Africa", "Uganda"),
         time >= "2020-03-08", type == "Time varying") %>% 
  ggplot(aes(x = time, y = contrib)) +
  geom_line(aes(col = var)) +
  facet_grid(country ~ type, scales = "free_x") +
  labs(col = "", x = "", y = "") +
  colorblindr::scale_color_OkabeIto(use_black = T) +
  scale_x_date(expand = c(0, 0)) +
  theme_bw() +
  theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5))

B1 <- contrib_ne %>% 
  na.omit() %>% 
  bind_rows(forlegend1) %>% 
  filter(country %in% c("Egypt", "Senegal", "South Africa", "Uganda"),
         type == "Time constant", time >= "2020-07-12") %>% 
  ggplot(aes(x = time, y = contrib)) +
  geom_line(aes(col = var)) +
  facet_grid(country ~ type, scales = "free_x") +
  labs(col = "", x = "", 
       y = expression("Contribution to"~log(phi[it]))) +
  ylim(c(-1, 2.55)) +
  colorblindr::scale_color_OkabeIto(use_black = T) +
  scale_x_date(expand = c(0, 0), 
               breaks = as.Date(c("2020-07-30")),
               labels = c("March to August")) +
  theme_bw() +
  theme(legend.position = "bottom",
        strip.background.y = element_blank(),
        strip.text.y = element_blank())

B2 <- contrib_ne %>% 
  na.omit() %>% 
  bind_rows(forlegend2) %>% 
  filter(country %in% c("Egypt", "Senegal", "South Africa", "Uganda"),
         time >= "2020-03-08", type == "Time varying") %>% 
  ggplot(aes(x = time, y = contrib)) +
  geom_line(aes(col = var)) +
  facet_grid(country ~ type, scales = "free_x") +
  labs(col = "", x = "", y = "") +
  colorblindr::scale_color_OkabeIto(use_black = T) +
  scale_x_date(expand = c(0, 0)) +
  ylim(c(-1, 2.55)) +
  theme_bw() +
  theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5))

ggpubr::ggarrange(A1, A2, B1, B2, labels = c("A", "", "B", ""),
                  widths = c(1, 2), nrow = 1,
                  common.legend = T)

ggsave("figs/figure2.pdf", width = 10, height = 7)  

# MAP FOR TIME CONSTANT VARIABLES ----------------------------------------------
pop <- map(x = "Pop2020", shape = africa, 
           palette = function(x) viridisLite::viridis(x, direction = -1),
           n = 8, digits = 0, legend_title = "Population", 
           panel_title = "Population in 2020")

africa$HDI <- africa$HDI_Level
HDI <- map(x = "HDI", shape = africa, breaks = c(0, 1, 2, 3),
           palette = "-Set3", labs = c("Low", "Medium", "High"),  
           legend_title = "HDI", 
           panel_title = "Human Development Index")

africa$landlock <- as.factor(africa$landlock)
LL <- map(x = "landlock", shape = africa, breaks = c(0, 1),
          labs = c("Coastal", "Landlocked"), palette = "Set3",  
          legend_title = "Access to Coast", 
          panel_title = "Access to Coast")

region <- map(x = "Region", shape = africa, breaks = c(0, 1),
              labs = c("Northen Africa", "Sub-Saharan Africa"), palette = "Set3",  
              legend_title = "UN Region", 
              panel_title = "Regions")

panel <- tmap_arrange(pop, HDI, LL, region, ncol = 2, nrow = 2)
tmap_save(panel, "figs/figureS14.pdf", width = 7 * 2, height = 8 * 2)

ggplot(africa, aes(x = Median_age, y = HDI_2018 * 100)) +
  geom_smooth(method = "lm") +
  geom_point() +
  labs(x = "Median Age (years)", y = "Human Development Index (%)") +
  theme_bw()
ggsave(filename = "figs/S15.pdf", width = 7, height = 5)

# WEIGHTS ----------------------------------------------------------------------
wmat <- getNEweights(fit) %>% 
  as_tibble()
wmat$from <- names(wmat)

wmat <- wmat %>% 
  tidyr::gather(key = "to", value = "wji", -from) 

wmat$from <- factor(wmat$from, levels = sort(unique(wmat$from), decreasing = T))  

wmat$wji[wmat$wji == 0] <- NA
ggplot(wmat, aes(x = to, y = from)) + 
  geom_tile(aes(fill = wji), colour = "white") +
  scale_fill_viridis_c(direction = -1, na.value = "white") +
  scale_x_discrete(position = "top") +
  labs(x = "", y = "") + 
  coord_equal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 0))

ggsave(filename = "figs/S16B.pdf", width = 10.5, height = 9.5)

rho <-  coef(fit)[["neweights.d"]]
rho_low <- confint(fit)["neweights.d", 1]
rho_up <- confint(fit)["neweights.d", 2]

what <- tibble(o = 1:9, w = (1:9) ^ -rho, wlow = (1:9) ^ -rho_low, wup = (1:9) ^ -rho_up)

ggplot(what, aes(x = o, y = w)) +
  geom_ribbon(aes(ymin = wlow, ymax = wup), alpha = .2) +
  geom_point() +
  geom_line() +
  scale_x_continuous(breaks = 1:9) +
  labs(x = "Adjacency order", y = "Non-normalized weights")

ggsave(filename = "figs/S16A.pdf", width = 6, height = 4)