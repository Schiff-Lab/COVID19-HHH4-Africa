# LOAD REQUIRED PACKAGES AND FUNCTIONS -----------------------------------------
if (!require("pacman")) install.packages("pacman")
pkgs = c("surveillance", "sf", "tmap", "ggplot2", "dplyr", "hhh4addon") # package names
pacman::p_load(pkgs, character.only = T)

source("R/functions.R")

# LOAD DATA --------------------------------------------------------------------

# Cases at country level
counts <- readRDS("data/processed/daily_cases.rds")

# Final model and auxiliary data 
fit <- readRDS("output/models/fitted_model_LAG7_RE.rds")

# Shapefile for ltlas
africa <- st_read("data/processed/geodata/africa.gpkg") 
africa_full <- st_read("data/original/geodata/africa.gpkg")

# DATA PROCESSING --------------------------------------------------------------
cases <- as_tibble(t(counts))
cases$COUNTRY <- colnames(counts)
cases_to_plot <- cases %>%
  tidyr::gather(key = "time", value = "observed", -COUNTRY) %>% 
  mutate(time = as.Date(time))

inference_days <- fit$control$subset 
inference_dates <- as.Date(rownames(counts)[inference_days])
start_day <- inference_days[1]
end_day <- inference_days[length(inference_days)]

# Predicted mean
fit$terms <- terms(fit)
fitted <- meanHHH(theta = unname(fit$coefficients), model = fit$terms, 
                  subset = start_day:end_day, total.only = F)
fitted$inference_dates = inference_dates
saveRDS(fitted,"output/models/fitted_mean.rds")



# RANDOM EFFECTS MAP -----------------------------------------------------------

reff_mean <- exp(ranef(fit))
africa$reff_mean_AR <- as.numeric(reff_mean[grep("ar.ri", names(reff_mean))])

africa_full <- africa_full %>% 
  full_join(st_drop_geometry(africa[c("name", "reff_mean_AR")]))

reff_mean_AR <- map(x = "reff_mean_AR", shape = africa_full, palette = "-RdYlBu", 
                    breaks = c(0.5, 0.7, 0.8, 0.9, 1, 1.5, 2, 2.5),
                    digits = 2, legend_title = "Mean relative risk", 
                    panel_title = "Random Effects - Within country")
    
tmap_save(reff_mean_AR, "figs/figure3AB.pdf", width = 7, height = 8)
  

# CONTRIBUTION OVER TIME FOR SPECIFIC COUNTRIES --------------------------------

units <- c("Egypt", "Senegal", "South Africa", "Uganda")
units_names <- units
  
ids <- match(units, colnames(fitted$mean))
results <- cbind(fitted$endemic[, ids[1]], fitted$epi.own[, ids[1]], fitted$epi.neighbours[, ids[1]]) %>% 
  rbind(cbind(fitted$endemic[, ids[2]], fitted$epi.own[, ids[2]], fitted$epi.neighbours[, ids[2]])) %>% 
  rbind(cbind(fitted$endemic[, ids[3]], fitted$epi.own[, ids[3]], fitted$epi.neighbours[, ids[3]])) %>% 
  rbind(cbind(fitted$endemic[, ids[4]], fitted$epi.own[, ids[4]], fitted$epi.neighbours[, ids[4]]))

results <- as.data.frame(results)
names(results) <- c("Endemic", "Within", "Between")
results$cname <- rep(units_names, each = length(inference_days))
results$COUNTRY <- rep(units, each = length(inference_days))
results$time <- inference_dates

results %>%
  inner_join(cases_to_plot, by = c("COUNTRY", "time")) %>% 
  select(-COUNTRY) %>% 
  tidyr::gather(key = "Contribution", value = "value", -cname, -time, -observed) %>% 
  mutate(Contribution = factor(Contribution, 
                               levels = c("Between", "Within", "Endemic"))) %>% 
  ggplot(aes(x = time)) +
  geom_area(aes(y = value, fill = Contribution), alpha = 1) +
  geom_line(aes(y = observed), linetype = 1, size = .3) +
  geom_point(aes(y = observed), size = 1) +
  facet_wrap(~ cname, scales = "free_y", ncol = 1) +
  scale_fill_brewer(type = "q", palette= 7)+
  scale_x_date(expand = c(0, 1)) +
  scale_y_continuous(expand = c(0.01, 0)) +
  labs(x = "", y = "Reported cases", fill = "") +
  theme_bw(base_size = 14) +
  theme(legend.position = "top", strip.text = element_text(face = "bold"))

ggsave("figs/figure4A.pdf", width = 6, height = 8, dpi = 320)

# CONTRIBUTION MAP -------------------------------------------------------------

# Model-components contributions over time
total <- apply(fitted$mean, 2, sum)
within <- apply(fitted$epi.own, 2, sum)
between <- apply(fitted$epi.neighbours, 2, sum)
endemic <- apply(fitted$endemic, 2, sum)

africa$within_prop <- within / total 
africa$between_prop <- between / total 
africa$endemic_prop <- endemic / total 

africa_full <- africa_full %>% 
  full_join(st_drop_geometry(africa[c("name", "within_prop", "between_prop",
                                      "endemic_prop")]))

within <- map(x = "within_prop", shape = africa_full, palette = "-RdYlBu",        
                    breaks = seq(0, 1, l = 11),
                    digits = 1, legend_title = "Contribution", 
                    panel_title = "Within cumulative contribution")

between <- map(x = "between_prop", shape = africa_full, palette = "-RdYlBu",         
                     breaks = seq(0, 1, l = 11),
                     digits = 1, legend_title = "Contribution", 
                     panel_title = "Between")

endemic <- map(x = "endemic_prop", shape = africa_full, palette = "-RdYlBu",         
                     breaks = seq(0, 1, l = 11),
                     digits = 1, legend_title = "Contribution", 
                     panel_title = "Endemic")

tmap_save(within, "figs/figure4B.pdf", width = 7, height = 8)


# PREDICTIONS FOR SPECIFIC COUNTRIES -------------------------------------------
pred <- oneStepAhead_hhh4lag(result = fit, tp = c(start_day, end_day) - 1, 
                             type = "final")

quants <- hhh4addon:::quantile.oneStepAhead(pred, 
                                            probs = c(0.025, 0.25, 0.5, 0.75, 0.975))

units_id <- match(units, africa$name)

quants <- quants[ , units_id, ]
quants <- rbind(quants[ , , 1], quants[ , , 2], quants[ , , 3], quants[ , , 4],
                quants[ , , 5]) 
quants <- as_tibble(quants)
quants$quantile <- rep(c("low95", "low50", "median", "up50", "up95"),
                       each = length(inference_days))
quants$time <- rep(inference_dates, time = 5)

code_to_name <- tibble(COUNTRY = units, cname = units_names)

fitted_mean <- predict(fit) %>% 
  as_tibble() %>% 
  select(units) %>%
  mutate(time = inference_dates) %>% 
  tidyr::gather(key = "COUNTRY", value = "mean", -time)
  

quants_insample <- quants %>% 
  tidyr::gather(key = "COUNTRY", value = "value", -quantile, -time) %>%
  tidyr::spread(key = quantile, value = value) %>% 
  inner_join(code_to_name) %>% 
  inner_join(cases_to_plot, by = c("COUNTRY", "time")) %>% 
  inner_join(fitted_mean, by = c("COUNTRY", "time")) %>% 
  select(-COUNTRY) 

nsims <- 10000
sims <- simulate(fit, y.start = observed(fit$stsObj)[(end_day - fit$max_lag + 1):end_day, ], 
                 subset = end_day:nrow(fit$stsObj), 
                 nsim = nsims, simplify = F)

sims <- do.call(rbind, lapply(sims, function(x) observed(x)[ ,units_id]))
sims <- as_tibble(sims)


pred_start <- as.Date(names(cases[,-1])[end_day])
pred_end <- pred_start + nrow(fit$stsObj) - end_day
sims$time <- rep(as.Date(pred_start:pred_end, origin = "1970/01/01"), 
                 time = nsims)

sims <- sims %>% 
  tidyr::gather(key = "COUNTRY", value = "value", -time) %>%
  group_by(COUNTRY, time) %>% 
  summarise(low95 = quantile(value, p = 0.025),
            low50 = quantile(value, p = 0.25),
            mean = mean(value), 
            median = quantile(value, p = 0.5),
            up50 = quantile(value, p = 0.75),
            up95 = quantile(value, p = 0.975)) %>% 
  ungroup() %>% 
  inner_join(code_to_name) %>%  
  inner_join(cases_to_plot, by = c("COUNTRY", "time")) %>% 
  select(-COUNTRY) 

sims <- bind_rows(quants_insample[quants_insample$time == pred_start -1, ], sims)
quants_insample$type <- "IS"
sims$type <- "OUS"


all_preds <- quants_insample %>% 
  bind_rows(sims)


ggplot(all_preds, aes(x = time)) +
  geom_ribbon(aes(ymin = low95, ymax = up95, fill = type), alpha = .3) +
  geom_ribbon(aes(ymin = low50, ymax = up50, fill = type), alpha = .3) +
  geom_line(aes(y = mean, col = type), size = .8) +
  geom_point(aes(y = observed, shape = type), size = 1.5) + 
  facet_wrap(~ cname, scales = "free_y") +
  labs(x = "", y = "Reported Cases") +
  scale_fill_manual("50% - 95% CI", values = c("dodgerblue3", "orange"),
                    labels = c("retrospective", "forecast")) +
  scale_color_manual("Predicted mean", values = c("dodgerblue3", "orange"),
                     labels = c("retrospective", "forecast")) +
  scale_shape_manual("Reported cases", values = c(19, 21), labels = c("", "")) +
  scale_x_date(date_breaks = "2 week", date_labels = "%b %d", 
               expand = c(0, 1)) +
  theme_gray(base_size = 21) +
  theme(strip.text = element_text(face = "bold"),
        legend.position = "top", legend.key=element_blank(),
        legend.text = element_text(size = 18),
        legend.key.size = unit(1, units = "cm")) +
  guides(color = guide_legend(override.aes = list(size = 3)),
         shape = guide_legend(override.aes = list(size = 4)))

ggsave("figs/figure5.pdf", width = 20, height = 10, dpi = 320)

# CONTRIBUTION OVER TIME FOR ALL COUNTRIES -------------------------------------

units <- colnames(fitted$mean)
units_names <- units

results <- cbind(fitted$endemic[, 1], fitted$epi.own[, 1], fitted$epi.neighbours[, 1])

for (i in 2:length(units_names)) {
    results <-  rbind(results, 
                      cbind(fitted$endemic[, i], 
                            fitted$epi.own[, i], 
                            fitted$epi.neighbours[, i])) 
}


results <- as.data.frame(results)
names(results) <- c("Endemic", "Within", "Between")
results$cname <- rep(units_names, each = length(inference_days))
results$COUNTRY <- rep(units, each = length(inference_days))
results$time <- inference_dates

results %>%
  inner_join(cases_to_plot, by = c("COUNTRY", "time")) %>% 
  select(-COUNTRY) %>% 
  tidyr::gather(key = "Contribution", value = "value", -cname, -time, -observed) %>% 
  mutate(Contribution = factor(Contribution, 
                               levels = c("Between", "Within", "Endemic"))) %>% 
  ggplot(aes(x = time)) +
  geom_area(aes(y = value, fill = Contribution), alpha = 1) +
  geom_line(aes(y = observed), linetype = 1, size = .3) +
  geom_point(aes(y = observed), size = 1) +
  facet_wrap(~ cname, scales = "free_y") +
  scale_fill_brewer(type = "q", palette= 7)+
  scale_x_date(expand = c(0, 1)) +
  scale_y_continuous(expand = c(0.01, 0)) +
  labs(x = "", y = "Reported cases", fill = "") +
  theme_bw(base_size = 14) +
  theme(legend.position = "top", strip.text = element_text(face = "bold"))

ggsave("figs/figureS11.pdf", 
       width = 6 * 5, height = 14, dpi = 320)



# PREDICTIONS FOR ALL COUNTRIES ------------------------------------------------
pred <- oneStepAhead_hhh4lag(result = fit, tp = c(start_day, end_day) - 1, 
                             type = "final")

quants <- hhh4addon:::quantile.oneStepAhead(pred, 
                                            probs = c(0.025, 0.25, 0.5, 0.75, 0.975))

units_id <- match(units, africa$name)

quants <- quants[ , units_id, ]
quants <- rbind(quants[ , , 1], quants[ , , 2], quants[ , , 3], quants[ , , 4],
                quants[ , , 5]) 
quants <- as_tibble(quants)
quants$quantile <- rep(c("low95", "low50", "median", "up50", "up95"),
                       each = length(inference_days))
quants$time <- rep(inference_dates, time = 5)

code_to_name <- tibble(COUNTRY = units, cname = units_names)

fitted_mean <- predict(fit) %>% 
  as_tibble() %>% 
  select(units) %>%
  mutate(time = inference_dates) %>% 
  tidyr::gather(key = "COUNTRY", value = "mean", -time)


quants_insample <- quants %>% 
  tidyr::gather(key = "COUNTRY", value = "value", -quantile, -time) %>%
  tidyr::spread(key = quantile, value = value) %>% 
  inner_join(code_to_name) %>% 
  inner_join(cases_to_plot, by = c("COUNTRY", "time")) %>% 
  inner_join(fitted_mean, by = c("COUNTRY", "time")) %>% 
  select(-COUNTRY) 

nsims <- 10000
sims <- simulate(fit, y.start = observed(fit$stsObj)[(end_day - fit$max_lag + 1):end_day, ], 
                 subset = end_day:nrow(fit$stsObj), 
                 nsim = nsims, simplify = F)

sims <- do.call(rbind, lapply(sims, function(x) observed(x)[ ,units_id]))
sims <- as_tibble(sims)


pred_start <- as.Date(names(cases[,-1])[end_day])
pred_end <- pred_start + nrow(fit$stsObj) - end_day
sims$time <- rep(as.Date(pred_start:pred_end, origin = "1970/01/01"), 
                 time = nsims)

sims <- sims %>% 
  tidyr::gather(key = "COUNTRY", value = "value", -time) %>%
  group_by(COUNTRY, time) %>% 
  summarise(low95 = quantile(value, p = 0.025),
            low50 = quantile(value, p = 0.25),
            mean = mean(value), 
            median = quantile(value, p = 0.5),
            up50 = quantile(value, p = 0.75),
            up95 = quantile(value, p = 0.975)) %>% 
  ungroup() %>% 
  inner_join(code_to_name) %>%  
  inner_join(cases_to_plot, by = c("COUNTRY", "time")) %>% 
  select(-COUNTRY) 

ggplot(sims, aes(x = time)) +
  geom_ribbon(aes(ymin = low95, ymax = up95), fill = "dodgerblue3", alpha = .5) +
  geom_line(aes(y = mean), col = "dodgerblue3", size = 1, linetype = 1) +
  geom_point(aes(y = observed)) + 
  facet_wrap(~ cname, scales = "free_y") +
  labs(x = "", y = "Reported cases") +
  scale_x_date(date_breaks = "1 day", date_labels = "%b %d")

sims <- bind_rows(quants_insample[quants_insample$time == pred_start -1, ], sims)
quants_insample$type <- "IS"
sims$type <- "OUS"


all_preds <- quants_insample %>% 
  bind_rows(sims)

saveRDS(all_preds, "output/models/predictions.rds")

ggplot(all_preds, aes(x = time)) +
  geom_ribbon(aes(ymin = low95, ymax = up95, fill = type), alpha = .3) +
  geom_ribbon(aes(ymin = low50, ymax = up50, fill = type), alpha = .3) +
  geom_line(aes(y = mean, col = type), size = .8) +
  geom_point(aes(y = observed, shape = type), size = 1.5) + 
  facet_wrap(~ cname, scales = "free_y") +
  labs(x = "", y = "Reported Cases") +
  scale_fill_manual("50% - 95% CI", values = c("dodgerblue3", "orange"),
                    labels = c("retrospective", "forecast")) +
  scale_color_manual("Predicted mean", values = c("dodgerblue3", "orange"),
                     labels = c("retrospective", "forecast")) +
  scale_shape_manual("Reported cases", values = c(19, 21), labels = c("", "")) +
  scale_x_date(expand = c(0, 1)) +
  theme_gray(base_size = 21) +
  theme(strip.text = element_text(face = "bold"),
        legend.position = "top", legend.key=element_blank(),
        legend.text = element_text(size = 18),
        legend.key.size = unit(1, units = "cm")) +
  guides(color = guide_legend(override.aes = list(size = 3)),
         shape = guide_legend(override.aes = list(size = 4)))

ggsave("figs/figureS12.pdf", 
       width = 20 * 2, height = 10 * 2, dpi = 320)
