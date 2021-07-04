# LOAD REQUIRED PACKAGES AND FUNCTIONS -----------------------------------------
if (!require("pacman")) install.packages("pacman")
pkgs = c("dplyr", "sf") # package names
pacman::p_load(pkgs, character.only = T)

# LOAD DATA --------------------------------------------------------------------

# Cases at country level  (JHU data)
cases <- readr::read_csv("https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_global.csv")
readr::write_csv(cases, "data/original/cases_raw.csv")


# Shapefile for Africa
africa <- st_read("data/original/geodata/africa.gpkg") 

# Policy index
policy <- readr::read_csv("https://raw.githubusercontent.com/OxCGRT/covid-policy-tracker/master/data/OxCGRT_latest.csv")
readr::write_csv(policy, "data/original/policy_raw.csv")

# DATA PREPARATION -------------------------------------------------------------

# Remove Countries with missing policy data/ not-reporting cases
countries_to_remove <-  c("Madagascar", "Western Sahara", "Guinea-Bissau", 
                          "Equatorial Guinea", "Tanzania")
africa <- africa[!(africa$name %in% countries_to_remove), ] 

# Clean JHU data :
#   - rename countries according to shapefile
#   - keep only countries that we are going to analyse
#   - calculate daily new cases
#   - replace the 4 negatives with 0

cases <- cases %>% 
  select(-`Province/State`, -Lat, -Long) %>% 
  rename(country = `Country/Region`) %>% 
  mutate(country = case_when(
    country == "Congo (Kinshasa)" ~ "Democratic Republic of the Congo",
    country == "Gambia" ~ "The Gambia",
    country == "Eswatini" ~ "Swaziland",
    country == "Congo (Brazzaville)" ~ "Republic of Congo",
    country == "Gambia" ~ "The Gambia",
    TRUE ~ country
  )) %>% 
  filter(country %in% africa$name)

# Check that the order of cases and countries in the shapefile are the same
cases <- cases[order(cases$country), ]
all(cases$country == africa$name)

# Reshape cases with times in the row and geographical units in the column
# so each column is the cases time series for the county
# the results should be of dimension T x I 

# First change the names of the columns to proper dates
dateNames <- function(x) {
  first <- substr(x, 1, nchar(x) - 2)
  last <- substr(x, nchar(x) - 1, nchar(x))
  paste0(first, 20, last)
}
names(cases)[-1] <- dateNames(names(cases)[-1]) %>% 
  as.Date(format = "%m/%d/%Y") %>% 
  as.character()
counts <- t(cases[, -1]) 
colnames(counts) <- cases$country
counts <- apply(counts, 2, diff)
counts[counts < 0] <- 0

# Clean policy data
policy_clean <- policy %>% 
  select(country = CountryName, date = Date, 
         testing = `H2_Testing policy`, sindex = StringencyIndex) %>% 
  mutate(date = as.Date(as.character(date), format = "%Y%m%d"),
         country = case_when(
           country == "Democratic Republic of Congo" ~ "Democratic Republic of the Congo",
           country == "Gambia" ~ "The Gambia",
           country == "Eswatini" ~ "Swaziland",
           country == "Congo" ~ "Republic of Congo",
           TRUE ~ country
         )) %>% 
  filter(country %in% africa$name)
  
testing <- policy_clean %>% 
  select(-sindex) %>% 
  tidyr::spread(key = country, value = testing) %>% 
  select(-date) %>% 
  as.matrix()

testing[is.na(testing)] <- 0

rownames(testing) <- unique(as.character(policy_clean$date))

sindex <- policy_clean %>% 
  select(-testing) %>% 
  tidyr::spread(key = country, value = sindex) %>% 
  select(-date) %>% 
  as.matrix()

sindex[is.na(sindex)] <- 0

rownames(sindex) <- unique(as.character(policy_clean$date))

# Weather data
# weather_clean <- readr::read_csv("data/original/AfricaCountries_2020-12-08_ALLEXTRACTEDDATA.csv")
  

# See what the common dates are for the time varying datasets and censor
# accordingly 
# final_dates <- Reduce(intersect, list(as.character(weather_clean$Date),
#                                       rownames(counts),
#                                       rownames(sindex)))
final_dates <- Reduce(intersect, list(rownames(counts), rownames(sindex)))


counts <- counts[rownames(counts) %in% final_dates, ]
sindex <- sindex[rownames(sindex) %in% final_dates, ]
testing <- testing[rownames(testing) %in% final_dates, ]

# Save the cleaned data 
saveRDS(sindex, "data/processed/stringency.rds")
saveRDS(testing, "data/processed/testing.rds")
saveRDS(counts, "data/processed/daily_cases.rds")
st_write(africa, "data/processed/geodata/africa.gpkg", delete_dsn = T)

