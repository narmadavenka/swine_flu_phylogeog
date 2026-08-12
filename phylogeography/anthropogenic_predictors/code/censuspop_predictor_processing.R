###### pop H3 Phylogeography- anthropogenic factors ######
## Predictor processing: Human populations from census bureau ##

## set working directory and load libraries
rm(list = ls())
wd <- "C:/Users/narma/OneDrive/Documents/PhD_backup/swine_flu/phylogeography/anthropogenic_predictors/"
setwd(wd)
options(scipen = 999)

## load packages
library(data.table)
library(readxl)
library(tidyverse)
library(tidycensus)
library(rnaturalearth)
library(rnaturalearthhires)
library(sf)
library(Biostrings)

## load state-level shapefile (to match up state-level identifiers)
us_states <- ne_states(country = "United States of America", returnclass = "sf")
us_states <- st_drop_geometry(us_states)

## load sequence reference file
ref <- read_tsv(paste0(wd,"../../sequences/tsvs/H3_1990-2026_metadata_final.tsv",sep=""))

## pull human population from census
  # API key already exists: 46ae53cac283ff1b9c2b3bb7c6d182c98d139364

  # use ACS to get estimates for 2005-2024

acs_years <- c(2005:2019, 2021:2024)  # 2020 not available, 2025 not released yet

get_acs_safe <- function(yr, retries = 5, wait = 10) {
  for (i in seq_len(retries)) {
    result <- tryCatch(
      get_acs(
        geography = "state",
        variables = "B01001_001",
        year      = yr,
        survey    = "acs1"
      ) %>% mutate(year = yr),
      error = function(e) {
        message("Year ", yr, " attempt ", i, " failed: ", e$message)
        Sys.sleep(wait * i)  # increasing wait: 10s, 20s, 30s...
        NULL
      }
    )
    if (!is.null(result)) return(result)
  }
  message("Year ", yr, " failed after ", retries, " attempts — skipping")
  return(NULL)
}

state_pop_acs <- map_dfr(acs_years, function(yr) {
  Sys.sleep(2)          # polite delay between every call- takes a while to run
  get_acs_safe(yr)
})

state_pop_acs %>% distinct(year) %>% arrange(year) 

  # pull in state intercensal populations for 1999-2004
intercensal <- read_xls(
  paste0(wd, "raw_data_sources/state_intercensal_populations_2000-2010.xls"),
  sheet     = 1,
  col_names = FALSE,
  skip      = 8,       # skip header rows
  n_max     = 51       # 50 states + DC, stops before Puerto Rico
) %>%
  select(1, 3:13) %>%   # drop apr2000 base col, jul2010 col
  setNames(c("geography", as.character(2000:2010))) %>%
  mutate(geography = str_remove(geography, "^\\.") %>% str_trim()) %>%
  pivot_longer(-geography, names_to = "year", values_to = "population") %>%
  mutate(year = as.integer(year), population = as.integer(population))

  # keep only 2000-2004 in intercensal and other artifact locations not wanted
intercensal <- intercensal %>% filter(year < 2005,
                                      !geography == "West")

  # rename variables in state_pop_acs
state_pop_acs <- state_pop_acs %>% dplyr::rename(geography = NAME,
                                                 population = estimate) %>%
  select(geography, year, population)

  # combine both together
census_pop <- rbind(intercensal, state_pop_acs)

  # filter to only include states in the sequences
states <- unique(ref$location)

census_pop <- census_pop %>% filter(geography %in% states)

# save dataframe
write.csv(census_pop, paste0(wd,"raw_data_sources/census_pop_2000-2024.csv"),row.names = FALSE)

##### MASCOT and MRM epoch predictor generator
  # generate separate csv matrices for each epoch
  # pre-2010 and post-2010

# average populations during the epochs

census_pop_epoch <- census_pop %>% 
  mutate(epoch = case_when(year < 2010 ~ "pre2010",
                           year >= 2010 ~ "post2010")) %>%
  group_by(epoch, geography) %>%
  mutate(avg_pop = mean(population)) %>%
  ungroup()

# log-transform and standardize
census_pop_epoch <- census_pop_epoch %>% mutate(avg_pop_logstand = log1p(avg_pop)) %>%
  mutate(avg_pop_logstand = scale(avg_pop_logstand))

# remove year and population, and turn into two matrices
  # also convert state names to match location names by putting underscores

census_pop_epoch <- census_pop_epoch %>% 
  select(-c(year,population,avg_pop)) %>% 
  distinct()

# create matrices- should be 4
  # censuspop_origin_pre_2010
  # censuspop_origin_post_2010
  # censuspop_destination_pre_2010
  # censuspop_destination_post_2010

# first create two dataframes for each epoch
epochs <- c("pre2010", "post2010")

census_pop_pre2010 <- census_pop_epoch %>%
  filter(epoch == "pre2010")
census_pop_post2010 <- census_pop_epoch %>%
  filter(epoch == "post2010")

epoch_census <- list(
  pre2010  = census_pop_pre2010,
  post2010 = census_pop_post2010
)

# create state conditions based on the sequence subsets

pre2010_states <- unique(ref$location)
  
  #ref %>% 
  #filter(name_cleaned %in% names(seq_1990.4_pre2010)) %>%
  #select(location) %>% distinct %>% pull()

post2010_states <- unique(ref$location)
  
  #ref %>% 
  #filter(name_cleaned %in% names(seq_1990.4_post2010)) %>%
  #select(location) %>% distinct %>% pull()

epoch_states <- list(
  pre2010  = pre2010_states,
  post2010 = post2010_states
)

pop_origin_matrices <- list()
pop_dest_matrices <- list()

# create origin and destination matrices for both datasets
for (ep in epochs) {
  
  states_ep <- epoch_states[[ep]]
  pop_ep    <- epoch_census[[ep]]
  
  # Build named vector over epoch-specific states (NA if state missing)
  pop_vector <- setNames(
    pop_ep$avg_pop_logstand[match(states_ep, pop_ep$geography)],
    states_ep
  )
  
  # ORIGIN matrix [i, j] = population of origin state i (repeated across columns)
  origin_mat <- matrix(
    rep(pop_vector, times = length(states_ep)),
    nrow = length(states_ep),
    ncol = length(states_ep),
    dimnames = list(states_ep, states_ep)
  )
  
  # DESTINATION matrix [i, j] = population of destination state j (repeated across rows)
  dest_mat <- matrix(
    rep(pop_vector, each = length(states_ep)),
    nrow = length(states_ep),
    ncol = length(states_ep),
    dimnames = list(states_ep, states_ep)
  )
  
  # Replace spaces with underscores in matrix dimension names after values are matched
  dimnames(origin_mat) <- list(gsub(" ", "_", states_ep), gsub(" ", "_", states_ep))
  dimnames(dest_mat)   <- list(gsub(" ", "_", states_ep), gsub(" ", "_", states_ep))
  
  pop_origin_matrices[[ep]] <- origin_mat
  pop_dest_matrices[[ep]]   <- dest_mat
}

# save the matrices 
matrix_list <- list(
  H3_pop_pre2010_origin  = pop_origin_matrices$pre2010,
  H3_pop_pre2010_dest    = pop_dest_matrices$pre2010,
  H3_pop_post2010_origin = pop_origin_matrices$post2010,
  H3_pop_post2010_dest   = pop_dest_matrices$post2010
)

lapply(names(matrix_list), function(x) {
  con <- file(paste0(wd, "processed_data_sources/", x, ".csv"), open = "wb")
  write.csv(
    matrix_list[[x]],
    file      = con,
    row.names = TRUE,
    quote     = FALSE
  )
  close(con)
})