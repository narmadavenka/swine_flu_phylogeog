###### Swine H3 Phylogeography- anthropogenic factors ######
## Predictor processing: Commercial swine populations at origin and destination states ##

## set working directory and load libraries
rm(list = ls())
wd <- "C:/Users/narma/OneDrive/Documents/PhD_backup/swine_flu/phylogeography/anthropogenic_predictors/"
setwd(wd)
options(scipen = 999)

## load packages
library(readxl)
library(readr)
library(tidyverse)
library(geosphere)
library(rnaturalearth)
library(rnaturalearthhires)
library(Biostrings)
library(sf)

## load state-level shapefile (to match up state-level identifiers)
us_states <- ne_states(country = "United States of America", returnclass = "sf")
us_states <- st_drop_geometry(us_states)

## create column FIPS which matches up with swine dataset

us_states$FIPS <- gsub("US", "", us_states$fips)

## load H3 sequence reference file
ref <- read_tsv(paste0(wd,"../../sequences/tsvs/H3_1990-2026_metadata_final.tsv",sep=""))

## load swine data 
swine <- read_csv(paste0(wd,"raw_data_sources/CoA_hog_inventory_1999-2025.csv",sep=""))

####################
#### Swine raw data processing ####

swine <- swine %>% 
  dplyr::select(Year, Period, State, state_code = `State ANSI`, Value) %>% 
  distinct()

# Make Value numeric
swine$Value <- gsub(",","",swine$Value)
swine$Value <- as.numeric(swine$Value) #the ones coming up as NA are those with (D)

## filter only to locations present in the fasta files
states <- unique(ref$location)
state_codes <- us_states %>% filter(name %in% states) %>% 
  select(FIPS) %>% distinct() %>% pull()
swine <- swine %>% filter(state_code %in% state_codes) %>%
  left_join(us_states[, c("FIPS","postal")], by = c("state_code" = "FIPS"))
## average across the different months per year
swine_fin <- swine %>% dplyr::group_by(Year, postal) %>% 
  mutate(swine_no = mean(Value)) %>% distinct()

swine_fin <- swine_fin %>% select(Year, State, postal, swine_no) %>% distinct()

##############################################
## convert into MASCOT epoch format- take average within each epoch
# needs to be log-transformed and standardized
# for each state pairing, have an origin and destination state predictor

# do two files here:
# 1) with epochs pre2010 and post 2010

# Define your epochs (pre/post-2010)
epochs <- c("pre2010", "post2010")

# average swine numbers over the epochs

swine_epoch <- swine_fin %>% 
  mutate(epoch = case_when(Year < 2010 ~ "pre2010",
                           Year >= 2010 ~ "post2010")) %>%
  group_by(epoch, postal) %>%
  mutate(avg_swine = mean(swine_no)) %>%
  ungroup()

# merge back the names onto the df
swine_epoch <- swine_epoch %>% left_join(us_states[ ,c("name", "postal")]) %>%
  select(-c(postal, State, Year)) %>% distinct()

# log-transform and standardize
swine_epoch <- swine_epoch %>% mutate(avg_swine_logstand = log1p(avg_swine)) %>%
  mutate(avg_swine_logstand = scale(avg_swine_logstand))

swine_pre2010 <- swine_epoch %>%
  filter(epoch == "pre2010")
swine_post2010 <- swine_epoch %>%
  filter(epoch == "post2010")

epoch_census <- list(pre2010  = swine_pre2010, post2010 = swine_post2010)

# create state conditions based on the sequence subsets

pre2010_states <- unique(ref$location)

post2010_states <- unique(ref$location)

epoch_states <- list(
  pre2010  = pre2010_states,
  post2010 = post2010_states
)

swine_origin_matrices <- list()
swine_dest_matrices <- list()

# create origin and destination matrices for both datasets
for (ep in epochs) {
  
  states_ep <- epoch_states[[ep]]
  pop_ep    <- epoch_census[[ep]]
  
  # Build named vector over epoch-specific states (NA if state missing)
  pop_vector <- setNames(
    pop_ep$avg_swine_logstand[match(states_ep, pop_ep$name)],
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
  
  swine_origin_matrices[[ep]] <- origin_mat
  swine_dest_matrices[[ep]]   <- dest_mat
}

# save the matrices 
matrix_list <- list(
  H3_swine_pre2010_origin  = swine_origin_matrices$pre2010,
  H3_swine_pre2010_dest    = swine_dest_matrices$pre2010,
  H3_swine_post2010_origin = swine_origin_matrices$post2010,
  H3_swine_post2010_dest   = swine_dest_matrices$post2010
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