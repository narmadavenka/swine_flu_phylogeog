###### Swine H3 Phylogeography- anthropogenic factors ######
## Predictor processing: Swine network (swine shipment flows between counties) ##

## set working directory and load libraries
rm(list = ls())
wd <- "C:/Users/narma/OneDrive/Documents/PhD_backup/swine_flu/phylogeography/anthropogenic_predictors/"
setwd(wd)
options(scipen = 999)

## load packages
library(data.table)
library(furrr)
library(future)
library(readxl)
library(tidyverse)
library(geosphere)
library(rnaturalearth)
library(rnaturalearthhires)

## load H3 sequence reference file
ref <- read_tsv(paste0(wd,"../../sequences/tsvs/H3_1990-2026_metadata_final.tsv",sep=""))

# load us_states for state metadata
us_states <- ne_states(country = "United States of America", returnclass = "sf")
us_states <- st_drop_geometry(us_states)

#### SKIP THIS ONCE YOU HAVE OUTPUT ####
## load network files
  # read in all files: 

plan(multisession, workers = parallel::detectCores() - 1)

folder_path <- paste0(wd,"/raw_data_sources/swine_networks/swine_networks") #change to directory location

files <- list.files(folder_path, pattern = "\\.network$", full.names = TRUE)

df_list <- future_map(files, ~ fread(
  .x,
  sep        = "\t",    # ← tab, not space
  header     = TRUE,    # ← header row already present in file
  data.table = FALSE
),
.progress = TRUE
)

## create combined dataframe with origin identifiers
networks_df <- bind_rows(df_list, .id = "source_df") %>%
  mutate(source_df = as.numeric(source_df))

## within each source_df, sum across all county shipments to just get state shipments
networks_avg <- networks_df %>% group_by(source_df, oStateAbbr, dStateAbbr) %>%
  mutate(sum_shipments = sum(volume))
networks_avg <- networks_avg %>% select(source_df, oStateAbbr, dStateAbbr, sum_shipments) %>% distinct()

rm(df_list, networks_df) # remove for space
gc()

## average across source_df now
network_final <- networks_avg %>% group_by(oStateAbbr, dStateAbbr) %>%
  mutate(avg_shipments = mean(sum_shipments)) %>% ungroup() %>%
  select(-c(source_df, sum_shipments)) %>%
  distinct()

# only keep states that are in the H3 sequence
states <- unique(ref$location)
state_abbr <- us_states %>% filter(name %in% states) %>% distinct() %>% pull(postal)

network_fin <- network_final %>% filter(oStateAbbr %in% state_abbr) %>% 
  filter(dStateAbbr %in% state_abbr)

## export final dataset (state shipments averaged across networks)
write.csv(network_fin,paste0(wd,"/raw_data_sources/swine_networks/state_swine_shipments_processed.csv",sep=""),row.names = FALSE)

############################
### Create MASCOT predictor 
  # create artificial epochs since this data isn't time varying

network_fin <- read_csv(paste0(wd,"/raw_data_sources/swine_networks/state_swine_shipments_processed.csv",sep=""))

# Define your epochs (pre/post-2010)
epochs <- c("pre2010", "post2010")

# log-transform and standardize values
network_fin <- network_fin %>% mutate(avg_shipments_logstand = log1p(avg_shipments),
                                      avg_shipments_logstand = scale(avg_shipments_logstand))

# create matrix
network_mat <- as.matrix(xtabs(avg_shipments_logstand ~ oStateAbbr + dStateAbbr, 
                         data = network_fin))

# change state abbr. names to full names
dimnames(network_mat) <- lapply(dimnames(network_mat), function(x) {
  # Replace unmatched names with their original names using ifelse
  new_names <- us_states$name[match(x, us_states$postal)]
  ifelse(is.na(new_names), x, new_names)
})

# define states in each
pre2010_states <- unique(ref$location)
  #ref %>% 
  #filter(name_cleaned %in% names(seq_1990.4_pre2010)) %>%
  #select(location) %>% distinct %>% pull()

post2010_states <- unique(ref$location)
  #ref %>% 
  #filter(name_cleaned %in% names(seq_1990.4_post2010)) %>%
  #select(location) %>% distinct %>% pull()

pre2010 <- network_mat[pre2010_states, pre2010_states]
post2010 <- network_mat[post2010_states, post2010_states]

# change these to have gsub with underscore 
rownames(pre2010) <- gsub(" ", "_", rownames(pre2010))
colnames(pre2010) <- gsub(" ", "_", colnames(pre2010))
rownames(post2010) <- gsub(" ", "_", rownames(post2010))
colnames(post2010) <- gsub(" ", "_", colnames(post2010))

matrix_list <- list(
  H3_swineshipment_pre2010  = pre2010,
  H3_swineshipment_post2010  = post2010
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