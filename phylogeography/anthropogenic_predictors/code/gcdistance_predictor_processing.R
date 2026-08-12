###### Swine H3 Phylogeography- anthropogenic factors ######
  ## Predictor processing: Great circle distance ##

## set working directory and load libraries

rm(list = ls())
wd <- "C:/Users/narma/OneDrive/Documents/PhD_backup/swine_flu/phylogeography/anthropogenic_predictors/"
setwd(wd)
options(scipen = 999)

## load packages
library(tidyverse)
library(geosphere)
library(sf)
library(rnaturalearth)
library(rnaturalearthhires)

## load H3 sequence reference file
ref <- read_tsv(paste0(wd,"../../sequences/tsvs/H3_1990-2026_metadata_final.tsv",sep=""))

## read in clade subsets you will be focusing on
#seq_1990.4_pre2010 <- readDNAStringSet("C:/Users/narma/OneDrive/Documents/PhD_backup/swine_flu/sequences/by_clade/06032026/big_downsample/1990_4/1990_4_pre_2010_rep1.fasta", format = "fasta")
#seq_1990.4_post2010 <- readDNAStringSet("C:/Users/narma/OneDrive/Documents/PhD_backup/swine_flu/sequences/by_clade/06032026/big_downsample/1990_4/1990_4_post_2010_rep1.fasta", format = "fasta")

## load state-level shapefile
us_states <- ne_states(country = "United States of America", returnclass = "sf")
us_df <- st_drop_geometry(us_states)

## extract state centroids

centroids <- us_states %>%
  st_transform(4326) %>%  # WGS84
  mutate(
    centroid = st_centroid(geometry),
    lon = st_coordinates(centroid)[, 1],
    lat = st_coordinates(centroid)[, 2],
    name = name,
    state = postal,
  ) %>%
  st_drop_geometry() %>%
  select(state, name, lat, lon)

# filter to states represented in sequences
states <- unique(ref$location)
centroids <- centroids %>% filter(name %in% states)

##########

## calculate great circle distances (km)
dist_mat <- distm(
  centroids[, c("lon", "lat")],
  fun = distGeo
) / 1000

rownames(dist_mat) <- centroids$name
colnames(dist_mat) <- centroids$name

## need to log-transform and standardize all of these
  # keep in matrix

dist_mat_final <- log1p(dist_mat) #log transform, handles zeroes
dist_mat_final <- scale(dist_mat_final)

###########

##### MASCOT epoch predictor generator
# generate separate csv matrices for each epoch
# pre-2010 and post-2010

# define states in each
pre2010_states <- unique(ref$location)
  #ref %>% 
  #filter(name_cleaned %in% names(seq_1990.4_pre2010)) %>%
  #select(location) %>% distinct %>% pull()

post2010_states <- unique(ref$location)
  #ref %>% 
  #filter(name_cleaned %in% names(seq_1990.4_post2010)) %>%
  #select(location) %>% distinct %>% pull()

pre2010 <- dist_mat_final[pre2010_states, pre2010_states]
post2010 <- dist_mat_final[post2010_states, post2010_states]

# change these to have gsub with underscore 
rownames(pre2010) <- gsub(" ", "_", rownames(pre2010))
colnames(pre2010) <- gsub(" ", "_", colnames(pre2010))
rownames(post2010) <- gsub(" ", "_", rownames(post2010))
colnames(post2010) <- gsub(" ", "_", colnames(post2010))

matrix_list <- list(
  H3_gcd_pre2010  = pre2010,
  H3_gcd_post2010  = post2010
  )

## generate subsetted matrices 

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