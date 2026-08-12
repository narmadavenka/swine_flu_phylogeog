###### Swine H3 Phylogeography- anthropogenic factors ######
## Predictor processing: Human mobility (state-to-state migration flows) ##

## set working directory and load libraries
rm(list = ls())
wd <- "C:/Users/narma/OneDrive/Documents/PhD_backup/swine_flu/phylogeography/anthropogenic_predictors/"
setwd(wd)
options(scipen = 999)

## load packages
library(readxl)
library(tidyverse)
library(geosphere)
library(rnaturalearth)
library(rnaturalearthhires)
library(Biostrings)

## load state-level shapefile (to match up state-level identifiers)
us_states <- ne_states(country = "United States of America", returnclass = "sf")

## load H3 sequence reference file
ref <- read_tsv(paste0(wd,"../../sequences/tsvs/H3_1990-2026_metadata_final.tsv",sep=""))

####################
#### Migration raw data processing ####

  ## Have files from 2000 (file labeled 2000 is 5 year flow between 1995 and 2000)-2024 with yearly flows, these are in XLS files right now
  ## put into one csv with year, origin, current_residence 

##### SKIP SECTION IF YOU GENERATED THIS INTERMEDIATE FILE #####  
#############################
### Raw migration data processing

# ── Valid state/territory names ───────────────────────────────────────────────
states <- c("Alabama", "Alaska", "Arizona", "Arkansas", "California", "Colorado",
  "Connecticut", "Delaware", "District of Columbia", "Florida", "Georgia",
  "Hawaii", "Idaho", "Illinois", "Indiana", "Iowa", "Kansas", "Kentucky",
  "Louisiana", "Maine", "Maryland", "Massachusetts", "Michigan", "Minnesota",
  "Mississippi", "Missouri", "Montana", "Nebraska", "Nevada", "New Hampshire",
  "New Jersey", "New Mexico", "New York", "North Carolina", "North Dakota",
  "Ohio", "Oklahoma", "Oregon", "Pennsylvania", "Puerto Rico", "Rhode Island",
  "South Carolina", "South Dakota", "Tennessee", "Texas", "Utah", "Vermont",
  "Virginia", "Washington", "West Virginia", "Wisconsin", "Wyoming")

# Extract years from file names 
extract_year <- function(filename) {
  years <- str_extract_all(filename, "20\\d{2}")[[1]]
  if (length(years) == 0) return(NA_integer_)
  as.integer(tail(years, 1))
}

########
#### Get functions to collect information from dataframes, and create long format column in one df year, origin, current_destination #### 
## Mixed formatting- each migration file is separated by year, and has different xls formats ##

# 1995-2000 data, wide

parse_2000 <- function(path, year) {
  df_raw <- read_excel(path, sheet = 1, col_names = FALSE, .name_repair = "minimal")
  
  # Section columns (0-indexed → +1 for R's 1-indexed)
  # Sections at cols: 3-10, 12-22, 24-33, 35-45, 47-57
  section_cols <- c(3:10, 12:22, 24:33, 35:45, 47:57) + 1
  
  # Origin names from row 9 (1-indexed) across section columns
  origins <- as.character(df_raw[9, section_cols])
  
  # Data rows: 13-65 (1-indexed), skipping rows 35 AND 36 (repeated header block)
  # Also skip row 11 (United States summary) 
  all_rows   <- 13:65
  skip_rows  <- c(11, 35, 36)   # United States + repeated header pair
  data_rows  <- setdiff(all_rows, skip_rows)
  
  # Destination state names from col 1
  dest_states <- as.character(df_raw[data_rows, 1][[1]])
  
  # Extract data matrix
  mat <- df_raw[data_rows, section_cols]
  colnames(mat) <- origins
  mat$current_residence <- dest_states
  
  # Pivot to long, clean values
  mat %>%
    pivot_longer(-current_residence, names_to = "origin", values_to = "estimate") %>%
    mutate(
      census_year      = year,
      estimate         = suppressWarnings(as.numeric(estimate)),  # '-' → NA
      reference_period = "5-year (1995-2000)"
    ) %>%
    filter(
      !is.na(estimate),
      !is.na(current_residence),
      !grepl("^\\s*$", current_residence)           # drop blank rows
    ) %>%
    select(census_year, origin, current_residence, estimate, reference_period)
}

# wide format (2005-2023)

parse_wide <- function(path, year) {
  raw <- if (str_ends(path, "\\.xls$")) read_xls(path, col_names = FALSE) else
    read_xlsx(path, col_names = FALSE)
  
  # Find subheader row: first row with >5 "Estimate" cells
  subheader_row <- which(
    map_int(seq_len(min(20, nrow(raw))), function(i) {
      sum(raw[i, ] == "Estimate", na.rm = TRUE)
    }) > 5
  )[1]
  if (is.na(subheader_row)) stop("Could not find Estimate header row")
  
  header_row  <- subheader_row - 1
  state_names <- as.character(raw[header_row, ])
  sub_labels  <- as.character(raw[subheader_row, ])
  
  # Estimate columns where the column header is a valid state/territory
  estimate_cols     <- which(sub_labels == "Estimate" & trimws(state_names) %in% states)
  current_residence <- trimws(state_names[estimate_cols])
  row_header_cols   <- setdiff(seq_len(ncol(raw)), c(estimate_cols, estimate_cols + 1))
  
  skip_patterns <- paste(
    c("^footnote", "^source", "^moe", "^note", "^current residence",
      "^residence 1 year", "^dataset", "^universe", "^table", "^\\*"),
    collapse = "|"
  )
  
  is_data_row <- function(i) {
    for (col in row_header_cols) {
      val <- trimws(as.character(raw[[col]][i]))
      if (!is.na(val) && val != "" && val != "NA")
        return(!str_detect(tolower(val), skip_patterns))
    }
    FALSE
  }
  
get_origin <- function(i) {
    for (col in row_header_cols) {
      val <- trimws(as.character(raw[[col]][i]))
      if (!is.na(val) && val != "" && val != "NA") return(val)
    }
    NA_character_
  }
  
data_rows <- which(map_lgl(seq(subheader_row + 2, nrow(raw)), is_data_row)
  ) + (subheader_row + 1)
  
map_dfr(data_rows, function(i) {
    origin <- get_origin(i)
    if (is.na(origin) || !origin %in% states) return(NULL)
    tibble(
      census_year = year,
      origin = origin,
      current_residence = current_residence,
      estimate = suppressWarnings(
        as.integer(as.numeric(unlist(raw[i, estimate_cols])))
      )
    )
  }) %>%
    filter(!is.na(estimate))   # drop same-state NAs (present in 2023+)
}

# long format (2024)

parse_long <- function(path, year) {
  raw <- read_xlsx(path, col_names = FALSE)
  
  data_start <- which(trimws(as.character(raw[[1]])) %in% states)[1]
  if (is.na(data_start)) stop("Could not find data rows")
  
  raw[data_start:nrow(raw), 1:3] %>%
    setNames(c("current_residence", "origin", "estimate")) |>
    mutate(
      census_year              = year,
      current_residence = trimws(current_residence),
      origin            = trimws(origin),
      estimate          = suppressWarnings(as.integer(estimate))
    ) |>
    filter(
      current_residence %in% states,
      origin %in% states,
      !is.na(estimate)
    ) %>%
    select(census_year, origin, current_residence, estimate)
}

#################
## Execute functions on files ##

parse_migration_file <- function(path, year) {
  if (year == 2024) parse_long(path, year) 
  else if (year == 2000) parse_2000(path,year)
  else parse_wide(path, year)
}

folder <- paste0(wd,"raw_data_sources/state_migration_flows")
files  <- list.files(folder, pattern = "\\.(xls|xlsx)$", full.names = TRUE)

migration <- map_dfr(files, function(path) {
  year <- extract_year(basename(path))
  if (is.na(year)) {
    message("SKIP ", basename(path), " (no year found)"); return(NULL)
  }
  message("Parsing ", basename(path), " (year = ", year, ") ...", appendLF = FALSE)
  df <- tryCatch(
    parse_migration_file(path, year),
    error = function(e) { message(" ERROR: ", e$message); NULL }
  )
  message(" ", nrow(df), " rows")
  df
}) %>%
  arrange(census_year, origin, current_residence)

## filter only to locations present in the fasta files
states <- unique(ref$location)

# filter to only states in the sequences and include state abbreviations needed for MASCOT

migration <- migration %>% filter(origin %in% states) %>% 
  filter(current_residence %in% states) %>% left_join(us_states[, c("name","postal")], by = c("origin" = "name")) %>%
  dplyr::rename(origin_abb = postal) %>%
  left_join(us_states[, c("name","postal")], by = c("current_residence" = "name")) %>%
  dplyr::rename(destination_abb = postal) %>%
  relocate(origin_abb, .after = origin) %>% 
  relocate(destination_abb, .after = current_residence) %>% select(-c(geometry.x, geometry.y))

write.csv(migration, file.path(folder, "statetostate_migration_combined.csv"), row.names = FALSE)

##############################################
## load in migration file  
migration <- read.csv(paste0(wd, "raw_data_sources/state_migration_flows/statetostate_migration_combined.csv"))

## convert into MASCOT epoch format- take average migration between places
  # needs to be log-transformed and standardized

  # do two files here with epochs pre2010 and post 2010

# average migration flows during the epochs

migration_epoch <- migration %>% 
  mutate(epoch = case_when(census_year < 2010 ~ "pre2010",
                           census_year >= 2010 ~ "post2010")) %>%
  group_by(epoch, origin, current_residence) %>%
  mutate(avg_migration = mean(estimate)) %>%
  ungroup()

# log transform and standardize
migration_epoch <- migration_epoch %>% mutate(avg_migration_logstand = log1p(avg_migration)) %>%
 mutate(avg_migration_logstand = scale(avg_migration_logstand))

##### Create matrices for each epoch
epochs <- c("pre2010", "post2010")

# define states in each
pre2010_states <- unique(ref$location)
  #ref %>% 
  #filter(name_cleaned %in% names(seq_1990.4_pre2010)) %>%
  #select(location) %>% distinct %>% pull()

post2010_states <- unique(ref$location)
  #ref %>% 
  #filter(name_cleaned %in% names(seq_1990.4_post2010)) %>%
  #select(location) %>% distinct %>% pull()

## create matrices

make_epoch_matrices <- function(df, epochs, pre_states, post_states) {
  lapply(setNames(epochs, epochs), function(ep) {
    
    # 1. Isolate the data for this specific epoch
    df_ep <- df %>% filter(epoch == ep)
    
    # 2. Pick the correct target states vector
    target_states <- if (ep == "pre2010") pre_states else post_states
    
    # --- DATA CLEANING PRE-CHECK ---
    # Strip spaces and make everything underscores BEFORE matching to prevent typos.
    target_states_clean  <- gsub(" ", "_", target_states)
    df_ep$origin            <- gsub(" ", "_", df_ep$origin)
    df_ep$current_residence <- gsub(" ", "_", df_ep$current_residence)
    
    # 3. Create a blank square matrix explicitly ordered by your state vector
    n_states <- length(target_states_clean)
    mat <- matrix(0, nrow = n_states, ncol = n_states,
                  dimnames = list(target_states_clean, target_states_clean))
    
    # 4. Only keep rows in the data frame that actually match your states
    df_filtered <- df_ep %>% 
      filter(origin %in% target_states_clean & current_residence %in% target_states_clean)
    
    # 5. MATRIX MATCHING FIX:
    # Use a 2-column character matrix to bind coordinates directly to values.
    # This acts like a GPS mapping coordinates directly onto the matrix cells.
    if(nrow(df_filtered) > 0) {
      coords <- as.matrix(df_filtered[, c("origin", "current_residence")])
      mat[coords] <- df_filtered$avg_migration_logstand
    }
    
    return(mat)
  })
}

matrix_list <- make_epoch_matrices(
  df = migration_epoch, 
  epochs = epochs, 
  pre_states = pre2010_states, 
  post_states = post2010_states
)

matrix_list <- list(
  H3_migration_pre2010  = matrix_list$pre2010,
  H3_migration_post2010  = matrix_list$post2010
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