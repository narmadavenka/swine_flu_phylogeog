############
#### H3 phylogeography: Set up spatial and temporal models that evaluate role of predictors in driving genetic events

# run separately for each epoch
  # 1. Univariate analysis
  # 2. Generalized Linear Model
      # Look at the role of great circle distance, human mobility (migration), 
      # human population, swine population, swine shipments on mean genetic 
      # distance

rm(list = ls())
wd <- "C:/Users/narma/OneDrive/Documents/PhD_backup/swine_flu/"
setwd(wd)

# load packages
library(ape)
library(phangorn)
library(ecodist)
library(tidyverse)
library(parallel)
library(lme4)
library(lmerTest)

# load files
##############
  # MCC tree 
tree <- read.nexus(paste0(wd,"trees/mcc_1990_v3.trees"))

  # trimmed posterior trees to build migration rates from (get phylogenetic uncertainty)
trees <- read.nexus(paste0(wd,"trees/H3_thinned_1000.trees"))

  # ref
ref <- read_tsv(paste0(wd,"sequences/tsvs/H3_1990-2026_metadata_final.tsv")) 
to_keep <- unique(ref$tree_label)

#change locations to have _ in spaces to match what is in predictors
ref$location <- gsub(" ", "_", ref$location)

# only keep sequences in ref
trees_clade <- lapply(trees, function(tr) {
  # Only keep tips that are both in to_keep AND present in this tree
  tips_present <- intersect(to_keep, tr$tip.label)
  keep.tip(tr, tips_present)
})
class(trees_clade) <- "multiPhylo"

tips_present <- intersect(to_keep, tree$tip.label)
tree <- keep.tip(tree, tips_present)

#### Univariate analyses
  # Association indices, parsimony score, and monophyletic clade size (MC)

### first analysis: do closely related sequences tend to come from the same state
  ## use mcc tree

  # define states, years, and clades
tip_states <- setNames(ref$location[match(tree$tip.label, ref$tree_label)],
                       tree$tip.label)
tip_years <- setNames(
  ref$year[match(tree$tip.label, ref$tree_label)],
  tree$tip.label
)
#### Operation
get_desc_tip_sets <- function(tree) {
  n_tips  <- length(tree$tip.label)
  n_nodes <- tree$Nnode
  lapply(seq_len(n_nodes), function(i) {
    desc <- phytools::getDescendants(tree, i + n_tips)
    tip_indices <- desc[desc <= n_tips]
    if (length(tip_indices) == 0) return(NULL)
    tree$tip.label[tip_indices]
  })
}

# ── Vectorised AI using pre-computed descendant sets ─────────────────────────
association_index_fast <- function(desc_sets, states) {
  ai_vals <- vapply(desc_sets, function(tips) {
    if (is.null(tips)) return(0)
    tips <- tips[tips %in% names(states)]
    n    <- length(tips)
    if (n == 0) return(0)
    max_count <- max(table(states[tips]))
    (1 - max_count / n) / (2^(n - 1))
  }, numeric(1))
  sum(ai_vals)
}

# ── Main loop ─────────────────────────────────────────────────────────────────
clades     <- unique(ref$clade_final[!is.na(ref$clade_final)])
results_AI <- list()
n_cores    <- max(1, detectCores() - 1)  # use parallel cores for permutations

for (cl in clades) {
  cat(sprintf("\n-- Clade: %s --\n", cl))
  
  clade_tips <- ref$tree_label[ref$clade_final == cl & 
                                 ref$tree_label %in% tree$tip.label &
                                 !is.na(ref$clade_final)]
  
  if (length(clade_tips) < 5) {
    cat("  Skipping — fewer than 5 tips\n"); next
  }
  
  clade_tree   <- keep.tip(tree, clade_tips)
  clade_states <- setNames(
    ref$location[match(clade_tree$tip.label, ref$tree_label)],
    clade_tree$tip.label
  )
  clade_states <- clade_states[!is.na(clade_states)]
  if (length(clade_states) < 5) {
    cat("  Skipping — too few tips with location\n"); next
  }
  clade_tree <- keep.tip(clade_tree, names(clade_states))
  
  n_states_obs <- length(unique(clade_states))
  cat(sprintf("  Tips: %d  States: %s\n",
              length(clade_states),
              paste(sort(unique(clade_states)), collapse=", ")))
  if (n_states_obs < 2) {
    cat("  Skipping — only one state\n"); next
  }
  
  # Pre-compute descendant sets ONCE for this clade tree
  desc_sets   <- get_desc_tip_sets(clade_tree)
  observed_AI <- association_index_fast(desc_sets, clade_states)
  
  # Parallel permutation test
  set.seed(42)
  null_AI <- replicate(999, {
    shuffled <- setNames(sample(clade_states), names(clade_states))
    association_index_fast(desc_sets, shuffled)
  })
  
  p_val <- mean(null_AI <= observed_AI)
  cat(sprintf("  AI = %.4f  p = %.4f\n", observed_AI, p_val))
  
  results_AI[[cl]] <- data.frame(
    clade        = cl,
    n_tips       = length(clade_states),
    n_states     = n_states_obs,
    observed_AI  = observed_AI,
    null_AI_mean = mean(null_AI),
    null_AI_sd   = sd(null_AI),
    p_value      = p_val,
    significant  = p_val < 0.05
  )
}

results_AI_df <- bind_rows(results_AI)
print(results_AI_df)

#### Next get the relationship between genetic distance and these other factors
###############
  # load predictors
    # get them all into matrices separated by epoch

  # set directory
pred_dir <- paste0(wd,"phylogeography/anthropogenic_predictors/processed_data_sources/")   # <-- change this

# ── Helper: read a predictor CSV into a clean named matrix ───────────────────
read_pred_matrix <- function(filename, pred_dir) {
  path <- file.path(pred_dir, filename)
  mat  <- read.csv(path, row.names = 1, check.names = FALSE)
  mat  <- as.matrix(mat)
  rownames(mat) <- trimws(rownames(mat))
  colnames(mat) <- trimws(colnames(mat))
  return(mat)
}

  predictors_pre <- list(
    gcd             = read_pred_matrix("H3_gcd_pre2010.csv",            pred_dir),
    human_migration = read_pred_matrix("H3_migration_pre2010.csv",      pred_dir),
    human_pop_orig  = read_pred_matrix("H3_pop_pre2010_origin.csv",     pred_dir),
    human_pop_dest  = read_pred_matrix("H3_pop_pre2010_dest.csv",       pred_dir),
    swine_pop_orig  = read_pred_matrix("H3_swine_pre2010_origin.csv",   pred_dir),
    swine_pop_dest  = read_pred_matrix("H3_swine_pre2010_dest.csv",     pred_dir),
    swine_shipment  = read_pred_matrix("H3_swineshipment_pre2010.csv",  pred_dir)
  )
  
  predictors_post = list(
    gcd             = read_pred_matrix("H3_gcd_post2010.csv",           pred_dir),
    human_migration = read_pred_matrix("H3_migration_post2010.csv",     pred_dir),
    human_pop_orig  = read_pred_matrix("H3_pop_post2010_origin.csv",    pred_dir),
    human_pop_dest  = read_pred_matrix("H3_pop_post2010_dest.csv",      pred_dir),
    swine_pop_orig  = read_pred_matrix("H3_swine_post2010_origin.csv",  pred_dir),
    swine_pop_dest  = read_pred_matrix("H3_swine_post2010_dest.csv",    pred_dir),
    swine_shipment  = read_pred_matrix("H3_swineshipment_post2010.csv", pred_dir)
  )

##############

## source function
source(paste0(wd,"phylogeography/glm/glm_functions.R"))

### Run GLM to find state-to-state pairs that drive these patterns in particular

glm_results <- list()

for (cl in clades) {
  for (ep in EPOCHS) {
    key <- paste(cl, ep, sep = "_")
    message("\n--- GLM: ", key, " ---")
    
    glm_results[[key]] <- run_clade_epoch_glm(
      trees              = trees_clade,
      ref                = ref,
      clade_id           = cl,
      epoch              = ep,
      year_cutoff        = YEAR_CUTOFF,
      predictors_pre     = predictors_pre,
      predictors_post    = predictors_post,
      min_states         = MIN_STATES,
      min_tips_per_state = MIN_TIPS_PER_STATE,
      model_type         = "multivariate"
    )
  }
}

# Combine outputs
coef_all  <- bind_rows(lapply(glm_results, `[[`, "coef_summary"))
pairs_all <- bind_rows(lapply(glm_results, `[[`, "pair_summary"))

# Save
write.csv(coef_all,  paste0(wd,"phylogeography/glm/glm_coef_summary.csv"),  row.names = FALSE)
write.csv(pairs_all, paste0(wd,"phylogeography/glm/glm_pair_summary.csv"),   row.names = FALSE)

# Quick look at outlier pairs
pairs_all %>%
  filter(prop_outlier > 0.5) %>%
  arrange(desc(prop_outlier)) %>%
  select(clade, epoch, origin, dest, mean_residual, prop_outlier, direction)