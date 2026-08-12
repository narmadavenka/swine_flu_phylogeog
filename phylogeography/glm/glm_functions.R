# =============================================================================
# GLM functions
# =============================================================================

# -----------------------------------------------------------------------------
# 1. STATE-LEVEL PATRISTIC DISTANCE
#    For a single tree, compute mean genetic distance between all pairs
#    of states for each clade and epoch, averaging over tip-to-tip distances.
# -----------------------------------------------------------------------------
state_patristic_dist <- function(tree, ref_epoch, focal_states) {
  
  # Keep only tips present in both tree and ref
  keep <- intersect(ref_epoch$tree_label, tree$tip.label)
  if (length(keep) < 2) return(NULL)
  
  tr      <- keep.tip(tree, keep)
  ref_sub <- ref_epoch[match(tr$tip.label, ref_epoch$tree_label), ]
  
  # Check every focal state has at least one tip
  present <- focal_states %in% ref_sub$location
  if (!all(present)) return(NULL)
  
  # All pairwise tip distances
  pd <- cophenetic(tr)
  
  n   <- length(focal_states)
  mat <- matrix(NA_real_, n, n,
                dimnames = list(focal_states, focal_states))
  
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      if (i == j) next
      si <- ref_sub$tree_label[ref_sub$location == focal_states[i]]
      sj <- ref_sub$tree_label[ref_sub$location == focal_states[j]]
      si <- intersect(si, rownames(pd))
      sj <- intersect(sj, colnames(pd))
      if (length(si) == 0 || length(sj) == 0) next
      mat[i, j] <- mean(pd[si, sj, drop = FALSE])
    }
  }
  
  # Reject if any off-diagonal cell is still NA
  if (any(is.na(mat[row(mat) != col(mat)]))) return(NULL)
  
  as.dist(mat)
}

# -----------------------------------------------------------------------------
# 2. Generalized Linear Model 
#    For a single tree, takes mean patristic distance between all pairs
#    of states per clade and epoch approach and regress with hypothesized drivers
# -----------------------------------------------------------------------------
run_clade_epoch_glm <- function(trees,
                                ref,
                                clade_id,
                                epoch,
                                year_cutoff        = 2010,
                                predictors_pre,
                                predictors_post,
                                min_states         = 3,
                                min_tips_per_state = 2,
                                model_type         = "univariate") {
  
  suppressPackageStartupMessages({
    library(lme4)
    library(lmerTest)
  })
  
  stopifnot(epoch %in% c("pre", "post"))
  stopifnot(model_type %in% c("univariate", "multivariate"))
  
  # --- 1. Filter reference (same as run_clade_epoch_mrm) ---
  ref_clade <- ref[ref$clade_final == clade_id, ]
  
  if (epoch == "pre") {
    ref_epoch <- ref_clade[ref_clade$year <  year_cutoff, ]
    pred_list <- predictors_pre
  } else {
    ref_epoch <- ref_clade[ref_clade$year >= year_cutoff, ]
    pred_list <- predictors_post
  }
  
  if (nrow(ref_epoch) == 0) {
    message("  [SKIP] ", clade_id, " ", epoch, ": no sequences in epoch")
    return(NULL)
  }
  
  # --- 2. Focal states (same intersection logic) ---
  state_counts <- table(ref_epoch$location)
  clade_states <- names(state_counts[state_counts >= min_tips_per_state])
  pred_states  <- rownames(pred_list[[1]])
  focal_states <- sort(intersect(clade_states, pred_states))
  
  if (length(focal_states) < min_states) {
    message("  [SKIP] ", clade_id, " ", epoch,
            ": only ", length(focal_states), " states")
    return(NULL)
  }
  
  dropped <- setdiff(clade_states, pred_states)
  if (length(dropped) > 0)
    message("  [NOTE] dropped states: ", paste(dropped, collapse = ", "))
  
  ref_epoch <- ref_epoch[ref_epoch$location %in% focal_states, ]
  
  message("  [GLM]  ", clade_id, " ", epoch,
          " | model: ", model_type,
          " | states: ", length(focal_states),
          " | tips: ", nrow(ref_epoch))
  
  # --- 3. Long-format predictor data frame (all directed pairs) ---
  pred_sub <- lapply(pred_list, function(mat) mat[focal_states, focal_states])
  
  pairs_base <- expand.grid(
    origin = focal_states,
    dest   = focal_states,
    stringsAsFactors = FALSE
  ) |> filter(origin != dest)
  
  for (pred_name in names(pred_sub)) {
    mat <- pred_sub[[pred_name]]
    pairs_base[[pred_name]] <- mapply(
      function(o, d) mat[o, d],
      pairs_base$origin, pairs_base$dest
    )
  }
  
  n_pairs    <- nrow(pairs_base)
  n_trees    <- length(trees)
  pred_names <- names(pred_list)
  
  # Accumulators for pair-level summaries
  pair_patristic_sum <- numeric(n_pairs)
  pair_fitted_sum    <- numeric(n_pairs)
  pair_resid_sum     <- numeric(n_pairs)
  pair_outlier_count <- numeric(n_pairs)
  pair_tree_count    <- 0
  
  coef_results <- vector("list", n_trees)
  
  # --- 4. Loop over trees ---
  for (i in seq_len(n_trees)) {
    
    pd <- state_patristic_dist(trees[[i]], ref_epoch, focal_states)
    if (is.null(pd)) next
    
    pd_mat <- as.matrix(pd)
    pairs_df <- pairs_base
    pairs_df$patristic <- mapply(
      function(o, d) pd_mat[o, d],
      pairs_df$origin, pairs_df$dest
    )
    if (any(is.na(pairs_df$patristic))) next
    
    pair_patristic_sum <- pair_patristic_sum + pairs_df$patristic
    pair_tree_count    <- pair_tree_count + 1
    
    if (model_type == "univariate") {
      
      tree_coefs <- lapply(pred_names, function(pred_name) {
        df_fit <- pairs_df
        df_fit[[pred_name]] <- scale(df_fit[[pred_name]])[, 1]
        formula_str <- paste("patristic ~", pred_name, "+ (1|origin) + (1|dest)")
        tryCatch({
          fit      <- lmer(as.formula(formula_str), data = df_fit,
                           REML = FALSE, control = lmerControl(optimizer = "bobyqa"))
          coef_tbl <- as.data.frame(coef(summary(fit)))
          data.frame(
            tree      = i,
            predictor = pred_name,
            coef      = coef_tbl[pred_name, "Estimate"],
            se        = coef_tbl[pred_name, "Std. Error"],
            t_val     = coef_tbl[pred_name, "t value"],
            pval      = coef_tbl[pred_name, "Pr(>|t|)"],
            stringsAsFactors = FALSE
          )
        }, error = function(e) {
          message("    lmer error tree ", i, " ", pred_name, ": ", e$message)
          NULL
        })
      })
      coef_results[[i]] <- bind_rows(Filter(Negate(is.null), tree_coefs))
      
    } else {
      
      # Multivariate: all predictors together, residuals identify outlier pairs
      df_fit <- pairs_df
      for (pred_name in pred_names)
        df_fit[[pred_name]] <- scale(df_fit[[pred_name]])[, 1]
      
      pred_terms  <- paste(pred_names, collapse = " + ")
      formula_str <- paste("patristic ~", pred_terms, "+ (1|origin) + (1|dest)")
      
      tryCatch({
        fit      <- lmer(as.formula(formula_str), data = df_fit,
                         REML = FALSE, control = lmerControl(optimizer = "bobyqa"))
        coef_tbl <- as.data.frame(coef(summary(fit)))
        pred_rows <- coef_tbl[pred_names, , drop = FALSE]
        
        coef_results[[i]] <- data.frame(
          tree      = i,
          predictor = rownames(pred_rows),
          coef      = pred_rows[, "Estimate"],
          se        = pred_rows[, "Std. Error"],
          t_val     = pred_rows[, "t value"],
          pval      = pred_rows[, "Pr(>|t|)"],
          stringsAsFactors = FALSE
        )
        
        fitted_vals <- fitted(fit)
        resids      <- residuals(fit)
        resid_sd    <- sd(resids)
        
        if (length(fitted_vals) == n_pairs) {
          pair_fitted_sum    <- pair_fitted_sum    + fitted_vals
          pair_resid_sum     <- pair_resid_sum     + resids
          pair_outlier_count <- pair_outlier_count +
            as.integer(abs(resids / resid_sd) > 1.96)
        }
        
      }, error = function(e) {
        message("    lmer error tree ", i, ": ", e$message)
        NULL
      })
    }
  }
  
  if (pair_tree_count == 0) return(NULL)
  
  # --- 5. Coefficient summary across trees ---
  coef_raw <- bind_rows(Filter(Negate(is.null), coef_results))
  if (nrow(coef_raw) == 0) return(NULL)
  
  coef_summary <- coef_raw |>
    group_by(predictor) |>
    summarise(
      n_trees   = n(),
      mean_coef = mean(coef,  na.rm = TRUE),
      sd_coef   = sd(coef,    na.rm = TRUE),
      ci_lower  = quantile(coef, 0.025, na.rm = TRUE),
      ci_upper  = quantile(coef, 0.975, na.rm = TRUE),
      prop_p05  = mean(pval < 0.05, na.rm = TRUE),
      prop_p10  = mean(pval < 0.10, na.rm = TRUE),
      .groups   = "drop"
    ) |>
    mutate(
      clade        = clade_id,
      epoch        = epoch,
      n_states     = length(focal_states),
      states_used  = paste(focal_states, collapse = "|"),
      model_type   = model_type,
      sig_robust   = prop_p05 >= 0.95,
      sig_moderate = prop_p05 >= 0.50
    )
  
  # --- 6. Per-pair summary ---
  # mean_residual > 0 : more genetically distant than predictors expect
  #                     (less viral gene flow than population structure predicts)
  # mean_residual < 0 : more genetically similar than predictors expect
  #                     (more viral gene flow than population structure predicts)
  # prop_outlier      : proportion of trees where |z-residual| > 1.96
  
  pair_summary <- pairs_base |>
    select(origin, dest) |>
    mutate(
      clade          = clade_id,
      epoch          = epoch,
      model_type     = model_type,
      n_trees        = pair_tree_count,
      mean_patristic = pair_patristic_sum / pair_tree_count
    )
  
  if (model_type == "multivariate") {
    pair_summary$mean_fitted   <- pair_fitted_sum   / pair_tree_count
    pair_summary$mean_residual <- pair_resid_sum    / pair_tree_count
    pair_summary$prop_outlier  <- pair_outlier_count / pair_tree_count
    pair_summary <- pair_summary |>
      mutate(
        direction = ifelse(mean_residual > 0,
                           "less_connected_than_expected",
                           "more_connected_than_expected")
      ) |>
      arrange(desc(prop_outlier), desc(abs(mean_residual)))
  }
  
  list(
    coef_summary = coef_summary,
    pair_summary = pair_summary,
    coef_raw     = coef_raw
  )
}
