#' Convert binary indicator vector to string key
#' @param ind Numeric vector of 0s and 1s.
#' @return A concatenated string representing the model structure.
#' @noRd
ind_to_key <- function(ind) {
  paste(ind, collapse="")
}

#' Look up ELPD in the evaluation cache
#' @param cache List containing keys, indicators, and logml values.
#' @param key String key of the model.
#' @return The ELPD value if found, otherwise NULL.
#' @noRd
cache_lookup <- function(cache, key) {
  idx <- which(cache$key == key)
  if(length(idx) == 0) return(NULL)
  return(cache$logml[idx])
}

#' Calculate Scott-Berger Log-Prior for model space
#' @param ind Numeric vector of 0s and 1s.
#' @param p Total number of predictors.
#' @return Numeric log-prior probability.
#' @noRd
log_prior_SB <- function(ind, p) {
  k <- sum(ind)
  if (k == 0 || k == p) return(-log(p + 1))
  return(-log(p + 1) - lchoose(p, k))
}

#' Audit and clean Standard Deviation columns to avoid singularities
#' @param data Data frame containing the predictors and SDs.
#' @param predictors Character vector of predictor names.
#' @param sd_cols Character vector of SD column names.
#' @param threshold_macro Numeric threshold to classify a variable as having prediction error.
#' @param threshold_micro Numeric threshold to replace near-zero SDs.
#' @param verbosity Integer for message printing.
#' @return A list with the cleaned data and a logical vector of clean variables.
#' @noRd
audit_and_clean_errors <- function(data, predictors, sd_cols, threshold_macro = 0.05, threshold_micro = 0.02, verbosity = 1) {
  p <- length(predictors)
  is_clean <- logical(p)
  names(is_clean) <- predictors

  micro_corrections <- 0

  for (j in 1:p) {
    sd_col <- sd_cols[j]
    media_sd <- mean(data[[sd_col]], na.rm = TRUE)

    if (media_sd < threshold_macro) {
      is_clean[j] <- TRUE # No error variable (OLS)
    } else {
      is_clean[j] <- FALSE # Variable with error (sigma)

      idx_reemplazo <- which(data[[sd_col]] < threshold_micro)
      if (length(idx_reemplazo) > 0) {
        data[[sd_col]][idx_reemplazo] <- threshold_micro
        micro_corrections <- micro_corrections + length(idx_reemplazo)
      }
    }
  }

  if (verbosity > 0) {
    cat(sprintf(">>> Error-free variables: %d | Variables with error (sigma): %d\n", sum(is_clean), sum(!is_clean)))
    if (micro_corrections > 0) {
      cat(sprintf(">>> %d error values were adjusted to the safe threshold of %.2f to avoid singularities.\n", micro_corrections, threshold_micro))
    }
  }

  return(list(data = data, is_clean = is_clean))
}

# ------------------------------------------------------------------------------
# 2. MODEL PROPOSAL FUNCTIONS (Tabu vs Classic)
# ------------------------------------------------------------------------------

#' Propose a new model structure (Classic MC3)
#' @param ind Current model indicator vector.
#' @return Proposed model indicator vector.
#' @noRd
propose_model <- function(ind) {
  new_ind <- ind
  move <- sample(c("flip","swap"), 1, prob=c(0.7,0.3))

  if(move == "flip") {
    j <- sample(1:length(ind), 1)
    new_ind[j] <- 1 - ind[j]
  } else {
    ones <- which(ind == 1)
    zeros <- which(ind == 0)
    if(length(ones) > 0 && length(zeros) > 0){
      j1 <- sample(ones, 1)
      j0 <- sample(zeros, 1)
      new_ind[j1] <- 0
      new_ind[j0] <- 1
    } else {
      j <- sample(1:length(ind), 1)
      new_ind[j] <- 1 - ind[j]
    }
  }
  return(new_ind)
}

#' Propose a new model structure with Tabu Search (Heuristic MC3)
#' @param ind_current Current model indicator vector.
#' @param cache List of evaluated models to avoid revisiting.
#' @param max_attempts Maximum number of proposal attempts.
#' @return Proposed model indicator vector.
#' @noRd
propose_model_tabu <- function(ind_current, cache, max_attempts = 50) {
  p <- length(ind_current)
  for (i in 1:max_attempts) {
    ind_prop <- ind_current
    if (runif(1) < 0.5 && p > 1) {
      idx0 <- which(ind_prop == 0)
      idx1 <- which(ind_prop == 1)
      if (length(idx0) > 0 && length(idx1) > 0) {
        s0 <- sample(idx0, 1); s1 <- sample(idx1, 1)
        ind_prop[s0] <- 1; ind_prop[s1] <- 0
      } else {
        s <- sample(1:p, 1)
        ind_prop[s] <- 1 - ind_prop[s]
      }
    } else {
      s <- sample(1:p, 1)
      ind_prop[s] <- 1 - ind_prop[s]
    }
    key_prop <- ind_to_key(ind_prop)
    if (is.na(match(key_prop, cache$key))) return(ind_prop)
  }
  return(ind_current)
}

# ------------------------------------------------------------------------------
# 3. MODEL EVALUATION
# ------------------------------------------------------------------------------

#' Evaluate a single model via LOO Pseudo-Likelihood
#' @param ind Model indicator vector.
#' @param data Prepared dataset.
#' @param y Target variable name.
#' @param predictors Predictor names.
#' @param sd_cols SD column names.
#' @param cache Evaluation cache.
#' @param base_models Compiled base stan templates.
#' @param is_clean Logical vector of error-free variables.
#' @param iter_stan Integer. Number of iterations for brms models.
#' @param cores_stan Integer. Number of cores for brms models.
#' @return List containing logml value and updated cache.
#' @noRd
evaluate_model_simple <- function(ind, data, y, predictors, sd_cols, cache, base_models, is_clean, iter_stan, cores_stan) {
  key <- ind_to_key(ind)
  cached <- cache_lookup(cache, key)
  if (!is.null(cached)) return(list(logml = cached, cache = cache))

  act_preds <- predictors[ind == 1]
  act_sds <- sd_cols[ind == 1]
  act_clean <- is_clean[ind == 1]

  c_preds <- act_preds[act_clean]
  n_preds <- act_preds[!act_clean]
  n_sds <- act_sds[!act_clean]

  c_count <- length(c_preds)
  n_count <- length(n_preds)

  df_temp <- data.frame(dummy_y = data[[y]])
  names(df_temp)[1] <- y

  if (c_count > 0) {
    for (i in 1:c_count) df_temp[[paste0("SlotC_", i)]] <- data[[c_preds[i]]]
  }
  if (n_count > 0) {
    for (i in 1:n_count) {
      df_temp[[paste0("SlotN_", i)]] <- data[[n_preds[i]]]
      df_temp[[paste0("sd_SlotN_", i)]] <- data[[n_sds[i]]]
    }
  }

  base_key <- paste0("C", c_count, "_N", n_count)

  resultado_seguro <- tryCatch({
    tmp_logml <- -Inf
    warmup_stan <- floor(iter_stan / 3) # Se calcula el warmup dinámicamente
    invisible(capture.output(suppressWarnings(suppressMessages({
      fit <- stats::update(base_models[[base_key]], newdata = df_temp, refresh = 0, silent = 2,
                           iter = iter_stan, warmup = warmup_stan, chains = 2, cores = cores_stan, control = list(adapt_delta = 0.95))

      # pBMA / Stacking: ELPD via LOO Pseudo-Likelihood
      loo_res <- suppressWarnings(brms::loo(fit))
      tmp_logml <- loo_res$estimates["elpd_loo", "Estimate"]
      rm(loo_res)

      rm(fit); gc(verbose = FALSE)
    }))))

    if (is.numeric(tmp_logml) && is.finite(tmp_logml)) {
      cache$key <- c(cache$key, key)
      cache$ind <- c(cache$ind, list(ind))
      cache$logml <- c(cache$logml, tmp_logml)
    }
    list(logml = tmp_logml, cache = cache)
  }, error = function(e) list(logml = -Inf, cache = cache))

  return(resultado_seguro)
}

#' Extract pointwise LOO (For Stacking optimization matrix)
#' @param ind Model indicator vector.
#' @param data Prepared dataset.
#' @param y Target variable name.
#' @param predictors Predictor names.
#' @param sd_cols SD column names.
#' @param base_models Compiled base stan templates.
#' @param is_clean Logical vector of error-free variables.
#' @param iter_stan Integer. Number of iterations for brms models.
#' @param cores_stan Integer. Number of cores for brms models.
#' @return Numeric vector of pointwise LOO values, or NULL on error.
#' @noRd
extraer_pointwise_loo <- function(ind, data, y, predictors, sd_cols, base_models, is_clean, iter_stan, cores_stan) {
  act_preds <- predictors[ind == 1]
  act_sds <- sd_cols[ind == 1]
  act_clean <- is_clean[ind == 1]

  c_preds <- act_preds[act_clean]
  n_preds <- act_preds[!act_clean]
  n_sds <- act_sds[!act_clean]

  c_count <- length(c_preds)
  n_count <- length(n_preds)

  df_temp <- data.frame(dummy_y = data[[y]])
  names(df_temp)[1] <- y

  if (c_count > 0) {
    for (i in 1:c_count) df_temp[[paste0("SlotC_", i)]] <- data[[c_preds[i]]]
  }
  if (n_count > 0) {
    for (i in 1:n_count) {
      df_temp[[paste0("SlotN_", i)]] <- data[[n_preds[i]]]
      df_temp[[paste0("sd_SlotN_", i)]] <- data[[n_sds[i]]]
    }
  }

  base_key <- paste0("C", c_count, "_N", n_count)

  pointwise <- tryCatch({
    res <- NULL
    warmup_stan <- floor(iter_stan / 3)
    invisible(capture.output(suppressWarnings(suppressMessages({
      fit <- stats::update(base_models[[base_key]], newdata = df_temp, refresh = 0, silent = 2,
                           iter = iter_stan, warmup = warmup_stan, chains = 2, cores = cores_stan, control = list(adapt_delta = 0.95))

      loo_res <- suppressWarnings(brms::loo(fit))
      res <- loo_res$pointwise[, "elpd_loo"]

      rm(fit, loo_res); gc(verbose = FALSE)
    }))))
    res
  }, error = function(e) NULL)

  return(pointwise)
}

# ------------------------------------------------------------------------------
# 4. EXPLORATION FUNCTION
# ------------------------------------------------------------------------------

#' Execute Model Space Search (Exhaustive or MC3)
#' @param data Prepared dataset.
#' @param y Target variable name.
#' @param predictors Predictor names.
#' @param sd_cols SD column names.
#' @param path_base_models File path to compiled base models.
#' @param p_max_exhaust Threshold for exhaustive search.
#' @param n_iter Number of iterations for MC3.
#' @param method Ensemble method ('stacking' or 'pbma').
#' @param verbosity Integer for message printing.
#' @param is_clean Logical vector of error-free variables.
#' @param iter_stan Integer. Number of iterations for brms models.
#' @param cores_stan Integer. Number of cores for brms models.
#' @return List containing final cache and target chain.
#' @noRd
run_discovery_phase <- function(data, y, predictors, sd_cols,
                                path_base_models,
                                p_max_exhaust = 9, n_iter = 1500,
                                method = c("stacking", "pbma"), verbosity = 1, is_clean,
                                iter_stan, cores_stan) {
  method <- match.arg(method)
  p <- length(predictors)
  cache <- list(key = character(), ind = list(), logml = numeric())

  if (verbosity > 0) cat(sprintf("--- Starting search phase (p=%d) ---\n", p))

  if (p <= p_max_exhaust) {
    # =======================================================
    # ROUTE A: EXHAUSTIVE SEARCH (p <= max_exhaust)
    # =======================================================

    n_total <- 2^p
    if (verbosity > 0) cat(sprintf("Route: Exhaustive. Evaluating %d models...\n", n_total))

    all_combs <- expand.grid(rep(list(0:1), p))
    lista_inds <- as.list(as.data.frame(t(all_combs)))
    progreso_pasos <- round(seq(1, n_total, length.out = 11))[-1]

    resultados <- list()
    for(i in 1:10) {
      inicio <- if(i==1) 1 else progreso_pasos[i-1] + 1
      fin <- progreso_pasos[i]
      if(inicio > n_total) break

      bloque_actual <- lista_inds[inicio:fin]

      res_bloque <- future.apply::future_lapply(bloque_actual, function(ind) {
        base_models_local <- readRDS(path_base_models)

        eval <- evaluate_model_simple(ind, data, y, predictors, sd_cols,
                                      list(key=character(), ind=list(), logml=numeric()),
                                      base_models_local, is_clean, iter_stan, cores_stan)

        rm(base_models_local)
        gc(verbose = FALSE)

        return(list(key = ind_to_key(ind), ind = ind, logml = eval$logml))
      }, future.seed = TRUE)

      resultados <- c(resultados, res_bloque)
      if (verbosity > 0) cat(sprintf("Exhaustive progress: %d%% (%d/%d models)\n", i*10, fin, n_total))
    }

    for (res in resultados) {
      if (is.finite(res$logml)) {
        cache$key <- c(cache$key, res$key); cache$ind <- c(cache$ind, list(res$ind)); cache$logml <- c(cache$logml, res$logml)
      }
    }
    return(list(cache = cache, type = "exhaustive"))

  } else {
    # =======================================================
    # ROUTE B: MCMC SEARCH (p > max_exhaust)
    # =======================================================

    if (method == "stacking") {
      if (verbosity > 0) cat(sprintf("Route: Heuristic MC3 (Tabu+Greedy) for Bayesian Stacking. Iterations: %d\n", n_iter))
      temperatures <- c(0.01, 2, 5, 10)
    } else {
      if (verbosity > 0) cat(sprintf("Route: Classic MC3 (No Tabu) for Pseudo-BMA. Iterations: %d\n", n_iter))
      temperatures <- c(1, 1.5, 2, 4)
    }

    K <- length(temperatures)

    init_res <- future.apply::future_lapply(1:K, function(k) {
      base_models_local <- readRDS(path_base_models)

      ind_init <- rbinom(p, 1, 0.5)
      eval <- evaluate_model_simple(ind_init, data, y, predictors, sd_cols,
                                    list(key=character(), ind=list(), logml=numeric()),
                                    base_models_local, is_clean, iter_stan, cores_stan)

      rm(base_models_local); gc(verbose = FALSE)
      return(list(ind = ind_init, logml = eval$logml, key = ind_to_key(ind_init)))
    }, future.seed = TRUE)

    inds_current <- lapply(init_res, `[[`, "ind")
    logmls_current <- sapply(init_res, `[[`, "logml")
    for(r in init_res) {
      cache$key <- c(cache$key, r$key); cache$ind <- c(cache$ind, list(r$ind)); cache$logml <- c(cache$logml, r$logml)
    }

    chain_target <- matrix(0, nrow=n_iter, ncol=p)
    progreso_pasos <- round(seq(1, n_iter, length.out = 11))[-1]

    for(i in 1:n_iter) {
      if (method == "stacking") {
        propu_modelos <- lapply(inds_current, propose_model_tabu, cache = cache)
      } else {
        propu_modelos <- lapply(inds_current, propose_model)
      }

      res_step <- future.apply::future_lapply(1:K, function(k) {
        base_models_local <- readRDS(path_base_models)

        eval <- evaluate_model_simple(propu_modelos[[k]], data, y, predictors, sd_cols,
                                      list(key=character(), ind=list(), logml=numeric()),
                                      base_models_local, is_clean, iter_stan, cores_stan)

        rm(base_models_local); gc(verbose = FALSE)
        return(list(logml = eval$logml, ind = propu_modelos[[k]], key = ind_to_key(propu_modelos[[k]])))
      }, future.seed = TRUE)

      # Acceptance check
      for(k in 1:K) {
        if(is.na(match(res_step[[k]]$key, cache$key))) {
          cache$key <- c(cache$key, res_step[[k]]$key)
          cache$ind <- c(cache$ind, list(res_step[[k]]$ind))
          cache$logml <- c(cache$logml, res_step[[k]]$logml)
        }

        log_prior_prop <- log_prior_SB(res_step[[k]]$ind, p)
        log_prior_curr <- log_prior_SB(inds_current[[k]], p)

        lp_prop <- res_step[[k]]$logml + log_prior_prop
        lp_curr <- logmls_current[k] + log_prior_curr

        if (is.infinite(lp_prop) && lp_prop < 0 && is.infinite(lp_curr) && lp_curr < 0) {
          log_alpha <- 0
        } else {
          log_alpha <- (lp_prop - lp_curr) / temperatures[k]
        }

        if(log(runif(1)) < min(0, log_alpha)) {
          inds_current[[k]] <- res_step[[k]]$ind
          logmls_current[k] <- res_step[[k]]$logml
        }
      }

      # Thermodynamic coupling swaps
      if (K > 1) {
        k1 <- sample(1:(K-1), 1)
        k2 <- k1 + 1

        lp1 <- logmls_current[k1] + log_prior_SB(inds_current[[k1]], p)
        lp2 <- logmls_current[k2] + log_prior_SB(inds_current[[k2]], p)

        if (is.infinite(lp1) && lp1 < 0 && is.infinite(lp2) && lp2 < 0) {
          log_alpha_swap <- 0
        } else {
          log_alpha_swap <- (lp2 - lp1) * ((1 / temperatures[k1]) - (1 / temperatures[k2]))
        }

        if(log(runif(1)) < min(0, log_alpha_swap)) {
          temp_ind <- inds_current[[k1]]
          inds_current[[k1]] <- inds_current[[k2]]
          inds_current[[k2]] <- temp_ind

          temp_log <- logmls_current[k1]
          logmls_current[k1] <- logmls_current[k2]
          logmls_current[k2] <- temp_log
        }
      }

      chain_target[i, ] <- inds_current[[1]]

      if (verbosity > 0 && i %in% progreso_pasos) {
        idx_prog <- which(progreso_pasos == i)
        cat(sprintf("Progress MC3: %d%% (%d/%d iterations) | Unique models found: %d\n",
                    idx_prog * 10, i, n_iter, length(cache$key)))
      }
    }

    tipo_retorno <- ifelse(method == "stacking", "mc3_heuristic", "mc3_classic")
    return(list(cache = cache, type = tipo_retorno, chain_target = chain_target))
  }
}
