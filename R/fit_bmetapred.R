#' Fit an Ensemble of bayesian linear meta-models with prediction errors
#'
#' @description
#' Fits a Pseudo-Bayesian model averaging (pBMA) or Stacking ensemble of models.
#' It accounts for prediction errors in the individual models.
#' Uses either an exhaustive search or a Markov Chain Monte Carlo (MC3) algorithm
#' to explore the meta-model space and calculate pBMA/Stacking weights.
#'
#' @param data_train A data frame containing the training data.
#' @param y A character string specifying the name of the target variable.
#' @param predictors A character vector specifying the names of the predictor variables (models).
#' @param sd_cols A character vector specifying the names of the columns containing
#'   the standard deviations (prediction errors) for each predictor (model). Must be in the
#'   same order as `predictors`.
#' @param method A character string specifying the ensemble method. Can be either
#'   \code{"stacking"} (default) or \code{"pbma"}.
#' @param p_max_exhaust Integer. The maximum number of predictors to perform an
#'   exhaustive search. If the number of predictors exceeds this value, the MC3
#'   algorithm is used. Default is 9.
#' @param n_iter_mc3 Integer. Number of iterations for the MC3 algorithm. Default is 1500.
#' @param burnin_pbma Integer. Number of burn-in iterations to discard when using
#'   the \code{"pbma"} method. If \code{NULL}, it defaults to 20% of \code{n_iter_mc3}.
#' @param min_K Integer. Minimum number of top models to evaluate before early stopping
#'   in the Stacking matrix compilation. Default is 10.
#' @param max_K Integer. Maximum number of top models to evaluate in the Stacking matrix.
#'   Default is 50.
#' @param pri_sd Numeric. Prior standard deviation for the regression coefficients.
#'   Default is 2.5.
#' @param pri_sd_int Numeric. Prior standard deviation for the intercepts.
#'   Default is 5.
#' @param pri_sd_int_sigma Numeric. Prior standard deviation for the sigma intercept. Default is 2.
#' @param pri_sd_sigma Numeric. Prior standard deviation for the sigma coefficients.
#'   Default is 0.5.
#' @param parallel_workers Integer. Number of workers for parallel evaluation of the model space.
#'   Default is 1 (sequential execution).
#' @param cores_stan Integer. Number of cores to use for each underlying brms model. Default is 1.
#' @param iter_stan Integer. Number of MCMC iterations for the underlying brms models. Default is 1500.
#' @param verbosity Integer. Level of printed messages. \code{0} suppresses output,
#'   \code{1} prints progress. Default is 1.
#' @param backend Character string specifying the package to use as the backend for \code{brms}. 
#'   Default is \code{"rstan"}. However, for intensive computational tasks and large model spaces, 
#'   we highly recommend installing \code{"cmdstanr"} (available from the Stan-dev repository) 
#'   and setting \code{backend = "cmdstanr"} to improve memory management and speed.
#'
#' @return An object of class \code{bmetapred_model} containing the optimized ensemble weights,
#'   the list of selected model indices, the compiled \code{brms} base templates,
#'   and internal parameters needed for prediction.
#'
#' @importFrom stats as.formula rnorm runif update var rbinom gaussian
#' @importFrom utils capture.output
#' @export
#'
#' @examples
#' \donttest{
#' set.seed(123)
#' df <- data.frame(
#'   Y = rnorm(30),
#'   M1 = rnorm(30), M2 = rnorm(30),
#'   sd_M1 = rep(0.1, 30), sd_M2 = rep(0.2, 30)
#' )
#'
#' my_ensemble <- fit_bmetapred(
#'   data_train = df,
#'   y = "Y",
#'   predictors = c("M1", "M2"),
#'   sd_cols = c("sd_M1", "sd_M2"),
#'   method = "stacking",
#'   parallel_workers = 1,
#'   cores_stan = 1,
#'   iter_stan = 300,
#'   verbosity = 0
#' )
#' summary(my_ensemble)
#' }

fit_bmetapred <- function(data_train, y, predictors, sd_cols,
                          method = c("stacking", "pbma"),
                          p_max_exhaust = 9, n_iter_mc3 = 1500, burnin_pbma = NULL,
                          min_K = 10, max_K = 50, pri_sd = 2.5, pri_sd_int = 5, pri_sd_sigma = 0.5, pri_sd_int_sigma = 5,
                          parallel_workers = 1, cores_stan = 1, iter_stan = 1500,
                          verbosity = 1, backend = "rstan") {
  
  method <- match.arg(method)
  p <- length(predictors)
  
  # Configure parallel plan politely (revert to user's plan on exit)
  if (parallel_workers > 1) {
    old_plan <- future::plan()
    future::plan(future::multisession, workers = parallel_workers)
    on.exit(future::plan(old_plan), add = TRUE)
  }
  
  # --- PHASE 0.1: Error cleaning and categorization ---
  auditoria <- audit_and_clean_errors(data_train, predictors, sd_cols, verbosity = verbosity)
  data_train <- auditoria$data
  is_clean <- auditoria$is_clean
  
  max_clean <- sum(is_clean)
  max_noisy <- sum(!is_clean)
  n_modelos_optimizados <- (max_clean + 1) * (max_noisy + 1)
  
  df_dummy <- data.frame(dummy_y = data_train[[y]])
  names(df_dummy)[1] <- y
  if (max_clean > 0) for(c in 1:max_clean) df_dummy[[paste0("SlotC_", c)]] <- rnorm(nrow(df_dummy))
  if (max_noisy > 0) {
    for(n in 1:max_noisy) {
      df_dummy[[paste0("SlotN_", n)]] <- rnorm(nrow(df_dummy))
      df_dummy[[paste0("sd_SlotN_", n)]] <- rep(0.1, nrow(df_dummy))
    }
  }
  
  grid_combs <- expand.grid(C = 0:max_clean, N = 0:max_noisy)
  
  backend_stan <- backend
  pars_stan <- brms::save_pars(all = FALSE)
  
  # Función helper para compilar el modelo base (evita repetir código)
  compile_fn <- function(i) {
    c <- grid_combs$C[i]
    n <- grid_combs$N[i]
    key <- paste0("C", c, "_N", n)
    
    if (c == 0 && n == 0) {
      form <- brms::bf(as.formula(paste(y, "~ 1")))
      priors_k <- c(brms::set_prior(paste0("normal(0, ", pri_sd_int, ")"), class = "Intercept"), brms::set_prior("exponential(1)", class = "sigma"))
    } else {
      terms_C <- if(c > 0) paste0("SlotC_", 1:c) else NULL
      terms_N <- if(n > 0) paste0("SlotN_", 1:n) else NULL
      form_mean <- as.formula(paste(y, "~", paste(c(terms_C, terms_N), collapse=" + ")))
      
      if (n > 0) {
        terms_sigma <- paste0("sd_SlotN_", 1:n)
        form_sigma <- as.formula(paste("sigma ~", paste(paste0("log(", terms_sigma, ")"), collapse=" + ")))
        form <- brms::bf(form_mean, form_sigma)
        
        priors_k <- c(brms::set_prior(paste0("normal(0, ", pri_sd, ")"), class = "b"),
                      brms::set_prior(paste0("normal(0, ", pri_sd_int, ")"), class = "Intercept"),
                      brms::set_prior(paste0("normal(0, ", pri_sd_int_sigma, ")"), class = "Intercept", dpar = "sigma"),
                      brms::set_prior(paste0("normal(0, ", pri_sd_sigma, ")"), class = "b", dpar = "sigma"))
      } else {
        form <- brms::bf(form_mean)
        priors_k <- c(brms::set_prior(paste0("normal(0, ", pri_sd, ")"), class = "b"),
                      brms::set_prior(paste0("normal(0, ", pri_sd_int, ")"), class = "Intercept"),
                      brms::set_prior("exponential(1)", class = "sigma"))
      }
    }
    
    fit <- suppressWarnings(suppressMessages(brms::brm(formula = form, data = df_dummy, prior = priors_k, family = gaussian(),
                                                       chains = 1, iter = 10, warmup=5,
                                                       backend = backend_stan,
                                                       cores = 1, # Base comp doesn't need internal multi-core
                                                       refresh = 0, silent = 2,
                                                       save_pars = pars_stan)))
    return(list(key = key, model = fit))
  }
  
  if (backend_stan == "cmdstanr" && parallel_workers > 1) {
    if (verbosity > 0) cat(sprintf("--- Compiling base models in parallel (%d models) ---\n", n_modelos_optimizados))
    base_models <- future.apply::future_lapply(1:nrow(grid_combs), compile_fn, future.seed = TRUE)
  } else {
    if (verbosity > 0) cat(sprintf("--- Compiling base models sequentially (%d models) ---\n", n_modelos_optimizados))
    base_models <- lapply(1:nrow(grid_combs), compile_fn)
  }
  
  base_models_named <- list()
  for (item in base_models) base_models_named[[item$key]] <- item$model
  
  path_base_models <- tempfile(pattern = "plantillas_stan_", fileext = ".rds")
  saveRDS(base_models_named, path_base_models)
  
  rm(base_models, base_models_named)
  gc(verbose = FALSE, full = TRUE)
  
  # --- PHASE 1: SEARCH ---
  search_res <- run_discovery_phase(data_train, y, predictors, sd_cols, path_base_models,
                                    p_max_exhaust, n_iter_mc3, method, verbosity, is_clean,
                                    iter_stan, cores_stan)

  base_models_named <- readRDS(path_base_models)

  # ============================================================================
  # PHASE 2 & 3: WEIGHT CALCULATION
  # ============================================================================
  if (method == "stacking") {
    if (verbosity > 0) cat("--- Phase 2: Building Stacking Matrix ---\n")
    df_cache <- data.frame(key = search_res$cache$key, ELPD = search_res$cache$logml, stringsAsFactors = FALSE)
    df_cache$ind <- search_res$cache$ind
    df_cache <- df_cache[order(df_cache$ELPD, decreasing = TRUE), ]

    N_rows <- nrow(data_train)
    matriz_lpd <- list()
    indices_finales <- list()
    k <- 1
    parar_bucle <- FALSE

    while(!parar_bucle && k <= nrow(df_cache)) {
      if (verbosity > 0) cat(sprintf("   Extracting pointwise LOO for model with Rank %d... ", k))
      ind_k <- df_cache$ind[[k]]

      lpd_k <- extraer_pointwise_loo(ind_k, data_train, y, predictors, sd_cols, base_models_named, is_clean, iter_stan, cores_stan)

      if (is.null(lpd_k)) {
        if (verbosity > 0) cat("Error. Skipping model.\n")
        k <- k + 1; next
      }

      matriz_lpd[[length(matriz_lpd) + 1]] <- lpd_k
      indices_finales[[length(indices_finales) + 1]] <- ind_k
      indice_matriz <- length(matriz_lpd)
      if (verbosity > 0) cat("OK.\n")

      if (indice_matriz > 1) {
        diffs <- matriz_lpd[[1]] - matriz_lpd[[indice_matriz]]
        elpd_diff <- sum(diffs)
        se_diff <- sqrt(N_rows * var(diffs))
        if (elpd_diff > (1.96 * se_diff) && indice_matriz >= min_K) {
          if (verbosity > 0) cat(sprintf(">> Early stop: Difference %.2f > 1.96*SE (%.2f).\n", elpd_diff, 1.96 * se_diff))
          parar_bucle <- TRUE
        }
      }
      if (indice_matriz >= max_K) {
        if (verbosity > 0) cat(">>  max_K Limit reached. Stopping.\n")
        parar_bucle <- TRUE
      }
      k <- k + 1
    }

    if (verbosity > 0) cat("--- Phase 3: Optimizing stacking weights ---\n")
    matriz_stacking <- do.call(cbind, matriz_lpd)
    pesos_finales <- loo::stacking_weights(matriz_stacking)

  } else if (method %in% c("pbma")) {
    if (verbosity > 0) cat("--- Phase 2/3: Calculating predictive p-BMA weights ---\n")

    if (search_res$type == "exhaustive") {
      df_cache <- data.frame(key = search_res$cache$key, ELPD = search_res$cache$logml, stringsAsFactors = FALSE)
      df_cache$ind <- search_res$cache$ind
      log_priors <- sapply(df_cache$ind, log_prior_SB, p = p)
      log_unnorm_post <- df_cache$ELPD + log_priors
      max_val <- max(log_unnorm_post)
      weights_unnorm <- exp(log_unnorm_post - max_val)
      pesos_totales <- weights_unnorm / sum(weights_unnorm)

      keep_idx <- which(pesos_totales > 1e-5)
      if(length(keep_idx) == 0) keep_idx <- order(pesos_totales, decreasing = TRUE)[1:min_K]

      pesos_filtrados <- pesos_totales[keep_idx]
      indices_filtrados <- df_cache$ind[keep_idx]
      ord <- order(pesos_filtrados, decreasing = TRUE)
      pesos_finales <- pesos_filtrados[ord]
      pesos_finales <- pesos_finales / sum(pesos_finales)
      indices_finales <- indices_filtrados[ord]
    } else {
      if (is.null(burnin_pbma)) burnin_pbma <- floor(n_iter_mc3 * 0.2)
      chain <- search_res$chain_target[(burnin_pbma + 1):n_iter_mc3, , drop = FALSE]
      chain_keys <- apply(chain, 1, ind_to_key)
      freq_table <- table(chain_keys)
      pesos_totales <- as.numeric(freq_table) / sum(freq_table)
      keys_unique <- names(freq_table)

      ord <- order(pesos_totales, decreasing = TRUE)
      pesos_finales <- pesos_totales[ord]
      keys_finales <- keys_unique[ord]
      indices_finales <- lapply(keys_finales, function(k) as.numeric(strsplit(k, "")[[1]]))
    }
  }

  # --- Constructing the result object ---
  ensamble_obj <- list(
    pesos = pesos_finales,
    top_indices = indices_finales,
    cache_completa = search_res$cache,
    base_models = base_models_named,
    params = list(y = y, predictors = predictors, sd_cols = sd_cols, is_clean = is_clean, threshold_micro = 0.02,
                  iter_stan = iter_stan, cores_stan = cores_stan),
    data_train = data_train,
    method = method
  )
  class(ensamble_obj) <- "bmetapred_model"

  if(file.exists(path_base_models)) file.remove(path_base_models)

  return(ensamble_obj)
}
