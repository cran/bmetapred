#' Predict Method for Ensemble Models
#'
#' @description
#' Generates posterior predictive samples and summary statistics for new data
#' based on a fitted ensemble model. It samples from
#' the underlying `brms` base models proportionally to their optimized weights.
#'
#' @param object An object of class \code{bmetapred_model} returned by \code{fit_bmetapred}.
#' @param newdata A data frame containing the new data to predict. Must include
#'   the predictor variables and their corresponding standard deviation columns
#'   as defined in the original training phase.
#' @param n_samples_totales Integer. The total number of posterior predictive
#'   samples to generate across the ensemble. Default is 2000.
#' @param estimator Character. For \code{"pbma"} models, dictates the prediction
#'   strategy. Use \code{"ensemble"} for weighted averaging or \code{"HPM"}
#'   (Highest Probability Model) to predict using only the top model.
#'   Ignored if the method is \code{"stacking"}.
#' @param occam_window Numeric. A value between 0 and 1 defining the Occam's window
#'   for predictive mass. Models falling outside this cumulative probability
#'   are discarded for predicting, and remaining weights are renormalized. Default is \code{NULL}.
#' @param iter_stan Integer. Number of MCMC iterations for underlying brms prediction.
#'   Defaults to the value used during training.
#' @param cores_stan Integer. Number of cores for underlying brms prediction.
#'   Defaults to the value used during training.
#' @param verbosity Integer. Level of printed messages. \code{0} suppresses output,
#'   \code{1} prints progress. Default is 1.
#' @param ... Additional arguments passed to methods.
#'
#' @return A data frame containing the predicted \code{estimate} (mean),
#'   \code{lower_95}, and \code{upper_95} credible intervals. If present in
#'   \code{newdata}, the \code{truth} and \code{datetime} columns are also appended.
#'   The raw predictive matrix is attached as an attribute named \code{"raw_samples"}.
#'
#' @importFrom stats update quantile
#' @method predict bmetapred_model
#' @export
#'
#' @examples
#' \donttest{
#' set.seed(123)
#' df_train <- data.frame(
#'   Y = rnorm(30),
#'   M1 = rnorm(30), M2 = rnorm(30),
#'   sd_M1 = rep(0.1, 30), sd_M2 = rep(0.2, 30)
#' )
#' df_test <- data.frame(
#'   M1 = rnorm(5), M2 = rnorm(5),
#'   sd_M1 = rep(0.1, 5), sd_M2 = rep(0.2, 5)
#' )
#'
#' my_ensemble <- fit_bmetapred(
#'   data_train = df_train,
#'   y = "Y",
#'   predictors = c("M1", "M2"),
#'   sd_cols = c("sd_M1", "sd_M2"),
#'   method = "stacking",
#'   parallel_workers = 1,
#'   cores_stan = 1,
#'   iter_stan = 300,
#'   verbosity = 0
#' )
#'
#' predictions <- predict(my_ensemble, newdata = df_test, n_samples_totales = 200, verbosity = 0)
#' head(predictions)
#' }

predict.bmetapred_model <- function(object, newdata, n_samples_totales = 2000,
                                    estimator = "ensemble", occam_window = NULL,
                                    iter_stan = NULL, cores_stan = NULL,
                                    verbosity = 1, ...) {

  if (verbosity > 0) cat(sprintf("--- Generating Predictions (%s) ---\n", toupper(object$method)))

  # Retrieve parameters from the trained object if not explicitly provided
  if (is.null(iter_stan)) iter_stan <- object$params$iter_stan
  if (is.null(cores_stan)) cores_stan <- object$params$cores_stan

  pesos <- object$pesos
  indices_finales <- object$top_indices
  base_models <- object$base_models
  y <- object$params$y
  predictors <- object$params$predictors
  sd_cols <- object$params$sd_cols

  is_clean <- object$params$is_clean
  data_train <- object$data_train
  t_micro <- object$params$threshold_micro

  if (object$method %in% c("pbma")) {
    if (estimator == "HPM") {
      idx_max <- which.max(pesos)
      pesos[] <- 0
      pesos[idx_max] <- 1
    } else if (!is.null(occam_window) && occam_window < 1) {
      masa_acumulada <- cumsum(pesos)
      top_occam <- which(masa_acumulada >= occam_window)[1]
      if (is.na(top_occam)) top_occam <- length(pesos)
      if (top_occam < length(pesos)) pesos[(top_occam + 1):length(pesos)] <- 0
      pesos <- pesos / sum(pesos)
    }
  }

  muestras_por_modelo <- round(as.numeric(pesos) * n_samples_totales)
  diff_samples <- n_samples_totales - sum(muestras_por_modelo)
  if (diff_samples != 0) {
    idx_max <- which.max(muestras_por_modelo)
    muestras_por_modelo[idx_max] <- muestras_por_modelo[idx_max] + diff_samples
  }

  lista_predicciones <- list()

  for (i in 1:length(indices_finales)) {
    if (muestras_por_modelo[i] > 0) {

      ind_w <- indices_finales[[i]]

      act_preds <- predictors[ind_w == 1]
      act_sds <- sd_cols[ind_w == 1]
      act_clean <- is_clean[ind_w == 1]

      c_preds <- act_preds[act_clean]
      n_preds <- act_preds[!act_clean]
      n_sds <- act_sds[!act_clean]

      c_count <- length(c_preds)
      n_count <- length(n_preds)

      df_temp_train <- data.frame(dummy_y = data_train[[y]])
      names(df_temp_train)[1] <- y

      df_temp_test <- data.frame(dummy_y = rep(NA, nrow(newdata)))
      if(y %in% names(newdata)) df_temp_test$dummy_y <- newdata[[y]]
      names(df_temp_test)[1] <- y

      if (c_count > 0) {
        for (j in 1:c_count) {
          df_temp_train[[paste0("SlotC_", j)]] <- data_train[[c_preds[j]]]
          df_temp_test[[paste0("SlotC_", j)]] <- newdata[[c_preds[j]]]
        }
      }

      if (n_count > 0) {
        for (j in 1:n_count) {
          df_temp_train[[paste0("SlotN_", j)]] <- data_train[[n_preds[j]]]
          df_temp_train[[paste0("sd_SlotN_", j)]] <- data_train[[n_sds[j]]]

          df_temp_test[[paste0("SlotN_", j)]] <- newdata[[n_preds[j]]]

          sd_test <- newdata[[n_sds[j]]]
          sd_test[sd_test < t_micro] <- t_micro
          df_temp_test[[paste0("sd_SlotN_", j)]] <- sd_test
        }
      }

      base_key <- paste0("C", c_count, "_N", n_count)

      invisible(capture.output(suppressWarnings(suppressMessages({
        warmup_stan <- floor(iter_stan / 2)
        fit_final <- stats::update(base_models[[base_key]], newdata = df_temp_train,
                                   refresh = 0, silent = 2, iter = iter_stan, warmup = warmup_stan, chains = 2, cores = cores_stan)

        preds_completas <- brms::posterior_predict(fit_final, newdata = df_temp_test)

        rm(fit_final); gc(verbose=FALSE)
      }))))

      idx_elegidos <- sample(1:nrow(preds_completas), muestras_por_modelo[i], replace = TRUE)
      lista_predicciones[[i]] <- preds_completas[idx_elegidos, , drop = FALSE]
    }
  }

  pred_final_samples <- do.call(rbind, lista_predicciones)
  res_final <- data.frame(
    estimate = colMeans(pred_final_samples),
    lower_95 = apply(pred_final_samples, 2, quantile, probs = 0.025),
    upper_95 = apply(pred_final_samples, 2, quantile, probs = 0.975)
  )

  if(y %in% names(newdata)) res_final$truth <- newdata[[y]]
  if("datetime" %in% names(newdata)) res_final$datetime <- newdata$datetime

  attr(res_final, "raw_samples") <- pred_final_samples

  return(res_final)
}

#' Extract Posterior Inclusion Probabilities (PIP)
#'
#' @description
#' Calculates the Posterior Inclusion Probability (PIP) for each predictor (model)
#' variable evaluated during the ensemble fitting process.
#'
#' @param object An object of class \code{bmetapred_model} returned by \code{fit_bmetapred}.
#'
#' @return A named numeric vector sorted in decreasing order, representing the
#'   cumulative probability (or weight) of inclusion for each base model.
#'
#' @export
#'
#' @examples
#' # Create a small bmetapred_model object for fast testing
#' mock_ensemble <- structure(
#'   list(
#'     pesos = c(0.6, 0.3, 0.1),
#'     top_indices = list(c(1, 1, 0), c(1, 0, 1), c(0, 1, 0)),
#'     params = list(predictors = c("M1", "M2", "M3")),
#'     method = "stacking"
#'   ),
#'   class = "bmetapred_model"
#' )
#' pips <- get_pip_bmetapred(mock_ensemble)

get_pip_bmetapred <- function(object) {
  pesos <- as.numeric(object$pesos)
  indices <- object$top_indices
  n_modelos <- length(indices)
  p <- length(indices[[1]])
  nombres_vars <- object$params$predictors

  pip <- numeric(p)
  names(pip) <- nombres_vars

  for (j in 1:p) {
    for (m in 1:n_modelos) {
      if (indices[[m]][j] == 1) {
        pip[j] <- pip[j] + pesos[m]
      }
    }
  }
  return(sort(pip, decreasing = TRUE))
}

#' Summarize a bmetapred Ensemble Model
#'
#' @description
#' Provides a tabular summary of the meta-models included in the final ensemble,
#' displaying their structural formula, rank, and assigned weight.
#'
#' @param object An object of class \code{bmetapred_model} returned by \code{fit_bmetapred}.
#' @param ... Additional arguments passed to methods.
#'
#' @return A data frame detailing the rank, weight, and structure of the
#'   meta-models that constitute the ensemble. The result is also printed to the console.
#'
#' @method summary bmetapred_model
#' @export
#'
#' @examples
#' # Create a lightweight mock bmetapred_model object for fast testing
#' mock_ensemble <- structure(
#'   list(
#'     pesos = c(0.6, 0.3, 0.1),
#'     top_indices = list(c(1, 1, 0), c(1, 0, 1), c(0, 1, 0)),
#'     params = list(predictors = c("M1", "M2", "M3")),
#'     method = "stacking"
#'   ),
#'   class = "bmetapred_model"
#' )
#' summary_table <- summary(mock_ensemble)

summary.bmetapred_model <- function(object, ...) {
  pesos <- as.numeric(object$pesos)
  indices <- object$top_indices
  nombres_vars <- object$params$predictors

  estructuras <- sapply(indices, function(ind) {
    vars_activas <- nombres_vars[ind == 1]
    if (length(vars_activas) == 0) return("Null Model (Intercept-only)")
    return(paste(vars_activas, collapse = " + "))
  })

  resumen <- data.frame(
    Rank = 1:length(pesos),
    Weight = round(pesos, 4),
    Structure = estructuras,
    stringsAsFactors = FALSE
  )

  resumen <- resumen[order(resumen$Weight, decreasing = TRUE), ]

  cat(sprintf("--- Ensemble summary (%s) ---\n", toupper(object$method)))
  cat(sprintf("Total number of models used: %d\n\n", nrow(resumen)))

  return(resumen)
}

#' Fit LOESS Error Models for Ensemble Predictors
#'
#' @description
#' Trains a Local Polynomial Regression (LOESS) model for each predictor to estimate
#' its prediction error (standard deviation) based on historical validation data.
#' This is particularly useful when base predictors do not output native uncertainty
#' estimates. The LOESS model predicts the squared error as a function of the
#' predicted value.
#'
#' @param data A data frame containing the historical or validation data.
#' @param y A character string specifying the name of the true target variable.
#' @param predictors A character vector specifying the names of the predictor variables.
#' @param span Numeric. The smoothing parameter for the LOESS model. Default is 0.4.
#'
#' @return An object of class \code{bmetapred_loess} containing the fitted LOESS models
#'   and the extrapolation boundaries for each predictor.
#'
#' @importFrom stats loess
#' @export
#'
#' @examples
#' set.seed(123)
#' df_val <- data.frame(Y = rnorm(50), M1 = rnorm(50), M2 = rnorm(50))
#' loess_estimators <- fit_LOESS(df_val, y = "Y", predictors = c("M1", "M2"))

fit_LOESS <- function(data, y, predictors, span = 0.4) {
  if (!y %in% names(data)) stop("Target variable 'y' not found in data.")

  loess_models <- list()

  for (pred in predictors) {
    if (!pred %in% names(data)) stop(sprintf("Predictor '%s' not found in data.", pred))

    Y_vec <- data[[y]]
    M_pred <- data[[pred]]

    # Calculate empirical variance (squared error)
    err_sq <- (Y_vec - M_pred)^2
    df_err <- data.frame(err_sq = err_sq, M_pred = M_pred)

    # Fit LOESS
    fit_l <- suppressWarnings(stats::loess(err_sq ~ M_pred, data = df_err, span = span))

    # Store the model and its training boundaries for robust extrapolation
    loess_models[[pred]] <- list(
      model = fit_l,
      min_M = min(M_pred, na.rm = TRUE),
      max_M = max(M_pred, na.rm = TRUE)
    )
  }

  obj <- list(
    models = loess_models,
    predictors = predictors,
    y = y
  )
  class(obj) <- "bmetapred_loess"

  return(obj)
}

#' Predict Standard Deviations using Fitted LOESS Error Models
#'
#' @description
#' Applies fitted LOESS error models to new data to estimate the standard
#' deviation (prediction error) of each predictor. It includes robust handling
#' of extrapolations by capping predictions at the boundary values of the
#' original training data to prevent \code{NA}s or collapsed variances.
#'
#' @param object An object of class \code{bmetapred_loess} returned by \code{fit_LOESS}.
#' @param newdata A data frame containing the new predictions to evaluate.
#' @param prefix Character. The prefix to add to the predictor names to create
#'   the new standard deviation columns. Default is \code{"sd_"}.
#' @param min_var Numeric. The minimum allowed variance to prevent zero or
#'   negative variances (singularities). Default is 0.0004 (SD = 0.02).
#' @param ... Additional arguments passed to methods.
#'
#' @return A data frame identical to \code{newdata} but with appended columns
#'   containing the estimated standard deviations for each predictor.
#'
#' @importFrom stats predict
#' @method predict bmetapred_loess
#' @export
#'
#' @examples
#' set.seed(123)
#' df_val <- data.frame(Y = rnorm(50), M1 = rnorm(50), M2 = rnorm(50))
#' loess_estimators <- fit_LOESS(df_val, y = "Y", predictors = c("M1", "M2"))
#'
#' df_new <- data.frame(M1 = rnorm(10), M2 = rnorm(10))
#' df_ready <- predict(loess_estimators, newdata = df_new)

predict.bmetapred_loess <- function(object, newdata, prefix = "sd_", min_var = 0.0004, ...) {
  out_data <- newdata
  predictors <- object$predictors

  for (pred in predictors) {
    if (!pred %in% names(newdata)) stop(sprintf("Predictor '%s' not found in newdata.", pred))

    nuevos_M <- newdata[[pred]]
    modelo_loess <- object$models[[pred]]$model
    min_M <- object$models[[pred]]$min_M
    max_M <- object$models[[pred]]$max_M

    df_new <- data.frame(M_pred = nuevos_M)
    preds_var <- suppressWarnings(stats::predict(modelo_loess, newdata = df_new))

    # 1. Handle extrapolations (NAs) by patching with the nearest known boundary
    if (any(is.na(preds_var))) {
      pred_min <- suppressWarnings(stats::predict(modelo_loess, newdata = data.frame(M_pred = min_M)))
      pred_max <- suppressWarnings(stats::predict(modelo_loess, newdata = data.frame(M_pred = max_M)))

      preds_var[is.na(preds_var) & nuevos_M < min_M] <- pred_min
      preds_var[is.na(preds_var) & nuevos_M > max_M] <- pred_max
    }

    # 2. Ensure variance doesn't collapse to 0 or negative values
    preds_var_safe <- pmax(preds_var, min_var)

    # 3. Convert Variance to Standard Deviation and append
    sd_name <- paste0(prefix, pred)
    out_data[[sd_name]] <- sqrt(preds_var_safe)
  }

  return(out_data)
}
