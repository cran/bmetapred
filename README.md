
# bmetapred: Error-Aware Bayesian Meta-Predictions and Ensembles

The bmetapred package provides a robust framework for Bayesian model
ensembling (Pseudo-BMA and Stacking) in M-Open settings, explicitly
incorporating the prediction errors of the base models.

When base models do not report native uncertainties, bmetapred includes
tools to estimate heteroscedastic errors from historical validation data
using Local Polynomial Regression (LOESS). The meta-model space is
efficiently explored using exhaustive search or a Markov Chain Monte
Carlo Model Composition (MC3) algorithm powered by brms.

## Installation

You can install the released version of bmetapred from CRAN with:

``` r
install.packages("bmetapred")
```

Alternatively, you can install the development version from GitHub with:

``` r
# install.packages("devtools")
devtools::install_github("JorgePiquerasMarques/bmetapred")
```

## Basic Workflow

This is a basic example showing the standard pipeline: from error
calibration to ensemble fitting and final prediction.

``` r
library(bmetapred)

# 1. Error Calibration (LOESS)
# Estimate prediction errors from historical data if base models lack native uncertainties
loess_estimators <- fit_LOESS(
  data = historical_data, 
  y = "Y", 
  predictors = c("M1", "M2", "M3")
)

# Inject the estimated standard deviations into the new dataset
data_with_sd <- predict(loess_estimators, newdata = current_data)

# 2. Ensemble Fitting
# Fit a Bayesian Stacking ensemble accounting for the estimated errors
my_ensemble <- fit_bmetapred(
  data_train = data_with_sd,
  y = "Y",
  predictors = c("M1", "M2", "M3"),
  sd_cols = c("sd_M1", "sd_M2", "sd_M3"),
  method = "stacking",
  parallel_workers = 2, # Explore model space in parallel
  cores_stan = 1,       # Cores for individual brms models
  iter_stan = 1500,     # MCMC iterations
  verbosity = 1
)

# 3. Model Inspection and Prediction
summary(my_ensemble)

# Generate robust predictive posterior samples
final_predictions <- predict(my_ensemble, newdata = test_data)
head(final_predictions)
```

## Computational Considerations

Adjusting Bayesian meta-models that account for heteroscedastic noise
via brms/Stan is computationally intensive. bmetapred implements a
dynamic caching system and Thermodynamic Integration (MC3 with heated
chains) to make the exploration of large model spaces computationally
viable, avoiding redundant compilations.

### Recommended Backend for brms

By default, `bmetapred` uses the standard `rstan` backend. However,
because exploring large meta-model spaces requires iterative C++
compilation and MCMC sampling, **we highly recommend using `cmdstanr`**
for production and heavy workloads.

You can install it via the Stan-dev repository and activate it in the
fit function:

``` r
# install.packages("cmdstanr", repos = c("[https://mc-stan.org/r-packages/](https://mc-stan.org/r-packages/)", getOption("repos")))
# cmdstanr::install_cmdstan()

my_ensemble <- fit_bmetapred(..., backend = "cmdstanr")
```

## FAQ

Q: \[Is the fit_LOESS required?\] A: \[No. bmetapred can be used
directly when prediction errors (standard deviations) are already
present in the dataset. fit_LOESS is an optional method in case errors
are not pre-estimated\].

Q: \[Can I fit and predict with LOESS in the same train data used to fit
the bmetapred method?\] A: \[No. The LOESS model should be pre-trained
in historical/external data to avoid overfitting, and estimate the
errors in the train/test data with the predict function. If historical
data is not available, the user may used an Out-of-Fold (OOF)
cross-validation strategy on the training data to prevent overfitting.
The process is defined as follows: Divide the training dataset into K
mutually exclusive folds. For each fold k (from 1 to K):Train the LOESS
model on the remaining K-1 folds. Predict the errors on the hold-out
fold k. After completing all K iterations, we obtain unbiased,
out-of-sample errors for the entire training set. The base model is
finally retrained on the full training dataset. When predicting on the
unseen test set, the OOF-calibrated LOESS model is applied to estimate
the predictive standard deviations robustly\].

Q: \[In which scenarios is bmetapred expected to outperform other
methods?\] A: \[Other methods do not take into account the prediction
error of the combined models. bmetapred shines when prediction errors
are heteroscedastic; for instance, if prediction models are worse at
predicting extreme values of Y.\]
