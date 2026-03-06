# -------------------------------------------------------------
# Spatiotemporal Bayesian analysis of COVID-19 incidence in Italy
#
# This script performs the full workflow for the hierarchical
# spatiotemporal model implemented in Stan:
#   1. Data acquisition and preprocessing
#   2. Preparation of Stan input data
#   3. Model fitting using cmdstanr
#   4. Posterior diagnostics and WAIC computation
#   5. Posterior predictive analysis and forecasting plots
#   6. Residual diagnostics for the temporal AR(1) component
# -------------------------------------------------------------


# -------------------------------------------------------------
# LIBRARIES
# Load required packages for data manipulation, Bayesian inference,
# diagnostics, and visualization
library(tidyverse)
library(lubridate)
library(cmdstanr)
library(posterior)
library(bayesplot)
library(loo)
library(tidybayes)
library(lmtest) # package used for residual diagnostics

# -------------------------------------------------------------
# 0. WORKING ENVIRONMENT
# Set the working directory containing the Stan model and output files
setwd("/Users/sguar/Desktop/DESKTOP ONE DRIVE/ELEONORA/STATISTICA COMPUTAZIONALE/MODELLO_LOG")

# -------------------------------------------------------------
# 1. DATA LOADING AND PREPROCESSING
# Import official COVID-19 regional data and prepare the dataset
dati_ufficiali <- read_csv(
  "https://raw.githubusercontent.com/pcm-dpc/COVID-19/master/dati-regioni/dpc-covid19-ita-regioni.csv",
  show_col_types = FALSE
)

dati <- dati_ufficiali %>%
  select(any_of(c("data","denominazione_regione",
                  "nuovi_positivi","lat","long"))) %>%
  mutate(data = as.Date(data)) %>%
  filter(!denominazione_regione %in% c("Sardegna","Sicilia"))

# Replace negative daily counts (reporting corrections) with zero
dati$nuovi_positivi[dati$nuovi_positivi < 0] <- 0

# Regional population data used as a demographic covariate in the model
popolazione_regioni <- tibble(
  denominazione_regione = c("Abruzzo","Basilicata","Calabria","Campania",
                            "Emilia-Romagna","Friuli Venezia Giulia","Lazio","Liguria","Lombardia",
                            "Marche","Molise","Piemonte","Puglia","Toscana","P.A. Trento",
                            "P.A. Bolzano","Umbria","Valle d'Aosta","Veneto"),
  popolazione = c(1293941,539999,1812871,5575025,4459477,1211357,5710272,
                  1508120,10035481,1501323,296547,4241875,3956430,3645691,
                  540958,532616,859300,125034,4851851)
)

data <- left_join(dati, popolazione_regioni,
                  by = "denominazione_regione")

# Small constant used to avoid log(0) in the response transformation
eps <- 0.5

# -------------------------------------------------------------
# 2. FUNCTIONS

## 2.1 Prepare data for a specific epidemic phase
prep_fase <- function(df, start, end){
  data_filtered <- df %>%
    filter(data >= as.Date(start),
           data <= as.Date(end))
  
  
  data_model <<- data_filtered %>%
    filter(!is.na(nuovi_positivi)) %>%
    mutate(
      regione_id = as.numeric(factor(denominazione_regione)),
      time_id    = as.numeric(factor(data)),
      pop        = popolazione,
      dow        = wday(data, week_start = 1),
      y_log = log(nuovi_positivi + eps)
      
    ) %>% arrange(regione_id, time_id)
  
  # Compute the spatial distance matrix between regional centroids
  coord_regioni <- data_model %>%
    select(regione_id, lat, long) %>% distinct() %>%
    arrange(regione_id)
  Dmat_space <<- as.matrix(dist(coord_regioni[,c("lat","long")]))
  
  # Extract the population vector ordered by region
  pop_vec <<- data_model %>% select(regione_id,popolazione) %>%
    distinct() %>% arrange(regione_id) %>% pull(popolazione)
  
  R <<- length(unique(data_model$regione_id))
  Time <<- length(unique(data_model$time_id))
  
  # Temporal 80–20 train/test split within the selected phase
  cut_idx <- floor(0.8*Time)
  data_train <<- data_model %>% filter(time_id <= cut_idx)
  data_test  <<- data_model %>% filter(time_id >  cut_idx)
  
  N_train <<- nrow(data_train)
  N_test  <<- nrow(data_test)
}

## 2.2 Build the list of inputs required by the Stan model
build_stan <- function(){
  stan_data <<- list(
    N_train = N_train, N_test = N_test,
    R = R, T = Time,
    y_log_train = data_train$y_log, 
    region_train = data_train$regione_id,
    time_train   = data_train$time_id,
    dow_train    = data_train$dow,
    region_test  = data_test$regione_id,
    time_test    = data_test$time_id,
    dow_test     = data_test$dow,
    pop = pop_vec,
    Dmat_space = Dmat_space
  )
}

## 2.3 Fit the Bayesian model using cmdstanr
fit_model <- function(tag){
  # Output file used to store the fitted object for the current phase
  fit_file <- paste0("fit_", tag, ".rds")
  
  mod1 <<- cmdstan_model("covid19_spatiotemporal_model.stan")
  fit1 <<- mod1$sample(
    data = stan_data,
    chains = 2, parallel_chains = 2,
    iter_sampling = 1000, iter_warmup = 1000,
    max_treedepth = 15, adapt_delta = 0.99,
    seed = 42, refresh = 200
  )
  
  # Save the fitted object for reproducibility and later inspection
  fit1$save_object(file = fit_file)
  
  message("→ Salvato ", fit_file)
}

## 2.4 Posterior diagnostics and WAIC computation
diagnostica <- function(tag){
  # Extract the full posterior sample
  posterior_draws <<- fit1$draws()
  
  # Parameters of primary interest for posterior summaries
  pars_interesse <- c(
    "mu",
    "sigma2_y",
    "gamma",
    "sigma2",
    "rho",
    "alpha",
    "phi_s",
    "beta"
  )
  summary_df <- fit1$summary(variables = pars_interesse)
  
  # Save a printed version of the posterior summary as a text file
  capture.output(
    print(summary_df),
    file = paste0("summary_", tag, ".txt")
  )
  
  # Open a PDF device for MCMC trace plots
  pdf(paste0("trace_", tag, ".pdf"), width = 10, height = 6)
  
  # Generate and export trace plots for convergence assessment
  p_trace <- mcmc_trace(
    posterior_draws,
    pars = c(
      "mu",
      "sigma2_y",
      paste0("gamma[", 1:6, "]"),
      "sigma2",
      "rho",
      "alpha",
      "phi_s",
      "beta"
    )
  )
  print(p_trace)    
  
  dev.off()
  
  # Compute WAIC from the pointwise log-likelihood
  loglik_1 <<- fit1$draws("log_lik") %>% as_draws_matrix()
  waic_val <<- waic(loglik_1)
  write_csv(
    as_tibble(waic_val$estimates, rownames = "metric"),
    paste0("waic_", tag, ".csv")
  )
}


## 2.5 Posterior predictive analysis and forecast plots
pred_plot <- function(tag){
  # Extract posterior predictive simulations on the log scale
  y_test_log <- fit1$draws("y_log_pred") %>% as_draws_matrix()
  
  # Posterior predictive means
  y_log_pred_mean <- colMeans(y_test_log)
  y_pred_mean     <- exp(y_log_pred_mean) - eps
  
  # 95% predictive intervals transformed back to the natural scale
  ci_nat <- apply(y_test_log, 2, function(v){
    qs <- quantile(v, probs = c(.025, .975))
    exp(qs) - eps
  })
  y_pred_lower <- ci_nat[1, ]
  y_pred_upper <- ci_nat[2, ]
  
  # Add posterior predictive quantities to the test set
  data_test <<- data_test %>%
    mutate(
      y_pred_mean     = y_pred_mean,
      y_pred_lower    = y_pred_lower,
      y_pred_upper    = y_pred_upper,
      y_log_pred_mean = y_log_pred_mean
    )
  
  # Add placeholder columns to the training set for plotting consistency
  data_train <<- data_train %>%
    mutate(
      y_pred_mean     = NA_real_,
      y_pred_lower    = NA_real_,
      y_pred_upper    = NA_real_,
      y_log_pred_mean = NA_real_
    )
  
  data_all <- bind_rows(data_train, data_test)
  
  # Export one forecast plot per region to a PDF file
  pdf(paste0("forecast_", tag, ".pdf"), width = 8, height = 6)
  for(reg in unique(data_all$denominazione_regione)){
    df_reg <- data_all %>% filter(denominazione_regione == reg)
    p <- ggplot(df_reg, aes(x = data)) +
      # Observed values over the full time series
      geom_line(aes(y = nuovi_positivi), color = "black") +
      # Posterior predictive mean on the test period only
      geom_line(aes(y = y_pred_mean),
                color    = "blue",
                linetype = "dashed",
                na.rm    = TRUE) +
      # Predictive interval on the test period
      geom_ribbon(aes(ymin = y_pred_lower, ymax = y_pred_upper),
                  alpha = .3,
                  na.rm = TRUE) +
      labs(
        title = paste0("Regione: ", reg, " (Fase ", tag, ")"),
        y     = "Nuovi positivi",
        x     = "Data"
      ) +
      theme_minimal()
    print(p)
  }
  dev.off()
  
  # Mean squared error on the natural scale
  mse_nat <<- mean((data_test$nuovi_positivi - data_test$y_pred_mean)^2, na.rm = TRUE)
  
  # Mean squared error on the log scale
  y_log_obs <- log(data_test$nuovi_positivi + eps)
  mse_log <<- mean((y_log_obs - data_test$y_log_pred_mean)^2, na.rm = TRUE)
  
  # Compute RMSE and normalized RMSE
  # Natural-scale RMSE and NRMSE normalized by the mean observed incidence
  rmse_nat   <- sqrt(mse_nat)
  mean_obs   <- mean(data_test$nuovi_positivi, na.rm = TRUE)
  nrmse_nat  <- rmse_nat / mean_obs
  
  # Log-scale RMSE and NRMSE normalized by the standard deviation of log-observed values
  rmse_log   <- sqrt(mse_log)
  sd_log_obs <- sd(y_log_obs, na.rm = TRUE)
  nrmse_log  <- rmse_log / sd_log_obs
  
  # Save predictive performance metrics
  write_csv(
    tibble(
      phase      = tag,
      mse_nat    = mse_nat,
      mse_log    = mse_log,
      rmse_nat   = rmse_nat,
      nrmse_nat  = nrmse_nat,
      rmse_log   = rmse_log,
      nrmse_log  = nrmse_log
    ),
    paste0("metrics_", tag, ".csv")
  )
  
  # Compute and save average predictions by region
  media_pred_mod1 <- data_test %>%
    group_by(denominazione_regione) %>%
    summarise(y_pred_mean_regione = mean(y_pred_mean, na.rm = TRUE),
              .groups = "drop") %>%
    arrange(denominazione_regione)
  
  media_obs <- data_test %>%
    group_by(denominazione_regione) %>%
    summarise(y_obs_mean_regione = mean(nuovi_positivi, na.rm = TRUE),
              .groups = "drop") %>%
    arrange(denominazione_regione)
  
  result_mod1 <- left_join(media_obs, media_pred_mod1, by = "denominazione_regione")
  write_csv(result_mod1, paste0("previsioni_", tag, ".csv"))
}


## 2.6 Residual diagnostics for the latent AR(1) temporal component
## Assess whether regional temporal residuals behave as white noise
## Durbin–Watson statistics are also computed to evaluate residual independence
residui_plot <- function(tag){
  # Extract all posterior draws of the temporal latent process w_t
  draws_df <- as_draws_df(fit1$draws(variables = "w_t"))
  
  n_iter <- nrow(draws_df)
  
  w_t_array <- array(NA, dim = c(n_iter, R, Time))
  
  for (r in 1:R) {
    for (t in 1:Time) {
      varname <- paste0("w_t[", r, ",", t, "]")
      w_t_array[, r, t] <- draws_df[[varname]]
    }
  }
  
  # Extract posterior draws of alpha
  alpha_vec <- as_draws_df(fit1$draws("alpha"))$alpha
  
  # Build the residual matrix for each region
  res_mat_list <- vector("list", R)
  
  for (r in 1:R) {
    res_r <- matrix(NA, nrow = n_iter, ncol = Time)
    
    for (iter in 1:n_iter) {
      # The first residual has no previous time point and is kept unchanged
      res_r[iter, 1] <- w_t_array[iter, r, 1]
      
      for (t in 2:Time) {
        res_r[iter, t] <- w_t_array[iter, r, t] - alpha_vec[iter] * w_t_array[iter, r, t - 1]
      }
    }
    
    res_mat_list[[r]] <- res_r
  }
  
  region_names <- data_model %>%
    distinct(regione_id, denominazione_regione) %>%
    arrange(regione_id) %>%
    pull(denominazione_regione)
  
  
  # Safety check to ensure consistency between residuals and region names
  stopifnot(length(res_mat_list) == length(region_names))
  
  # Export ACF, PACF, and density plots of residual means to PDF
  pdf(paste0("residui_", tag, ".pdf"), width = 10, height = 5)
  
  for (r in seq_along(res_mat_list)) {
    # One row and three panels per region
    par(mfrow = c(1, 3))
    
    # Posterior mean residual at each time point
    res_mean_r <- colMeans(res_mat_list[[r]], na.rm = TRUE)
    
    # Residual autocorrelation function
    acf(res_mean_r, main = paste("ACF residui -", region_names[r]))
    
    # Partial autocorrelation function of residuals
    pacf(res_mean_r, main = paste("PACF residui -", region_names[r]))
    
    # Kernel density estimate of residual means
    plot(density(res_mean_r), main = paste("Densità residui -", region_names[r]))
  }
  
  # Close the graphics device
  dev.off()
  
  
  # Durbin–Watson test for residual independence
  
  # Tibble used to store region-specific DW statistics
  dw_stats <- tibble(
    regione = region_names,
    D_w     = NA_real_
  )
  
  # Loop over regions
  for (r in seq_along(res_mat_list)) {
    res_mean_r <- colMeans(res_mat_list[[r]], na.rm = TRUE)
    
    # Lagged residuals
    r_tm1 <- res_mean_r[1:(Time-1)]
    
    # Current residuals
    r_t   <- res_mean_r[2:Time]
    
    # Linear regression without intercept
    fit_ar1 <- lm(r_t ~ r_tm1 - 1)
    
    # Durbin–Watson test
    dw_res  <- dwtest(fit_ar1)
    
    # Store the test statistic for each region
    dw_stats$D_w[r] <- as.numeric(dw_res$statistic)
  }
  
  # Lower and upper bounds of the acceptable DW interval
  dw_lower <- 1.7
  dw_upper <- 2.3
  
  # Reorder regions by DW value for clearer visualization
  dw_stats <- dw_stats %>%
    arrange(desc(D_w)) %>%
    mutate(regione = factor(regione, levels = regione))
  
  # Plot the Durbin–Watson statistic by region
  pdf(paste0("test_", tag, ".pdf"), width = 10, height = 5)
  
  p_dw <- ggplot(dw_stats, aes(x = regione, y = D_w)) +
    # Highlight the acceptable DW interval
    annotate(
      "rect",
      xmin = -Inf,
      xmax = Inf,
      ymin = dw_lower,
      ymax = dw_upper,
      alpha = 0.2,
      fill = "lightgreen"
    ) +
    # Regional DW statistics
    geom_col(fill = "steelblue", width = 0.7) +
    # Interval boundaries
    geom_hline(yintercept = dw_lower, color = "darkgreen", linetype = "dashed") +
    geom_hline(yintercept = dw_upper, color = "darkgreen", linetype = "dashed") +
    # Reference line at DW = 2
    geom_hline(yintercept = 2, color = "red", linetype = "solid", size = 0.8) +
    coord_flip() +
    labs(
      title = paste0("Durbin–Watson statistic on mean residuals (interval [", 
                     dw_lower, ", ", dw_upper, "])"),
      x = "Regione",
      y = expression(D[w])
    ) +
    theme_minimal(base_size = 12)
  
  print(p_dw)
  dev.off()
} 




# -------------------------------------------------------------
# 3. EPIDEMIC PHASES
# Define the epidemic phases corresponding to different stages
# of COVID-19 spread and public health interventions in Italy
fasi <- tribble(
  ~fase,  ~start,        ~end,
  "F1",   "2020-03-09",  "2020-05-03",  # national lockdown
  "F2",  "2020-05-04",  "2020-12-31",   # pre-vaccination color-coded restrictions
  "F3",  "2021-01-01",  "2021-03-31",   # vaccination campaign start
  "F4",  "2021-04-01",  "2021-11-30",   # Alpha–Delta variant dominance
  #  "F5",  "2021-12-01",  "2022-03-31",   # Omicron wave
  #  "F6",   "2022-04-01",  "2023-05-05"   # late emergency period
)

# -------------------------------------------------------------
# 4. MAIN LOOP
# For each epidemic phase:
#   - prepare the dataset
#   - build Stan inputs
#   - fit the model
#   - perform posterior diagnostics
#   - generate predictions and forecast plots
#   - run residual diagnostics
for(i in 1:nrow(fasi)){
  tag <- fasi$fase[i]
  prep_fase(data, fasi$start[i], fasi$end[i])
  build_stan()
  fit_model(tag)
  diagnostica(tag)
  pred_plot(tag)
  residui_plot(tag)
  cat(">>> completata", tag, "\n")
}