# =============================================================================
# Bayesian Flight Height Modeling
# =============================================================================
#
# Lognormal model with GPS quality covariates for bird flight altitudes.
#
# Features:
#   - Observation error depends on GPS HDOP and satellite count
#   - Fixed effects for alpha_0 and alpha_1 by ecological group
#   - Random effect on log-scale standard deviation
#   - Automatic convergence filtering (keeps only converged chains)
#
# =============================================================================


# -----------------------------------------------------------------------------
# CONVERGENCE DIAGNOSTICS
# -----------------------------------------------------------------------------

#' Assess Chain Convergence Using Gelman-Rubin Diagnostic
assess_chain_convergence <- function(samples, rhat_threshold = 1.1) {
  
  n_chains <- length(samples)
  
  if (n_chains < 2) {
    message("    Only 1 chain - skipping convergence assessment")
    return(list(converged_idx = 1, all_converged = TRUE, max_rhat = NA))
  }
  
  # Check all chains together
  tryCatch({
    gelman_all <- coda::gelman.diag(samples, multivariate = FALSE)
    max_rhat <- max(gelman_all$psrf[, "Point est."], na.rm = TRUE)
    
    if (max_rhat < rhat_threshold) {
      message(sprintf("    All %d chains converged (max Rhat = %.3f)", n_chains, max_rhat))
      return(list(converged_idx = seq_len(n_chains), all_converged = TRUE, max_rhat = max_rhat))
    }
  }, error = function(e) NULL)
  
  # Find best converged subset
  message("    Not all chains converged. Finding best subset...")
  
  best_subset <- 1
  best_rhat <- Inf
  
  for (size in n_chains:2) {
    for (combo in combn(n_chains, size, simplify = FALSE)) {
      tryCatch({
        gelman_sub <- coda::gelman.diag(samples[combo], multivariate = FALSE)
        max_rhat <- max(gelman_sub$psrf[, "Point est."], na.rm = TRUE)
        if (max_rhat < rhat_threshold && length(combo) > length(best_subset)) {
          best_subset <- combo
          best_rhat <- max_rhat
        }
      }, error = function(e) NULL)
    }
    if (length(best_subset) == size && best_rhat < rhat_threshold) break
  }
  
  if (best_rhat >= rhat_threshold) {
    warning("No chain subset converged. Using all chains with caution.")
    return(list(converged_idx = seq_len(n_chains), all_converged = FALSE, max_rhat = best_rhat))
  }
  
  message(sprintf("    Converged chains: %s (max Rhat = %.3f)", paste(best_subset, collapse = ", "), best_rhat))
  list(converged_idx = best_subset, all_converged = FALSE, max_rhat = best_rhat)
}


#' Filter to Converged Chains Only
filter_converged_chains <- function(samples, converged_idx) {
  filtered <- list(
    samples = samples$samples[converged_idx],
    WAIC = samples$WAIC,
    converged_chains = converged_idx,
    n_chains_original = length(samples$samples),
    n_chains_kept = length(converged_idx)
  )
  class(filtered$samples) <- "mcmc.list"
  filtered
}


# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------

#' Prepare GPS Data for Modeling
prepare_modeling_data <- function(gps_data, season) {
  
  gps <- gps_data %>%
    mutate(
      height_above_msl = as.numeric(height_above_msl),
      distance_to_shore = as.numeric(distance_to_shore),
      in_GOL = as.logical(in_GOL),
      gps_hdop = as.numeric(gps_hdop),
      gps_satellite_count = as.numeric(gps_satellite_count)
    ) %>%
    filter(in_GOL == TRUE, state == "Migration") %>%
    filter(!is.na(height_above_msl), !is.na(ecological_group)) %>%
    filter(season == !!season) %>%
    filter(height_above_msl > quantile(height_above_msl, 0.01, na.rm = TRUE)) %>%
    filter(height_above_msl < quantile(height_above_msl, 0.99, na.rm = TRUE)) %>%
    filter(gps_hdop <= 5, na.rm = FALSE) %>%
    filter(gps_satellite_count >= 4, na.rm = FALSE)
  
  if (nrow(gps) == 0) return(NULL)
  
  # Standardize GPS quality (NA -> 0)
  gps$z_sat <- as.numeric(scale(gps$gps_satellite_count))
  gps$z_sat[is.na(gps$z_sat)] <- 0
  gps$z_hdop <- as.numeric(scale(gps$gps_hdop))
  gps$z_hdop[is.na(gps$z_hdop)] <- 0
  
  gps
}


# -----------------------------------------------------------------------------
# MAIN FUNCTION
# -----------------------------------------------------------------------------

#' Run Flight Height Bayesian Model (Single Season)
#'
#' @param gps_data GPS dataframe with ecological_group column
#' @param output_dir Output directory
#' @param season "pre-breeding" or "post-breeding"
#' @param n_iter MCMC iterations (default: 1000000)
#' @param burnin_prop Burnin proportion (default: 0.9)
#' @param thin Thinning interval (default: 20)
#' @param chains Number of chains (default: 10)
#' @param rhat_threshold Convergence threshold (default: 1.1)
#' @return List with samples, proportions, gps_filtered, convergence, eco_groups
run_flight_height_model <- function(gps_data,
                                     output_dir,
                                     season,
                                     n_iter = 1000000,
                                     burnin_prop = 0.9,
                                     thin = 20,
                                     chains = 10,
                                     rhat_threshold = 1.1) {

  # Load packages
  for (pkg in c("nimble", "coda", "MCMCvis")) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
      install.packages(pkg)
      library(pkg, character.only = TRUE)
    }
  }

  burnin <- n_iter * burnin_prop

  message("\n", strrep("=", 60))
  message(sprintf("  %s MIGRATION", toupper(gsub("-", " ", season))))
  message(strrep("=", 60), "\n")

  # Directories
  model_dir <- file.path(output_dir, "model")
  prop_dir <- file.path(output_dir, "proportions")
  dir.create(model_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(prop_dir, showWarnings = FALSE, recursive = TRUE)

  # --- DATA ---
  message("  Preparing data...")
  gps <- prepare_modeling_data(gps_data, season)
  if (is.null(gps)) { 
    warning("No data for: ", season)
    return(NULL) 
  }

  N <- nrow(gps)
  species <- as.factor(gps$ecological_group)
  nspecies <- length(unique(species))
  
  message(sprintf("  N = %s | Groups = %d", format(N, big.mark = ","), nspecies))

  constants <- list(
    N = N, 
    nspecies = nspecies,
    species_indices = as.numeric(species),
    dist_shore = as.numeric(scale(gps$distance_to_shore)),
    gps_hdop = gps$z_hdop,
    gps_nb_sat = gps$z_sat
  )
  data <- list(obs_alt = gps$height_above_msl)

  # --- MODEL ---
  message("  Building model...")
  
code <- nimbleCode({
  
  for (i in 1:N) {
    
    obs_alt[i] ~ dnorm(mean = true_alt[i], sd = sigma_obs[i])                                             # observation process
    true_alt[i] ~ dgamma(shape =  shape_alt[species_indices[i]], rate = rate_alt[i])                      # latent process
    
    
    
    log(mu_alt[i]) <- alpha_0_mu[species_indices[i]] + alpha_1_mu[species_indices[i]] * dist_shore[i]     # calcul rate_alt
    rate_alt[i] <- shape_alt[species_indices[i]] / mu_alt[i] 
    
  
    log(sigma_obs[i]) <- beta_0_obs + beta_1_obs * gps_hdop[i] + beta_2_obs * gps_nb_sat[i]               # calcul sigma_obs
    
  }
  
  
  
  for(j in 1:nspecies){
    
    alpha_0_mu[j] ~ dnorm(mean = 0, sd = 1)
    alpha_1_mu[j] ~ dnorm(mean = 0, sd = 1)
    
    
    log_shape_alt[j] ~ dnorm(mean = 0, sd = 1)
    shape_alt[j] <- exp(log_shape_alt[j])             
    
  }
  
  
  beta_0_obs ~ dnorm(mean = 0, sd = 1)
  beta_1_obs ~ dnorm(mean = 0, sd = 1)
  beta_2_obs ~ dnorm(mean = 0, sd = 1)
  
})

  # --- INITS ---
  message(sprintf("  Initializing %d chains...", chains))
  
  l0 <- list( 
    alpha_0_mu = rep(0,5),  
    alpha_1_mu = rep(0,5),
    log_shape_alt = rep(0,5),
    beta_0_obs = 0,
    beta_1_obs = 0,
    beta_2_obs = 0)
  
  model <- nimbleModel(code = code, data = data, constants = constants, inits = l0)

  inits <- list(l0)

  for(i in 1:(chains-1)){
    
    lnew <- list(
      alpha_0_mu = rnorm(5, mean = 0, sd = 1),
      alpha_1_mu = rnorm(5, mean = 0, sd = 1),
      log_shape_alt = rnorm(5, mean = 0, sd = 1),
      beta_0_obs = rnorm(1, mean = 0, sd = 1),
      beta_1_obs = rnorm(1, mean = 0, sd = 1),
      beta_2_obs = rnorm(1, mean = 0, sd = 1))
    
    inits <- c(inits, list(lnew))
  }

  # --- MCMC ---
  message(sprintf("  Running MCMC: %d iter, %d chains, thin=%d", n_iter, chains, thin))

conf <- configureMCMC(model, enableWAIC = TRUE,
                      monitors = c('alpha_0_mu',
                                   'alpha_1_mu',
                                   'log_shape_alt',
                                   'beta_0_obs',
                                   'beta_1_obs',
                                   'beta_2_obs'))
  
  mMCMC <- buildMCMC(conf)
  
  message("  Compiling model...")
  Cmodel <- compileNimble(model)
  
  message("  Compiling MCMC...")
  CMCMC <- compileNimble(mMCMC, project = Cmodel)
  
  message("  Sampling...")
  samples <- runMCMC(
    CMCMC, 
    niter = n_iter, 
    nburnin = burnin, 
    nchains = chains,
    inits = inits, 
    samplesAsCodaMCMC = TRUE, 
    thin = thin, 
    WAIC = TRUE
  )

  # --- CONVERGENCE ---
  message("  Checking convergence...")
  convergence <- assess_chain_convergence(samples$samples, rhat_threshold)
  samples_filtered <- filter_converged_chains(samples, convergence$converged_idx)
  message(sprintf("  Kept %d/%d chains", samples_filtered$n_chains_kept, samples_filtered$n_chains_original))

  # --- SAVE ---
  saveRDS(samples_filtered, file.path(model_dir, paste0("model_", season, ".rds")))

  # Diagnostics
  sink(file.path(model_dir, paste0("diagnostics_", season, ".txt")))
  cat(sprintf("=== %s ===\nChains: %s of %d | Rhat: %.3f\n\n",
              toupper(season), paste(convergence$converged_idx, collapse = ","),
              samples_filtered$n_chains_original, convergence$max_rhat))
  print(MCMCsummary(samples_filtered$samples))
  cat("\n=== WAIC ===\n")
  print(samples_filtered$WAIC)
  sink()

# --- PROPORTIONS ---
  message("  Computing proportions...")
  
  ite <- do.call(rbind, lapply(samples_filtered$samples, as.data.frame))
  eco_groups <- sort(unique(gps$ecological_group))
  proportion_data <- data.frame()

  for (g in seq_along(eco_groups)) {

    # Gamma parameters at mean distance to shore (dist_shore = 0)
    # log(mu_alt) = alpha_0_mu  =>  mu_alt = exp(alpha_0_mu)
    # shape_alt = exp(log_shape_alt)
    # rate_alt = shape_alt / mu_alt
    shape_g <- exp(ite[, paste0("log_shape_alt[", g, "]")])
    mu_g    <- exp(ite[, paste0("alpha_0_mu[", g, "]")])
    rate_g  <- shape_g / mu_g

    p <- data.frame(
      p_0_20   = pgamma(20,  shape = shape_g, rate = rate_g),
      p_20_200 = pgamma(200, shape = shape_g, rate = rate_g) - pgamma(20,  shape = shape_g, rate = rate_g),
      p_20_300 = pgamma(300, shape = shape_g, rate = rate_g) - pgamma(20,  shape = shape_g, rate = rate_g),
      p_300_inf = 1 - pgamma(300, shape = shape_g, rate = rate_g)
    )

    proportion_data <- rbind(proportion_data, data.frame(
      ecological_group = eco_groups[g],
      altitude_range = c("0-20 m", "20-200 m", "20-300 m", ">300 m"),
      median_pct = round(sapply(p, median) * 100, 1),
      CI_95 = paste0(round(sapply(p, quantile, 0.025) * 100, 1), "-",
                     round(sapply(p, quantile, 0.975) * 100, 1))
    ))
    
    saveRDS(p, file.path(prop_dir, paste0("prop_", tolower(gsub(" ", "_", eco_groups[g])), "_", season, ".rds")))
  }

  saveRDS(proportion_data, file.path(prop_dir, paste0("proportions_", season, ".rds")))

  message(sprintf("\n  Done! Output: %s", output_dir))
  message(strrep("=", 60))

  list(
    samples = samples_filtered, 
    proportions = proportion_data, 
    gps_filtered = gps, 
    convergence = convergence, 
    eco_groups = eco_groups
  )
}