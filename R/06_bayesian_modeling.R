# =============================================================================
# Bayesian Flight Height Modeling
# =============================================================================
#
# Hierarchical Bayesian analysis of bird flight altitudes during migration
# using a lognormal distribution model. Supports both pre-breeding and
# post-breeding migration seasons.
#
# =============================================================================


# -----------------------------------------------------------------------------
# FLIGHT HEIGHT DISTRIBUTION PLOTTING
# -----------------------------------------------------------------------------

#' Plot Flight Height Distribution for Single Group
#'
#' Creates a lognormal density plot showing flight height distribution with
#' altitude bands color-coded and proportion estimates in legend.
#'
#' @param meanlog Lognormal meanlog parameter (median estimate)
#' @param sdlog Lognormal sdlog parameter (median estimate)
#' @param p_recap Dataframe with proportion summaries (median, CI_025, CI_975)
#' @param group_name Name of ecological group
#' @param season_label Season label for subtitle
#' @param logo_grob Optional logo grob to add to plot
#' @return ggplot object
plot_flight_height_distribution_group <- function(meanlog,
                                                   sdlog,
                                                   p_recap,
                                                   group_name,
                                                   season_label = "Post-breeding migration",
                                                   logo_grob = NULL) {

  # Generate density curve data
  x_vals <- seq(0, 1700, length.out = 10000)
  dens_vals <- dlnorm(x_vals, meanlog = meanlog, sdlog = sdlog)

  # Create data frame with altitude band classification
  pg_data <- data.frame(x = x_vals, y = dens_vals) %>%
    mutate(
      fill_group = case_when(
        x > 300 ~ "300_inf",
        x <= 20 ~ "0_20",
        x > 20 & x <= 200 ~ "20_200",
        x > 200 & x <= 300 ~ "200_300"
      )
    )

  # Create legend labels with proportion estimates and CIs
  label_0_20 <- sprintf("between 0-20 m = %.1f%% (%.1f-%.1f%%)",
                        p_recap$p_0_20[1] * 100,
                        p_recap$p_0_20[2] * 100,
                        p_recap$p_0_20[3] * 100)

  label_20_200 <- sprintf("between 20-200 m = %.1f%% (%.1f-%.1f%%)",
                          p_recap$p_20_200[1] * 100,
                          p_recap$p_20_200[2] * 100,
                          p_recap$p_20_200[3] * 100)

  label_200_300 <- sprintf("between 20-300 m = %.1f%% (%.1f-%.1f%%)",
                           p_recap$p_20_300[1] * 100,
                           p_recap$p_20_300[2] * 100,
                           p_recap$p_20_300[3] * 100)

  label_300_inf <- sprintf("above 300 m = %.1f%% (%.1f-%.1f%%)",
                           p_recap$p_300_inf[1] * 100,
                           p_recap$p_300_inf[2] * 100,
                           p_recap$p_300_inf[3] * 100)

  # Define fill colors
  fill_colors <- c(
    "0_20" = "#cc4778",
    "20_200" = "#fde725",
    "200_300" = "#ed7953",
    "300_inf" = "#9c179e"
  )

  # Define labels for legend
  fill_labels <- c(
    "0_20" = label_0_20,
    "20_200" = label_20_200,
    "200_300" = label_200_300,
    "300_inf" = label_300_inf
  )

  # Create base plot
  plot_base <- ggplot(pg_data, aes(x = x, y = y, fill = fill_group)) +
    geom_area(alpha = 0.6) +
    geom_line(color = "black") +
    geom_vline(xintercept = 20, linetype = "longdash", col = "grey") +
    geom_vline(xintercept = 200, linetype = "longdash", col = "grey") +
    geom_vline(xintercept = 300, linetype = "longdash") +
    coord_flip(xlim = c(0, 1600)) +
    scale_fill_manual(
      name = "Proportion of points (95% CI)",
      values = fill_colors,
      labels = fill_labels
    ) +
    labs(
      title = "Distribution of flight heights in the Gulf of Lion",
      subtitle = paste(season_label, group_name, sep = ", "),
      x = "Flight height (m)",
      y = "Frequency"
    ) +
    theme_classic() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.position = c(0.65, 0.4)
    )

  # Add logo if provided
  if (!is.null(logo_grob)) {
    if (requireNamespace("cowplot", quietly = TRUE)) {
      plot_base <- cowplot::ggdraw() +
        cowplot::draw_plot(plot_base) +
        cowplot::draw_grob(logo_grob, y = 0.45, width = 0.07)
    }
  }

  plot_base
}


#' Create All Flight Height Distribution Plots
#'
#' Generates flight height distribution plots for all ecological groups
#' and the population level.
#'
#' @param samples MCMC samples object
#' @param gps_data Filtered GPS data
#' @param plots_dir Output directory for plots
#' @param season_filter Season label
#' @param logo_path Optional path to logo image
create_flight_height_plots <- function(samples,
                                        gps_data,
                                        plots_dir,
                                        season_filter = "post-breeding",
                                        logo_path = NULL) {

  message("  Generating flight height distribution plots...")

  # Load logo if provided
  logo_grob <- NULL
  if (!is.null(logo_path) && file.exists(logo_path)) {
    if (requireNamespace("png", quietly = TRUE) && requireNamespace("grid", quietly = TRUE)) {
      logo_img <- png::readPNG(logo_path)
      logo_grob <- grid::rasterGrob(logo_img, interpolate = TRUE)
    }
  }

  # Get MCMC summary
  sum_v <- MCMCsummary(samples$samples)

  # Get ecological groups
  eco_groups <- gps_data %>%
    select(ecological_group) %>%
    distinct() %>%
    arrange(ecological_group) %>%
    pull(ecological_group)

  nspecies <- length(eco_groups)

  # Rename population-level parameters in summary
  rownames(sum_v)[rownames(sum_v) == "alpha_0_mu"] <- paste0("alpha_0_re_mu[", nspecies + 1, "]")
  rownames(sum_v)[rownames(sum_v) == "alpha_1_mu"] <- paste0("alpha_1_re_mu[", nspecies + 1, "]")
  rownames(sum_v)[rownames(sum_v) == "beta_0_sigma"] <- paste0("beta_0_re_sigma[", nspecies + 1, "]")

  # Combine MCMC chains for iteration-level calculations
  ite <- rbind(as.data.frame(samples$samples[[1]]), as.data.frame(samples$samples[[2]]))
  names(ite)[names(ite) == "alpha_0_mu"] <- paste0("alpha_0_re_mu[", nspecies + 1, "]")
  names(ite)[names(ite) == "alpha_1_mu"] <- paste0("alpha_1_re_mu[", nspecies + 1, "]")
  names(ite)[names(ite) == "beta_0_sigma"] <- paste0("beta_0_re_sigma[", nspecies + 1, "]")

  # Define group names including population
  group_names <- c(eco_groups, "Population")

  # Format season label
  season_label <- paste(tools::toTitleCase(gsub("-", " ", season_filter)), "migration")

  # Generate plot for each group
  for (g in seq_along(group_names)) {
    group_name <- group_names[g]
    group_idx <- if (g <= nspecies) g else nspecies + 1

    message(sprintf("    - %s...", group_name))

    # Calculate proportions for all iterations
    p_ite <- data.frame(
      p_0_20 = numeric(nrow(ite)),
      p_20_200 = numeric(nrow(ite)),
      p_20_300 = numeric(nrow(ite)),
      p_300_inf = numeric(nrow(ite))
    )

    for (i in 1:nrow(ite)) {
      meanlog <- ite[i, paste0("alpha_0_re_mu[", group_idx, "]")]
      sdlog <- ite[i, paste0("beta_0_re_sigma[", group_idx, "]")]

      p_ite$p_0_20[i] <- plnorm(20, meanlog, sdlog)
      p_ite$p_20_200[i] <- plnorm(200, meanlog, sdlog) - plnorm(20, meanlog, sdlog)
      p_ite$p_20_300[i] <- plnorm(300, meanlog, sdlog) - plnorm(20, meanlog, sdlog)
      p_ite$p_300_inf[i] <- 1 - plnorm(300, meanlog, sdlog)
    }

    # Create summary table for proportions
    p_recap <- data.frame(
      p_0_20 = c(
        quantile(p_ite$p_0_20, 0.5),
        quantile(p_ite$p_0_20, 0.025),
        quantile(p_ite$p_0_20, 0.975)
      ),
      p_20_200 = c(
        quantile(p_ite$p_20_200, 0.5),
        quantile(p_ite$p_20_200, 0.025),
        quantile(p_ite$p_20_200, 0.975)
      ),
      p_20_300 = c(
        quantile(p_ite$p_20_300, 0.5),
        quantile(p_ite$p_20_300, 0.025),
        quantile(p_ite$p_20_300, 0.975)
      ),
      p_300_inf = c(
        quantile(p_ite$p_300_inf, 0.5),
        quantile(p_ite$p_300_inf, 0.025),
        quantile(p_ite$p_300_inf, 0.975)
      )
    )
    rownames(p_recap) <- c("median", "CI_025", "CI_975")

    # Get median parameter estimates for density curve
    meanlog <- sum_v[paste0("alpha_0_re_mu[", group_idx, "]"), "50%"]
    sdlog <- sum_v[paste0("beta_0_re_sigma[", group_idx, "]"), "50%"]

    # Create plot
    plot <- plot_flight_height_distribution_group(
      meanlog = meanlog,
      sdlog = sdlog,
      p_recap = p_recap,
      group_name = group_name,
      season_label = season_label,
      logo_grob = logo_grob
    )

    # Save plot
    filename <- sprintf("flight_height_%s_%s.png",
                        tolower(gsub(" ", "_", group_name)),
                        gsub("-", "_", season_filter))

    ggsave(
      filename = file.path(plots_dir, filename),
      plot = plot,
      width = 8,
      height = 6,
      dpi = 300
    )
  }

  message("  Flight height plots saved.")
  invisible(NULL)
}


# -----------------------------------------------------------------------------
# CORE MODELING FUNCTION (SINGLE SEASON)
# -----------------------------------------------------------------------------

#' Run Flight Height Model for Single Season
#'
#' Internal function that fits the Bayesian model for one season.
#'
#' @param gps_data Processed GPS dataframe with ecological_group column
#' @param output_dir Output directory for model results
#' @param plots_dir Directory for plots
#' @param season_filter Season to analyze
#' @param n_iter MCMC iterations
#' @param burnin_prop Burnin proportion
#' @param thin Thinning interval
#' @param chains Number of chains
#' @param logo_path Optional path to logo image
#' @param create_plots Whether to generate plots
#' @return List with samples, proportions, and filtered GPS data
run_single_season_model <- function(gps_data,
                                     output_dir,
                                     plots_dir,
                                     season_filter,
                                     n_iter,
                                     burnin_prop,
                                     thin,
                                     chains,
                                     logo_path,
                                     create_plots) {

  burnin <- n_iter * burnin_prop

  # Create output directories
  model_dir <- file.path(output_dir, "model")
  prop_dir <- file.path(output_dir, "proportions")
  dir.create(model_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(prop_dir, showWarnings = FALSE, recursive = TRUE)


  # --- DATA PREPARATION ---

  message("  Preparing data for modeling...")

  gps <- gps_data %>%
    filter(in_GOL == TRUE, state == "Migration") %>%
    filter(!is.na(height_above_msl), !is.na(ecological_group)) %>%
    filter(season == season_filter) %>%
    filter(height_above_msl > quantile(height_above_msl, 0.02, na.rm = TRUE)) %>%
    filter(height_above_msl < 5000)

  message(sprintf("  Filtered data: %s observations", format(nrow(gps), big.mark = ",")))

  if (nrow(gps) == 0) {
    warning(sprintf("No data available for season: %s", season_filter))
    return(NULL)
  }

  # Prepare model inputs
  obs_alt <- gps$height_above_msl
  N <- length(obs_alt)

  species <- as.numeric(as.factor(gps$ecological_group))
  nspecies <- length(unique(species))

  dist_shore <- scale(gps$distance_to_shore)[, 1]

  constants <- list(N = N, dist_shore = dist_shore, nspecies = nspecies, species_indices = species)
  data <- list(obs_alt = obs_alt)


  # --- MODEL DEFINITION ---

  message("  Building model...")

  code <- nimbleCode({
    for (i in 1:N) {
      obs_alt[i] ~ dnorm(mean = true_alt[i], sd = sigma_obs)
      true_alt[i] ~ dlnorm(meanlog = mu[i], sdlog = sigma_lognorm[i])
      mu[i] <- alpha_0_re_mu[species_indices[i]] + alpha_1_re_mu[species_indices[i]] * dist_shore[i]
      sigma_lognorm[i] <- beta_0_re_sigma[species_indices[i]]
    }

    for (j in 1:nspecies) {
      alpha_0_re_mu[j] ~ dnorm(alpha_0_mu, sd = alpha_0_sigma)
      beta_0_re_sigma[j] ~ dnorm(beta_0_sigma, sd = beta_0_sigma_sigma)
      alpha_1_re_mu[j] ~ dnorm(alpha_1_mu, sd = alpha_1_sigma)
    }

    sigma_obs ~ dunif(min = 0, max = 50)
    alpha_0_sigma ~ dlnorm(meanlog = 0, sdlog = 2)
    alpha_1_sigma ~ dlnorm(meanlog = 0, sdlog = 2)
    beta_0_sigma_sigma ~ dlnorm(meanlog = 0, sdlog = 2)
    alpha_0_mu ~ dnorm(mean = 0, sd = 2)
    alpha_1_mu ~ dnorm(mean = 0, sd = 2)
    beta_0_sigma ~ dlnorm(meanlog = 0, sdlog = 2)
  })


  # --- INITIAL VALUES ---

  set.seed(42)

  inits <- lapply(1:chains, function(x) {
    list(
      sigma_obs = runif(1, 0, 30),
      alpha_0_mu = rnorm(1, 0, 2),
      alpha_1_mu = rnorm(1, 0, 2),
      beta_0_sigma_sigma = rlnorm(1, 0, 2),
      alpha_0_sigma = rlnorm(1, 0, 2),
      alpha_1_sigma = rlnorm(1, 0, 2),
      beta_0_sigma = rlnorm(1, 0, 2)
    )
  })


  # --- RUN MCMC ---

  message(sprintf("  Running MCMC: %d iterations, %d chains...", n_iter, chains))
  message("  (This may take several minutes)")

  model <- nimbleModel(code = code, data = data, constants = constants, inits = inits)

  conf <- configureMCMC(
    model,
    enableWAIC = TRUE,
    monitors = c('alpha_0_mu', 'alpha_1_mu', 'beta_0_sigma',
                 'alpha_0_sigma', 'alpha_1_sigma', 'beta_0_sigma_sigma',
                 'sigma_obs', 'alpha_0_re_mu', 'alpha_1_re_mu', 'beta_0_re_sigma')
  )

  mMCMC <- buildMCMC(conf)
  Cmodel <- compileNimble(model)
  CMCMC <- compileNimble(mMCMC, project = Cmodel)

  samples <- runMCMC(
    CMCMC,
    niter = n_iter,
    nburnin = burnin,
    nchains = chains,
    samplesAsCodaMCMC = TRUE,
    thin = thin,
    WAIC = TRUE,
    setSeed = 1:chains
  )

  # Save model
  model_file <- file.path(model_dir, paste0("model_lognorm_", season_filter, ".rds"))
  saveRDS(samples, file = model_file)
  message(sprintf("  Model saved: %s", model_file))


  # --- DIAGNOSTICS ---

  message("  Saving diagnostics...")

  diag_file <- file.path(model_dir, paste0("convergence_diagnostics_", season_filter, ".txt"))
  sink(diag_file)
  cat(sprintf("=== MCMC Summary Statistics (%s) ===\n\n", season_filter))
  print(MCMCsummary(samples$samples))
  cat("\n\n=== WAIC ===\n")
  print(samples$WAIC)
  sink()


  # --- PROPORTION CALCULATIONS ---

  message("  Calculating altitude proportions...")

  sum_v <- MCMCsummary(samples$samples)

  # Rename population-level parameters
  rownames(sum_v)[rownames(sum_v) == "alpha_0_mu"] <- paste0("alpha_0_re_mu[", nspecies + 1, "]")
  rownames(sum_v)[rownames(sum_v) == "alpha_1_mu"] <- paste0("alpha_1_re_mu[", nspecies + 1, "]")
  rownames(sum_v)[rownames(sum_v) == "beta_0_sigma"] <- paste0("beta_0_re_sigma[", nspecies + 1, "]")

  ite <- rbind(as.data.frame(samples$samples[[1]]), as.data.frame(samples$samples[[2]]))
  names(ite)[names(ite) == "alpha_0_mu"] <- paste0("alpha_0_re_mu[", nspecies + 1, "]")
  names(ite)[names(ite) == "alpha_1_mu"] <- paste0("alpha_1_re_mu[", nspecies + 1, "]")
  names(ite)[names(ite) == "beta_0_sigma"] <- paste0("beta_0_re_sigma[", nspecies + 1, "]")

  # Get group names
  eco_groups <- gps %>%
    select(ecological_group) %>%
    distinct() %>%
    arrange(ecological_group) %>%
    pull(ecological_group)

  group_names <- c(eco_groups, "Population")

  # Calculate proportions for each group
  proportion_data <- data.frame()

  for (g in seq_along(group_names)) {
    group_name <- group_names[g]
    group_idx <- if (g <= nspecies) g else nspecies + 1

    p_ite <- data.frame(
      p_0_20 = numeric(nrow(ite)),
      p_20_200 = numeric(nrow(ite)),
      p_20_300 = numeric(nrow(ite)),
      p_300_inf = numeric(nrow(ite))
    )

    for (i in 1:nrow(ite)) {
      meanlog <- ite[i, paste0("alpha_0_re_mu[", group_idx, "]")]
      sdlog <- ite[i, paste0("beta_0_re_sigma[", group_idx, "]")]

      p_ite$p_0_20[i] <- plnorm(20, meanlog, sdlog)
      p_ite$p_20_200[i] <- plnorm(200, meanlog, sdlog) - plnorm(20, meanlog, sdlog)
      p_ite$p_20_300[i] <- plnorm(300, meanlog, sdlog) - plnorm(20, meanlog, sdlog)
      p_ite$p_300_inf[i] <- 1 - plnorm(300, meanlog, sdlog)
    }

    group_data <- data.frame(
      ecological_group = group_name,
      altitude_range = c("0-20 m", "20-200 m", "20-300 m", "Above 300 m"),
      proportion_median = round(c(
        median(p_ite$p_0_20),
        median(p_ite$p_20_200),
        median(p_ite$p_20_300),
        median(p_ite$p_300_inf)
      ) * 100, 1),
      CI_95 = paste0(
        round(c(
          quantile(p_ite$p_0_20, 0.025),
          quantile(p_ite$p_20_200, 0.025),
          quantile(p_ite$p_20_300, 0.025),
          quantile(p_ite$p_300_inf, 0.025)
        ) * 100, 1),
        "-",
        round(c(
          quantile(p_ite$p_0_20, 0.975),
          quantile(p_ite$p_20_200, 0.975),
          quantile(p_ite$p_20_300, 0.975),
          quantile(p_ite$p_300_inf, 0.975)
        ) * 100, 1)
      )
    )

    proportion_data <- rbind(proportion_data, group_data)

    # Save individual group results
    saveRDS(p_ite, file = file.path(prop_dir, paste0("proportions_", 
                                                      tolower(gsub(" ", "_", group_name)), 
                                                      "_", season_filter, ".rds")))
  }

  # Save combined proportion data
  saveRDS(proportion_data, file = file.path(prop_dir, paste0("proportion_data_", season_filter, ".rds")))

  # Create summary table
  proportion_table <- proportion_data %>%
    gt(groupname_col = "ecological_group") %>%
    cols_label(
      altitude_range = "Altitude Range",
      proportion_median = "Proportion (%)",
      CI_95 = "95% CI"
    ) %>%
    tab_header(
      title = "Flight height proportions by altitude range",
      subtitle = paste(tools::toTitleCase(gsub("-", " ", season_filter)), "migration")
    )

  gtsave(proportion_table, file.path(plots_dir, paste0("proportions_table_", season_filter, ".html")))
  gtsave(proportion_table, file.path(plots_dir, paste0("proportions_table_", season_filter, ".png")))


  # --- FLIGHT HEIGHT DISTRIBUTION PLOTS ---

  if (create_plots) {
    create_flight_height_plots(
      samples = samples,
      gps_data = gps,
      plots_dir = plots_dir,
      season_filter = season_filter,
      logo_path = logo_path
    )
  }

  list(samples = samples, proportions = proportion_data, gps_filtered = gps)
}


# -----------------------------------------------------------------------------
# MAIN MODELING FUNCTIONS
# -----------------------------------------------------------------------------

#' Run Flight Height Bayesian Model
#'
#' Fits a hierarchical Bayesian model to flight height data and generates
#' diagnostic plots, proportion estimates, and distribution visualizations.
#'
#' @param gps_data Processed GPS dataframe with ecological_group column
#' @param output_dir Output directory for model results
#' @param plots_dir Directory for plots
#' @param season_filter Season to analyze: "post-breeding", "pre-breeding", or "both"
#' @param n_iter MCMC iterations (default: 50000)
#' @param burnin_prop Burnin proportion (default: 0.2)
#' @param thin Thinning interval (default: 5)
#' @param chains Number of chains (default: 2)
#' @param logo_path Optional path to logo image for plots
#' @param create_plots Whether to generate flight height distribution plots (default: TRUE)
#' @return List with results for each season analyzed
run_flight_height_model <- function(gps_data,
                                     output_dir,
                                     plots_dir,
                                     season_filter = "both",
                                     n_iter = 50000,
                                     burnin_prop = 0.2,
                                     thin = 5,
                                     chains = 2,
                                     logo_path = NULL,
                                     create_plots = TRUE) {

  # Load modeling packages
  packages <- c("nimble", "coda", "MCMCvis", "gt")
  for (pkg in packages) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
      install.packages(pkg, dependencies = TRUE)
      library(pkg, character.only = TRUE)
    }
  }

  # Determine which seasons to analyze
  if (season_filter == "both") {
    seasons <- c("pre-breeding", "post-breeding")
  } else {
    seasons <- season_filter
  }

  results <- list()

  for (season in seasons) {
    message("\n", strrep("=", 60))
    message(sprintf("ANALYZING %s MIGRATION", toupper(gsub("-", " ", season))))
    message(strrep("=", 60), "\n")

    results[[season]] <- run_single_season_model(
      gps_data = gps_data,
      output_dir = output_dir,
      plots_dir = plots_dir,
      season_filter = season,
      n_iter = n_iter,
      burnin_prop = burnin_prop,
      thin = thin,
      chains = chains,
      logo_path = logo_path,
      create_plots = create_plots
    )

    message(sprintf("\n  %s modeling complete!", tools::toTitleCase(gsub("-", " ", season))))
  }

  message("\n", strrep("=", 60))
  message("ALL MODELING COMPLETE")
  message(strrep("=", 60))
  message(sprintf("\nResults saved to: %s", output_dir))

  invisible(results)
}


# -----------------------------------------------------------------------------
# STANDALONE PLOTTING FROM SAVED MODEL
# -----------------------------------------------------------------------------

#' Generate Plots from Saved Model Results
#'
#' Loads a previously saved model and generates flight height distribution plots.
#' Useful for re-generating plots without re-running the MCMC.
#'
#' @param model_path Path to saved model RDS file
#' @param gps_data GPS dataframe (same filtering will be applied)
#' @param plots_dir Output directory for plots
#' @param season_filter Season filter to apply to GPS data
#' @param logo_path Optional path to logo image
plot_from_saved_model <- function(model_path,
                                   gps_data,
                                   plots_dir,
                                   season_filter = "post-breeding",
                                   logo_path = NULL) {

  # Load required packages
  if (!require("MCMCvis", quietly = TRUE)) {
    stop("Package 'MCMCvis' is required. Install with: install.packages('MCMCvis')")
  }

  message("Loading saved model...")
  samples <- readRDS(model_path)

  # Filter GPS data (same as in model)
  gps <- gps_data %>%
    filter(in_GOL == TRUE, state == "Migration") %>%
    filter(!is.na(height_above_msl), !is.na(ecological_group)) %>%
    filter(season == season_filter) %>%
    filter(height_above_msl > quantile(height_above_msl, 0.02, na.rm = TRUE)) %>%
    filter(height_above_msl < 5000)

  message(sprintf("Filtered data: %s observations", format(nrow(gps), big.mark = ",")))

  # Generate plots
  create_flight_height_plots(
    samples = samples,
    gps_data = gps,
    plots_dir = plots_dir,
    season_filter = season_filter,
    logo_path = logo_path
  )

  message("Plotting complete!")
  invisible(NULL)
}