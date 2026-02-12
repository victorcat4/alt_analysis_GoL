# =============================================================================
# Plot Bayesian Model Results
# =============================================================================
#
# Visualization functions for Gamma flight height model outputs.
#
# =============================================================================


# -----------------------------------------------------------------------------
# FLIGHT HEIGHT DISTRIBUTION PLOTTING
# -----------------------------------------------------------------------------

#' Plot Flight Height Distribution for Single Group
#'
#' Creates a Gamma density plot showing flight height distribution with
#' altitude bands color-coded and proportion estimates in legend.
#'
#' @param shape Gamma shape parameter (median estimate)
#' @param rate Gamma rate parameter (median estimate)
#' @param p_recap Dataframe with proportion summaries (median, CI_025, CI_975)
#' @param group_name Name of ecological group
#' @param season_label Season label for subtitle
#' @param logo_grob Optional logo grob to add to plot
#' @return ggplot object
plot_flight_height_distribution_group <- function(shape,
                                                   rate,
                                                   p_recap,
                                                   group_name,
                                                   season_label = "Post-breeding migration",
                                                   logo_grob = NULL) {

  x_vals <- seq(0.01, 1700, length.out = 10000)
  dens_vals <- dgamma(x_vals, shape = shape, rate = rate)

  pg_data <- data.frame(x = x_vals, y = dens_vals) %>%
    mutate(
      fill_group = case_when(
        x > 300 ~ "300_inf",
        x <= 20 ~ "0_20",
        x > 20 & x <= 200 ~ "20_200",
        x > 200 & x <= 300 ~ "200_300"
      )
    )

  label_0_20 <- sprintf("between 0-20 m = %.1f%% (%.1f-%.1f%%)",
                        p_recap$p_0_20[1] * 100,
                        p_recap$p_0_20[2] * 100,
                        p_recap$p_0_20[3] * 100)

  label_20_200 <- sprintf("between 20-200 m = %.1f%% (%.1f-%.1f%%)",
                          p_recap$p_20_200[1] * 100,
                          p_recap$p_20_200[2] * 100,
                          p_recap$p_20_200[3] * 100)

  label_200_300 <- sprintf("between 200-300 m = %.1f%% (%.1f-%.1f%%)",
                           p_recap$p_200_300[1] * 100,
                           p_recap$p_200_300[2] * 100,
                           p_recap$p_200_300[3] * 100)

  label_300_inf <- sprintf("above 300 m = %.1f%% (%.1f-%.1f%%)",
                           p_recap$p_300_inf[1] * 100,
                           p_recap$p_300_inf[2] * 100,
                           p_recap$p_300_inf[3] * 100)

  fill_colors <- c(
    "0_20" = "#cc4778",
    "20_200" = "#fde725",
    "200_300" = "#ed7953",
    "300_inf" = "#9c179e"
  )

  fill_labels <- c(
    "0_20" = label_0_20,
    "20_200" = label_20_200,
    "200_300" = label_200_300,
    "300_inf" = label_300_inf
  )

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
    ) +
    guides(fill = guide_legend(reverse = TRUE))

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
#' Generates flight height distribution plots for all ecological groups.
#'
#' @param samples MCMC samples object (filtered to converged chains)
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

  logo_grob <- NULL
  if (!is.null(logo_path) && file.exists(logo_path)) {
    if (requireNamespace("png", quietly = TRUE) && requireNamespace("grid", quietly = TRUE)) {
      logo_img <- png::readPNG(logo_path)
      logo_grob <- grid::rasterGrob(logo_img, interpolate = TRUE)
    }
  }

  sum_v <- MCMCsummary(samples$samples)

  eco_groups <- gps_data %>%
    select(ecological_group) %>%
    distinct() %>%
    arrange(ecological_group) %>%
    pull(ecological_group)

  nspecies <- length(eco_groups)

  ite <- do.call(rbind, lapply(samples$samples, as.data.frame))

  season_label <- paste(tools::toTitleCase(gsub("-", " ", season_filter)), "migration")

  for (g in seq_along(eco_groups)) {
    group_name <- eco_groups[g]

    message(sprintf("    - %s...", group_name))

    # Calculate proportions for all iterations using Gamma distribution
    # Gamma parameterization: shape = exp(log_shape_alt), mu = exp(alpha_0_mu), rate = shape / mu
    shape_g <- exp(ite[, paste0("log_shape_alt[", g, "]")])
    mu_g    <- exp(ite[, paste0("alpha_0_mu[", g, "]")])
    rate_g  <- shape_g / mu_g

    p_ite <- data.frame(
      p_0_20    = pgamma(20,  shape = shape_g, rate = rate_g),
      p_20_200  = pgamma(200, shape = shape_g, rate = rate_g) - pgamma(20,  shape = shape_g, rate = rate_g),
      p_200_300 = pgamma(300, shape = shape_g, rate = rate_g) - pgamma(200, shape = shape_g, rate = rate_g),
      p_300_inf = 1 - pgamma(300, shape = shape_g, rate = rate_g)
    )

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
      p_200_300 = c(
        quantile(p_ite$p_200_300, 0.5),
        quantile(p_ite$p_200_300, 0.025),
        quantile(p_ite$p_200_300, 0.975)
      ),
      p_300_inf = c(
        quantile(p_ite$p_300_inf, 0.5),
        quantile(p_ite$p_300_inf, 0.025),
        quantile(p_ite$p_300_inf, 0.975)
      )
    )
    rownames(p_recap) <- c("median", "CI_025", "CI_975")

    # Median parameters for the density curve
    shape_med <- exp(sum_v[paste0("log_shape_alt[", g, "]"), "50%"])
    mu_med    <- exp(sum_v[paste0("alpha_0_mu[", g, "]"), "50%"])
    rate_med  <- shape_med / mu_med

    plot <- plot_flight_height_distribution_group(
      shape = shape_med,
      rate = rate_med,
      p_recap = p_recap,
      group_name = group_name,
      season_label = season_label,
      logo_grob = logo_grob
    )

    filename <- sprintf("flight_height_%s_%s.pdf",
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


#' Create Proportion Summary Table
#'
#' Creates and saves a formatted table of altitude proportions by group.
#'
#' @param proportion_data Dataframe with proportion summaries
#' @param plots_dir Output directory for table
#' @param season_filter Season label
create_proportion_table <- function(proportion_data, plots_dir, season_filter) {
  
  if (!require("gt", quietly = TRUE)) {
    warning("Package 'gt' not available. Skipping table creation.")
    return(invisible(NULL))
  }
  
  proportion_table <- proportion_data %>%
    gt(groupname_col = "ecological_group") %>%
    cols_label(
      altitude_range = "Altitude Range",
      median_pct = "Proportion (%)",
      CI_95 = "95% CI"
    ) %>%
    tab_header(
      title = "Flight height proportions by altitude range",
      subtitle = paste(tools::toTitleCase(gsub("-", " ", season_filter)), "migration")
    )

  gtsave(proportion_table, file.path(plots_dir, paste0("proportions_table_", season_filter, ".html")))
  
  message(sprintf("  Proportion table saved: %s", file.path(plots_dir, paste0("proportions_table_", season_filter, ".pdf"))))
  
  invisible(proportion_table)
}

#' Create Quantile Summary Table
#'
#' Computes estimated Q25, Q50 (median), and Q75 of flight height from the
#' Gamma posterior, with 95% credible intervals, per ecological group.
#' Also reports sample size (N) per group.
#'
#' @param samples MCMC samples object (filtered to converged chains)
#' @param gps_data Filtered GPS data (same as used for model fitting)
#' @param plots_dir Output directory for table
#' @param season_filter Season label
create_quantile_summary_table <- function(samples,
                                           gps_data,
                                           plots_dir,
                                           season_filter) {

  if (!require("gt", quietly = TRUE)) {
    warning("Package 'gt' not available. Skipping quantile table creation.")
    return(invisible(NULL))
  }

  message("  Computing quantile summary table...")

  eco_groups <- sort(unique(gps_data$ecological_group))
  ite <- do.call(rbind, lapply(samples$samples, as.data.frame))

  # Sample sizes per group
  n_per_group <- gps_data %>%
    count(ecological_group, name = "N") %>%
    arrange(ecological_group)

  results <- list()

  for (g in seq_along(eco_groups)) {
    group_name <- eco_groups[g]

    # Posterior Gamma parameters at mean distance (dist_shore = 0)
    shape_g <- exp(ite[, paste0("log_shape_alt[", g, "]")])
    mu_g    <- exp(ite[, paste0("alpha_0_mu[", g, "]")])
    rate_g  <- shape_g / mu_g

    # Compute quantiles of the Gamma distribution for each posterior draw
    q25_draws  <- qgamma(0.25, shape = shape_g, rate = rate_g)
    q50_draws  <- qgamma(0.50, shape = shape_g, rate = rate_g)
    q75_draws  <- qgamma(0.75, shape = shape_g, rate = rate_g)

    N_g <- n_per_group$N[n_per_group$ecological_group == group_name]

    results[[g]] <- tibble(
      ecological_group = group_name,
      N = N_g,
      Q25_median = round(median(q25_draws), 1),
      Q25_CI95   = paste0(round(quantile(q25_draws, 0.025), 1), "-",
                          round(quantile(q25_draws, 0.975), 1)),
      Q50_median = round(median(q50_draws), 1),
      Q50_CI95   = paste0(round(quantile(q50_draws, 0.025), 1), "-",
                          round(quantile(q50_draws, 0.975), 1)),
      Q75_median = round(median(q75_draws), 1),
      Q75_CI95   = paste0(round(quantile(q75_draws, 0.025), 1), "-",
                          round(quantile(q75_draws, 0.975), 1))
    )
  }

  quantile_df <- bind_rows(results)

  season_label <- tools::toTitleCase(gsub("-", " ", season_filter))

  quantile_table <- quantile_df %>%
    gt() %>%
    cols_label(
      ecological_group = "Ecological Group",
      N = "N",
      Q25_median = "Q25 (m)",
      Q25_CI95   = "95% CI",
      Q50_median = "Q50 (m)",
      Q50_CI95   = "95% CI",
      Q75_median = "Q75 (m)",
      Q75_CI95   = "95% CI"
    ) %>%
    tab_spanner(label = "25th percentile", columns = c(Q25_median, Q25_CI95)) %>%
    tab_spanner(label = "Median", columns = c(Q50_median, Q50_CI95)) %>%
    tab_spanner(label = "75th percentile", columns = c(Q75_median, Q75_CI95)) %>%
    tab_header(
      title = "Estimated flight height quantiles (Gamma model)",
      subtitle = paste(season_label, "migration")
    ) %>%
    cols_align(align = "center", columns = -ecological_group)

  out_file <- file.path(plots_dir, paste0("quantile_summary_", season_filter, ".html"))
  gtsave(quantile_table, out_file)

  # Also save the underlying data as RDS
  saveRDS(quantile_df, file.path(plots_dir, paste0("quantile_summary_", season_filter, ".rds")))

  message(sprintf("  Quantile summary table saved: %s", out_file))
  invisible(quantile_table)
}

# -----------------------------------------------------------------------------
# PROPORTION BY DISTANCE TO SHORE
# -----------------------------------------------------------------------------

#' Plot Predicted Proportion Below 300m by Distance to Shore
#'
#' Computes the predicted proportion of flight heights below 300m as a function
#' of distance to shore from the Gamma model posterior, per ecological group.
#' Uses exp(alpha_0_mu + alpha_1_mu * dist_z) as the mean and exp(log_shape_alt)
#' as the shape parameter.
#'
#' @param samples MCMC samples object (filtered to converged chains)
#' @param gps_data Filtered GPS data (same as used for model fitting)
#' @param plots_dir Output directory for plots
#' @param season_filter Season label
#' @param n_dist Number of distance points to predict at (default: 500)
#' @param max_iter Maximum posterior iterations to use (default: 5000, for speed)
create_prop_by_distance_plot <- function(samples,
                                          gps_data,
                                          plots_dir,
                                          season_filter,
                                          n_dist = 500,
                                          max_iter = 5000) {

  message("  Generating proportion by distance to shore plot...")

  # Recover distance scaling parameters (same scale() as in model fitting)
  dist_mean <- mean(gps_data$distance_to_shore, na.rm = TRUE)
  dist_sd   <- sd(gps_data$distance_to_shore, na.rm = TRUE)

  # Use observed data range (don't extrapolate beyond what the model saw)
  dist_range_m <- range(gps_data$distance_to_shore, na.rm = TRUE)
  message(sprintf("    Data distance range: %.0f - %.0f m", dist_range_m[1], dist_range_m[2]))

  eco_groups <- sort(unique(gps_data$ecological_group))
  nspecies   <- length(eco_groups)

  # Combine posterior chains
  ite <- do.call(rbind, lapply(samples$samples, as.data.frame))

  # Subsample for speed
  n_ite <- nrow(ite)
  if (n_ite > max_iter) {
    ite <- ite[sample(n_ite, max_iter), ]
  }

  # Distance sequence
  dist_seq_m <- seq(dist_range_m[1], dist_range_m[2], length.out = n_dist)
  dist_seq_z <- (dist_seq_m - dist_mean) / dist_sd

  # Compute predictions per group
  results_list <- vector("list", nspecies)

  for (g in seq_along(eco_groups)) {
    alpha_0  <- ite[[paste0("alpha_0_mu[", g, "]")]]
    alpha_1  <- ite[[paste0("alpha_1_mu[", g, "]")]]
    shape_g  <- exp(ite[[paste0("log_shape_alt[", g, "]")]])

    # log(mu(d)) = alpha_0 + alpha_1 * d_z  =>  mu(d) = exp(alpha_0 + alpha_1 * d_z)
    # rate(d) = shape / mu(d)
    # Compute matrix [n_ite x n_dist]
    log_mu_mat <- outer(alpha_0, rep(1, length(dist_seq_z))) +
                  outer(alpha_1, dist_seq_z)
    mu_mat   <- exp(log_mu_mat)
    rate_mat <- shape_g / mu_mat  # shape_g is recycled across columns

    p_mat <- pgamma(300, shape = shape_g, rate = rate_mat)

    results_list[[g]] <- tibble(
      ecological_group = eco_groups[g],
      dist_m           = dist_seq_m,
      prop_med         = apply(p_mat, 2, median) * 100,
      prop_q025        = apply(p_mat, 2, quantile, 0.025) * 100,
      prop_q975        = apply(p_mat, 2, quantile, 0.975) * 100
    )
  }

  data_all <- bind_rows(results_list)

  # Convert distance to km
  data_all <- data_all %>%
    mutate(dist_km = dist_m / 1000)

  # Build plot
  season_label <- tools::toTitleCase(gsub("-", " ", season_filter))

  plot <- ggplot(
    data = data_all,
    aes(x = dist_km, fill = ecological_group, color = ecological_group)
  ) +
    geom_line(aes(y = prop_med)) +
    geom_ribbon(aes(y = prop_med, ymin = prop_q025, ymax = prop_q975), alpha = 0.2, color = NA) +
    labs(
      y = "Proportion of points below 300 m (%)",
      x = "Distance to shore (km)",
      title = paste(season_label, "migration"),
      fill = "Ecological group",
      color = "Ecological group"
    ) +
    coord_cartesian(ylim = c(0, 100)) +
    scale_x_continuous(breaks = seq(0, 80, by = 20), limits = c(0, 80)) +
    theme_classic(base_size = 12) +
    theme(
      legend.position = "bottom",
      panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
      panel.grid.minor = element_line(color = "grey95", linewidth = 0.15)
    )

  # Save
  base_name <- paste0("gamma_", gsub("-", "_", season_filter), "_prop_under_300")

  ggsave(plot = plot, filename = file.path(plots_dir, paste0(base_name, ".pdf")),
         device = "pdf", width = 7, height = 7)

  message(sprintf("  Distance plot saved: %s.pdf", base_name))
  invisible(plot)
}


# -----------------------------------------------------------------------------
# MCMC DIAGNOSTICS PLOTS
# -----------------------------------------------------------------------------

#' Create MCMC Trace Plots
#'
#' Generates trace plots for convergence assessment.
#'
#' @param samples MCMC samples object
#' @param plots_dir Output directory
#' @param season_filter Season label
create_trace_plots <- function(samples, plots_dir, season_filter) {
  
  message("  Generating trace plots...")
  
  filename <- file.path(plots_dir, paste0("mcmc_trace_", season_filter, ".pdf"))
  
  pdf(filename, width = 12, height = 8)
  MCMCtrace(samples$samples, pdf = FALSE)
  dev.off()
  
  message(sprintf("  Trace plots saved: %s", filename))
  invisible(NULL)
}


# -----------------------------------------------------------------------------
# MAIN PLOTTING FUNCTION
# -----------------------------------------------------------------------------

#' Generate All Model Plots
#'
#' Creates all visualization outputs for a fitted model.
#'
#' @param model_results Results from run_flight_height_model() or run_single_season_model()
#' @param plots_dir Output directory for plots
#' @param season_filter Season label
#' @param logo_path Optional path to logo image
#' @param create_trace Whether to create trace plots (default: TRUE)
generate_model_plots <- function(model_results,
                                  plots_dir,
                                  season_filter,
                                  logo_path = NULL,
                                  create_trace = TRUE) {
  
  if (!require("MCMCvis", quietly = TRUE)) {
    stop("Package 'MCMCvis' is required. Install with: install.packages('MCMCvis')")
  }
  
  dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)
  
  message(sprintf("\nGenerating plots for %s...", season_filter))
  
  # Flight height distribution plots
  create_flight_height_plots(
    samples = model_results$samples,
    gps_data = model_results$gps_filtered,
    plots_dir = plots_dir,
    season_filter = season_filter,
    logo_path = logo_path
  )
  
  # Proportion table
  create_proportion_table(
    proportion_data = model_results$proportions,
    plots_dir = plots_dir,
    season_filter = season_filter
  )

  # Quantile summary table (Q25, Q50, Q75 with CI)
  create_quantile_summary_table(
    samples = model_results$samples,
    gps_data = model_results$gps_filtered,
    plots_dir = plots_dir,
    season_filter = season_filter
  )

  # Proportion by distance to shore
  create_prop_by_distance_plot(
    samples = model_results$samples,
    gps_data = model_results$gps_filtered,
    plots_dir = plots_dir,
    season_filter = season_filter
  )
  
  # Trace plots
  if (create_trace) {
    create_trace_plots(
      samples = model_results$samples,
      plots_dir = plots_dir,
      season_filter = season_filter
    )
  }
  
  message("  All plots generated.")
  invisible(NULL)
}


# -----------------------------------------------------------------------------
# STANDALONE PLOTTING FROM SAVED MODEL
# -----------------------------------------------------------------------------

#' Generate Plots from Saved Model Results
#'
#' Loads a previously saved model and generates all plots.
#'
#' @param model_path Path to saved model RDS file
#' @param gps_data GPS dataframe (same filtering will be applied)
#' @param plots_dir Output directory for plots
#' @param season_filter Season filter to apply to GPS data
#' @param logo_path Optional path to logo image
#' @param create_trace Whether to create trace plots (default: TRUE)
plot_from_saved_model <- function(model_path,
                                   gps_data,
                                   plots_dir,
                                   season_filter = "post-breeding",
                                   logo_path = NULL,
                                   create_trace = TRUE) {

  if (!require("MCMCvis", quietly = TRUE)) {
    stop("Package 'MCMCvis' is required. Install with: install.packages('MCMCvis')")
  }

  message("Loading saved model...")
  samples <- readRDS(model_path)

  # Filter GPS data (same as in model)
  gps <- gps_data %>%
    mutate(
      location_long = as.numeric(location_long),
      location_lat = as.numeric(location_lat),
      height_above_msl = as.numeric(height_above_msl),
      distance_to_shore = as.numeric(distance_to_shore),
      in_GOL = as.logical(in_GOL)
    ) %>%
    filter(in_GOL == TRUE, state == "Migration") %>%
    filter(!is.na(height_above_msl), !is.na(ecological_group)) %>%
    filter(season == season_filter) %>%
    filter(height_above_msl > quantile(height_above_msl, 0.02, na.rm = TRUE)) %>%
    filter(height_above_msl < quantile(height_above_msl, 0.98, na.rm = TRUE))

  message(sprintf("Filtered data: %s observations", format(nrow(gps), big.mark = ",")))
  message(sprintf("Model used %d converged chains", samples$n_chains_kept))

  # Load proportion data if available
  prop_dir <- dirname(dirname(model_path))
  prop_file <- file.path(prop_dir, "proportions", paste0("proportions_", season_filter, ".rds"))  

  if (file.exists(prop_file)) {
    proportions <- readRDS(prop_file)
  } else {
    message("  Proportion data not found. Recalculating...")
    proportions <- NULL
  }

  # Create mock results object
  model_results <- list(
    samples = samples,
    gps_filtered = gps,
    proportions = proportions
  )

  dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)

  # Generate plots
  create_flight_height_plots(
    samples = samples,
    gps_data = gps,
    plots_dir = plots_dir,
    season_filter = season_filter,
    logo_path = logo_path
  )
  
  if (!is.null(proportions)) {
    create_proportion_table(
      proportion_data = proportions,
      plots_dir = plots_dir,
      season_filter = season_filter
    )
  }

  # Quantile summary table (Q25, Q50, Q75 with CI)
  create_quantile_summary_table(
    samples = model_results$samples,
    gps_data = model_results$gps_filtered,
    plots_dir = plots_dir,
    season_filter = season_filter
  )

  # Proportion by distance to shore
  create_prop_by_distance_plot(
    samples = samples,
    gps_data = gps,
    plots_dir = plots_dir,
    season_filter = season_filter
  )
  
  if (create_trace) {
    create_trace_plots(
      samples = samples,
      plots_dir = plots_dir,
      season_filter = season_filter
    )
  }

  message("\nPlotting complete!")
  invisible(NULL)
}