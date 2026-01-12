# ==============================================================================
# POST-BREEDING MIGRATION FLIGHT HEIGHT ANALYSIS
# ==============================================================================
# Description: Hierarchical Bayesian analysis of bird flight altitudes during
#              post-breeding migration using a lognormal distribution model
# Authors: Victor Cat, Maud Queroué, Nicolas Courbin, Coline Canonne
# ==============================================================================

# SETUP ------------------------------------------------------------------------

# Load required packages
library(nimble)
library(coda)
library(dplyr)
library(ggplot2)
library(MCMCvis)
library(patchwork)
library(transport)
library(sf)
library(readxl)
library(gt)
library(grid)
library(png)
library(cowplot)

theme_set(theme_bw())

# MCMC parameters
n_iter <- 50000 
burnin <- n_iter * 0.2 
thin <- 5
chains <- 2
n_sim <- 2000

# Ecological group names
ecological_group_names <- c(
  "1" = "Ardéidé",
  "2" = "Grand échassier", 
  "3" = "Laridé",
  "4" = "Limicole",
  "5" = "Terrestre",
  "6" = "Population"
)

# DATA PREPARATION -------------------------------------------------------------

# Load ecological classification
EC <- read_excel("data/supp/EC.xlsx")

# Load and filter GPS data for post-breeding migration
gps <- readRDS("data/gps_migration_4.rds") %>% 
  filter(in_GOL == TRUE & state == "Migration") %>% 
  filter(!is.na(height_above_msl)) %>%
  filter(season == "pre-breeding") %>%
  left_join(EC, by = c("individual_taxon_canonical_name" = "specie_scientific")) %>% 
  filter(height_above_msl > quantile(height_above_msl, 0.02)) %>%
  filter(height_above_msl < 5000)

# Prepare model data
obs_alt <- gps %>% pull(height_above_msl)
N <- length(obs_alt)

species <- gps %>% pull(ecological_group) %>% as.factor()
nspecies <- length(unique(species))
species_indices <- as.numeric(species)

dist_shore <- gps %>% 
  pull(distance_to_shore) %>% 
  scale(center = TRUE, scale = TRUE) %>% 
  as.vector()

constants <- list(
  N = N, 
  dist_shore = dist_shore,
  nspecies = nspecies, 
  species_indices = species_indices
)

data <- list(obs_alt = obs_alt)

# MODEL DEFINITION -------------------------------------------------------------

code <- nimbleCode({
  # Observation and latent process models
  for (i in 1:N) {
    obs_alt[i] ~ dnorm(mean = true_alt[i], sd = sigma_obs)
    true_alt[i] ~ dlnorm(meanlog = mu[i], sdlog = sigma_lognorm[i])
    mu[i] <- alpha_0_re_mu[species_indices[i]] + alpha_1_re_mu[species_indices[i]] * dist_shore[i]
    sigma_lognorm[i] <- beta_0_re_sigma[species_indices[i]]
  }
  
  # Random effects by ecological group
  for (j in 1:nspecies) {
    alpha_0_re_mu[j] ~ dnorm(alpha_0_mu, sd = alpha_0_sigma)         
    beta_0_re_sigma[j] ~ dnorm(beta_0_sigma, sd = beta_0_sigma_sigma)          
    alpha_1_re_mu[j] ~ dnorm(alpha_1_mu, sd = alpha_1_sigma)
  }

  # Priors
  sigma_obs ~ dunif(min = 0, max = 50)
  alpha_0_sigma ~ dlnorm(meanlog = 0, sdlog = 2)
  alpha_1_sigma ~ dlnorm(meanlog = 0, sdlog = 2)
  beta_0_sigma_sigma ~ dlnorm(meanlog = 0, sdlog = 2)
  alpha_0_mu ~ dnorm(mean = 0, sd = 2)
  alpha_1_mu ~ dnorm(mean = 0, sd = 2)
  beta_0_sigma ~ dlnorm(meanlog = 0, sdlog = 2)
})

# INITIAL VALUES ---------------------------------------------------------------

set.seed(10)

inits <- list(
  list(
    sigma_obs = runif(1, min = 0, max = 30),
    alpha_0_mu = rnorm(1, mean = 0, sd = 2),
    alpha_1_mu = rnorm(1, mean = 0, sd = 2),
    beta_0_sigma_sigma = rlnorm(1, meanlog = 0, sdlog = 2),
    alpha_0_sigma = rlnorm(1, meanlog = 0, sdlog = 2),
    alpha_1_sigma = rlnorm(1, meanlog = 0, sdlog = 2),
    beta_0_sigma = rlnorm(1, meanlog = 0, sdlog = 2)
  ),
  list(
    sigma_obs = runif(1, min = 0, max = 30),
    alpha_0_mu = rnorm(1, mean = 0, sd = 2),
    alpha_1_mu = rnorm(1, mean = 0, sd = 2),
    beta_0_sigma_sigma = rlnorm(1, meanlog = 0, sdlog = 2),
    alpha_0_sigma = rlnorm(1, meanlog = 0, sdlog = 2),
    alpha_1_sigma = rlnorm(1, meanlog = 0, sdlog = 2),
    beta_0_sigma = rlnorm(1, meanlog = 0, sdlog = 2)
  )
)

print(inits)

# RUN MCMC ---------------------------------------------------------------------

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
  setSeed = c(1, 2)
)

saveRDS(samples, file = "output/model_lognorm_preB.rds")

# CONVERGENCE DIAGNOSTICS ------------------------------------------------------

MCMCsummary(samples$samples)
MCMCtrace(samples$samples, pdf = FALSE)
MCMCplot(samples$samples)

# MODEL VALIDATION -------------------------------------------------------------

# Simulate data from fitted model
est_df <- data.frame(
  param = rownames(MCMCsummary(samples$samples)), 
  estimate = MCMCsummary(samples$samples)$mean, 
  sd = MCMCsummary(samples$samples)$sd
)

sim <- tibble()

for (i in 1:length(obs_alt)) {
  idx <- species_indices[i]
  
  alpha_0_re_mu <- est_df$estimate[est_df$param == paste0("alpha_0_re_mu[", idx, "]")]
  alpha_1_re_mu <- est_df$estimate[est_df$param == paste0("alpha_1_re_mu[", idx, "]")]
  sigma <- est_df$estimate[est_df$param == paste0("beta_0_re_sigma[", idx, "]")]
  
  alt_theo <- rlnorm(100, meanlog = alpha_0_re_mu + alpha_1_re_mu * dist_shore[i], sdlog = sigma)
  
  temp <- tibble(
    alt_theo = alt_theo,
    ecological_group = rep(gps$ecological_group[i], times = 100),
    dist = rep(dist_shore[i], times = 100),
    alt_obs = rnorm(100, mean = alt_theo, sd = sigma)
  )
  
  sim <- bind_rows(sim, temp)
}

# Compare simulated vs observed distributions
colors <- c("Simulated theoric alt" = "blue", 
            "Simulated observed alt" = "red", 
            "Observed alt" = "black")

plot <- ggplot() +
  geom_histogram(
    data = sim, 
    aes(x = alt_theo, y = after_stat(density), fill = "Simulated theoric alt"), 
    binwidth = 50, color = "black", position = "identity", size = 0.2
  ) +
  geom_histogram(
    data = gps, 
    aes(x = height_above_msl, y = after_stat(density), fill = "Observed alt"), 
    binwidth = 50, color = "black", alpha = 0.3, size = 0.2
  ) +
  labs(x = "Hauteur de vol (m)", y = "Nombre de points", title = "") +
  scale_fill_manual(values = colors) +
  theme_minimal() +
  facet_wrap(~ecological_group, nrow = 2) +
  coord_flip(xlim = c(-200, 2000)) + 
  theme(
    legend.position = "none", 
    axis.text.x = element_blank(),
    axis.title = element_text(size = 9), 
    axis.text.y = element_text(size = 7),
    panel.grid.major = element_line(size = 0.2),
    panel.grid.minor = element_line(size = 0.1)
  )

print(plot)
ggsave(plot = plot, filename = "output/plots/compare_out_lognorm_preB.png")

# PROPORTION CALCULATIONS ------------------------------------------------------

# Load logo
logo_img <- readPNG("data/supp/north_logo.png")
logo_grob <- rasterGrob(logo_img, interpolate = TRUE)

out <- readRDS(file = paste0(path, "output/model_lognorm_preB.rds"))
sum_v <- MCMCsummary(out$samples)

# Add this renaming:
rownames(sum_v)[rownames(sum_v) == "alpha_0_mu"] <- "alpha_0_re_mu[6]"
rownames(sum_v)[rownames(sum_v) == "alpha_1_mu"] <- "alpha_1_re_mu[6]"
rownames(sum_v)[rownames(sum_v) == "beta_0_sigma"] <- "beta_0_re_sigma[6]"


ite <- rbind(as.data.frame(out$samples[[1]]), as.data.frame(out$samples[[2]])) %>%
   rename(
    `alpha_0_re_mu[6]` = alpha_0_mu,
    `alpha_1_re_mu[6]` = alpha_1_mu,
    `beta_0_re_sigma[6]` = beta_0_sigma
  ) 

# Define groups to analyze
groups <- c(1, 2, 3, 4, 5, 6)

# Define group names (adjust as needed)
group_names <- c(
  "1" = "Ardeid",
  "2" = "Large wader",
  "3" = "Larid",
  "4" = "Shorebird",
  "5" = "Terrestrial",
  "6" = "Population"  # Adjust this name
)

# Initialize data frame to store all proportions
proportion_data_preB <- data.frame()

# Loop through each group
for (group_sp in groups) {
  
  group_name <- group_names[as.character(group_sp)]
  
  # Proportions in each altitude class by iterations 
  cols <- c("p_0_20", "p_20_200", "p_20_300", "p_300_inf")
  p_ite <- data.frame(matrix(NA, nrow = nrow(ite), ncol = length(cols)))
  colnames(p_ite) <- cols
  
  for (i in 1:nrow(ite)) {
    
    meanlog <- ite[i, paste0("alpha_0_re_mu[", group_sp, "]")]
    sdlog <- ite[i, paste0("beta_0_re_sigma[", group_sp, "]")]
    
    # Calculate proportions 
    p_ite$p_0_20[i] <- plnorm(20, meanlog = meanlog, sdlog = sdlog) 
    p_ite$p_20_200[i] <- plnorm(200, meanlog = meanlog, sdlog = sdlog) - plnorm(20, meanlog = meanlog, sdlog = sdlog) 
    p_ite$p_20_300[i] <- plnorm(300, meanlog = meanlog, sdlog = sdlog) - plnorm(20, meanlog = meanlog, sdlog = sdlog) 
    p_ite$p_300_inf[i] <- 1 - plnorm(300, meanlog = meanlog, sdlog = sdlog)
    
  }
  
  # Summary table
  p_recap <- data.frame(matrix(NA, nrow = 3, ncol = length(cols)))
  colnames(p_recap) <- cols
  rownames(p_recap) <- c("median", "CI_025", "CI_975")
  
  p_recap$p_0_20[1] <- quantile(p_ite$p_0_20, probs = 0.5)
  p_recap$p_0_20[2] <- quantile(p_ite$p_0_20, probs = 0.025)
  p_recap$p_0_20[3] <- quantile(p_ite$p_0_20, probs = 0.975)
  
  p_recap$p_20_200[1] <- quantile(p_ite$p_20_200, probs = 0.5)
  p_recap$p_20_200[2] <- quantile(p_ite$p_20_200, probs = 0.025)
  p_recap$p_20_200[3] <- quantile(p_ite$p_20_200, probs = 0.975)
  
  p_recap$p_20_300[1] <- quantile(p_ite$p_20_300, probs = 0.5)
  p_recap$p_20_300[2] <- quantile(p_ite$p_20_300, probs = 0.025)
  p_recap$p_20_300[3] <- quantile(p_ite$p_20_300, probs = 0.975)
  
  p_recap$p_300_inf[1] <- quantile(p_ite$p_300_inf, probs = 0.5) 
  p_recap$p_300_inf[2] <- quantile(p_ite$p_300_inf, probs = 0.025)
  p_recap$p_300_inf[3] <- quantile(p_ite$p_300_inf, probs = 0.975)
  
  # Print results for this group
  cat("\n=== Group", group_sp, "-", group_name, "===\n")
  print(p_recap)
  
  # Add to proportion data frame
  group_data <- data.frame(
    ecological_group = group_name,
    altitude_range = c("0-20 m", "20-200 m", "20-300 m", "Above 300 m"),
    proportion_median = round(c(p_recap$p_0_20[1], p_recap$p_20_200[1], 
                                 p_recap$p_20_300[1], p_recap$p_300_inf[1]) * 100, 1),
    CI_95 = paste0(round(c(p_recap$p_0_20[2], p_recap$p_20_200[2], 
                           p_recap$p_20_300[2], p_recap$p_300_inf[2]) * 100, 1),
                   "-",
                   round(c(p_recap$p_0_20[3], p_recap$p_20_200[3], 
                           p_recap$p_20_300[3], p_recap$p_300_inf[3]) * 100, 1))
  )
  
  proportion_data_preB <- rbind(proportion_data_preB, group_data)
  
  # Plot 
  meanlog <- sum_v[paste0("alpha_0_re_mu[", group_sp, "]"), 4] 
  sdlog <- sum_v[paste0("beta_0_re_sigma[", group_sp, "]"), 4]
  
  # Generate a sequence of x values
  x_vals <- seq(0, 1700, length.out = 10000)
  
  # Calculate lognormal density for each x
  dens_vals <- dlnorm(x_vals, meanlog = meanlog, sdlog = sdlog)
  
  # Put into a data frame
  pg_data <- data.frame(x = x_vals, y = dens_vals) %>%
    dplyr::mutate(
      fill_group = case_when(
        x_vals > 300 ~ "300_inf",
        x_vals <= 20 ~ "0_20",
        x_vals > 20 & x_vals <= 200 ~ "20_200",
        x_vals > 200 & x_vals <= 300 ~ "200_300"
      )
    )
  
  # Create labels for legend
  label_0_20 <- paste0("between 0–20 m = ", round(p_recap$p_0_20[1] * 100, digits = 1), 
                       "% (", round(p_recap$p_0_20[2] * 100, digits = 1), 
                       "-", round(p_recap$p_0_20[3] * 100, digits = 1), "%)")
  
  label_20_200 <- paste0("between 20–200 m = ", round(p_recap$p_20_200[1] * 100, digits = 1), 
                         "% (", round(p_recap$p_20_200[2] * 100, digits = 1), 
                         "-", round(p_recap$p_20_200[3] * 100, digits = 1), "%)")
  
  label_200_300 <- paste0("between 20–300 m = ", round(p_recap$p_20_300[1] * 100, digits = 1), 
                          "% (", round(p_recap$p_20_300[2] * 100, digits = 1), 
                          "-", round(p_recap$p_20_300[3] * 100, digits = 1), "%)")
  
  label_300_inf <- paste0("above 300 m = ", round(p_recap$p_300_inf[1] * 100, digits = 1), 
                          "% (", round(p_recap$p_300_inf[2] * 100, digits = 1), 
                          "-", round(p_recap$p_300_inf[3] * 100, digits = 1), "%)")
  
  # Define fill colors with correct factor levels as names
  fill_colors <- c(
    "0_20" = "#cc4778",
    "20_200" = "#fde725",
    "200_300" = "#ed7953",
    "300_inf" = "#9c179e"
  )
  
  # Define labels separately
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
    scale_fill_manual(name = "Proportion of points (95% CI)", 
                      values = fill_colors,
                      labels = fill_labels) +
    labs(
      title = "Distribution of flight heights in the Gulf of Lion",
      subtitle = paste("Post-breeding migration,", group_name, sep = " "),
      x = "Flight height (m)", 
      y = "Frequency"
    ) +
    theme_classic() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.position = c(0.65, 0.4)
    )
  
  # Add logo
  plot <- ggdraw() +
    draw_plot(plot_base) +
    draw_grob(logo_grob, y = 0.45, width = 0.07)
  
  # Save plot
  ggsave(
    filename = paste0(path, "output/plots/group_", group_sp, "_preB.png"),
    plot = plot,
    width = 8,
    height = 6,
    dpi = 300
  )
  
  # Save results table
  saveRDS(p_recap, file = paste0(path, "output/proportions/proportions_group_", group_sp, "_preB.rds"))
  
}

# Create and save summary table
proportion_table <- proportion_data_preB %>%
  gt(groupname_col = "ecological_group") %>%
  cols_label(
    altitude_range = "Altitude Range",
    proportion_median = "Proportion (%)",
    CI_95 = "95% CI"
  ) %>%
  tab_header(
    title = "Flight height proportions by altitude range",
    subtitle = "Post-breeding migration"
  )

# Display the table
print(proportion_table)

# Save the table as HTML
gtsave(proportion_table, 
       filename = paste0(path, "output/plots/proportions_table_preB.html"))

# Save the table as PNG
gtsave(proportion_table, 
       filename = paste0(path, "output/plots/proportions_table_preB.png"))

# Save the data frame as RDS for later use
saveRDS(proportion_data_preB, 
        file = paste0(path, "output/proportions/proportion_data_preB.rds"))

cat("\n========================================\n")
cat("POST-BREEDING ANALYSIS COMPLETE\n")
cat("========================================\n")