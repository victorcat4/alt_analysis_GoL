# Gulf of Lion Bird Migration Analysis

Processing and analysis of GPS tracking data from multiple bird migration studies over the Gulf of Lion (Mediterranean Sea).

**Author:** Victor Cat

---

## Project Structure

```
project/
├── Analysis.R                          # Main analysis script
├── 00_setup.R                          # Package loading
├── 01_formatting_functions.R           # Data cleaning, column standards, trajectory assignment
├── 02_spatial_functions.R              # Spatial filtering and distance calculations
├── 03_data_readers.R                   # Readers for each data source
├── 04_visualization_functions.R        # Summary plots and tables
├── 05_processing_pipeline.R            # Processing orchestration (parallel + sequential)
├── 06_bayesian_modeling.R              # Gamma flight height model (nimble)
├── 07_plot_modeling.R                  # Model result visualizations
├── 08_process_honey_buzzard_data2.R    # (placeholder)
├── 09_process_Red_Kit_data.R           # (placeholder)
├── data/supp/
│   ├── europe.rds                      # European coastline
│   ├── GOL_extended.rds                # Extended study area polygon
│   ├── GOL_restricted.rds              # Restricted Gulf of Lion polygon
│   └── EC.xlsx                         # Species → ecological group mapping
└── output/
    ├── gps_migration.csv               # Final processed dataset
    ├── model/                          # Saved MCMC results
    ├── proportions/                    # Altitude proportion estimates
    └── plots/                          # All generated figures and tables
```

---

## Workflow

The analysis runs in five phases, controlled by flags in `Analysis.R`:

| Flag | Phase | Description |
|------|-------|-------------|
| `RUN_PROCESSING` | 1 | Process all data sources → `gps_migration.csv` |
| — | 2 | Summary statistics (always runs) |
| — | 3 | Load data, add ecological groups (always runs) |
| — | 4 | Generate visualizations (always runs) |
| `RUN_MODELING` | 5a | Bayesian model fitting (computationally intensive) |
| `RUN_MODEL_PLOTS` | 5b | Generate plots from model results |

### Phase 1: Data Processing

Processes five data sources, each individual going through an 11-step pipeline:

1. `clean_row()` — Remove invalid coordinates/timestamps, deduplicate
2. `first_loc()` — Record equipment deployment location
3. `select_pts()` — Keep only individuals crossing the restricted Gulf of Lion area; retain points within the extended area
4. `downsample_tracks()` — Reduce to 10-minute intervals
5. `assign_trajectories()` — Split at July 1st boundaries and gaps > 5 days
6. `remove_short_trajectories()` — Drop trajectories with < 5 points
7. `add_season()` — Classify as "pre-breeding" (Jan–Jun) or "post-breeding" (Jul–Dec)
8. `check_restricted_area()` — Flag points inside Gulf of Lion (`in_GOL`)
9. `classify_staging_migration()` — DBSCAN clustering (eps=0.1°, minPts=10): clustered → "Staging", unclustered → "Migration"
10. `calculate_distance_to_shore()` — Distance to nearest coastline (meters, EPSG:3035)
11. Final coordinate validation

Movebank data is processed in parallel (future/furrr) with temporary files to avoid write conflicts. Non-Movebank sources are processed sequentially.

### Data Sources

| Source | Species | Reader |
|--------|---------|--------|
| Movebank (individual CSVs) | Multiple | `read_movebank_individual()` |
| Turtle Dove & Stone Curlew | *Streptopelia turtur*, *Burhinus oedicnemus* | `read_turtle_dove_stone_curlew()` |
| Honey Buzzard (Byholm) | *Pernis apivorus* | `read_honey_buzzard_byholm()` |
| BIRDMAP | Multiple | `read_birdmap()` |
| Montagu's Harrier (Schlaich et al.) | *Circus pygargus* | `read_montagu_harrier()` |

### Excluded Species

*Puffinus yelkouan*, *Calonectris diomedea*, *Thalasseus sandvicensis*, *Larus michahellis*, *Aquila clanga*, *Aquila chrysaetos*, *Ciconia ciconia*.

### Ecological Groups

Ardeid, Large wader, Larid, Shorebird, Terrestrial bird — mapped via `EC.xlsx`.

---

## Phase 4: Visualizations

Generated in `output/plots/`:

- Migration maps (both seasons, by ecological group)
- Migration vs Staging spatial map
- Flight height histograms (overall and by group)
- Day vs Night migration proportions (using suncalc sunrise/sunset)
- Species summary table (HTML, with per-season and per-state breakdowns)

---

## Phase 5: Bayesian Flight Height Model

### Model (Gamma with observation error)

```
obs_alt[i]  ~ Normal(true_alt[i], sigma_obs[i])
true_alt[i] ~ Gamma(shape[group], rate[i])

log(mu[i])        = alpha_0[group] + alpha_1[group] * dist_shore[i]
rate[i]           = shape[group] / mu[i]
log(sigma_obs[i]) = beta_0 + beta_1 * hdop[i] + beta_2 * sat_count[i]
```

Fixed effects per ecological group for intercept (`alpha_0`), distance slope (`alpha_1`), and shape. Observation error depends on GPS quality (HDOP, satellite count).

### MCMC Settings (defaults in `Analysis.R`)

- 1,000,000 iterations, 90% burn-in, thin=20, 4 chains
- Automatic convergence filtering via Gelman-Rubin diagnostic

### Model Outputs

- Flight height density plots per group (color-coded altitude bands with proportion CIs)
- Proportion summary tables (0–20 m, 20–200 m, 20–300 m, >300 m)
- Quantile summary tables (Q25, Q50, Q75 with 95% CIs)
- Predicted proportion below 300 m as a function of distance to shore
- MCMC trace plots

---

## Output Columns (`gps_migration.csv`)

| Column | Description |
|--------|-------------|
| `individual_id` | Unique bird identifier |
| `study_id` | Source study |
| `individual_taxon_canonical_name` | Scientific name |
| `location_long`, `location_lat` | Coordinates (WGS84) |
| `timestamp` | Date/time (UTC) |
| `height_above_msl` | Flight height (m) |
| `gps_hdop`, `gps_vdop`, `gps_satellite_count` | GPS quality metrics |
| `equip_loc` | First recorded location |
| `trajectory` | Trajectory ID (`{id}_traj_{n}`) |
| `season` | "pre-breeding" or "post-breeding" |
| `in_GOL` | Point inside restricted Gulf of Lion area |
| `state` | "Migration" or "Staging" |
| `db_cluster` | DBSCAN cluster ID (0 = migration) |
| `distance_to_shore` | Distance to coastline (meters) |

---

## Required Packages

**Core:** tidyverse, sf, readxl, lubridate, dbscan, scales, suncalc, future, furrr, progressr

**Modeling (Phase 5):** nimble, coda, MCMCvis

**Tables/Plots:** gt, cowplot, png, grid

---

## Usage

```r
# Full run (edit flags in Analysis.R first)
source("Analysis.R")
```

Set `RUN_PROCESSING = TRUE` for first run (slow). On subsequent runs, set it to `FALSE` to skip processing and work from the saved `gps_migration.csv`.